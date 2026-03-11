# =============================================================================
# External helpers
# =============================================================================

make_dgp_params <- function(
    Tm = 102, burn = 100, draws = 500, seed = 1712,
    p1 = 5,
    n_m_flow = 10, n_m_stock = 10,
    n_q_flow = 3,  n_q_stock = 2
) {
  
  Nm <- n_m_flow + n_m_stock
  Nq <- n_q_flow + n_q_stock
  p2 <- Nm + Nq
  
  list(
    Tm = Tm,
    burn = burn,
    draws = draws,
    seed = seed,
    
    p1 = p1,
    p2 = p2,
    Nm = Nm,
    Nq = Nq,
    
    r1 = 1,
    r2 = 1,
    v1 = 2,
    v2 = 2,
    
    rhoF = 0.9,
    rhoG = 0.3,
    rhoE = 0.0,
    
    R2_BI = 0.5,
    sigma_y = NULL,
    sigma_z = 0.25,
    
    d_row = 0,
    d_col = 0,
    
    w_ea = NULL,
    
    var_type = c(
      rep("monthly_flow", n_m_flow),
      rep("monthly_stock", n_m_stock),
      rep("quarterly_flow", n_q_flow),
      rep("quarterly_stock", n_q_stock)
    ),
    
    monthly_release_delay = c(
      rep(0, ceiling(Nm / 3)),
      rep(1, ceiling(Nm / 3)),
      rep(2, Nm - 2 * ceiling(Nm / 3))
    ),
    
    quarterly_release_delay = rep(0, Nq),
    
    quarterly_obs_month = 3,
    ragged_edge = TRUE,
    target_quarter = NULL
  )
}

make_banded_cov <- function(n, offdiag = 0.3, diag_scale = 1) {
  S <- diag(diag_scale, n)
  if (n >= 2) {
    for (i in 1:(n - 1)) {
      S[i, i + 1] <- offdiag
      S[i + 1, i] <- offdiag
    }
  }
  S + 1e-8 * diag(n)
}

rmatrixnorm <- function(n = 1, mean = NULL, Sigma_r, Sigma_c) {
  n_rows <- nrow(Sigma_r)
  n_cols <- nrow(Sigma_c)
  
  if (is.null(mean)) {
    mean <- matrix(0, n_rows, n_cols)
  }
  
  chol_r <- chol(Sigma_r)
  chol_c <- chol(Sigma_c)
  
  out <- vector("list", n)
  for (k in 1:n) {
    Z <- matrix(rnorm(n_rows * n_cols), n_rows, n_cols)
    out[[k]] <- mean + chol_r %*% Z %*% t(chol_c)
  }
  
  if (n == 1) return(out[[1]])
  out
}

mariano_murasawa_agg <- function(y_monthly) {
  # y_monthly: T_months x n_units
  T_months <- nrow(y_monthly)
  n_units  <- ncol(y_monthly)
  T_quarters <- floor(T_months / 3)
  
  y_quarterly <- matrix(NA_real_, nrow = T_quarters, ncol = n_units)
  
  for (q in 1:T_quarters) {
    t <- 3 * q
    if (t >= 5) {
      y_quarterly[q, ] <- (
        y_monthly[t, ] +
          2 * y_monthly[t - 1, ] +
          3 * y_monthly[t - 2, ] +
          2 * y_monthly[t - 3, ] +
          y_monthly[t - 4, ]
      ) / 3
    }
  }
  
  y_quarterly
}

vectorize_panel <- function(X_array, col_prefix = "unit_var") {
  # X_array: T x p1 x p2
  Tt <- dim(X_array)[1]
  p1 <- dim(X_array)[2]
  p2 <- dim(X_array)[3]
  
  X_vec <- matrix(NA_real_, nrow = Tt, ncol = p1 * p2)
  col_names <- character(p1 * p2)
  
  k <- 1
  for (i in 1:p1) {
    for (j in 1:p2) {
      X_vec[, k] <- X_array[, i, j]
      col_names[k] <- paste0("unit", i, "_var", j)
      k <- k + 1
    }
  }
  
  colnames(X_vec) <- col_names
  X_vec
}

aggregate_lf_panel <- function(X_monthly_full, var_type) {
  # X_monthly_full: Tm x p1 x p2
  Tm <- dim(X_monthly_full)[1]
  p1 <- dim(X_monthly_full)[2]
  p2 <- dim(X_monthly_full)[3]
  Tq <- floor(Tm / 3)
  
  stopifnot(length(var_type) == p2)
  
  X_lf <- array(NA_real_, dim = c(Tq, p1, p2))
  
  for (q in 1:Tq) {
    idx <- (3 * q - 2):(3 * q)
    
    for (j in 1:p2) {
      block <- X_monthly_full[idx, , j, drop = FALSE]
      block <- block[, , 1]
      
      if (var_type[j] %in% c("monthly_stock", "quarterly_stock")) {
        X_lf[q, , j] <- colMeans(block)
      } else if (var_type[j] %in% c("monthly_flow", "quarterly_flow")) {
        X_lf[q, , j] <- colSums(block)
      } else {
        stop("Unknown variable type in var_type.")
      }
    }
  }
  
  X_lf
}

get_target_quarter_months <- function(q_star) {
  c(3 * q_star - 2, 3 * q_star - 1, 3 * q_star)
}

build_base_mf_panel <- function(X_monthly_full, params) {
  # Creates a baseline mixed-frequency panel:
  # - monthly variables observed every month
  # - quarterly variables observed only in quarter release month
  # then ragged-edge adjustments are applied later in the target quarter
  
  Tm <- dim(X_monthly_full)[1]
  p1 <- dim(X_monthly_full)[2]
  p2 <- dim(X_monthly_full)[3]
  Tq <- floor(Tm / 3)
  
  Nm <- params$Nm
  Nq <- params$Nq
  quarterly_obs_month <- if (is.null(params$quarterly_obs_month)) 3 else params$quarterly_obs_month
  
  stopifnot(p2 == Nm + Nq)
  stopifnot(quarterly_obs_month %in% c(1, 2, 3))
  
  X_mf <- X_monthly_full
  
  if (Nq > 0) {
    quarterly_cols <- (Nm + 1):p2
    X_mf[, , quarterly_cols] <- NA_real_
    
    for (q in 1:Tq) {
      rel_m <- 3 * q - (3 - quarterly_obs_month)
      X_mf[rel_m, , quarterly_cols] <- X_monthly_full[rel_m, , quarterly_cols]
    }
  }
  
  X_mf
}

apply_vintage_mask <- function(X_base_mf, params, q_star, vintage = c("M1", "M2", "M3")) {
  vintage <- match.arg(vintage)
  
  Tm <- dim(X_base_mf)[1]
  Nm <- params$Nm
  Nq <- params$Nq
  
  m_target <- get_target_quarter_months(q_star)
  m1 <- m_target[1]
  m2 <- m_target[2]
  m3 <- m_target[3]
  
  current_month_observed <- switch(
    vintage,
    M1 = m1,
    M2 = m2,
    M3 = m3
  )
  
  # Truncate sample to the current vintage month
  X_obs <- X_base_mf[1:current_month_observed, , , drop = FALSE]
  
  # Target-quarter months that are still inside the truncated sample
  m_target_in_sample <- m_target[m_target <= current_month_observed]
  
  # Apply ragged-edge to monthly variables in the target quarter
  if (isTRUE(params$ragged_edge) && Nm > 0 && length(m_target_in_sample) > 0) {
    delays_m <- params$monthly_release_delay
    if (is.null(delays_m)) delays_m <- rep(0, Nm)
    stopifnot(length(delays_m) == Nm)
    
    for (j in 1:Nm) {
      delay_j <- delays_m[j]
      latest_available_month_j <- current_month_observed - delay_j
      
      for (mt in m_target_in_sample) {
        if (mt > latest_available_month_j) {
          X_obs[mt, , j] <- NA_real_
        }
      }
    }
  }
  
  # Apply release delay to quarterly variables in target quarter
  if (Nq > 0 && length(m_target_in_sample) > 0) {
    delays_q <- params$quarterly_release_delay
    if (is.null(delays_q)) delays_q <- rep(0, Nq)
    stopifnot(length(delays_q) == Nq)
    
    quarterly_obs_month <- if (is.null(params$quarterly_obs_month)) 3 else params$quarterly_obs_month
    quarterly_cols <- (Nm + 1):(Nm + Nq)
    
    release_month_target_q <- m1 + (quarterly_obs_month - 1)
    
    for (k in 1:Nq) {
      j <- quarterly_cols[k]
      release_month_j <- release_month_target_q + delays_q[k]
      
      # Remove target-quarter quarterly observations inside truncated sample
      for (mt in m_target_in_sample) {
        X_obs[mt, , j] <- NA_real_
      }
      
      # Reinsert only if already released by the current vintage
      if (release_month_j <= current_month_observed) {
        X_obs[release_month_j, , j] <- X_base_mf[release_month_j, , j]
      }
    }
  }
  
  obs_mask <- array(1L, dim = dim(X_obs))
  obs_mask[is.na(X_obs)] <- 0L
  
  list(
    X_mf = X_obs,
    obs_mask = obs_mask,
    current_month = current_month_observed
  )
}

simulate_matrix_dgp_truth <- function(params) {
  stopifnot(is.list(params))
  if (!is.null(params$seed)) set.seed(params$seed)
  
  # ===========================================================================
  # 1. Parameters
  # ===========================================================================
  
  T_months <- params$Tm
  burn_in  <- if (is.null(params$burn)) 100 else params$burn
  T_total  <- T_months + burn_in
  
  n_units <- params$p1
  n_vars  <- params$p2
  
  n_monthly_vars   <- params$Nm
  n_quarterly_vars <- params$Nq
  
  stopifnot(n_vars == n_monthly_vars + n_quarterly_vars)
  
  r_row <- params$r1
  r_col <- params$r2
  
  g_row <- params$v1
  g_col <- params$v2
  
  rho_F <- params$rhoF
  rho_G <- params$rhoG
  rho_E <- params$rhoE
  
  stopifnot(abs(rho_F) < 1, abs(rho_G) < 1, abs(rho_E) < 1)
  
  target_R2 <- params$R2_BI
  stopifnot(target_R2 > 0, target_R2 < 1)
  
  sigma_z <- if (is.null(params$sigma_z)) 0.25 else params$sigma_z
  d_row   <- if (is.null(params$d_row)) 0.3 else params$d_row
  d_col   <- if (is.null(params$d_col)) 0.3 else params$d_col
  
  proxy_weights <- params$w_ea
  if (is.null(proxy_weights)) {
    proxy_weights <- rep(1 / n_units, n_units)
  }
  stopifnot(length(proxy_weights) == n_units)
  stopifnot(abs(sum(proxy_weights) - 1) < 1e-8)
  
  # ===========================================================================
  # 2. Loadings
  # ===========================================================================
  
  R_signal <- matrix(rnorm(n_units * r_row), n_units, r_row)
  C_signal <- matrix(rnorm(n_vars  * r_col), n_vars,  r_col)
  
  R_noise  <- matrix(rnorm(n_units * g_row), n_units, g_row)
  C_noise  <- matrix(rnorm(n_vars  * g_col), n_vars,  g_col)
  
  R_signal <- qr.Q(qr(R_signal))
  C_signal <- qr.Q(qr(C_signal))
  R_noise  <- qr.Q(qr(R_noise))
  C_noise  <- qr.Q(qr(C_noise))
  
  # ===========================================================================
  # 3. Factor dynamics
  # ===========================================================================
  
  A_F <- sqrt(rho_F) * diag(r_row)
  B_F <- sqrt(rho_F) * diag(r_col)
  
  A_G <- sqrt(rho_G) * diag(g_row)
  B_G <- sqrt(rho_G) * diag(g_col)
  
  g_dim <- g_row * g_col
  if (g_dim == 4) {
    Sigma_G <- diag(c(1.25, 1.75, 2.25, 2.75))
  } else {
    Sigma_G <- diag(seq(1.25, by = 0.5, length.out = g_dim))
  }
  
  F_signal <- array(0, dim = c(T_total, r_row, r_col))
  G_noise  <- array(0, dim = c(T_total, g_row, g_col))
  
  F_signal[1, , ] <- matrix(0, r_row, r_col)
  G_noise[1, , ]  <- matrix(0, g_row, g_col)
  
  for (t in 2:T_total) {
    F_prev <- matrix(F_signal[t - 1, , ], r_row, r_col)
    G_prev <- matrix(G_noise[t - 1, , ],  g_row, g_col)
    
    U_F_t <- matrix(rnorm(r_row * r_col), r_row, r_col)
    U_G_t <- matrix(
      MASS::mvrnorm(1, mu = rep(0, g_dim), Sigma = Sigma_G),
      g_row, g_col
    )
    
    F_signal[t, , ] <- A_F %*% F_prev %*% t(B_F) + U_F_t
    G_noise[t, , ]  <- A_G %*% G_prev %*% t(B_G) + U_G_t
  }
  
  # ===========================================================================
  # 4. Idiosyncratic component
  # ===========================================================================
  
  A_E <- sqrt(rho_E) * diag(n_units)
  B_E <- sqrt(rho_E) * diag(n_vars)
  
  Sigma_row <- make_banded_cov(n_units, offdiag = d_row, diag_scale = 1)
  Sigma_col <- make_banded_cov(n_vars,  offdiag = d_col, diag_scale = 1)
  
  E_idio <- array(0, dim = c(T_total, n_units, n_vars))
  E_idio[1, , ] <- matrix(0, n_units, n_vars)
  
  for (t in 2:T_total) {
    E_prev <- matrix(E_idio[t - 1, , ], n_units, n_vars)
    E_innov <- rmatrixnorm(1, Sigma_r = Sigma_row, Sigma_c = Sigma_col)
    E_idio[t, , ] <- A_E %*% E_prev %*% t(B_E) + E_innov
  }
  
  # ===========================================================================
  # 5. Predictor panel X_t
  # ===========================================================================
  
  X_all <- array(0, dim = c(T_total, n_units, n_vars))
  
  for (t in 1:T_total) {
    Ft <- matrix(F_signal[t, , ], r_row, r_col)
    Gt <- matrix(G_noise[t, , ],  g_row, g_col)
    Et <- matrix(E_idio[t, , ],   n_units, n_vars)
    
    signal_part <- R_signal %*% Ft %*% t(C_signal)
    noise_part  <- R_noise  %*% Gt %*% t(C_noise)
    
    X_all[t, , ] <- signal_part + noise_part + Et
  }
  
  # ===========================================================================
  # 6. Monthly target
  # ===========================================================================
  
  stopifnot(r_row == 1, r_col == 1)
  
  common_signal <- matrix(0, nrow = T_total, ncol = n_units)
  for (t in 1:T_total) {
    f_t <- F_signal[t, 1, 1]
    common_signal[t, ] <- rep(f_t, n_units)
  }
  
  sigma_y <- params$sigma_y
  if (is.null(sigma_y)) {
    signal_var <- var(common_signal[, 1])
    sigma_y <- sqrt(signal_var * (1 - target_R2) / target_R2)
  }
  
  y_monthly_signal_all <- common_signal
  y_monthly_all <- common_signal + matrix(
    rnorm(T_total * n_units, sd = sigma_y),
    T_total, n_units
  )
  
  # ===========================================================================
  # 7. Quarterly target and proxy
  # ===========================================================================
  
  y_quarterly_all        <- mariano_murasawa_agg(y_monthly_all)
  y_quarterly_signal_all <- mariano_murasawa_agg(y_monthly_signal_all)
  
  T_quarters_total <- nrow(y_quarterly_all)
  
  z_quarterly_all <- as.numeric(
    y_quarterly_all %*% proxy_weights + rnorm(T_quarters_total, sd = sigma_z)
  )
  
  # ===========================================================================
  # 8. Drop burn-in
  # ===========================================================================
  
  keep_months <- (burn_in + 1):T_total
  
  X_monthly_full   <- X_all[keep_months, , , drop = FALSE]
  y_monthly        <- y_monthly_all[keep_months, , drop = FALSE]
  y_monthly_signal <- y_monthly_signal_all[keep_months, , drop = FALSE]
  F_signal_kept    <- F_signal[keep_months, , , drop = FALSE]
  G_noise_kept     <- G_noise[keep_months, , , drop = FALSE]
  
  T_quarters <- floor(T_months / 3)
  keep_quarters <- (T_quarters_total - T_quarters + 1):T_quarters_total
  
  y_quarterly        <- y_quarterly_all[keep_quarters, , drop = FALSE]
  y_quarterly_signal <- y_quarterly_signal_all[keep_quarters, , drop = FALSE]
  z_quarterly        <- z_quarterly_all[keep_quarters]
  
  X_lf_full <- aggregate_lf_panel(X_monthly_full, params$var_type)
  X_lf_full_vec <- vectorize_panel(X_lf_full)
  
  list(
    params = params,
    dims = list(
      T_months = T_months,
      T_quarters = T_quarters,
      n_units = n_units,
      n_vars = n_vars,
      n_monthly_vars = n_monthly_vars,
      n_quarterly_vars = n_quarterly_vars
    ),
    loadings = list(
      R_signal = R_signal,
      C_signal = C_signal,
      R_noise  = R_noise,
      C_noise  = C_noise
    ),
    factors = list(
      F_signal = F_signal_kept,
      G_noise  = G_noise_kept
    ),
    covariances = list(
      Sigma_G   = Sigma_G,
      Sigma_row = Sigma_row,
      Sigma_col = Sigma_col
    ),
    truth = list(
      X_monthly          = X_monthly_full,
      X_quarterly_lf     = X_lf_full,
      X_quarterly_lf_vec = X_lf_full_vec,
      y_monthly          = y_monthly,
      y_quarterly        = y_quarterly,
      z_quarterly        = z_quarterly,
      y_monthly_signal   = y_monthly_signal,
      y_quarterly_signal = y_quarterly_signal
    )
  )
}

build_info_sets <- function(truth_obj, params) {
  X_monthly_full <- truth_obj$truth$X_monthly
  X_lf_full      <- truth_obj$truth$X_quarterly_lf
  y_q            <- truth_obj$truth$y_quarterly
  z_q            <- truth_obj$truth$z_quarterly
  
  Tm <- dim(X_monthly_full)[1]
  Tq <- dim(X_lf_full)[1]
  
  q_star <- params$target_quarter
  if (is.null(q_star)) q_star <- Tq
  
  stopifnot(q_star >= 2, q_star <= Tq)
  
  target_months <- get_target_quarter_months(q_star)
  stopifnot(max(target_months) <= Tm)
  
  # Base MF panel: quarterlies only in release month
  X_base_mf <- build_base_mf_panel(X_monthly_full, params)
  
  # Forecast set (quarterly balanced, no NA, train up to q_star - 1)
  forecast_set <- list(
    matrix = list(
      X_lf_train = X_lf_full[1:(q_star - 1), , , drop = FALSE],
      X_lf_target = X_lf_full[q_star, , , drop = FALSE]
    ),
    vector = list(
      X_lf_train_vec  = vectorize_panel(X_lf_full[1:(q_star - 1), , , drop = FALSE]),
      X_lf_target_vec = vectorize_panel(X_lf_full[q_star, , , drop = FALSE])
    ),
    target = list(
      y_train = y_q[1:(q_star - 1), , drop = FALSE],
      y_target = y_q[q_star, , drop = FALSE],
      z_train = z_q[1:(q_star - 1)],
      z_target = z_q[q_star]
    ),
    meta = list(
      q_star = q_star,
      vintage = "forecast"
    )
  )
  
  # Nowcast vintages M1, M2, M3
  make_nowcast_set <- function(vintage_name) {
    tmp <- apply_vintage_mask(X_base_mf, params, q_star = q_star, vintage = vintage_name)
    
    list(
      matrix = list(
        X_mf = tmp$X_mf,
        obs_mask = tmp$obs_mask
      ),
      vector = list(
        X_mf_vec = vectorize_panel(tmp$X_mf)
      ),
      target = list(
        y_train = y_q[1:(q_star - 1), , drop = FALSE],
        y_target = y_q[q_star, , drop = FALSE],
        z_train = z_q[1:(q_star - 1)],
        z_target = z_q[q_star]
      ),
      meta = list(
        q_star = q_star,
        current_month = tmp$current_month,
        target_months = target_months,
        vintage = vintage_name
      )
    )
  }
  
  list(
    forecast = forecast_set,
    M1 = make_nowcast_set("M1"),
    M2 = make_nowcast_set("M2"),
    M3 = make_nowcast_set("M3")
  )
}

simulate_matrix_dgp <- function(params) {
  # Basic checks
  stopifnot(is.list(params))
  stopifnot(params$p2 == params$Nm + params$Nq)
  stopifnot(length(params$var_type) == params$p2)
  stopifnot(length(params$monthly_release_delay) == params$Nm)
  stopifnot(length(params$quarterly_release_delay) == params$Nq)
  
  truth_obj <- simulate_matrix_dgp_truth(params)
  info_sets <- build_info_sets(truth_obj, params)
  
  list(
    params = params,
    dims = truth_obj$dims,
    loadings = truth_obj$loadings,
    factors = truth_obj$factors,
    covariances = truth_obj$covariances,
    truth = truth_obj$truth,
    info_sets = info_sets
  )
}