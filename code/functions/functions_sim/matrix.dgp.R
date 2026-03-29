# =============================================================================
# Matrix DGP for Monte Carlo nowcasting experiments
# -----------------------------------------------------------------------------
# Goal:
#   Compare nowcasts at M1, M2, M3 of the current quarter between:
#   (i)  Matrix MF-TPRF
#   (ii) Vector MF-TPRF
#
# Main changes relative to the previous version:
#   1. No full-sample standardization anymore
#   2. Vintage-safe standardization only, using information available up to
#      the current vintage
#   3. Explicit storage of the true predictive factor f_true
#   4. Oracle objects for the infeasible best benchmark (BI) at M1/M2/M3
#   5. Slightly safer loading generation with rank checks
#
# Baseline design:
#   Tm in {100, 200, 400}
#   p1 in {10, 20}
#   p2 in {25, 50, 100}
#   r1 = r2 = 1
#   v1 = 1, v2 = 4
#   Sigma_G fixed = diag(1.25, 1.75, 2.25, 2.75)
# =============================================================================


# =============================================================================
# Function: make_dgp_params
# Purpose : Create a clean parameter list for the baseline DGP
# =============================================================================
make_dgp_params <- function(
    Tm = 100,
    burn = 400,
    draws = 500,
    seed = 1712,
    p1 = 10,
    p2 = 25,
    r1 = 1,
    r2 = 1,
    v1 = 1,
    v2 = 4,
    rhoF = 0.9,
    rhoG = 0.3,
    rhoE = 0.0,
    R2_BI = 0.5,
    sigma_z = 0.25,
    share_monthly = 0.8,
    share_m_flow = 0.6,
    share_q_flow = 0.5,
    delta_r = 0.3,
    delta_c = 0.3,
    quarterly_obs_month = 3,
    ragged_edge = TRUE,
    target_quarter = NULL,
    bi_num_month_lags = 5
) {
  stopifnot(Tm > 10, burn >= 0, draws >= 1)
  stopifnot(p1 >= 1, p2 >= 5)
  stopifnot(r1 == 1, r2 == 1)
  stopifnot(v1 == 1, v2 == 4)
  stopifnot(abs(rhoF) < 1, abs(rhoG) < 1, abs(rhoE) < 1)
  stopifnot(R2_BI > 0, R2_BI < 1)
  stopifnot(share_monthly >= 0, share_monthly <= 1)
  stopifnot(share_m_flow >= 0, share_m_flow <= 1)
  stopifnot(share_q_flow >= 0, share_q_flow <= 1)
  stopifnot(quarterly_obs_month %in% c(1, 2, 3))
  stopifnot(bi_num_month_lags >= 1)
  
  Nm <- round(p2 * share_monthly)
  Nm <- max(1, min(Nm, p2 - 1))
  Nq <- p2 - Nm
  
  n_m_flow  <- round(Nm * share_m_flow)
  n_m_flow  <- max(0, min(n_m_flow, Nm))
  n_m_stock <- Nm - n_m_flow
  
  n_q_flow  <- round(Nq * share_q_flow)
  n_q_flow  <- max(0, min(n_q_flow, Nq))
  n_q_stock <- Nq - n_q_flow
  
  stopifnot(n_m_flow + n_m_stock == Nm)
  stopifnot(n_q_flow + n_q_stock == Nq)
  stopifnot(Nm + Nq == p2)
  
  var_type <- c(
    rep("monthly_flow",    n_m_flow),
    rep("monthly_stock",   n_m_stock),
    rep("quarterly_flow",  n_q_flow),
    rep("quarterly_stock", n_q_stock)
  )
  
  # ---------------------------------------------------------------------------
  # Monthly release delays:
  #   split monthly variables into 3 groups with delays 0,1,2
  # ---------------------------------------------------------------------------
  monthly_release_delay <- c(
    rep(0, ceiling(Nm / 3)),
    rep(1, ceiling(Nm / 3)),
    rep(2, Nm - 2 * ceiling(Nm / 3))
  )
  
  # ---------------------------------------------------------------------------
  # Quarterly variables:
  #   by default released in the month quarterly_obs_month of each quarter
  #   with additional variable-specific delay equal to 0
  # ---------------------------------------------------------------------------
  quarterly_release_delay <- rep(0, Nq)
  
  list(
    Tm = Tm,
    burn = burn,
    draws = draws,
    seed = seed,
    
    p1 = p1,
    p2 = p2,
    Nm = Nm,
    Nq = Nq,
    
    n_m_flow = n_m_flow,
    n_m_stock = n_m_stock,
    n_q_flow = n_q_flow,
    n_q_stock = n_q_stock,
    var_type = var_type,
    
    r1 = r1,
    r2 = r2,
    v1 = v1,
    v2 = v2,
    
    rhoF = rhoF,
    rhoG = rhoG,
    rhoE = rhoE,
    
    R2_BI = R2_BI,
    sigma_y = NULL,
    sigma_z = sigma_z,
    w_ea = NULL,
    
    delta_r = delta_r,
    delta_c = delta_c,
    
    monthly_release_delay = monthly_release_delay,
    quarterly_release_delay = quarterly_release_delay,
    quarterly_obs_month = quarterly_obs_month,
    ragged_edge = ragged_edge,
    target_quarter = target_quarter,
    
    bi_num_month_lags = bi_num_month_lags
  )
}


# =============================================================================
# Function: make_tridiag_cov
# Purpose : Build a symmetric tridiagonal covariance matrix
#           with unit diagonal and constant first off-diagonal
# =============================================================================
make_tridiag_cov <- function(n, offdiag = 0.3, diag_scale = 1) {
  stopifnot(n >= 1)
  stopifnot(abs(offdiag) < diag_scale)
  
  S <- diag(diag_scale, n)
  
  if (n >= 2) {
    for (i in 1:(n - 1)) {
      S[i, i + 1] <- offdiag
      S[i + 1, i] <- offdiag
    }
  }
  
  S + 1e-8 * diag(n)
}


# =============================================================================
# Function: rmatrixnorm
# Purpose : Draw from a matrix normal distribution with separable covariance
#           vec(X) ~ N(vec(mean), Sigma_c kron Sigma_r)
# =============================================================================
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


# =============================================================================
# Function: orthonormal_scaled
# Purpose : Build a matrix with orthogonal columns scaled so that
#           Q'Q = scale * I
# Notes   : Includes basic rank checks
# =============================================================================
orthonormal_scaled <- function(M, scale, tol = 1e-10) {
  qrM <- qr(M, tol = tol)
  if (qrM$rank < ncol(M)) {
    stop("Input matrix is rank deficient in orthonormal_scaled().")
  }
  Q <- qr.Q(qrM)
  Q <- Q[, 1:ncol(M), drop = FALSE]
  sqrt(scale) * Q
}


# =============================================================================
# Function: orthogonalize_against
# Purpose : Remove from M the projection on the column space of B
# =============================================================================
orthogonalize_against <- function(M, B) {
  if (ncol(B) == 0) return(M)
  P_B <- B %*% solve(t(B) %*% B) %*% t(B)
  M - P_B %*% M
}


# =============================================================================
# Function: build_random_loading_block
# Purpose : Generate a Gaussian random loading block, optionally orthogonalized
#           against a reference loading space, and scaled so that Q'Q = scale I
# =============================================================================
build_random_loading_block <- function(n_rows, n_cols, scale, B_ref = NULL,
                                       max_tries = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  for (iter in 1:max_tries) {
    M_raw <- matrix(rnorm(n_rows * n_cols), n_rows, n_cols)
    
    if (!is.null(B_ref)) {
      M_raw <- orthogonalize_against(M_raw, B_ref)
    }
    
    qrM <- qr(M_raw)
    if (qrM$rank == n_cols) {
      return(orthonormal_scaled(M_raw, scale = scale))
    }
  }
  
  stop("Unable to generate a full-rank loading block after max_tries attempts.")
}


# =============================================================================
# Function: mariano_murasawa_agg
# Purpose : Aggregate a monthly latent target into a quarterly series using
#           the Mariano-Murasawa approximation
# =============================================================================
mariano_murasawa_agg <- function(y_monthly) {
  T_months   <- nrow(y_monthly)
  n_units    <- ncol(y_monthly)
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


# =============================================================================
# Function: vectorize_panel
# Purpose : Vectorize a T x p1 x p2 array into a T x (p1*p2) matrix
# =============================================================================
vectorize_panel <- function(X_array) {
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


# =============================================================================
# Function: aggregate_lf_panel
# Purpose : Convert monthly predictors into a low-frequency quarterly panel
#           - stock variables: quarterly average
#           - flow variables : quarterly sum
# Notes   : Can work with NA values
# =============================================================================
aggregate_lf_panel <- function(X_monthly_full, var_type) {
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
        X_lf[q, , j] <- colMeans(block, na.rm = TRUE)
      } else if (var_type[j] %in% c("monthly_flow", "quarterly_flow")) {
        X_lf[q, , j] <- colSums(block, na.rm = TRUE)
      } else {
        stop("Unknown variable type in var_type.")
      }
      
      # -----------------------------------------------------------------------
      # If all observations are NA in the block, colMeans/colSums may return
      # NaN. Convert those back to NA.
      # -----------------------------------------------------------------------
      X_lf[q, , j][is.nan(X_lf[q, , j])] <- NA_real_
    }
  }
  
  X_lf
}


# =============================================================================
# Function: get_target_quarter_months
# Purpose : Return the three monthly indices corresponding to target quarter q*
# =============================================================================
get_target_quarter_months <- function(q_star) {
  c(3 * q_star - 2, 3 * q_star - 1, 3 * q_star)
}


# =============================================================================
# Function: build_base_mf_panel
# Purpose : Create a baseline mixed-frequency monthly panel
#           - monthly variables observed every month
#           - quarterly variables observed only in their release month
# Notes   : This function works on the raw monthly panel
# =============================================================================
build_base_mf_panel <- function(X_monthly_full_raw, params) {
  Tm <- dim(X_monthly_full_raw)[1]
  p2 <- dim(X_monthly_full_raw)[3]
  Tq <- floor(Tm / 3)
  
  Nm <- params$Nm
  Nq <- params$Nq
  quarterly_obs_month <- params$quarterly_obs_month
  
  stopifnot(p2 == Nm + Nq)
  stopifnot(quarterly_obs_month %in% c(1, 2, 3))
  
  X_mf <- X_monthly_full_raw
  
  # ---------------------------------------------------------------------------
  # Quarterly variables are only observed in their release month within each
  # quarter; outside that month they are set to NA.
  # ---------------------------------------------------------------------------
  if (Nq > 0) {
    quarterly_cols <- (Nm + 1):p2
    X_mf[, , quarterly_cols] <- NA_real_
    
    for (q in 1:Tq) {
      rel_m <- 3 * q - (3 - quarterly_obs_month)
      X_mf[rel_m, , quarterly_cols] <- X_monthly_full_raw[rel_m, , quarterly_cols]
    }
  }
  
  X_mf
}


# =============================================================================
# Function: apply_vintage_mask
# Purpose : Build the information set available at a given nowcast vintage
#           M1, M2, M3 of target quarter q*
# Notes   : This applies release-delay masking only within the target quarter
# =============================================================================
apply_vintage_mask <- function(X_base_mf_raw, params, q_star,
                               vintage = c("M1", "M2", "M3")) {
  vintage <- match.arg(vintage)
  
  Nm <- params$Nm
  Nq <- params$Nq
  
  m_target <- get_target_quarter_months(q_star)
  m1 <- m_target[1]
  
  current_month_observed <- switch(
    vintage,
    M1 = m_target[1],
    M2 = m_target[2],
    M3 = m_target[3]
  )
  
  X_obs <- X_base_mf_raw[1:current_month_observed, , , drop = FALSE]
  m_target_in_sample <- m_target[m_target <= current_month_observed]
  
  # ---------------------------------------------------------------------------
  # Monthly variables:
  #   some of the target-quarter months are not yet available because of
  #   publication delays
  # ---------------------------------------------------------------------------
  if (isTRUE(params$ragged_edge) && Nm > 0 && length(m_target_in_sample) > 0) {
    delays_m <- params$monthly_release_delay
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
  
  # ---------------------------------------------------------------------------
  # Quarterly variables:
  #   even within the target quarter, they are observed only in their release
  #   month if that month has already occurred at the vintage
  # ---------------------------------------------------------------------------
  if (Nq > 0 && length(m_target_in_sample) > 0) {
    delays_q <- params$quarterly_release_delay
    stopifnot(length(delays_q) == Nq)
    
    quarterly_obs_month <- params$quarterly_obs_month
    quarterly_cols <- (Nm + 1):(Nm + Nq)
    
    release_month_target_q <- m1 + (quarterly_obs_month - 1)
    
    for (k in 1:Nq) {
      j <- quarterly_cols[k]
      release_month_j <- release_month_target_q + delays_q[k]
      
      # remove target-quarter entries first
      for (mt in m_target_in_sample) {
        X_obs[mt, , j] <- NA_real_
      }
      
      # reveal only if release has occurred
      if (release_month_j <= current_month_observed) {
        X_obs[release_month_j, , j] <- X_base_mf_raw[release_month_j, , j]
      }
    }
  }
  
  obs_mask <- array(1L, dim = dim(X_obs))
  obs_mask[is.na(X_obs)] <- 0L
  
  list(
    X_mf_raw = X_obs,
    obs_mask = obs_mask,
    current_month = current_month_observed
  )
}


# =============================================================================
# Function: fit_standardization_tensor
# Purpose : Fit cell-by-cell standardization parameters on a tensor
# Notes   : Uses only the data available in the input tensor, with NA support
# =============================================================================
fit_standardization_tensor <- function(X_tens, center = TRUE, scale = TRUE) {
  stopifnot(length(dim(X_tens)) == 3)
  
  p1 <- dim(X_tens)[2]
  p2 <- dim(X_tens)[3]
  
  mu   <- matrix(0, nrow = p1, ncol = p2)
  sdev <- matrix(1, nrow = p1, ncol = p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      x_ij <- X_tens[, i, j]
      x_obs <- x_ij[!is.na(x_ij)]
      
      mu_ij <- if (center && length(x_obs) > 0) mean(x_obs) else 0
      sd_ij <- if (scale  && length(x_obs) > 1) stats::sd(x_obs) else 1
      
      if (is.na(sd_ij) || sd_ij <= 1e-12) sd_ij <- 1
      
      mu[i, j]   <- mu_ij
      sdev[i, j] <- sd_ij
    }
  }
  
  list(
    center = mu,
    scale  = sdev
  )
}


# =============================================================================
# Function: apply_standardization_tensor
# Purpose : Apply cell-by-cell standardization to a tensor using pre-fitted
#           means and standard deviations
# =============================================================================
apply_standardization_tensor <- function(X_tens, center_mat, scale_mat) {
  stopifnot(length(dim(X_tens)) == 3)
  stopifnot(all(dim(center_mat) == dim(scale_mat)))
  stopifnot(dim(X_tens)[2] == nrow(center_mat))
  stopifnot(dim(X_tens)[3] == ncol(center_mat))
  
  Tt <- dim(X_tens)[1]
  p1 <- dim(X_tens)[2]
  p2 <- dim(X_tens)[3]
  
  X_std <- array(NA_real_, dim = dim(X_tens))
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      X_std[, i, j] <- (X_tens[, i, j] - center_mat[i, j]) / scale_mat[i, j]
    }
  }
  
  X_std
}


# =============================================================================
# Function: build_bi_monthly_design
# Purpose : Build a quarterly design matrix for the infeasible best benchmark
#           based on the true monthly predictive factor f_true
# Notes   : This is the oracle counterpart of a mixed-frequency regression:
#           y_q = alpha + beta' * [f_t, f_{t-1}, ..., f_{t-L+1}] + error
# =============================================================================
build_bi_monthly_design <- function(f_true, Tq, num_lags = 5) {
  stopifnot(length(f_true) >= 3 * Tq)
  stopifnot(num_lags >= 1)
  
  X_bi <- matrix(NA_real_, nrow = Tq, ncol = num_lags)
  colnames(X_bi) <- paste0("f_lag", 0:(num_lags - 1))
  
  for (q in 1:Tq) {
    t_end_q <- 3 * q
    idx <- t_end_q:(t_end_q - num_lags + 1)
    
    if (min(idx) >= 1) {
      X_bi[q, ] <- f_true[idx]
    }
  }
  
  X_bi
}


# =============================================================================
# Function: build_bi_vintage_row
# Purpose : Build the oracle BI regressor row available at a given vintage
#           for the target quarter q_star
# Notes   : Uses the true monthly factor only up to the current vintage month
#           and pads unavailable target-quarter months with NA
# =============================================================================
build_bi_vintage_row <- function(f_true, q_star, current_month, num_lags = 5) {
  stopifnot(num_lags >= 1)
  
  t_ref <- current_month
  idx <- t_ref:(t_ref - num_lags + 1)
  
  out <- rep(NA_real_, num_lags)
  valid <- idx >= 1
  out[valid] <- f_true[idx[valid]]
  
  names(out) <- paste0("f_lag", 0:(num_lags - 1))
  out
}


# =============================================================================
# Function: simulate_matrix_dgp_truth
# Purpose : Simulate the full DGP and return the latent truth
# Notes   : No standardization is done here
# =============================================================================
simulate_matrix_dgp_truth <- function(params) {
  stopifnot(is.list(params))
  if (!is.null(params$seed)) set.seed(params$seed)
  
  # ---------------------------------------------------------------------------
  # 1. Read parameters
  # ---------------------------------------------------------------------------
  T_months <- params$Tm
  burn_in  <- params$burn
  T_total  <- T_months + burn_in
  
  n_units <- params$p1
  n_vars  <- params$p2
  
  Nm <- params$Nm
  Nq <- params$Nq
  stopifnot(n_vars == Nm + Nq)
  
  r_row <- params$r1
  r_col <- params$r2
  g_row <- params$v1
  g_col <- params$v2
  
  stopifnot(r_row == 1, r_col == 1)
  stopifnot(g_row == 1, g_col == 4)
  
  rho_F <- params$rhoF
  rho_G <- params$rhoG
  rho_E <- params$rhoE
  
  target_R2 <- params$R2_BI
  sigma_z   <- params$sigma_z
  delta_r   <- params$delta_r
  delta_c   <- params$delta_c
  
  stopifnot(abs(rho_F) < 1, abs(rho_G) < 1, abs(rho_E) < 1)
  stopifnot(target_R2 > 0, target_R2 < 1)
  
  proxy_weights <- params$w_ea
  if (is.null(proxy_weights)) {
    proxy_weights <- rep(1 / n_units, n_units)
  }
  stopifnot(length(proxy_weights) == n_units)
  stopifnot(abs(sum(proxy_weights) - 1) < 1e-8)
  
  # ---------------------------------------------------------------------------
  # 2. Generate matrix loadings
  #    These are the matrix-model analogues of vector loadings Phi
  # ---------------------------------------------------------------------------
  R_signal <- build_random_loading_block(
    n_rows = n_units,
    n_cols = r_row,
    scale  = n_units
  )
  
  C_signal <- build_random_loading_block(
    n_rows = n_vars,
    n_cols = r_col,
    scale  = n_vars
  )
  
  R_noise <- build_random_loading_block(
    n_rows = n_units,
    n_cols = g_row,
    scale  = n_units,
    B_ref  = R_signal
  )
  
  C_noise <- build_random_loading_block(
    n_rows = n_vars,
    n_cols = g_col,
    scale  = n_vars,
    B_ref  = C_signal
  )
  
  # ---------------------------------------------------------------------------
  # 3. Factor dynamics
  #    Predictive factor block F_t and non-predictive block G_t
  # ---------------------------------------------------------------------------
  A_F <- sqrt(rho_F) * diag(r_row)
  B_F <- sqrt(rho_F) * diag(r_col)
  
  A_G <- sqrt(rho_G) * diag(g_row)
  B_G <- sqrt(rho_G) * diag(g_col)
  
  g_dim <- g_row * g_col
  stopifnot(g_dim == 4)
  
  Sigma_G <- diag(c(1.25, 1.75, 2.25, 2.75))
  
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
  
  # ---------------------------------------------------------------------------
  # 4. Idiosyncratic matrix disturbance
  # ---------------------------------------------------------------------------
  A_E <- sqrt(rho_E) * diag(n_units)
  B_E <- sqrt(rho_E) * diag(n_vars)
  
  Sigma_row <- make_tridiag_cov(n_units, offdiag = delta_r, diag_scale = 1)
  Sigma_col <- make_tridiag_cov(n_vars,  offdiag = delta_c, diag_scale = 1)
  
  E_idio <- array(0, dim = c(T_total, n_units, n_vars))
  E_idio[1, , ] <- matrix(0, n_units, n_vars)
  
  for (t in 2:T_total) {
    E_prev  <- matrix(E_idio[t - 1, , ], n_units, n_vars)
    E_innov <- rmatrixnorm(1, Sigma_r = Sigma_row, Sigma_c = Sigma_col)
    
    E_idio[t, , ] <- A_E %*% E_prev %*% t(B_E) + E_innov
  }
  
  # ---------------------------------------------------------------------------
  # 5. Raw monthly predictor tensor
  # ---------------------------------------------------------------------------
  X_all_raw <- array(0, dim = c(T_total, n_units, n_vars))
  
  for (t in 1:T_total) {
    Ft <- matrix(F_signal[t, , ], r_row, r_col)
    Gt <- matrix(G_noise[t, , ],  g_row, g_col)
    Et <- matrix(E_idio[t, , ],   n_units, n_vars)
    
    signal_part <- R_signal %*% Ft %*% t(C_signal)
    noise_part  <- R_noise  %*% Gt %*% t(C_noise)
    
    X_all_raw[t, , ] <- signal_part + noise_part + Et
  }
  
  # ---------------------------------------------------------------------------
  # 6. True scalar predictive factor f_true
  #    Since r1 = r2 = 1, the predictive matrix factor reduces to a scalar
  # ---------------------------------------------------------------------------
  f_true_all <- as.numeric(F_signal[, 1, 1])
  
  # ---------------------------------------------------------------------------
  # 7. Monthly latent target
  #    Each unit shares the same common predictive signal f_true_t
  # ---------------------------------------------------------------------------
  common_signal <- matrix(f_true_all, nrow = T_total, ncol = n_units)
  
  sigma_y <- params$sigma_y
  if (is.null(sigma_y)) {
    signal_var <- stats::var(f_true_all)
    sigma_y <- sqrt(signal_var * (1 - target_R2) / target_R2)
  }
  
  y_monthly_signal_all <- common_signal
  y_monthly_all <- common_signal + matrix(
    rnorm(T_total * n_units, sd = sigma_y),
    T_total, n_units
  )
  
  # ---------------------------------------------------------------------------
  # 8. Quarterly target and proxy
  # ---------------------------------------------------------------------------
  y_quarterly_all        <- mariano_murasawa_agg(y_monthly_all)
  y_quarterly_signal_all <- mariano_murasawa_agg(y_monthly_signal_all)
  
  T_quarters_total <- nrow(y_quarterly_all)
  
  z_quarterly_signal_all <- as.numeric(y_quarterly_signal_all %*% proxy_weights)
  z_quarterly_all <- z_quarterly_signal_all + rnorm(T_quarters_total, sd = sigma_z)
  
  # ---------------------------------------------------------------------------
  # 9. Drop burn-in
  # ---------------------------------------------------------------------------
  keep_months <- (burn_in + 1):T_total
  
  X_monthly_raw    <- X_all_raw[keep_months, , , drop = FALSE]
  y_monthly        <- y_monthly_all[keep_months, , drop = FALSE]
  y_monthly_signal <- y_monthly_signal_all[keep_months, , drop = FALSE]
  f_true           <- f_true_all[keep_months]
  
  F_signal_kept <- F_signal[keep_months, , , drop = FALSE]
  G_noise_kept  <- G_noise[keep_months, , , drop = FALSE]
  
  T_quarters <- floor(T_months / 3)
  keep_quarters <- (T_quarters_total - T_quarters + 1):T_quarters_total
  
  y_quarterly        <- y_quarterly_all[keep_quarters, , drop = FALSE]
  y_quarterly_signal <- y_quarterly_signal_all[keep_quarters, , drop = FALSE]
  z_quarterly        <- z_quarterly_all[keep_quarters]
  z_quarterly_signal <- z_quarterly_signal_all[keep_quarters]
  
  # ---------------------------------------------------------------------------
  # 10. Full raw low-frequency panel (not standardized yet)
  # ---------------------------------------------------------------------------
  X_quarterly_lf_raw <- aggregate_lf_panel(X_monthly_raw, params$var_type)
  
  list(
    params = modifyList(params, list(sigma_y = sigma_y)),
    dims = list(
      T_months = T_months,
      T_quarters = T_quarters,
      n_units = n_units,
      n_vars = n_vars,
      n_monthly_vars = Nm,
      n_quarterly_vars = Nq
    ),
    loadings = list(
      R_signal = R_signal,
      C_signal = C_signal,
      R_noise  = R_noise,
      C_noise  = C_noise
    ),
    factors = list(
      F_signal = F_signal_kept,
      G_noise  = G_noise_kept,
      f_true   = f_true
    ),
    covariances = list(
      Sigma_G   = Sigma_G,
      Sigma_row = Sigma_row,
      Sigma_col = Sigma_col
    ),
    truth = list(
      X_monthly_raw      = X_monthly_raw,
      X_quarterly_lf_raw = X_quarterly_lf_raw,
      y_monthly          = y_monthly,
      y_quarterly        = y_quarterly,
      z_quarterly        = z_quarterly,
      y_monthly_signal   = y_monthly_signal,
      y_quarterly_signal = y_quarterly_signal,
      z_quarterly_signal = z_quarterly_signal
    )
  )
}


# =============================================================================
# Function: build_info_sets
# Purpose : Build nowcast information sets M1, M2, M3 for the target quarter
#           using vintage-safe standardization and storing BI oracle objects
# =============================================================================
build_info_sets <- function(truth_obj, params) {
  X_monthly_full_raw <- truth_obj$truth$X_monthly_raw
  X_lf_full_raw      <- truth_obj$truth$X_quarterly_lf_raw
  y_q                <- truth_obj$truth$y_quarterly
  z_q                <- truth_obj$truth$z_quarterly
  f_true             <- truth_obj$factors$f_true
  
  Tm <- dim(X_monthly_full_raw)[1]
  Tq <- dim(X_lf_full_raw)[1]
  
  q_star <- params$target_quarter
  if (is.null(q_star)) q_star <- Tq
  
  stopifnot(q_star >= 2, q_star <= Tq)
  
  target_months <- get_target_quarter_months(q_star)
  stopifnot(max(target_months) <= Tm)
  
  # ---------------------------------------------------------------------------
  # Raw mixed-frequency monthly panel with quarterly variables observed only at
  # release dates
  # ---------------------------------------------------------------------------
  X_base_mf_raw <- build_base_mf_panel(X_monthly_full_raw, params)
  
  make_nowcast_set <- function(vintage_name) {
    # -------------------------------------------------------------------------
    # 1. Apply the vintage mask to the raw mixed-frequency panel
    # -------------------------------------------------------------------------
    tmp <- apply_vintage_mask(X_base_mf_raw, params, q_star = q_star, vintage = vintage_name)
    
    current_month <- tmp$current_month
    current_quarter_train <- q_star - 1
    
    # -------------------------------------------------------------------------
    # 2. Training samples available at this vintage
    #    - High-frequency sample: all months up to current_month
    #    - Low-frequency sample : all completed quarters up to q_star - 1
    # -------------------------------------------------------------------------
    X_hf_train_raw <- tmp$X_mf_raw
    X_lf_train_raw <- X_lf_full_raw[1:current_quarter_train, , , drop = FALSE]
    
    # -------------------------------------------------------------------------
    # 3. Vintage-safe standardization:
    #    fit means and sds only on the raw information available at this vintage
    # -------------------------------------------------------------------------
    std_obj_hf <- fit_standardization_tensor(X_hf_train_raw)
    X_hf_train_std <- apply_standardization_tensor(
      X_hf_train_raw,
      center_mat = std_obj_hf$center,
      scale_mat  = std_obj_hf$scale
    )
    
    # -------------------------------------------------------------------------
    # 4. Apply the SAME monthly standardization parameters to the low-frequency
    #    quarterly panel constructed from raw monthly data
    #    This keeps the scale coherent across HF and LF objects
    # -------------------------------------------------------------------------
    X_lf_train_std <- apply_standardization_tensor(
      X_lf_train_raw,
      center_mat = std_obj_hf$center,
      scale_mat  = std_obj_hf$scale
    )
    
    # -------------------------------------------------------------------------
    # 5. Vectorized versions for the vector MF-TPRF benchmark
    # -------------------------------------------------------------------------
    X_hf_train_vec_std <- vectorize_panel(X_hf_train_std)
    X_lf_train_vec_std <- vectorize_panel(X_lf_train_std)
    
    # -------------------------------------------------------------------------
    # 6. Oracle objects:
    #    - true raw monthly tensor up to the vintage
    #    - true factor up to the vintage
    #    - BI quarterly design based on the true factor
    # -------------------------------------------------------------------------
    X_hf_true_raw <- X_monthly_full_raw[1:current_month, , , drop = FALSE]
    X_hf_true_vec_raw <- vectorize_panel(X_hf_true_raw)
    
    X_lf_true_raw <- X_lf_full_raw[1:current_quarter_train, , , drop = FALSE]
    X_lf_true_vec_raw <- vectorize_panel(X_lf_true_raw)
    
    f_true_up_to_vintage <- f_true[1:current_month]
    
    # BI training design based on the true monthly factor for completed quarters
    X_bi_train <- build_bi_monthly_design(
      f_true = f_true,
      Tq = current_quarter_train,
      num_lags = params$bi_num_month_lags
    )
    
    # BI target row for quarter q_star at this vintage
    x_bi_target <- build_bi_vintage_row(
      f_true = f_true,
      q_star = q_star,
      current_month = current_month,
      num_lags = params$bi_num_month_lags
    )
    
    # -------------------------------------------------------------------------
    # 7. Observed target and proxy objects
    # -------------------------------------------------------------------------
    y_train  <- y_q[1:current_quarter_train, , drop = FALSE]
    y_target <- y_q[q_star, , drop = FALSE]
    
    z_train  <- z_q[1:current_quarter_train]
    z_target <- z_q[q_star]
    
    list(
      matrix = list(
        X_mf_std = X_hf_train_std,
        X_lf_std = X_lf_train_std,
        obs_mask = tmp$obs_mask
      ),
      vector = list(
        X_mf_vec_std = X_hf_train_vec_std,
        X_lf_vec_std = X_lf_train_vec_std
      ),
      oracle = list(
        X_hf_true_raw     = X_hf_true_raw,
        X_hf_true_vec_raw = X_hf_true_vec_raw,
        X_lf_true_raw     = X_lf_true_raw,
        X_lf_true_vec_raw = X_lf_true_vec_raw,
        f_true_up_to_vintage = f_true_up_to_vintage,
        X_bi_train = X_bi_train,
        x_bi_target = x_bi_target
      ),
      target = list(
        y_train  = y_train,
        y_target = y_target,
        z_train  = z_train,
        z_target = z_target
      ),
      standardization = list(
        center = std_obj_hf$center,
        scale  = std_obj_hf$scale
      ),
      meta = list(
        q_star = q_star,
        current_month = current_month,
        current_quarter_train = current_quarter_train,
        target_months = target_months,
        vintage = vintage_name
      )
    )
  }
  
  list(
    M1 = make_nowcast_set("M1"),
    M2 = make_nowcast_set("M2"),
    M3 = make_nowcast_set("M3")
  )
}


# =============================================================================
# Function: simulate_matrix_dgp
# Purpose : Full wrapper returning truth + nowcast information sets
# =============================================================================
simulate_matrix_dgp <- function(params) {
  stopifnot(is.list(params))
  stopifnot(params$p2 == params$Nm + params$Nq)
  stopifnot(length(params$var_type) == params$p2)
  stopifnot(length(params$monthly_release_delay) == params$Nm)
  stopifnot(length(params$quarterly_release_delay) == params$Nq)
  
  truth_obj <- simulate_matrix_dgp_truth(params)
  info_sets <- build_info_sets(truth_obj, truth_obj$params)
  
  list(
    params = truth_obj$params,
    dims = truth_obj$dims,
    loadings = truth_obj$loadings,
    factors = truth_obj$factors,
    covariances = truth_obj$covariances,
    truth = truth_obj$truth,
    info_sets = info_sets
  )
}