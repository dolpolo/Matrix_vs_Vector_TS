# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

unbalancedness_tensor <- function(X_full, Unb, current_t) {
  # X_full: monthly tensor [T_m x P x K]
  # Unb   : matrix [P x K] of publication delays in months
  # current_t: current monthly vintage index
  
  stopifnot(length(dim(X_full)) == 3)
  
  T_m <- dim(X_full)[1]
  P   <- dim(X_full)[2]
  K   <- dim(X_full)[3]
  
  stopifnot(T_m >= current_t)
  stopifnot(all(dim(Unb) == c(P, K)))
  
  X_cut <- X_full[1:current_t, , , drop = FALSE]
  
  for (p in seq_len(P)) {
    for (k in seq_len(K)) {
      delay_pk <- Unb[p, k]
      
      # Missing or negative delay means: keep the series unchanged
      if (is.na(delay_pk) || delay_pk < 0) next
      
      delay_pk  <- as.integer(delay_pk)
      release_t <- current_t - delay_pk
      
      if (release_t < 1) {
        # The series is not yet available at the current vintage
        X_cut[, p, k] <- NA_real_
      } else if (release_t < current_t) {
        # The series is available only up to release_t
        X_cut[(release_t + 1):current_t, p, k] <- NA_real_
      }
    }
  }
  
  X_cut
}


# ==============================================================================
# MONTH-IN-QUARTER INDEX
# ==============================================================================

compute_m_tr <- function(date_t, dates_q) {
  # Return:
  #   1 = first month after the latest published quarter
  #   2 = second month after the latest published quarter
  #   3 = third month after the latest published quarter
  
  published_quarters <- dates_q[dates_q < date_t]
  if (length(published_quarters) == 0) return(NA_integer_)
  
  last_q_date <- max(published_quarters)
  last_q_idx  <- which(dates_q == last_q_date)
  
  M1 <- dates_q[last_q_idx] %m+% months(1)
  M2 <- dates_q[last_q_idx] %m+% months(2)
  M3 <- dates_q[last_q_idx] %m+% months(3)
  
  if (date_t == M1) return(1L)
  if (date_t == M2) return(2L)
  if (date_t == M3) return(3L)
  
  return(NA_integer_)
}


# ==============================================================================
# VERBOSE PRINT HELPER
# ==============================================================================

print_if_verbose <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) cat(...)
}


# ==============================================================================
# BUILD REAL-TIME DATASET WITH FIXED IMPUTATION RANK
# ==============================================================================

build_realtime_dataset_fixed_rank <- function(
    current_t,
    X_full,
    Y_q_all_full,
    dates_m,
    dates_q,
    Unb,
    agg,
    N_m,
    N_q,
    r_impute
) {
  # Build the real-time dataset at a given vintage using a FIXED imputation rank.
  # Rank selection is not performed here.
  
  X_cut <- unbalancedness_tensor(
    X_full    = X_full,
    Unb       = Unb,
    current_t = current_t
  )
  
  W_cut <- array(
    as.integer(!is.na(X_cut)),
    dim      = dim(X_cut),
    dimnames = dimnames(X_cut)
  )
  
  out_std <- standardize_mat_with_na(X_cut)
  X_std   <- out_std$X_scaled
  
  # Fixed-rank Cen & Lam imputation
  out_cl <- mfm.cl(Y = X_std, W = W_cut, r = r_impute)
  
  X_imputed <- X_std
  X_imputed[W_cut == 0] <- out_cl$Y_hat[W_cut == 0]
  
  X_quarterly <- aggregate_tensor_to_quarterly(
    X_tens = X_imputed,
    agg    = agg,
    N_m    = N_m,
    N_q    = N_q
  )
  
  current_date <- dates_m[current_t]
  idx_pub <- which(dates_q < current_date)
  
  if (length(idx_pub) < 2) return(NULL)
  
  T_q_use <- min(max(idx_pub), dim(X_quarterly)[1], nrow(Y_q_all_full))
  if (T_q_use < 2) return(NULL)
  
  Y_cut <- Y_q_all_full[1:T_q_use, , drop = FALSE]
  X_lf  <- X_quarterly[1:T_q_use, , , drop = FALSE]
  X_hf  <- X_imputed
  
  list(
    current_t     = current_t,
    current_date  = current_date,
    T_q_use       = T_q_use,
    X_lf          = X_lf,
    X_hf          = X_hf,
    Y_cut         = Y_cut,
    W_cut         = W_cut,
    X_std         = X_std,
    X_imputed     = X_imputed
  )
}


# ==============================================================================
# SELECT OR ACCEPT HYPERPARAMETERS FOR ONE REGIME
# ==============================================================================

select_regime_hyperparameters <- function(
    calibration_t,
    regime_label,
    X_full,
    Y_q_all_full,
    proxy_name,
    proxy_mode = c("scalar", "multivariate"),
    params,
    dates_m,
    dates_q,
    Unb,
    agg,
    N_m,
    N_q,
    user_hyper = list(
      r_impute   = NULL,
      Lproxy     = NULL,
      p_AR       = NULL,
      L_midas    = NULL,
      r_targeted = NULL
    ),
    ils_maxit = 100,
    ils_tol   = 1e-8,
    verbose   = TRUE
) {
  
  proxy_mode <- match.arg(proxy_mode)
  
  print_if_verbose(
    "\n============================================================\n",
    "Selecting hyperparameters for regime: ", regime_label, "\n",
    "Calibration vintage: ", as.character(dates_m[calibration_t]), "\n",
    "Proxy mode: ", proxy_mode, "\n",
    "============================================================\n",
    verbose = verbose
  )
  
  X_cut <- unbalancedness_tensor(
    X_full    = X_full,
    Unb       = Unb,
    current_t = calibration_t
  )
  
  W_cut <- array(
    as.integer(!is.na(X_cut)),
    dim      = dim(X_cut),
    dimnames = dimnames(X_cut)
  )
  
  out_std <- standardize_mat_with_na(X_cut)
  X_std   <- out_std$X_scaled
  
  if (!is.null(user_hyper$r_impute)) {
    
    r_impute <- as.integer(user_hyper$r_impute)
    
    out_cl <- mfm.cl(
      Y = X_std,
      W = W_cut,
      r = r_impute
    )
    
    X_imputed <- X_std
    X_imputed[W_cut == 0] <- out_cl$Y_hat[W_cut == 0]
    
  } else {
    
    out_init <- init_CL_Yu(
      Y_std = X_std,
      W     = W_cut,
      kmax  = params$Kmax
    )
    
    r_impute  <- out_init$r
    X_imputed <- out_init$Y_init
  }
  
  X_quarterly <- aggregate_tensor_to_quarterly(
    X_tens = X_imputed,
    agg    = agg,
    N_m    = N_m,
    N_q    = N_q
  )
  
  calibration_date <- dates_m[calibration_t]
  idx_pub <- which(dates_q < calibration_date)
  
  if (length(idx_pub) < 2L) {
    stop("Not enough published quarterly observations at calibration date.")
  }
  
  T_q_use <- min(
    max(idx_pub),
    dim(X_quarterly)[1],
    nrow(Y_q_all_full)
  )
  
  if (T_q_use < 2L) {
    stop("Quarterly sample too short at calibration date.")
  }
  
  X_lf  <- X_quarterly[1:T_q_use, , , drop = FALSE]
  X_hf  <- X_imputed
  Y_cut <- Y_q_all_full[1:T_q_use, , drop = FALSE]
  
  country_order <- dimnames(X_lf)[[2]]
  
  if (is.null(country_order) || length(country_order) == 0L) {
    stop("Country names are missing from the row dimension of X_lf.")
  }
  
  if (proxy_name %in% country_order) {
    stop("The aggregate target cannot be part of the row dimension of X_lf.")
  }
  
  if (!setequal(
    country_order,
    setdiff(colnames(Y_cut), proxy_name)
  )) {
    stop("Country rows of X_lf and national GDP columns of Y_cut are not aligned.")
  }
  
  Y_cut <- Y_cut[
    ,
    c(proxy_name, country_order),
    drop = FALSE
  ]
  
  y_proxy <- as.numeric(Y_cut[, proxy_name])
  
  proxy_selection <- NULL
  Z_q_grid <- NULL
  
  X_lf_grid <- X_lf
  X_hf_grid <- X_hf
  y_q_grid  <- y_proxy
  
  target_proxy_window <- NULL
  
  if (proxy_mode == "scalar") {
    
    X_lf_vec <- tensor_to_vector(
      X_tens = X_lf,
      N_m    = N_m,
      N_q    = N_q
    )
    
    if (!is.null(user_hyper$Lproxy)) {
      
      Lproxy <- as.integer(user_hyper$Lproxy)
      
      if (Lproxy < 1L) {
        stop("user_hyper$Lproxy must be at least one.")
      }
      
      proxy_selection <- NULL
      
    } else {
      
      proxy_selection <- select_L_autoproxy_3prf(
        X_lf = X_lf_vec,
        y_q  = y_proxy,
        Zmax = params$Zmax
      )
      
      Lproxy <- proxy_selection$L_opt
    }
    
    Z_q_grid <- NULL
    
    X_lf_grid <- X_lf
    X_hf_grid <- X_hf
    y_q_grid  <- y_proxy
    
    use_autoproxy_grid <- TRUE
    pass1_orthonormalize_grid <- TRUE
    pass1_center_grid <- TRUE
    
  } else {
    
    target_proxy_window <- prepare_multivariate_proxy_window(
      X_lf        = X_lf,
      X_hf        = X_hf,
      Y_q_all     = Y_cut,
      country_cols = country_order
    )
    
    Z_q_grid <- target_proxy_window$Z_q
    
    X_lf_grid <- target_proxy_window$X_lf
    X_hf_grid <- target_proxy_window$X_hf
    y_q_grid  <- y_proxy[target_proxy_window$q_idx]
    
    Lproxy <- ncol(Z_q_grid)
    
    if (!is.null(user_hyper$Lproxy) &&
        as.integer(user_hyper$Lproxy) != Lproxy) {
      stop(
        "In multivariate mode, user_hyper$Lproxy must equal the number of countries."
      )
    }
    
    use_autoproxy_grid <- FALSE
    pass1_orthonormalize_grid <- FALSE
    pass1_center_grid <- FALSE
  }
  
  need_grid_selection <- is.null(user_hyper$p_AR) ||
    is.null(user_hyper$L_midas) ||
    is.null(user_hyper$r_targeted)
  
  if (need_grid_selection) {
    
    lag_sel <- choose_UMIDAS_grid_tensor_MF(
      X_lf = X_lf_grid,
      X_hf = X_hf_grid,
      y_q  = y_q_grid,
      
      rmax     = r_impute,
      Lproxy   = Lproxy,
      Lmax     = params$Lmax,
      p_AR_max = params$p_AR_max,
      
      use_autoproxy = use_autoproxy_grid,
      y_proxy       = y_proxy,
      Z_q           = Z_q_grid,
      
      standardize_proxy        = TRUE,
      orthonormalize_each_iter = TRUE,
      orthonormalize_final_Z   = TRUE,
      
      pass1_orthonormalize_Z = pass1_orthonormalize_grid,
      pass1_center_Z         = pass1_center_grid,
      
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol
    )
    
  } else {
    
    lag_sel <- NULL
  }
  
  p_AR <- if (!is.null(user_hyper$p_AR)) {
    as.integer(user_hyper$p_AR)
  } else {
    as.integer(lag_sel$best_BIC$p_AR)
  }
  
  L_midas <- if (!is.null(user_hyper$L_midas)) {
    as.integer(user_hyper$L_midas)
  } else {
    as.integer(lag_sel$best_BIC$L)
  }
  
  r_targeted <- if (!is.null(user_hyper$r_targeted)) {
    as.integer(user_hyper$r_targeted)
  } else {
    as.integer(lag_sel$r_selected)
  }
  
  list(
    regime_label     = regime_label,
    calibration_t    = calibration_t,
    calibration_date = calibration_date,
    T_q_use          = T_q_use,
    
    proxy_mode       = proxy_mode,
    proxy_columns    = if (proxy_mode == "scalar") {
      proxy_name
    } else {
      country_order
    },
    
    r_impute         = r_impute,
    Lproxy           = Lproxy,
    p_AR             = p_AR,
    L_midas          = L_midas,
    r_targeted       = r_targeted,
    
    proxy_selection  = proxy_selection
  )
}


# ==============================================================================
# PSEUDO REAL-TIME MATRIX MF-TPRF
# ==============================================================================

pseudo_realtime_tensor_mf_tprf_fixed_hyper <- function(
    X_full,
    Y_q_all_full,
    proxy_name = "EA",
    proxy_mode = c("scalar", "multivariate"),
    forecast_mode = c("countries", "aggregate", "both"),
    params,
    dates_m,
    dates_q,
    W_full,
    Unb,
    agg,
    N_m,
    N_q,
    do_post_covid_recalibration = TRUE,
    regime_data_pre  = NULL,
    regime_data_post = NULL,
    post_start_date  = params$covid_end %m+% months(1),
    user_hyper_pre  = list(
      r_impute   = NULL,
      Lproxy     = NULL,
      p_AR       = NULL,
      L_midas    = NULL,
      r_targeted = NULL
    ),
    user_hyper_post = list(
      r_impute   = NULL,
      Lproxy     = NULL,
      p_AR       = NULL,
      L_midas    = NULL,
      r_targeted = NULL
    ),
    ils_maxit = 100,
    ils_tol   = 1e-8,
    verbose   = TRUE
) {
  
  proxy_mode <- match.arg(proxy_mode)
  forecast_mode <- match.arg(forecast_mode)
  
  stopifnot(length(dim(X_full)) == 3)
  
  if (!all(dim(W_full) == dim(X_full))) {
    stop("W_full must have the same dimensions as X_full.")
  }
  
  if (!proxy_name %in% colnames(Y_q_all_full)) {
    stop("proxy_name not found in Y_q_all_full.")
  }
  
  check_regime_data <- function(regime_data, name) {
    
    required <- c("X_full", "W_full", "Unb", "agg", "N_m", "N_q")
    
    if (!all(required %in% names(regime_data))) {
      stop(
        name,
        " is missing: ",
        paste(setdiff(required, names(regime_data)), collapse = ", ")
      )
    }
    
    if (length(dim(regime_data$X_full)) != 3L) {
      stop(name, "$X_full must be a three-dimensional array.")
    }
    
    if (!all(dim(regime_data$W_full) == dim(regime_data$X_full))) {
      stop(name, "$W_full must have the same dimensions as ", name, "$X_full.")
    }
    
    if (dim(regime_data$X_full)[3] != regime_data$N_m + regime_data$N_q) {
      stop(name, ": dim(X_full)[3] must equal N_m + N_q.")
    }
    
    if (!all(dim(regime_data$Unb) == c(
      dim(regime_data$X_full)[2],
      dim(regime_data$X_full)[3]
    ))) {
      stop(name, "$Unb has incompatible dimensions.")
    }
    
    if (!all(dim(regime_data$agg) == c(
      dim(regime_data$X_full)[2],
      dim(regime_data$X_full)[3]
    ))) {
      stop(name, "$agg has incompatible dimensions.")
    }
    
    if (is.null(regime_data$label)) {
      regime_data$label <- name
    }
    
    if (is.null(regime_data$selection_end)) {
      regime_data$selection_end <- NA
    }
    
    regime_data
  }
  
  if (is.null(regime_data_pre)) {
    
    regime_data_pre <- list(
      label         = "pre",
      selection_end = NA,
      X_full        = X_full,
      W_full        = W_full,
      Unb           = Unb,
      agg           = agg,
      N_m           = N_m,
      N_q           = N_q
    )
  }
  
  if (is.null(regime_data_post)) {
    regime_data_post <- regime_data_pre
  }
  
  regime_data_pre <- check_regime_data(
    regime_data = regime_data_pre,
    name = "regime_data_pre"
  )
  
  regime_data_post <- check_regime_data(
    regime_data = regime_data_post,
    name = "regime_data_post"
  )
  
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  
  if (length(t_start) != 1L || length(t_end) != 1L) {
    stop("start_eval and end_eval must each appear exactly once in dates_m.")
  }
  
  if (t_start > t_end) {
    stop("start_eval must be weakly before end_eval.")
  }
  
  t_est_end <- max(which(dates_m < params$start_eval))
  
  if (!is.finite(t_est_end) || t_est_end < 24L) {
    stop("Estimation sample is too short before start_eval.")
  }
  
  t_recalib <- which(dates_m >= params$covid_end)[1]
  t_post_start <- which(dates_m >= post_start_date)[1]
  
  do_recal <- isTRUE(do_post_covid_recalibration) &&
    !is.na(t_recalib) &&
    !is.na(t_post_start) &&
    t_recalib <= t_end &&
    t_post_start <= t_end
  
  countries_eval <- setdiff(colnames(Y_q_all_full), proxy_name)
  
  if (length(countries_eval) == 0L) {
    stop("No national GDP targets found in Y_q_all_full.")
  }
  
  hyper_pre <- select_regime_hyperparameters(
    calibration_t = t_est_end,
    regime_label  = "PRE-COVID",
    
    X_full        = regime_data_pre$X_full,
    Y_q_all_full  = Y_q_all_full,
    
    proxy_name = proxy_name,
    proxy_mode = proxy_mode,
    
    params  = params,
    dates_m = dates_m,
    dates_q = dates_q,
    
    Unb = regime_data_pre$Unb,
    agg = regime_data_pre$agg,
    N_m = regime_data_pre$N_m,
    N_q = regime_data_pre$N_q,
    
    user_hyper = user_hyper_pre,
    
    ils_maxit = ils_maxit,
    ils_tol   = ils_tol,
    verbose   = verbose
  )
  
  hyper_post <- hyper_pre
  
  if (do_recal) {
    
    hyper_post <- select_regime_hyperparameters(
      calibration_t = t_recalib,
      regime_label  = "POST-COVID",
      
      X_full        = regime_data_post$X_full,
      Y_q_all_full  = Y_q_all_full,
      
      proxy_name = proxy_name,
      proxy_mode = proxy_mode,
      
      params  = params,
      dates_m = dates_m,
      dates_q = dates_q,
      
      Unb = regime_data_post$Unb,
      agg = regime_data_post$agg,
      N_m = regime_data_post$N_m,
      N_q = regime_data_post$N_q,
      
      user_hyper = user_hyper_post,
      
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol,
      verbose   = verbose
    )
  }
  
  n_vintages <- t_end - t_start + 1L
  
  print_if_verbose(
    "\n============================================================\n",
    "PSEUDO REAL-TIME MATRIX MF-TPRF\n",
    "============================================================\n",
    "Proxy mode: ", proxy_mode, "\n",
    "Evaluation vintages: ", n_vintages, "\n",
    "Pre-COVID hyperparameters:\n",
    "  r_impute = (", hyper_pre$r_impute[1], ", ", hyper_pre$r_impute[2], ")\n",
    "  Lproxy   = ", hyper_pre$Lproxy, "\n",
    "  p_AR     = ", hyper_pre$p_AR, "\n",
    "  L_midas  = ", hyper_pre$L_midas, "\n",
    "  r_cap    = (", hyper_pre$r_targeted[1], ", ", hyper_pre$r_targeted[2], ")\n",
    "Post-COVID hyperparameters:\n",
    "  r_impute = (", hyper_post$r_impute[1], ", ", hyper_post$r_impute[2], ")\n",
    "  Lproxy   = ", hyper_post$Lproxy, "\n",
    "  p_AR     = ", hyper_post$p_AR, "\n",
    "  L_midas  = ", hyper_post$L_midas, "\n",
    "  r_cap    = (", hyper_post$r_targeted[1], ", ", hyper_post$r_targeted[2], ")\n",
    "============================================================\n\n",
    verbose = verbose
  )
  
  now_M1 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M2 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M3 <- setNames(vector("list", length(countries_eval)), countries_eval)
  
  now_A_M1 <- list()
  now_A_M2 <- list()
  now_A_M3 <- list()
  
  vintage_log <- vector("list", t_end - t_start + 1L)
  log_idx <- 0L
  
  for (tt in seq.int(t_start, t_end)) {
    
    use_post <- do_recal && tt >= t_post_start
    
    active_hyper <- if (use_post) hyper_post else hyper_pre
    active_data  <- if (use_post) regime_data_post else regime_data_pre
    active_regime <- if (use_post) "POST-COVID" else "PRE-COVID"
    
    rt <- build_realtime_dataset_fixed_rank(
      current_t    = tt,
      X_full       = active_data$X_full,
      Y_q_all_full = Y_q_all_full,
      dates_m      = dates_m,
      dates_q      = dates_q,
      Unb          = active_data$Unb,
      agg          = active_data$agg,
      N_m          = active_data$N_m,
      N_q          = active_data$N_q,
      r_impute     = active_hyper$r_impute
    )
    
    if (is.null(rt)) {
      next
    }
    
    country_order_rt <- dimnames(rt$X_lf)[[2]]
    
    if (is.null(country_order_rt) ||
        !setequal(
          country_order_rt,
          setdiff(colnames(rt$Y_cut), proxy_name)
        )) {
      stop("Country rows of real-time X_lf and columns of Y_cut are not aligned.")
    }
    
    rt$Y_cut <- rt$Y_cut[
      ,
      c(proxy_name, country_order_rt),
      drop = FALSE
    ]
    
    out_rt <- Tensor_MF_TPRF(
      X_lf    = rt$X_lf,
      X_hf    = rt$X_hf,
      Y_q_all = rt$Y_cut,
      
      proxy_name   = proxy_name,
      proxy_mode   = proxy_mode,
      forecast_mode = forecast_mode,
      
      Lproxy  = active_hyper$Lproxy,
      L_midas = active_hyper$L_midas,
      p_AR    = active_hyper$p_AR,
      rmax    = active_hyper$r_targeted,
      
      standardize_proxy        = TRUE,
      orthonormalize_each_iter = TRUE,
      orthonormalize_final_Z   = TRUE,
      
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol
    )
    
    m_tr <- compute_m_tr(
      date_t  = dates_m[tt],
      dates_q = dates_q
    )
    
    if (is.na(m_tr)) {
      next
    }
    
    key <- as.character(dates_m[tt])
    
    if (forecast_mode %in% c("aggregate", "both")) {
      
      y_nowcast_A <- out_rt$aggregate$y_nowcast
      value_A <- tail(y_nowcast_A, 1L)
      
      if (m_tr == 1L) now_A_M1[[key]] <- value_A
      if (m_tr == 2L) now_A_M2[[key]] <- value_A
      if (m_tr == 3L) now_A_M3[[key]] <- value_A
    }
    
    if (forecast_mode %in% c("countries", "both")) {
      
      for (cc in countries_eval) {
        
        y_nowcast_cc <- out_rt$by_country[[cc]]$y_nowcast
        value_cc <- tail(y_nowcast_cc, 1L)
        
        if (m_tr == 1L) now_M1[[cc]][[key]] <- value_cc
        if (m_tr == 2L) now_M2[[cc]][[key]] <- value_cc
        if (m_tr == 3L) now_M3[[cc]][[key]] <- value_cc
      }
    }
    
    log_idx <- log_idx + 1L
    
    vintage_log[[log_idx]] <- list(
      vintage_date    = dates_m[tt],
      active_regime   = active_regime,
      proxy_mode      = proxy_mode,
      proxy_columns   = active_hyper$proxy_columns,
      selection_label = active_data$label,
      selection_end   = active_data$selection_end,
      T_q_use         = rt$T_q_use,
      month_position  = m_tr,
      r_impute        = active_hyper$r_impute,
      Lproxy          = active_hyper$Lproxy,
      p_AR            = active_hyper$p_AR,
      L_midas         = active_hyper$L_midas,
      r_targeted      = active_hyper$r_targeted
    )
    
    vintage_number <- tt - t_start + 1L
    
    show_progress <- vintage_number == 1L ||
      vintage_number == n_vintages ||
      vintage_number %% 12L == 0L ||
      (do_recal && tt == t_post_start)
    
    if (show_progress) {
      print_if_verbose(
        "Vintage ", vintage_number, "/", n_vintages,
        " | ", as.character(dates_m[tt]),
        " | ", active_regime,
        " | Tq = ", rt$T_q_use,
        " | M", m_tr,
        "\n",
        verbose = verbose
      )
    }
  }
  
  vintage_log <- vintage_log[seq_len(log_idx)]
  
  list(
    proxy_mode = proxy_mode,
    
    hyper_pre  = hyper_pre,
    hyper_post = hyper_post,
    
    regime_info = list(
      pre = list(
        label         = regime_data_pre$label,
        selection_end = regime_data_pre$selection_end,
        N_m           = regime_data_pre$N_m,
        N_q           = regime_data_pre$N_q,
        K             = dim(regime_data_pre$X_full)[3]
      ),
      post = list(
        label         = regime_data_post$label,
        selection_end = regime_data_post$selection_end,
        N_m           = regime_data_post$N_m,
        N_q           = regime_data_post$N_q,
        K             = dim(regime_data_post$X_full)[3]
      )
    ),
    
    forecast_mode = forecast_mode,
    countries = if (forecast_mode %in% c("countries", "both")) {
      countries_eval
    } else {
      character(0)
    },
    
    M1 = if (forecast_mode %in% c("countries", "both")) {
      lapply(now_M1, unlist)
    } else {
      NULL
    },
    
    M2 = if (forecast_mode %in% c("countries", "both")) {
      lapply(now_M2, unlist)
    } else {
      NULL
    },
    
    M3 = if (forecast_mode %in% c("countries", "both")) {
      lapply(now_M3, unlist)
    } else {
      NULL
    },
    
    aggregate = list(
      target = proxy_name,
      M1 = unlist(now_A_M1),
      M2 = unlist(now_A_M2),
      M3 = unlist(now_A_M3)
    ),
    
    vintage_log = vintage_log
  )
}