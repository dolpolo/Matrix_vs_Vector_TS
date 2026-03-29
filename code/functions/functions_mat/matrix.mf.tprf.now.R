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
    params,
    dates_m,
    dates_q,
    Unb,
    agg,
    N_m,
    N_q,
    user_hyper = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL, r_targeted = NULL),
    ils_maxit = 100,
    ils_tol   = 1e-8,
    verbose   = TRUE
) {
  # Select regime-specific hyperparameters once, or accept user-fixed values.
  # User-supplied values always take precedence over internal selection.
  
  print_if_verbose(
    "\n============================================================\n",
    "Selecting hyperparameters for regime: ", regime_label, "\n",
    "Calibration vintage: ", as.character(dates_m[calibration_t]), "\n",
    "============================================================\n",
    verbose = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 1. Build calibration dataset with real-time publication delays
  # --------------------------------------------------------------------------
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
  
  # --------------------------------------------------------------------------
  # Step 2. Determine imputation rank
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$r_impute)) {
    r_impute <- user_hyper$r_impute
    
    print_if_verbose(
      "Imputation rank provided by user: (",
      r_impute[1], ", ", r_impute[2], ")\n",
      verbose = verbose
    )
    
    out_cl <- mfm.cl(Y = X_std, W = W_cut, r = r_impute)
    
    X_imputed <- X_std
    X_imputed[W_cut == 0] <- out_cl$Y_hat[W_cut == 0]
    
  } else {
    print_if_verbose(
      "Imputation rank not provided. Selecting it internally with init_CL_Yu().\n",
      verbose = verbose
    )
    
    out_init <- init_CL_Yu(
      Y_std = X_std,
      W     = W_cut,
      kmax  = params$Kmax
    )
    
    r_impute <- out_init$r
    X_imputed <- out_init$Y_init
  }
  
  # --------------------------------------------------------------------------
  # Step 3. Aggregate monthly tensor to quarterly tensor
  # --------------------------------------------------------------------------
  X_quarterly <- aggregate_tensor_to_quarterly(
    X_tens = X_imputed,
    agg    = agg,
    N_m    = N_m,
    N_q    = N_q
  )
  
  calibration_date <- dates_m[calibration_t]
  idx_pub <- which(dates_q < calibration_date)
  
  if (length(idx_pub) < 2) {
    stop("Not enough published quarterly observations at calibration date.")
  }
  
  T_q_use <- min(max(idx_pub), dim(X_quarterly)[1], nrow(Y_q_all_full))
  if (T_q_use < 2) {
    stop("Quarterly sample too short at calibration date.")
  }
  
  Y_cut <- Y_q_all_full[1:T_q_use, , drop = FALSE]
  X_lf  <- X_quarterly[1:T_q_use, , , drop = FALSE]
  X_hf  <- X_imputed
  
  y_proxy  <- as.numeric(Y_cut[, proxy_name])
  X_lf_vec <- tensor_to_vector(X_tens = X_lf, N_m = N_m, N_q = N_q)
  
  # --------------------------------------------------------------------------
  # Step 4. Determine number of proxy factors
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$Lproxy)) {
    Lproxy <- user_hyper$Lproxy
    
    print_if_verbose(
      "Number of proxy factors provided by user: ", Lproxy, "\n",
      verbose = verbose
    )
    
  } else {
    print_if_verbose(
      "Number of proxy factors not provided. Selecting it internally.\n",
      verbose = verbose
    )
    
    proxy_sel <- select_L_autoproxy_3prf(
      X    = X_lf_vec,
      y    = y_proxy,
      Zmax = params$Zmax
    )
    
    Lproxy <- proxy_sel$L_opt
  }
  
  # --------------------------------------------------------------------------
  # Step 5. Determine p_AR, L_midas, and targeted rank
  # --------------------------------------------------------------------------
  need_grid_selection <- is.null(user_hyper$p_AR) ||
    is.null(user_hyper$L_midas) ||
    is.null(user_hyper$r_targeted)
  
  if (need_grid_selection) {
    print_if_verbose(
      "Some dynamic hyperparameters are missing. Running internal grid selection.\n",
      verbose = verbose
    )
    
    lag_sel <- choose_UMIDAS_grid_tensor_MF(
      X_lf = X_lf,
      X_hf = X_hf,
      y_q  = y_proxy,
      rmax = r_impute,
      Lproxy   = Lproxy,
      Lmax     = params$Lmax,
      p_AR_max = params$p_AR_max,
      use_autoproxy = TRUE,
      y_proxy = y_proxy,
      standardize_proxy = TRUE,
      orthonormalize_each_iter = TRUE,
      orthonormalize_final_Z   = TRUE,
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol
    )
  } else {
    lag_sel <- NULL
  }
  
  if (!is.null(user_hyper$p_AR)) {
    p_AR <- user_hyper$p_AR
    print_if_verbose("AR lag order provided by user: ", p_AR, "\n", verbose = verbose)
  } else {
    p_AR <- lag_sel$best_BIC$p_AR
    print_if_verbose("AR lag order selected internally: ", p_AR, "\n", verbose = verbose)
  }
  
  if (!is.null(user_hyper$L_midas)) {
    L_midas <- user_hyper$L_midas
    print_if_verbose("MIDAS lag order provided by user: ", L_midas, "\n", verbose = verbose)
  } else {
    L_midas <- lag_sel$best_BIC$L
    print_if_verbose("MIDAS lag order selected internally: ", L_midas, "\n", verbose = verbose)
  }
  
  if (!is.null(user_hyper$r_targeted)) {
    r_targeted <- user_hyper$r_targeted
    print_if_verbose(
      "Targeted factor rank provided by user: (",
      r_targeted[1], ", ", r_targeted[2], ")\n",
      verbose = verbose
    )
  } else {
    r_targeted <- lag_sel$r_selected
    print_if_verbose(
      "Targeted factor rank selected internally: (",
      r_targeted[1], ", ", r_targeted[2], ")\n",
      verbose = verbose
    )
  }
  
  print_if_verbose(
    "\nFrozen hyperparameters for regime ", regime_label, ":\n",
    "  - Calibration date     : ", as.character(calibration_date), "\n",
    "  - Available quarters   : ", T_q_use, "\n",
    "  - Imputation rank      : (", r_impute[1], ", ", r_impute[2], ")\n",
    "  - Number of proxies    : ", Lproxy, "\n",
    "  - AR lag order         : ", p_AR, "\n",
    "  - MIDAS lag order      : ", L_midas, "\n",
    "  - Targeted factor rank : (", r_targeted[1], ", ", r_targeted[2], ")\n\n",
    verbose = verbose
  )
  
  list(
    regime_label     = regime_label,
    calibration_t    = calibration_t,
    calibration_date = calibration_date,
    T_q_use          = T_q_use,
    r_impute         = r_impute,
    Lproxy           = Lproxy,
    p_AR             = p_AR,
    L_midas          = L_midas,
    r_targeted       = r_targeted
  )
}


# ==============================================================================
# PSEUDO REAL-TIME (EXPANDING) — MATRIX MF-TPRF WITH FIXED HYPERPARAMETERS
# ==============================================================================

pseudo_realtime_tensor_mf_tprf_fixed_hyper <- function(
    X_full,
    Y_q_all_full,
    proxy_name = "EA",
    params,
    dates_m,
    dates_q,
    W_full,
    Unb,
    agg,
    N_m,
    N_q,
    do_post_covid_recalibration = TRUE,
    user_hyper_pre  = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL, r_targeted = NULL),
    user_hyper_post = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL, r_targeted = NULL),
    ils_maxit = 100,
    ils_tol   = 1e-8,
    verbose   = TRUE
) {
  # Expanding pseudo real-time nowcast with fixed complexity within each regime.
  # Hyperparameters are selected only once per regime unless explicitly provided
  # by the user.
  
  # --------------------------------------------------------------------------
  # Step 0. Basic checks
  # --------------------------------------------------------------------------
  stopifnot(length(dim(X_full)) == 3)
  
  T_m <- dim(X_full)[1]
  P   <- dim(X_full)[2]
  K   <- dim(X_full)[3]
  
  if (K != (N_m + N_q)) stop("K must equal N_m + N_q (GDP excluded).")
  if (!all(dim(W_full) == dim(X_full))) stop("W_full must have same dim as X_full.")
  if (!(proxy_name %in% colnames(Y_q_all_full))) stop("proxy_name not found in Y_q_all_full.")
  
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("start_eval/end_eval not found in dates_m.")
  }
  if (t_start > t_end) stop("start_eval > end_eval.")
  
  t_est_end <- max(which(dates_m < params$start_eval))
  if (!is.finite(t_est_end) || t_est_end < 24) {
    stop("Estimation sample too short before evaluation.")
  }
  
  t_recalib <- which(dates_m >= params$covid_end)[1]
  do_recal  <- isTRUE(do_post_covid_recalibration) &&
    !is.na(t_recalib) &&
    (t_recalib <= t_end)
  
  countries_eval <- setdiff(colnames(Y_q_all_full), proxy_name)
  if (length(countries_eval) < 1) stop("No countries beyond proxy in Y_q_all_full.")
  
  print_if_verbose(
    "\n############################################################\n",
    "EXPANDING PSEUDO REAL-TIME MATRIX MF-TPRF\n",
    "Model complexity is fixed within each regime\n",
    "############################################################\n",
    "Evaluation window:\n",
    "  - Start: ", as.character(dates_m[t_start]), "\n",
    "  - End  : ", as.character(dates_m[t_end]), "\n",
    "Pre-COVID calibration date : ", as.character(dates_m[t_est_end]), "\n",
    "Post-COVID recalibration   : ", if (do_recal) as.character(dates_m[t_recalib]) else "disabled", "\n\n",
    verbose = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 1. Select or accept PRE-COVID hyperparameters
  # --------------------------------------------------------------------------
  hyper_pre <- select_regime_hyperparameters(
    calibration_t = t_est_end,
    regime_label  = "PRE-COVID",
    X_full        = X_full,
    Y_q_all_full  = Y_q_all_full,
    proxy_name    = proxy_name,
    params        = params,
    dates_m       = dates_m,
    dates_q       = dates_q,
    Unb           = Unb,
    agg           = agg,
    N_m           = N_m,
    N_q           = N_q,
    user_hyper    = user_hyper_pre,
    ils_maxit     = ils_maxit,
    ils_tol       = ils_tol,
    verbose       = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 2. Select or accept POST-COVID hyperparameters
  # --------------------------------------------------------------------------
  hyper_post <- hyper_pre
  
  if (do_recal) {
    hyper_post <- select_regime_hyperparameters(
      calibration_t = t_recalib,
      regime_label  = "POST-COVID",
      X_full        = X_full,
      Y_q_all_full  = Y_q_all_full,
      proxy_name    = proxy_name,
      params        = params,
      dates_m       = dates_m,
      dates_q       = dates_q,
      Unb           = Unb,
      agg           = agg,
      N_m           = N_m,
      N_q           = N_q,
      user_hyper    = user_hyper_post,
      ils_maxit     = ils_maxit,
      ils_tol       = ils_tol,
      verbose       = verbose
    )
  }
  
  # --------------------------------------------------------------------------
  # Step 3. Create output containers
  # --------------------------------------------------------------------------
  now_M1 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M2 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M3 <- setNames(vector("list", length(countries_eval)), countries_eval)
  
  vintage_log <- vector("list", length = t_end - t_start + 1)
  log_idx <- 0L
  
  # --------------------------------------------------------------------------
  # Step 4. Expanding real-time loop
  # --------------------------------------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    
    use_post <- do_recal && (tt >= t_recalib)
    active_hyper <- if (use_post) hyper_post else hyper_pre
    active_regime <- if (use_post) "POST-COVID" else "PRE-COVID"
    
    print_if_verbose(
      "\n------------------------------------------------------------\n",
      "Real-time vintage: ", as.character(date_t), "\n",
      "Active regime    : ", active_regime, "\n",
      "Frozen hyperparameters used:\n",
      "  - Imputation rank      : (", active_hyper$r_impute[1], ", ", active_hyper$r_impute[2], ")\n",
      "  - Number of proxies    : ", active_hyper$Lproxy, "\n",
      "  - AR lag order         : ", active_hyper$p_AR, "\n",
      "  - MIDAS lag order      : ", active_hyper$L_midas, "\n",
      "  - Targeted factor rank : (", active_hyper$r_targeted[1], ", ", active_hyper$r_targeted[2], ")\n",
      verbose = verbose
    )
    
    rt <- build_realtime_dataset_fixed_rank(
      current_t    = tt,
      X_full       = X_full,
      Y_q_all_full = Y_q_all_full,
      dates_m      = dates_m,
      dates_q      = dates_q,
      Unb          = Unb,
      agg          = agg,
      N_m          = N_m,
      N_q          = N_q,
      r_impute     = active_hyper$r_impute
    )
    
    if (is.null(rt)) {
      print_if_verbose(
        "Skipped vintage: not enough published quarterly information.\n",
        verbose = verbose
      )
      next
    }
    
    out_rt <- Tensor_MF_TPRF(
      X_lf       = rt$X_lf,
      X_hf       = rt$X_hf,
      Y_q_all    = rt$Y_cut,
      proxy_name = proxy_name,
      Lproxy     = active_hyper$Lproxy,
      L_midas    = active_hyper$L_midas,
      p_AR       = active_hyper$p_AR,
      rmax       = active_hyper$r_targeted,
      standardize_proxy = TRUE,
      orthonormalize_each_iter = TRUE,
      orthonormalize_final_Z   = TRUE,
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol
    )
    
    m_tr <- compute_m_tr(date_t, dates_q)
    if (is.na(m_tr)) {
      print_if_verbose(
        "Skipped vintage: invalid month-within-quarter mapping.\n",
        verbose = verbose
      )
      next
    }
    
    key <- as.character(date_t)
    
    for (cc in countries_eval) {
      y_cc <- out_rt$by_country[[cc]]$y_nowcast
      if (is.null(y_cc) || length(y_cc) == 0) next
      
      y_last <- tail(y_cc, 1)
      
      if (m_tr == 1L) now_M1[[cc]][[key]] <- y_last
      if (m_tr == 2L) now_M2[[cc]][[key]] <- y_last
      if (m_tr == 3L) now_M3[[cc]][[key]] <- y_last
    }
    
    log_idx <- log_idx + 1L
    vintage_log[[log_idx]] <- list(
      vintage_date   = date_t,
      active_regime  = active_regime,
      T_q_use        = rt$T_q_use,
      month_position = m_tr,
      r_impute       = active_hyper$r_impute,
      Lproxy         = active_hyper$Lproxy,
      p_AR           = active_hyper$p_AR,
      L_midas        = active_hyper$L_midas,
      r_targeted     = active_hyper$r_targeted
    )
    
    print_if_verbose(
      "Vintage completed successfully.\n",
      "  - Published quarters used : ", rt$T_q_use, "\n",
      "  - Month position          : M", m_tr, "\n",
      verbose = verbose
    )
  }
  
  vintage_log <- vintage_log[seq_len(log_idx)]
  
  # --------------------------------------------------------------------------
  # Step 5. Return output using the original naming convention
  # --------------------------------------------------------------------------
  list(
    hyper_pre  = hyper_pre,
    hyper_post = hyper_post,
    countries  = countries_eval,
    M1         = lapply(now_M1, function(x) unlist(x)),
    M2         = lapply(now_M2, function(x) unlist(x)),
    M3         = lapply(now_M3, function(x) unlist(x)),
    vintage_log = vintage_log
  )
}