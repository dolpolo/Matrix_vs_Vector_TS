# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

unbalancedness <- function(X_full, dates, Freq, Unb, current_t) {
  # Construct the pseudo-real-time dataset up to current_t by applying
  # publication delays to each variable.
  
  stopifnot(nrow(X_full) >= current_t)
  
  X_cut <- X_full[1:current_t, , drop = FALSE]
  
  for (j in seq_len(ncol(X_cut))) {
    delay_j   <- Unb[j]
    release_t <- current_t - delay_j
    
    if (release_t < 1) {
      # The variable has not been released yet at the current vintage
      X_cut[, j] <- NA
    } else {
      # The variable is available only up to release_t
      if (release_t + 1 <= current_t) {
        X_cut[(release_t + 1):current_t, j] <- NA
      }
    }
  }
  
  return(X_cut)
}


# ==============================================================================
# MONTH-IN-QUARTER INDEX
# ==============================================================================

compute_m_tr <- function(date_t, dates_q) {
  # Return:
  #   1 = first month after the latest published quarter
  #   2 = second month after the latest published quarter
  #   3 = third month after the latest published quarter
  
  last_q_date <- max(dates_q[dates_q < date_t])
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
# BUILD XP-IMPUTED AND AGGREGATED DATASET UP TO A GIVEN VINTAGE
# ==============================================================================

build_XP_agg_up_to <- function(
    X_full,
    y_q,
    dates,
    dates_q,
    Freq,
    Unb,
    current_t,
    agg_m,
    agg_q,
    Kmax = 10,
    r_impute = NULL
) {
  # Build the real-time dataset up to current_t.
  #
  # If r_impute is NULL:
  #   - select the number of factors internally with init_XP_ER()
  # If r_impute is fixed:
  #   - estimate factors with estimate_factors_XP()
  #   - replace missing values with the common component C_hat
  
  # --------------------------------------------------------------------------
  # Step 1. Apply real-time publication delays
  # --------------------------------------------------------------------------
  X_cut <- unbalancedness(
    X_full    = X_full,
    dates     = dates,
    Freq      = Freq,
    Unb       = Unb,
    current_t = current_t
  )
  
  # --------------------------------------------------------------------------
  # Step 2. Standardize with missing values preserved
  # --------------------------------------------------------------------------
  out_std <- standardize_with_na(X_cut)
  X_std   <- out_std$X_std
  
  # --------------------------------------------------------------------------
  # Step 3. Impute missing values
  # --------------------------------------------------------------------------
  if (is.null(r_impute)) {
    # Automatic rank selection + XP imputation
    imp_xp <- init_XP_ER(X_std, Kmax = Kmax)
    X_xp   <- imp_xp$X_init
    r_used <- imp_xp$r
  } else {
    # Fixed-rank XP factor estimation + imputation
    fac_xp <- estimate_factors_XP(X_std, r = r_impute)
    C_hat  <- fac_xp$C_hat
    
    X_xp <- X_std
    na_idx <- is.na(X_std)
    X_xp[na_idx] <- C_hat[na_idx]
    
    r_used <- r_impute
  }
  
  # --------------------------------------------------------------------------
  # Step 4. Split monthly and quarterly variables and aggregate
  # --------------------------------------------------------------------------
  X_m_xp <- X_xp[, Freq == "M", drop = FALSE]
  X_q_xp <- X_xp[, Freq == "Q", drop = FALSE]
  
  X_mq_xp <- agg_mq(X_m_xp, agg_m)
  X_qq_xp <- agg_qq(X_q_xp, agg_q)
  
  X_xp_agg <- cbind(X_mq_xp, X_qq_xp)
  
  # --------------------------------------------------------------------------
  # Step 5. Keep only quarters published at current_t
  # --------------------------------------------------------------------------
  idx_pub <- which(dates_q < dates[current_t])
  if (length(idx_pub) < 2) return(NULL)
  
  T_q_current <- tail(idx_pub, 1)
  T_q_current <- min(T_q_current, nrow(X_xp_agg), length(y_q))
  
  if (T_q_current < 2) return(NULL)
  
  list(
    X_hf = X_xp,
    X_lf = X_xp_agg,
    y_q  = as.numeric(y_q)[1:T_q_current],
    T_q  = T_q_current,
    r    = r_used
  )
}


# ==============================================================================
# SELECT OR ACCEPT HYPERPARAMETERS FOR ONE REGIME
# ==============================================================================

select_regime_hyperparameters_XP <- function(
    calibration_t,
    regime_label,
    X_full,
    y_q,
    params,
    dates,
    dates_q,
    Freq,
    Unb,
    agg_m,
    agg_q,
    user_hyper = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL),
    verbose = TRUE
) {
  # Select regime-specific hyperparameters once, or accept user-fixed values.
  # User-supplied values always take precedence over internal selection.
  
  print_if_verbose(
    "\n============================================================\n",
    "Selecting hyperparameters for regime: ", regime_label, "\n",
    "Calibration vintage: ", as.character(dates[calibration_t]), "\n",
    "============================================================\n",
    verbose = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 1. Build calibration dataset
  # --------------------------------------------------------------------------
  obj_est <- build_XP_agg_up_to(
    X_full    = X_full,
    y_q       = y_q,
    dates     = dates,
    dates_q   = dates_q,
    Freq      = Freq,
    Unb       = Unb,
    current_t = calibration_t,
    agg_m     = agg_m,
    agg_q     = agg_q,
    Kmax      = params$Kmax,
    r_impute  = user_hyper$r_impute
  )
  
  if (is.null(obj_est)) {
    stop("Not enough published GDP at calibration date.")
  }
  
  r_use <- obj_est$r
  
  if (!is.null(user_hyper$r_impute)) {
    print_if_verbose(
      "Imputation rank provided by user: ", r_use, "\n",
      verbose = verbose
    )
  } else {
    print_if_verbose(
      "Imputation rank selected internally: ", r_use, "\n",
      verbose = verbose
    )
  }
  
  X_lf_cut <- obj_est$X_lf[1:obj_est$T_q, , drop = FALSE]
  y_q_cut  <- obj_est$y_q
  
  # --------------------------------------------------------------------------
  # Step 2. Determine number of proxy factors
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$Lproxy)) {
    Lproxy <- user_hyper$Lproxy
    
    print_if_verbose(
      "Number of proxy factors provided by user: ", Lproxy, "\n",
      verbose = verbose
    )
    
  } else {
    print_if_verbose(
      "Number of proxy factors not provided. Selecting internally.\n",
      verbose = verbose
    )
    
    pls <- select_L_autoproxy_3prf(
      X    = X_lf_cut,
      y    = y_q_cut,
      Zmax = params$Zmax
    )
    
    Lproxy <- pls$L_opt
  }
  
  # --------------------------------------------------------------------------
  # Step 3. Determine AR and MIDAS lag orders
  # --------------------------------------------------------------------------
  need_lag_selection <- is.null(user_hyper$p_AR) || is.null(user_hyper$L_midas)
  
  if (need_lag_selection) {
    print_if_verbose(
      "Some dynamic hyperparameters are missing. Running internal lag selection.\n",
      verbose = verbose
    )
    
    lag_sel <- choose_UMIDAS_lag(
      X_lf        = X_lf_cut,
      X_hf        = obj_est$X_hf,
      y_q         = y_q_cut,
      Lmax        = params$Lmax,
      Lproxy      = Lproxy,
      p_AR_max    = params$p_AR_max,
      Robust_F    = params$Robust_F,
      alpha       = params$alpha,
      robust_type = params$robust_type,
      nw_lag      = params$nw_lag
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
  
  print_if_verbose(
    "\nFrozen hyperparameters for regime ", regime_label, ":\n",
    "  - Calibration date   : ", as.character(dates[calibration_t]), "\n",
    "  - Available quarters : ", obj_est$T_q, "\n",
    "  - Imputation rank    : ", r_use, "\n",
    "  - Number of proxies  : ", Lproxy, "\n",
    "  - AR lag order       : ", p_AR, "\n",
    "  - MIDAS lag order    : ", L_midas, "\n\n",
    verbose = verbose
  )
  
  list(
    regime_label     = regime_label,
    calibration_t    = calibration_t,
    calibration_date = dates[calibration_t],
    T_q              = obj_est$T_q,
    r                = r_use,
    Lproxy           = Lproxy,
    p_AR             = p_AR,
    L_midas          = L_midas
  )
}


# ==============================================================================
# EXPANDING NOWCAST WITH FIXED HYPERPARAMETERS
# ==============================================================================

pseudo_realtime_MF_TPRF_XP <- function(
    X_full,
    y_q,
    params,
    dates,
    dates_q,
    Freq,
    Unb,
    agg_m,
    agg_q,
    do_post_covid_recalibration = TRUE,
    user_hyper_pre  = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL),
    user_hyper_post = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL),
    verbose = TRUE
) {
  # Expanding pseudo real-time nowcast with fixed complexity within each regime.
  # Hyperparameters are selected only once per regime unless explicitly provided
  # by the user.
  
  # --------------------------------------------------------------------------
  # Step 0. Evaluation window
  # --------------------------------------------------------------------------
  t_start <- which(dates == params$start_eval)
  t_end   <- which(dates == params$end_eval)
  
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("start_eval or end_eval not found in 'dates'.")
  }
  if (t_start > t_end) {
    stop("start_eval must be <= end_eval.")
  }
  
  t_est_end <- max(which(dates < params$start_eval))
  if (length(t_est_end) == 0 || t_est_end < 24) {
    stop("Estimation sample too short or not defined.")
  }
  
  t_recalib <- which(dates >= params$covid_end)[1]
  do_recalib <- isTRUE(do_post_covid_recalibration) &&
    !is.na(t_recalib) &&
    t_recalib <= t_end
  
  print_if_verbose(
    "\n############################################################\n",
    "EXPANDING PSEUDO REAL-TIME VECTOR MF-TPRF\n",
    "Model complexity is fixed within each regime\n",
    "############################################################\n",
    "Evaluation window:\n",
    "  - Start: ", as.character(dates[t_start]), "\n",
    "  - End  : ", as.character(dates[t_end]), "\n",
    "Pre-COVID calibration date : ", as.character(dates[t_est_end]), "\n",
    "Post-COVID recalibration   : ", if (do_recalib) as.character(dates[t_recalib]) else "disabled", "\n\n",
    verbose = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 1. Select or accept PRE-COVID hyperparameters
  # --------------------------------------------------------------------------
  hyper_pre <- select_regime_hyperparameters_XP(
    calibration_t = t_est_end,
    regime_label  = "PRE-COVID",
    X_full        = X_full,
    y_q           = y_q,
    params        = params,
    dates         = dates,
    dates_q       = dates_q,
    Freq          = Freq,
    Unb           = Unb,
    agg_m         = agg_m,
    agg_q         = agg_q,
    user_hyper    = user_hyper_pre,
    verbose       = verbose
  )
  
  # --------------------------------------------------------------------------
  # Step 2. Select or accept POST-COVID hyperparameters
  # --------------------------------------------------------------------------
  hyper_post <- hyper_pre
  
  if (do_recalib) {
    hyper_post <- select_regime_hyperparameters_XP(
      calibration_t = t_recalib,
      regime_label  = "POST-COVID",
      X_full        = X_full,
      y_q           = y_q,
      params        = params,
      dates         = dates,
      dates_q       = dates_q,
      Freq          = Freq,
      Unb           = Unb,
      agg_m         = agg_m,
      agg_q         = agg_q,
      user_hyper    = user_hyper_post,
      verbose       = verbose
    )
  }
  
  # --------------------------------------------------------------------------
  # Step 3. Containers for monthly nowcasts
  # --------------------------------------------------------------------------
  now_M1 <- list()
  now_M2 <- list()
  now_M3 <- list()
  vintage_log <- list()
  
  # --------------------------------------------------------------------------
  # Step 4. Real-time expanding loop with fixed hyperparameters
  # --------------------------------------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates[tt]
    
    use_post <- do_recalib && tt >= t_recalib
    active_hyper  <- if (use_post) hyper_post else hyper_pre
    active_regime <- if (use_post) "POST-COVID" else "PRE-COVID"
    
    print_if_verbose(
      "\n------------------------------------------------------------\n",
      "Real-time vintage: ", as.character(date_t), "\n",
      "Active regime    : ", active_regime, "\n",
      "Frozen hyperparameters used:\n",
      "  - Imputation rank   : ", active_hyper$r, "\n",
      "  - Number of proxies : ", active_hyper$Lproxy, "\n",
      "  - AR lag order      : ", active_hyper$p_AR, "\n",
      "  - MIDAS lag order   : ", active_hyper$L_midas, "\n",
      verbose = verbose
    )
    
    obj_rt <- build_XP_agg_up_to(
      X_full    = X_full,
      y_q       = y_q,
      dates     = dates,
      dates_q   = dates_q,
      Freq      = Freq,
      Unb       = Unb,
      current_t = tt,
      agg_m     = agg_m,
      agg_q     = agg_q,
      Kmax      = params$Kmax,
      r_impute  = active_hyper$r
    )
    
    if (is.null(obj_rt) || length(obj_rt$y_q) < 2) {
      print_if_verbose(
        "Skipped vintage: not enough published quarterly information.\n",
        verbose = verbose
      )
      next
    }
    
    X_lf_cut <- obj_rt$X_lf[1:obj_rt$T_q, , drop = FALSE]
    y_q_cut  <- obj_rt$y_q
    X_hf_cut <- obj_rt$X_hf
    
    out <- MF_TPRF(
      X_lf        = X_lf_cut,
      X_hf        = X_hf_cut,
      y_q         = y_q_cut,
      Lproxy      = active_hyper$Lproxy,
      L_midas     = active_hyper$L_midas,
      p_AR        = active_hyper$p_AR,
      Robust_F    = params$Robust_F,
      alpha       = params$alpha,
      robust_type = params$robust_type,
      nw_lag      = params$nw_lag
    )
    
    y_rt_last <- tail(out$y_nowcast, 1)
    
    m_tr <- compute_m_tr(date_t, dates_q)
    if (is.na(m_tr)) {
      print_if_verbose(
        "Skipped vintage: invalid month-within-quarter mapping.\n",
        verbose = verbose
      )
      next
    }
    
    key <- as.character(date_t)
    
    if (m_tr == 1L) now_M1[[key]] <- y_rt_last
    if (m_tr == 2L) now_M2[[key]] <- y_rt_last
    if (m_tr == 3L) now_M3[[key]] <- y_rt_last
    
    vintage_log[[length(vintage_log) + 1L]] <- list(
      vintage_date   = date_t,
      active_regime  = active_regime,
      T_q            = obj_rt$T_q,
      month_position = m_tr,
      r              = active_hyper$r,
      Lproxy         = active_hyper$Lproxy,
      p_AR           = active_hyper$p_AR,
      L_midas        = active_hyper$L_midas
    )
    
    print_if_verbose(
      "Vintage completed successfully.\n",
      "  - Published quarters used : ", obj_rt$T_q, "\n",
      "  - Month position          : M", m_tr, "\n",
      verbose = verbose
    )
  }
  
  # --------------------------------------------------------------------------
  # Step 5. Return output using the original naming convention
  # --------------------------------------------------------------------------
  list(
    hyper_pre  = list(
      Lproxy = hyper_pre$Lproxy,
      L_midas = hyper_pre$L_midas,
      p_AR = hyper_pre$p_AR,
      r = hyper_pre$r
    ),
    hyper_post = list(
      Lproxy = hyper_post$Lproxy,
      L_midas = hyper_post$L_midas,
      p_AR = hyper_post$p_AR,
      r = hyper_post$r,
      t_recalib = if (do_recalib) dates[t_recalib] else NA
    ),
    M1 = now_M1,
    M2 = now_M2,
    M3 = now_M3,
    vintage_log = vintage_log
  )
}