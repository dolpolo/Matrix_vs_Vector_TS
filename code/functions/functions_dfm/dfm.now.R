# ==============================================================================
# VECTOR DFM - PSEUDO REAL-TIME NOWCASTING WITH REGIME-SPECIFIC HYPERPARAMETERS
# ------------------------------------------------------------------------------
# This file implements:
#   1. Real-time publication delays (unbalancedness)
#   2. Month-in-quarter identification
#   3. Regime-specific DFM hyperparameter selection
#   4. Expanding pseudo real-time nowcasting with:
#        - one PRE-COVID hyperparameter calibration
#        - one POST-COVID hyperparameter recalibration
#   5. Monthly nowcasts stored as M1 / M2 / M3 on:
#        - standardized scale
#        - original scale
# ==============================================================================


# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

unbalancedness <- function(X_full, dates, Freq, Unb, current_t) {
  # Build the pseudo-real-time panel up to current_t by applying publication
  # delays to each variable.
  #
  # X_full    : T x N raw matrix (may already contain structural NA)
  # dates     : monthly dates vector of length T
  # Freq      : variable frequency label ("M"/"Q"), length N
  # Unb       : publication delay in months, length N
  # current_t : current monthly vintage index
  
  stopifnot(nrow(X_full) >= current_t)
  
  X_cut <- X_full[1:current_t, , drop = FALSE]
  
  for (j in seq_len(ncol(X_cut))) {
    
    delay_j   <- Unb[j]
    release_t <- current_t - delay_j
    
    if (release_t < 1) {
      # Variable not yet released at this vintage
      X_cut[, j] <- NA
    } else {
      # Variable available only up to release_t
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
  # Returns:
  #   1 = first month of the target quarter
  #   2 = second month of the target quarter
  #   3 = third month of the target quarter
  #
  # The target quarter is the quarter immediately following the last published
  # quarter before date_t.
  
  idx_past <- which(dates_q < date_t)
  if (length(idx_past) == 0) return(NA_integer_)
  
  last_q_idx <- max(idx_past)
  
  # If the last published quarter is already the last available one,
  # there is no future target quarter to nowcast.
  if (last_q_idx == length(dates_q)) return(NA_integer_)
  
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
# SELECT OR ACCEPT DFM HYPERPARAMETERS FOR ONE REGIME
# ==============================================================================

select_regime_hyperparameters_dfm <- function(
    calibration_t,
    regime_label,
    X_full,
    params,
    dates_m,
    Freq,
    Unb,
    NQ,
    user_hyper = list(r = NULL, p = NULL, q = NULL),
    pmax = 9,
    qmax = 9,
    min_est_T = 24,
    verbose = TRUE
) {
  # Select DFM hyperparameters once for a given regime, using the same logic
  # as in the full-sample DFM:
  #   - r via covariance proxy + Eigenvalue Ratio
  #   - p via BIC on preliminary XP factors
  #   - q via BIC on idiosyncratic residuals
  
  print_if_verbose(
    "\n============================================================\n",
    "Selecting DFM hyperparameters for regime: ", regime_label, "\n",
    "Calibration vintage: ", as.character(dates_m[calibration_t]), "\n",
    "============================================================\n",
    verbose = verbose
  )
  
  if (calibration_t < min_est_T) {
    stop("Calibration sample too short.")
  }
  
  # --------------------------------------------------------------------------
  # Step 1. Build pseudo-real-time calibration panel
  # --------------------------------------------------------------------------
  X_cut <- unbalancedness(
    X_full    = X_full,
    dates     = dates_m,
    Freq      = Freq,
    Unb       = Unb,
    current_t = calibration_t
  )
  
  X_est <- X_cut[1:calibration_t, , drop = FALSE]
  
  # --------------------------------------------------------------------------
  # Step 2. Standardize with NA preserved
  # --------------------------------------------------------------------------
  out_std <- standardize_with_na(X_est)
  X_std   <- out_std$X_std
  
  if (nrow(X_std) < 6) {
    stop("Too few observations after pseudo-real-time cut.")
  }
  
  N  <- ncol(X_std)
  NM <- N - NQ
  
  if (NM < 1 || NQ < 1) {
    stop("Invalid dimensions: both monthly and quarterly variables are required.")
  }
  
  # --------------------------------------------------------------------------
  # Step 3. Select or accept r
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$r)) {
    
    r_use <- user_hyper$r
    cov_proxy_out <- NULL
    ER_proxy_out  <- NULL
    
    print_if_verbose(
      "Number of factors r provided by user: ", r_use, "\n",
      verbose = verbose
    )
    
  } else {
    
    cov_proxy_out <- all_purpose_covariance(X_std)
    ER_proxy_out  <- select_num_factors_ER(
      Sigma       = cov_proxy_out$Sigma_tilde,
      Kmax        = params$Kmax
    )
    
    r_use <- ER_proxy_out$r
    
    if (is.null(r_use) || r_use < 1) {
      stop("Could not identify a valid r at calibration date.")
    }
    
    print_if_verbose(
      "Number of factors r selected internally: ", r_use, "\n",
      verbose = verbose
    )
  }
  
  # --------------------------------------------------------------------------
  # Step 4. Preliminary XP factor extraction
  # --------------------------------------------------------------------------
  xp <- estimate_factors_XP(X_std, r = r_use)
  
  # --------------------------------------------------------------------------
  # Step 5. Select or accept p
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$p)) {
    
    p_use <- user_hyper$p
    factor_var_ic <- NULL
    
    print_if_verbose(
      "Factor VAR order p provided by user: ", p_use, "\n",
      verbose = verbose
    )
    
  } else {
    
    factor_var_ic <- select_p_var_ic(
      F_hat         = xp$F_hat,
      pmax          = pmax,
      include_const = TRUE
    )
    
    p_use <- factor_var_ic$p_BIC
    
    print_if_verbose(
      "Factor VAR order p selected internally: ", p_use, "\n",
      verbose = verbose
    )
  }
  
  # --------------------------------------------------------------------------
  # Step 6. Select or accept q
  # --------------------------------------------------------------------------
  if (!is.null(user_hyper$q)) {
    
    q_use <- user_hyper$q
    idio_ar_ic <- NULL
    
    print_if_verbose(
      "Idiosyncratic AR order q provided by user: ", q_use, "\n",
      verbose = verbose
    )
    
  } else {
    
    Ehat <- X_std - xp$C_hat
    
    idio_ar_ic <- select_p_ar_ic_global(
      Ehat[, 1:NM, drop = FALSE],
      pmax          = qmax,
      include_const = FALSE,
      min_T         = 30
    )
    
    q_use <- idio_ar_ic$p_BIC
    
    print_if_verbose(
      "Idiosyncratic AR order q selected internally: ", q_use, "\n",
      verbose = verbose
    )
  }
  
  print_if_verbose(
    "\nFrozen DFM hyperparameters for regime ", regime_label, ":\n",
    "  - Calibration date : ", as.character(dates_m[calibration_t]), "\n",
    "  - r                : ", r_use, "\n",
    "  - p                : ", p_use, "\n",
    "  - q                : ", q_use, "\n\n",
    verbose = verbose
  )
  
  return(list(
    regime_label     = regime_label,
    calibration_t    = calibration_t,
    calibration_date = dates_m[calibration_t],
    r                = r_use,
    p                = p_use,
    q                = q_use,
    out_std          = out_std,
    xp               = xp,
    cov_proxy_out    = cov_proxy_out,
    ER_proxy_out     = ER_proxy_out,
    factor_var_ic    = factor_var_ic,
    idio_ar_ic       = idio_ar_ic
  ))
}


# ==============================================================================
# EXPANDING PSEUDO REAL-TIME DFM WITH FIXED HYPERPARAMETERS WITHIN REGIME
# ==============================================================================

pseudo_realtime_DFM_EM_reestimate <- function(
    X_full,        # T x N raw panel (NOT standardized), includes GDP
    NQ,            # number of quarterly variables (at the end of X_full)
    params,        # must contain at least start_eval, end_eval, covid_end, restr, Kmax, kappa
    dates_m,       # monthly Date vector, length T
    dates_q,       # quarter-end Date vector
    Freq,          # length N, "M"/"Q"
    Unb,           # length N, publication delays
    gdp_col,       # GDP column index in X_full / C
    agg,           # length N or length NQ
    do_post_covid_recalibration = TRUE,
    user_hyper_pre  = list(r = NULL, p = NULL, q = NULL),
    user_hyper_post = list(r = NULL, p = NULL, q = NULL),
    max_iter_em = 50,
    tol_em      = 1e-4,
    pmax        = 9,
    qmax        = 9,
    min_est_T   = 24,
    verbose     = TRUE
) {
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  make_aggN <- function(agg, NM, NQ) {
    N <- NM + NQ
    
    if (is.null(agg)) return(NULL)
    if (length(agg) == N)  return(as.integer(agg))
    if (length(agg) == NQ) return(as.integer(c(rep(NA_integer_, NM), agg)))
    
    stop("agg must have length N or length NQ.")
  }
  
  A_power_h <- function(A, h) {
    k <- nrow(A)
    if (h <= 0) return(diag(k))
    
    Ap <- diag(k)
    for (j in seq_len(h)) {
      Ap <- Ap %*% A
    }
    Ap
  }
  
  # ---------------------------------------------------------------------------
  # 0. Evaluation window
  # ---------------------------------------------------------------------------
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("start_eval / end_eval not found in dates_m.")
  }
  if (t_start > t_end) {
    stop("start_eval must be <= end_eval.")
  }
  
  # ---------------------------------------------------------------------------
  # 0b. PRE-COVID calibration date
  # ---------------------------------------------------------------------------
  t_est_end <- max(which(dates_m < params$start_eval))
  
  if (length(t_est_end) == 0 || t_est_end < min_est_T) {
    stop("Estimation sample too short or not defined.")
  }
  
  # ---------------------------------------------------------------------------
  # 0c. POST-COVID recalibration date
  # ---------------------------------------------------------------------------
  t_recalib <- which(dates_m >= params$covid_end)[1]
  
  do_recalib <- isTRUE(do_post_covid_recalibration) &&
    !is.na(t_recalib) &&
    t_recalib <= t_end
  
  print_if_verbose(
    "\n############################################################\n",
    "EXPANDING PSEUDO REAL-TIME VECTOR DFM\n",
    "Model complexity is fixed within each regime\n",
    "############################################################\n",
    "Evaluation window:\n",
    "  - Start: ", as.character(dates_m[t_start]), "\n",
    "  - End  : ", as.character(dates_m[t_end]), "\n",
    "Pre-COVID calibration date : ", as.character(dates_m[t_est_end]), "\n",
    "Post-COVID recalibration   : ", if (do_recalib) as.character(dates_m[t_recalib]) else "disabled", "\n\n",
    verbose = verbose
  )
  
  # ---------------------------------------------------------------------------
  # 1. Select or accept PRE-COVID hyperparameters
  # ---------------------------------------------------------------------------
  hyper_pre <- select_regime_hyperparameters_dfm(
    calibration_t = t_est_end,
    regime_label  = "PRE-COVID",
    X_full        = X_full,
    params        = params,
    dates_m       = dates_m,
    Freq          = Freq,
    Unb           = Unb,
    NQ            = NQ,
    user_hyper    = user_hyper_pre,
    pmax          = pmax,
    qmax          = qmax,
    min_est_T     = min_est_T,
    verbose       = verbose
  )
  
  # ---------------------------------------------------------------------------
  # 2. Select or accept POST-COVID hyperparameters
  # ---------------------------------------------------------------------------
  hyper_post <- hyper_pre
  
  if (do_recalib) {
    hyper_post <- select_regime_hyperparameters_dfm(
      calibration_t = t_recalib,
      regime_label  = "POST-COVID",
      X_full        = X_full,
      params        = params,
      dates_m       = dates_m,
      Freq          = Freq,
      Unb           = Unb,
      NQ            = NQ,
      user_hyper    = user_hyper_post,
      pmax          = pmax,
      qmax          = qmax,
      min_est_T     = min_est_T,
      verbose       = verbose
    )
  }
  
  # ---------------------------------------------------------------------------
  # 3. Storage containers
  # ---------------------------------------------------------------------------
  now_M1_std  <- list()
  now_M2_std  <- list()
  now_M3_std  <- list()
  
  now_M1_orig <- list()
  now_M2_orig <- list()
  now_M3_orig <- list()
  
  vintage_log <- list()
  
  # ---------------------------------------------------------------------------
  # 4. Expanding real-time loop with fixed hyperparameters within regime
  # ---------------------------------------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    
    use_post <- do_recalib && tt >= t_recalib
    active_hyper  <- if (use_post) hyper_post else hyper_pre
    active_regime <- if (use_post) "POST-COVID" else "PRE-COVID"
    
    print_if_verbose(
      "\n------------------------------------------------------------\n",
      "DFM real-time vintage: ", as.character(date_t), "\n",
      "Active regime        : ", active_regime, "\n",
      "Frozen hyperparameters used:\n",
      "  - r : ", active_hyper$r, "\n",
      "  - p : ", active_hyper$p, "\n",
      "  - q : ", active_hyper$q, "\n",
      verbose = verbose
    )
    
    # ------------------------------------------------------------------------
    # Step 4.1. Apply publication delays
    # ------------------------------------------------------------------------
    X_cut <- unbalancedness(
      X_full    = X_full,
      dates     = dates_m,
      Freq      = Freq,
      Unb       = Unb,
      current_t = tt
    )
    
    X_tt <- X_cut[1:tt, , drop = FALSE]
    
    # ------------------------------------------------------------------------
    # Step 4.2. Standardize with NA preserved
    # ------------------------------------------------------------------------
    out_std_tt <- standardize_with_na(X_tt)
    X_std_tt   <- out_std_tt$X_std
    
    T_tt <- nrow(X_std_tt)
    N    <- ncol(X_std_tt)
    
    if (T_tt < 6) {
      print_if_verbose("Skipped vintage: too few observations.\n", verbose = verbose)
      next
    }
    
    # ------------------------------------------------------------------------
    # Step 4.3. Dimensions
    # ------------------------------------------------------------------------
    NQ_tt <- NQ
    NM_tt <- N - NQ_tt
    
    if (NM_tt < 1 || NQ_tt < 1) {
      stop("Check NQ: both monthly and quarterly blocks are required.")
    }
    
    # ------------------------------------------------------------------------
    # Step 4.4. Normalize aggregation vector
    # ------------------------------------------------------------------------
    aggN_tt <- make_aggN(agg, NM = NM_tt, NQ = NQ_tt)
    
    # ------------------------------------------------------------------------
    # Step 4.5. Initialize DFM with regime-specific hyperparameters
    # ------------------------------------------------------------------------
    Init_tt <- InitialCond(
      X_std    = X_std_tt,
      r        = active_hyper$r,
      p_factor = active_hyper$p,
      q_idio   = active_hyper$q,
      NM       = NM_tt,
      NQ       = NQ_tt,
      restr    = params$restr,
      agg      = aggN_tt,
      kappa    = params$kappa
    )
    
    # ------------------------------------------------------------------------
    # Step 4.6. EM estimation on pseudo-real-time panel
    # ------------------------------------------------------------------------
    res_tt <- DFM_EM(
      X        = X_std_tt,
      Init     = Init_tt,
      max_iter = max_iter_em,
      tol      = tol_em
    )
    
    A  <- res_tt$A
    C  <- res_tt$C
    Qm <- res_tt$Q
    Rk <- res_tt$R
    Z0 <- res_tt$Z0
    V0 <- res_tt$V0
    
    k_state <- nrow(A)
    
    # ------------------------------------------------------------------------
    # Step 4.7. Kalman filtering at vintage t
    # ------------------------------------------------------------------------
    kf <- kalman.na(
      par   = list(t = A, l = C, q = Qm, h = Rk),
      y     = t(X_std_tt),
      k     = k_state,
      start = list(f = Z0, P = V0),
      O     = lapply(seq_len(T_tt), function(t) which(!is.na(X_std_tt[t, ])))
    )
    
    # Prefer filtered state if available, otherwise use smoothed final state
    if (!is.null(kf$fu)) {
      a_T <- kf$fu[, T_tt, drop = FALSE]
    } else {
      a_T <- kf$fs[, T_tt, drop = FALSE]
    }
    
    # ------------------------------------------------------------------------
    # Step 4.8. Identify target-quarter month position
    # ------------------------------------------------------------------------
    m_tr <- compute_m_tr(date_t, dates_q)
    if (is.na(m_tr)) {
      print_if_verbose("Skipped vintage: invalid month-within-quarter mapping.\n", verbose = verbose)
      next
    }
    
    h <- 3 - m_tr
    if (h < 0) next
    
    # ------------------------------------------------------------------------
    # Step 4.9. h-step state projection
    # ------------------------------------------------------------------------
    Ap     <- A_power_h(A, h)
    z_pred <- Ap %*% a_T
    
    # ------------------------------------------------------------------------
    # Step 4.10. GDP nowcast on standardized scale
    # ------------------------------------------------------------------------
    yhat_gdp_std <- as.numeric(C[gdp_col, , drop = FALSE] %*% z_pred)
    
    # ------------------------------------------------------------------------
    # Step 4.11. GDP nowcast on original scale
    # ------------------------------------------------------------------------
    mu_gdp <- out_std_tt$mean[gdp_col]
    sd_gdp <- out_std_tt$sd[gdp_col]
    yhat_gdp_orig <- mu_gdp + sd_gdp * yhat_gdp_std
    
    # ------------------------------------------------------------------------
    # Step 4.12. Store monthly nowcast
    # ------------------------------------------------------------------------
    key <- as.character(date_t)
    
    if (m_tr == 1L) {
      now_M1_std[[key]]  <- yhat_gdp_std
      now_M1_orig[[key]] <- yhat_gdp_orig
    }
    if (m_tr == 2L) {
      now_M2_std[[key]]  <- yhat_gdp_std
      now_M2_orig[[key]] <- yhat_gdp_orig
    }
    if (m_tr == 3L) {
      now_M3_std[[key]]  <- yhat_gdp_std
      now_M3_orig[[key]] <- yhat_gdp_orig
    }
    
    vintage_log[[length(vintage_log) + 1L]] <- list(
      vintage_date   = date_t,
      active_regime  = active_regime,
      month_position = m_tr,
      r              = active_hyper$r,
      p              = active_hyper$p,
      q              = active_hyper$q
    )
    
    print_if_verbose(
      "Vintage completed successfully.\n",
      "  - Month position : M", m_tr, "\n",
      "  - GDP nowcast    : ", round(yhat_gdp_orig, 6), " (orig scale)\n",
      verbose = verbose
    )
  }
  
  # ---------------------------------------------------------------------------
  # 5. Return output in MF-TPRF-like structure
  # ---------------------------------------------------------------------------
  return(list(
    hyper_pre = list(
      r = hyper_pre$r,
      p = hyper_pre$p,
      q = hyper_pre$q
    ),
    
    hyper_post = list(
      r         = hyper_post$r,
      p         = hyper_post$p,
      q         = hyper_post$q,
      t_recalib = if (do_recalib) dates_m[t_recalib] else NA
    ),
    
    M1_std  = unlist(now_M1_std),
    M2_std  = unlist(now_M2_std),
    M3_std  = unlist(now_M3_std),
    
    M1_orig = unlist(now_M1_orig),
    M2_orig = unlist(now_M2_orig),
    M3_orig = unlist(now_M3_orig),
    
    vintage_log = vintage_log
  ))
}