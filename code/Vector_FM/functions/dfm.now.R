# current_t <- 238 M1(q0)
# current_t <- 239 M2(q0)
# current_t <- 240 M3(q0)
# current_t <- 241 M1(q1)

unbalancedness <- function(X_full, dates, Freq, Unb, current_t) {
  
  stopifnot(nrow(X_full) >= current_t)
  X_cut <- X_full[1:current_t, , drop = FALSE]
  
  for (j in seq_len(ncol(X_cut))) {
    
    delay_j  <- Unb[j]          # ritardo di pubblicazione (in mesi)
    release_t <- current_t - delay_j
    
    if (release_t < 1) {
      # la variabile NON è mai stata pubblicata entro current_t
      X_cut[, j] <- NA
    } else {
      # la variabile è disponibile solo fino a release_t
      if (release_t + 1 <= current_t) {
        X_cut[(release_t + 1):current_t, j] <- NA
      }
    }
  }
  
  return(X_cut)
}


compute_m_tr <- function(date_t, dates_q) {
  
  # trimestri chiusi prima del mese corrente
  idx_past <- which(dates_q < date_t)
  if (length(idx_past) == 0) return(NA_integer_)
  
  last_q_idx <- max(idx_past)
  
  # se l'ultimo trimestre chiuso è anche l'ultimo disponibile, non c'è un target quarter
  if (last_q_idx == length(dates_q)) return(NA_integer_)
  
  # trimestre target = successivo
  # mesi del trimestre target
  M1 <- dates_q[last_q_idx] %m+% months(1)
  M2 <- dates_q[last_q_idx] %m+% months(2)
  M3 <- dates_q[last_q_idx] %m+% months(3)
  
  if (date_t == M1) return(1L)
  if (date_t == M2) return(2L)
  if (date_t == M3) return(3L)
  
  NA_integer_
}

pseudo_realtime_DFM_EM_reestimate <- function(
    X_full,        # T x N raw (NOT standardized), includes GDP
    NQ,            # number of quarterly variables (at the end of X_full)
    params,        # list with at least start_eval, end_eval, restr, Kmax, kappa
    dates_m,       # monthly Date vector (length T)
    dates_q,       # quarter-end Date vector
    Freq,          # "M"/"Q" per series (length N)
    Unb,           # publication delays (in months), length N
    gdp_col,       # GDP column index in X_full (matches row in C)
    agg,           # length N (1 stock / 2 flow) OR length NQ (only quarterly)
    max_iter_em = 50,
    tol_em      = 1e-4,
    pmax        = 9,
    qmax        = 9,
    min_est_T   = 24
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
    for (j in 1:h) Ap <- Ap %*% A
    Ap
  }
  
  # ---------------------------------------------------------------------------
  # 0) Evaluation window
  # ---------------------------------------------------------------------------
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("start_eval / end_eval not found in dates_m.")
  }
  
  # ---------------------------------------------------------------------------
  # 0b) Pre-start_eval sample for hyper-selection (r_fix, p_fix, q_fix)
  # ---------------------------------------------------------------------------
  t_est_end <- max(which(dates_m < params$start_eval))
  if (length(t_est_end) == 0 || t_est_end < min_est_T) {
    stop("Estimation sample too short or not defined.")
  }
  
  cat("\n>>> DFM HYPER-SELECTION using data up to", as.character(dates_m[t_est_end]), "\n")
  
  # ragged-edge cut at t_est_end
  X_cut_est <- unbalancedness(
    X_full    = X_full,
    dates     = dates_m,
    Freq      = Freq,
    Unb       = Unb,
    current_t = t_est_end
  )
  
  # standardize with NA (full-sample pipeline)
  out_std_est <- standardize_with_na(X_cut_est[1:t_est_end, , drop = FALSE])
  X_std_est   <- out_std_est$X_std
  if (nrow(X_std_est) < 6) stop("Too few observations in estimation sample after cut.")
  
  # r via all-purpose covariance + ER (full-sample pipeline)
  cov_proxy_out <- all_purpose_covariance(X_std_est)
  ER_proxy_out  <- select_num_factors_ER(cov_proxy_out$Sigma_tilde, Kmax = params$Kmax)
  r_fix         <- ER_proxy_out$r
  if (is.null(r_fix) || r_fix < 1) stop("r_fix not identified in estimation sample.")
  
  # p via IC on XP factors (full-sample pipeline)
  xp_est <- estimate_factors_XP(X_std_est, r = r_fix)
  ic_est <- select_p_var_ic(xp_est$F_hat, pmax = pmax, include_const = TRUE)
  p_fix  <- ic_est$p_BIC
  
  # q via global AR IC on idio residuals (full-sample pipeline)
  Ehat_est <- X_std_est - xp_est$C_hat  # keeps NA where X_std_est is NA
  ic_q_est <- select_p_ar_ic_global(Ehat_est, pmax = qmax, include_const = FALSE, min_T = 30)
  q_fix    <- ic_q_est$p_BIC
  
  cat("# [Hyper] r_fix =", r_fix, "| p_fix =", p_fix, "| q_fix =", q_fix, "\n")
  
  # ---------------------------------------------------------------------------
  # Storage
  # ---------------------------------------------------------------------------
  now_M1_std  <- list(); now_M2_std  <- list(); now_M3_std  <- list()
  now_M1_orig <- list(); now_M2_orig <- list(); now_M3_orig <- list()
  
  # ---------------------------------------------------------------------------
  # 1) Rolling loop
  # ---------------------------------------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    cat("\n>>> DFM REAL-TIME at", as.character(date_t),
        "| r_fix =", r_fix, "| p_fix =", p_fix, "| q_fix =", q_fix, "\n")
    
    # 1) Ragged-edge (publication delays)
    # NOTE: keep X_full's own "structural NA" (e.g., quarterly only at quarter-end).
    X_cut <- unbalancedness(
      X_full    = X_full,
      dates     = dates_m,
      Freq      = Freq,
      Unb       = Unb,
      current_t = tt
    )
    X_tt <- X_cut[1:tt, , drop = FALSE]
    
    # 2) Standardize (with NA) + keep mean/sd for destandardization (pseudo real-time)
    out_std_tt <- standardize_with_na(X_tt)
    X_std_tt   <- out_std_tt$X_std
    
    T_tt <- nrow(X_std_tt)
    N    <- ncol(X_std_tt)
    if (T_tt < 6) next
    
    # 3) Dimensions
    NQ_tt <- NQ
    NM_tt <- N - NQ_tt
    if (NM_tt < 1 || NQ_tt < 1) stop("Check NQ: must have both monthly and quarterly series.")
    
    # 4) agg normalization (allow agg length N or NQ)
    aggN_tt <- make_aggN(agg, NM = NM_tt, NQ = NQ_tt)
    
    # 5) Initialization (must match full-sample InitialCond signature/logic)
    #    IMPORTANT: InitialCond expects X_std with NA (NOT imputed by user).
    Init_tt <- InitialCond(
      X_std    = X_std_tt,
      r        = r_fix,
      p_factor = p_fix,
      q_idio   = q_fix,
      NM       = NM_tt,
      NQ       = NQ_tt,
      restr    = params$restr,
      agg      = aggN_tt,
      kappa    = params$kappa
    )
    
    # 6) EM estimation on TRUE ragged-edge standardized data (same as full-sample)
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
    
    # 7) Kalman FILTER (real-time): we want a_{t|t} at the last available month tt
    #    Best practice: use filtered state (fu). If your kalman.na does not return it,
    #    either (i) add fu/Pu to kalman.na, or (ii) use fs at last time (usually equals).
    kf <- kalman.na(
      par   = list(t = A, l = C, q = Qm, h = Rk),
      y     = t(X_std_tt),   # N x T
      k     = k_state,
      start = list(f = Z0, P = V0),
      O     = lapply(1:T_tt, function(t) which(!is.na(X_std_tt[t, ])))
    )
    
    # Prefer filtered if available, else fallback to smoothed at final time
    if (!is.null(kf$fu)) {
      a_T <- kf$fu[, T_tt, drop = FALSE]     # a_{t|t}
    } else {
      a_T <- kf$fs[, T_tt, drop = FALSE]     # fallback (at final time should match)
    }
    
    # 8) Identify month-in-quarter (M1/M2/M3) of the TARGET quarter
    m_tr <- compute_m_tr(date_t, dates_q)
    if (is.na(m_tr)) next
    
    h <- 3 - m_tr
    if (h < 0) next
    
    # 9) h-step projection z_{t+h|t} = A^h a_{t|t}
    Ap <- A_power_h(A, h)
    z_pred <- Ap %*% a_T
    
    # 10) GDP nowcast (STANDARDIZED): yhat = c_gdp' z_pred
    yhat_gdp_std <- as.numeric(C[gdp_col, , drop = FALSE] %*% z_pred)
    
    # 11) GDP nowcast (ORIGINAL SCALE) using rolling mean/sd up to tt
    mu_gdp <- out_std_tt$mean[gdp_col]
    sd_gdp <- out_std_tt$sd[gdp_col]
    yhat_gdp_orig <- mu_gdp + sd_gdp * yhat_gdp_std
    
    # 12) Store by month (key = current month date)
    key <- as.character(date_t)
    if (m_tr == 1) { now_M1_std[[key]]  <- yhat_gdp_std;  now_M1_orig[[key]] <- yhat_gdp_orig }
    if (m_tr == 2) { now_M2_std[[key]]  <- yhat_gdp_std;  now_M2_orig[[key]] <- yhat_gdp_orig }
    if (m_tr == 3) { now_M3_std[[key]]  <- yhat_gdp_std;  now_M3_orig[[key]] <- yhat_gdp_orig }
  }
  
  return(list(
    r_fix = r_fix,
    p_fix = p_fix,
    q_fix = q_fix,
    M1_std  = unlist(now_M1_std),
    M2_std  = unlist(now_M2_std),
    M3_std  = unlist(now_M3_std),
    M1_orig = unlist(now_M1_orig),
    M2_orig = unlist(now_M2_orig),
    M3_orig = unlist(now_M3_orig)
  ))
}
