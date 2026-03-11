# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

# X_full <- X
# current_t <- 303
# current_t <- 304
# current_t <- 305

unbalancedness_tensor <- function(X_full, Unb, current_t) {
  # X_full: array [T_m x P x K]
  # Unb   : matrice [P x K] con ritardi (mesi) per ciascun (paese, variabile)
  # current_t: indice del mese corrente (1..T_m)
  
  stopifnot(length(dim(X_full)) == 3)
  T_m <- dim(X_full)[1]
  P   <- dim(X_full)[2]
  K   <- dim(X_full)[3]
  
  stopifnot(T_m >= current_t)
  stopifnot(all(dim(Unb) == c(P, K)))
  
  # taglio nel tempo: 1..current_t
  X_cut <- X_full[1:current_t, , , drop = FALSE]   # [current_t x P x K]
  
  for (p in seq_len(P)) {
    for (k in seq_len(K)) {
      delay_pk <- Unb[p, k]   # ritardo in mesi
      
      # se NA o negativo → non tocco quella serie
      if (is.na(delay_pk) || delay_pk < 0) next
      
      delay_pk  <- as.integer(delay_pk)
      release_t <- current_t - delay_pk
      
      if (release_t < 1) {
        # mai osservata entro current_t → tutta NA
        X_cut[ , p, k] <- NA_real_
      } else if (release_t < current_t) {
        # disponibile solo fino a release_t
        X_cut[(release_t + 1):current_t, p, k] <- NA_real_
      } else {
        # release_t >= current_t → nessun taglio
      }
    }
  }
  
  X_cut
}



compute_m_tr <- function(date_t, dates_q) {
  # ultimo trimestre rilasciato PRIMA del mese corrente
  last_q_date <- max(dates_q[dates_q < date_t])
  last_q_idx  <- which(dates_q == last_q_date)
  
  M1 <- dates_q[last_q_idx] %m+% months(1)
  M2 <- dates_q[last_q_idx] %m+% months(2)
  M3 <- dates_q[last_q_idx] %m+% months(3)
  
  if (date_t == M1) return(1)
  if (date_t == M2) return(2)
  if (date_t == M3) return(3)
  return(NA_integer_)
}
# ==============================================================================
# PSEUDO REAL-TIME (EXPANDING) — Matrix MF-TPRF
#   Coherent with full-sample pipeline:
#     standardize_mat_with_na -> init_CL_Yu -> aggregate_tensor_to_quarterly
#     -> Lproxy (Kraemer) -> (p_AR, L_midas) grid BIC -> Tensor_MF_TPRF
#   Two regimes:
#     Hyper #1: up to start_eval-1
#     Hyper #2: once at first date >= covid_end
# ==============================================================================

pseudo_realtime_Tensor_MF_TPRF <- function(
    X_full,            # monthly predictors tensor: [T_m x P x K] (GDP excluded)
    Y_q_all_full,      # quarterly targets matrix: [T_q_full x (1+P_eval)] (EA + countries)
    proxy_name = "EA",
    params,
    dates_m, dates_q,  # monthly & quarterly date vectors
    W_full,            # missingness mask for predictors (same dim as X_full): 1 obs, 0 miss
    agg,               # aggregation rule for predictors (same shape used by aggregate_tensor_to_quarterly)
    N_m, N_q,           # number of monthly & quarterly predictors (GDP excluded)
    # controls
    do_recalib = TRUE,
    ils_maxit = 100, ils_tol = 1e-8,
    verbose = TRUE
) {
  
  # -----------------------------
  # 0) checks + eval window
  # -----------------------------
  stopifnot(length(dim(X_full)) == 3)
  T_m <- dim(X_full)[1]
  P   <- dim(X_full)[2]
  K   <- dim(X_full)[3]
  if (K != (N_m + N_q)) stop("K must equal N_m + N_q (GDP excluded).")
  if (!all(dim(W_full) == dim(X_full))) stop("W_full must have same dim as X_full.")
  if (!(proxy_name %in% colnames(Y_q_all_full))) stop("proxy_name not in Y_q_all_full colnames.")
  
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  if (length(t_start) == 0 || length(t_end) == 0) stop("start_eval/end_eval not found in dates_m.")
  if (t_start > t_end) stop("start_eval > end_eval.")
  
  t_est_end <- max(which(dates_m < params$start_eval))
  if (!is.finite(t_est_end) || t_est_end < 24) stop("Estimation sample too short (pre-eval).")
  
  # first date >= covid_end (for recalibration)
  t_recalib <- which(dates_m >= params$covid_end)[1]
  do_recal <- isTRUE(do_recalib) && !is.na(t_recalib) && (t_recalib <= t_end)
  
  countries_eval <- setdiff(colnames(Y_q_all_full), proxy_name)
  if (length(countries_eval) < 1) stop("No countries beyond proxy in Y_q_all_full.")
  
  # helper to print
  say <- function(...) if (isTRUE(verbose)) cat(...)
  
  # --------------------------------------------------
  # CORE BUILDER: given current_t (month index), build
  #   X_cl (monthly imputed), X_cl_q (quarterly agg), r_hat, and proxy y_EA_q_cut
  # --------------------------------------------------
  build_data_up_to_t <- function(current_t) {
    
    # cut predictors to current_t (expanding)
    X_cut <- X_full[1:current_t, , , drop = FALSE]
    W_cut <- W_full[1:current_t, , , drop = FALSE]
    
    # standardize with NA (uses mask)
    out_std <- standardize_mat_with_na(X_cut)
    X_std   <- out_std$X_scaled
    
    # Cen&Lam + Yu init (imputation + r selection)
    imp_cl <- init_CL_Yu(X_std, W_cut, params$Kmax)
    X_cl   <- imp_cl$Y_init
    r_hat  <- imp_cl$r
    
    # aggregate to quarterly
    X_cl_q <- aggregate_tensor_to_quarterly(
      X_tens = X_cl,
      agg    = agg,
      N_m    = N_m,
      N_q    = N_q
    )
    
    # quarters published up to date_t: use dates_q < dates_m[current_t]
    date_t <- dates_m[current_t]
    idx_pub <- which(dates_q < date_t)
    if (length(idx_pub) < 2) return(NULL)
    
    T_q_use <- min(max(idx_pub), dim(X_cl_q)[1], nrow(Y_q_all_full))
    if (T_q_use < 2) return(NULL)
    
    # cut targets and predictors to published quarters
    Y_cut <- Y_q_all_full[1:T_q_use, , drop = FALSE]
    X_lf  <- X_cl_q[1:T_q_use, , , drop = FALSE]
    X_hf  <- X_cl       # monthly up to current_t already
    
    list(
      out_std = out_std,
      imp_cl  = imp_cl,
      r_hat   = r_hat,
      T_q_use = T_q_use,
      X_lf    = X_lf,
      X_hf    = X_hf,
      Y_cut   = Y_cut
    )
  }
  
  # ==================================================
  # 1) HYPER #1 selection (pre-eval): up to start_eval-1
  # ==================================================
  say("\n>>> HYPER #1 selection using data up to ", as.character(dates_m[t_est_end]), "\n")
  est1 <- build_data_up_to_t(t_est_end)
  if (is.null(est1)) stop("Hyper #1: not enough quarterly GDP published.")
  
  # Lproxy via Kraemer on vectorized X_lf (quarterly)
  X_lf_vec1 <- tensor_to_vector(X_tens = est1$X_lf, N_m = N_m, N_q = N_q)
  y_proxy1  <- as.numeric(est1$Y_cut[, proxy_name])
  
  pls1   <- select_L_autoproxy_3prf(X_lf_vec1, y_proxy1, Zmax = params$Zmax)
  Lproxy_1 <- pls1$L_opt
  
  # (p_AR, L_midas) via BIC grid (matrix)
  lag1 <- choose_UMIDAS_grid_tensor_MF(
    X_lf = est1$X_lf,
    X_hf = est1$X_hf,
    y_q  = y_proxy1,
    r    = est1$r_hat,
    Lproxy   = Lproxy_1,
    Lmax     = params$Lmax,
    p_AR_max = params$p_AR_max,
    use_autoproxy = TRUE,
    y_proxy = y_proxy1,
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z   = TRUE,
    ils_maxit = ils_maxit,
    ils_tol   = ils_tol
  )
  
  L_midas_1 <- lag1$best_BIC$L
  p_ar_1    <- lag1$best_BIC$p_AR
  r_1       <- est1$r_hat
  
  say("# Hyper #1: Lproxy=", Lproxy_1,
      " | L_midas=", L_midas_1,
      " | p_AR=", p_ar_1,
      " | r=(", r_1[1], ",", r_1[2], ")\n")
  
  # ==================================================
  # 2) HYPER #2 selection (post-covid) once at first date >= covid_end
  # ==================================================
  Lproxy_2 <- Lproxy_1; L_midas_2 <- L_midas_1; p_ar_2 <- p_ar_1; r_2 <- r_1
  
  if (do_recal) {
    say("\n>>> HYPER #2 selection using data up to ", as.character(dates_m[t_recalib]), "\n")
    est2 <- build_data_up_to_t(t_recalib)
    
    if (!is.null(est2)) {
      X_lf_vec2 <- tensor_to_vector(X_tens = est2$X_lf, N_m = N_m, N_q = N_q)
      y_proxy2  <- as.numeric(est2$Y_cut[, proxy_name])
      
      pls2 <- select_L_autoproxy_3prf(X_lf_vec2, y_proxy2, Zmax = params$Zmax)
      Lproxy_2 <- pls2$L_opt
      
      lag2 <- choose_UMIDAS_grid_tensor_MF(
        X_lf = est2$X_lf,
        X_hf = est2$X_hf,
        y_q  = y_proxy2,
        r    = est2$r_hat,
        Lproxy   = Lproxy_2,
        Lmax     = params$Lmax,
        p_AR_max = params$p_AR_max,
        use_autoproxy = TRUE,
        y_proxy = y_proxy2,
        standardize_proxy = TRUE,
        orthonormalize_each_iter = TRUE,
        orthonormalize_final_Z   = TRUE,
        ils_maxit = ils_maxit,
        ils_tol   = ils_tol
      )
      
      L_midas_2 <- lag2$best_BIC$L
      p_ar_2    <- lag2$best_BIC$p_AR
      r_2       <- est2$r_hat
      
      say("# Hyper #2: Lproxy=", Lproxy_2,
          " | L_midas=", L_midas_2,
          " | p_AR=", p_ar_2,
          " | r=(", r_2[1], ",", r_2[2], ")\n")
    } else {
      say("(!) Hyper #2 skipped: not enough GDP published at recalib date.\n")
    }
  }
  
  # ==================================================
  # 3) REAL-TIME LOOP (expanding) with fixed hypers by regime
  # ==================================================
  now_M1 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M2 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M3 <- setNames(vector("list", length(countries_eval)), countries_eval)
  
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    
    use_post <- do_recal && (tt >= t_recalib)
    Lproxy_use <- if (use_post) Lproxy_2 else Lproxy_1
    L_midas_use <- if (use_post) L_midas_2 else L_midas_1
    p_ar_use <- if (use_post) p_ar_2 else p_ar_1
    r_use <- if (use_post) r_2 else r_1
    
    say("\n>>> REAL-TIME at ", as.character(date_t),
        " | regime=", if (use_post) "post-COVID" else "pre-COVID",
        " | Lproxy=", Lproxy_use,
        " | L_midas=", L_midas_use,
        " | p_AR=", p_ar_use,
        " | r=(", r_use[1], ",", r_use[2], ")\n")
    
    rt <- build_data_up_to_t(tt)
    if (is.null(rt)) next
    
    # run MF-TPRF (full-sample estimator on the cut dataset)
    out_rt <- Tensor_MF_TPRF(
      X_lf        = rt$X_lf,
      X_hf        = rt$X_hf,
      Y_q_all     = rt$Y_cut,
      proxy_name  = proxy_name,
      Lproxy      = Lproxy_use,
      L_midas     = L_midas_use,
      p_AR        = p_ar_use,
      r           = r_use,
      standardize_proxy = TRUE,
      orthonormalize_each_iter = TRUE,
      orthonormalize_final_Z   = TRUE,
      ils_maxit = ils_maxit,
      ils_tol   = ils_tol
    )
    
    # last nowcast month in this vintage
    m_tr <- compute_m_tr(date_t, dates_q)  # must return 1/2/3
    if (is.na(m_tr)) next
    key <- as.character(date_t)
    
    for (cc in countries_eval) {
      y_cc <- out_rt$by_country[[cc]]$y_nowcast
      if (is.null(y_cc) || length(y_cc) == 0) next
      y_last <- tail(y_cc, 1)
      
      if (m_tr == 1) now_M1[[cc]][[key]] <- y_last
      if (m_tr == 2) now_M2[[cc]][[key]] <- y_last
      if (m_tr == 3) now_M3[[cc]][[key]] <- y_last
    }
  }
  
  list(
    hyper_pre  = list(Lproxy = Lproxy_1, L_midas = L_midas_1, p_AR = p_ar_1, r = r_1, t_est_end = dates_m[t_est_end]),
    hyper_post = list(Lproxy = Lproxy_2, L_midas = L_midas_2, p_AR = p_ar_2, r = r_2,
                      t_recalib = if (do_recal) dates_m[t_recalib] else NA),
    countries = countries_eval,
    M1 = lapply(now_M1, function(x) unlist(x)),
    M2 = lapply(now_M2, function(x) unlist(x)),
    M3 = lapply(now_M3, function(x) unlist(x))
  )
}