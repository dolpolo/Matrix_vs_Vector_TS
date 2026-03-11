# ==============================================================================
# UNBALANCEDNESS
# ==============================================================================

# X_full <- X
# current_t <- 303

unbalancedness <- function(X_full, dates, Freq, Unb, current_t) {
  
  stopifnot(nrow(X_full) >= current_t)
  X_cut <- X_full[1:current_t, , drop = FALSE]
  
  for (j in seq_len(ncol(X_cut))) {
    
    delay_j <- Unb[j]          # ritardo di pubblicazione (in mesi)
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
  
  # ultimo trimestre rilasciato PRIMA del mese corrente
  last_q_date <- max(dates_q[dates_q < date_t])
  
  # indice del trimestre
  last_q_idx <- which(dates_q == last_q_date)
  
  # mesi del trimestre successivo (target quarter)
  M1 <- dates_q[last_q_idx] %m+% months(1)
  M2 <- dates_q[last_q_idx] %m+% months(2)
  M3 <- dates_q[last_q_idx] %m+% months(3)
  
  if (date_t == M1) return(1)
  if (date_t == M2) return(2)
  if (date_t == M3) return(3)
  
  return(NA)
}


# ==============================================================================
# prepara dataset “XP-imputed + aggregated” fino a una data tt
# ==============================================================================

build_XP_agg_up_to <- function(X_full, y_q, dates, dates_q, Freq, Unb,
                               current_t, agg_m, agg_q,
                               Kmax = 10) {
  
  # 1) ragged edge
  X_cut <- unbalancedness(
    X_full    = X_full,
    dates     = dates,
    Freq      = Freq,
    Unb       = Unb,
    current_t = current_t
  )
  
  # 2) standardize
  out_std <- standardize_with_na(X_cut)
  X_std   <- out_std$X_std
  
  # 3) XP imputation + r selection
  imp_xp <- init_XP_ER(X_std)   # se hai Kmax interno, passa Kmax = Kmax
  X_xp   <- imp_xp$X_init
  r_est  <- imp_xp$r
  
  # 4) split M/Q and aggregate
  X_m_xp <- X_xp[, Freq == "M", drop = FALSE]
  X_q_xp <- X_xp[, Freq == "Q", drop = FALSE]
  
  X_mq_xp <- agg_mq(X_m_xp, agg_m)
  X_qq_xp <- agg_qq(X_q_xp, agg_q)
  
  X_xp_agg <- cbind(X_mq_xp, X_qq_xp)
  
  # 5) GDP quarters available at current_t
  idx_pub <- which(dates_q < dates[current_t])
  if (length(idx_pub) < 2) return(NULL)
  
  T_q_current <- tail(idx_pub, 1)
  T_q_current <- min(T_q_current, nrow(X_xp_agg))
  
  list(
    X_hf = X_xp,                       # monthly completed
    X_lf = X_xp_agg,                   # quarterly aggregated predictors
    y_q  = as.numeric(y_q)[1:T_q_current],
    T_q  = T_q_current,
    r    = r_est
  )
}

# ==============================================================================
# EXPANDING NOWCAST
# ==============================================================================

# X_full <- X
# dates <- dates_m
# tt <- 304 # [M1]
# tt <- 305 # [M2]
# tt <- 306 # [M3]

pseudo_realtime_MF_TPRF_XP <- function(
    X_full, y_q, params,
    dates, dates_q,
    Freq, Unb,
    agg_m, agg_q
) {
  
  # --------------------------------------------
  # 0) Evaluation window
  # --------------------------------------------
  t_start <- which(dates == params$start_eval)
  t_end   <- which(dates == params$end_eval)
  if (length(t_start) == 0 || length(t_end) == 0)
    stop("start_eval o end_eval non trovati in 'dates'.")
  
  # --------------------------------------------
  # 1) Hyper selection #1 (pre-eval): up to start_eval-1
  # --------------------------------------------
  t_est_end <- max(which(dates < params$start_eval))
  if (length(t_est_end) == 0 || t_est_end < 24)
    stop("Estimation sample troppo corto o non definito.")
  
  cat("\n>>> HYPER #1 selection using data up to",
      as.character(dates[t_est_end]), "\n")
  
  obj_est1 <- build_XP_agg_up_to(
    X_full   = X_full,
    y_q      = y_q,
    dates    = dates,
    dates_q  = dates_q,
    Freq     = Freq,
    Unb      = Unb,
    current_t = t_est_end,
    agg_m    = agg_m,
    agg_q    = agg_q,
    Kmax     = params$Kmax
  )
  if (is.null(obj_est1)) stop("Non abbastanza GDP pubblicato per hyper #1.")
  
  # Lproxy
  pls1 <- select_L_autoproxy_3prf(obj_est1$X_lf[1:obj_est1$T_q, , drop=FALSE],
                                  obj_est1$y_q,
                                  Zmax = params$Zmax)
  Lproxy_1 <- pls1$L_opt
  
  # (p_AR, L) via BIC
  lag1 <- choose_UMIDAS_lag(
    X_lf        = obj_est1$X_lf[1:obj_est1$T_q, , drop=FALSE],
    X_hf        = obj_est1$X_hf,
    y_q         = obj_est1$y_q,
    Lmax        = params$Lmax,
    Lproxy      = Lproxy_1,
    p_AR_max    = params$p_AR_max,
    Robust_F    = params$Robust_F,
    alpha       = params$alpha,
    robust_type = params$robust_type,
    nw_lag      = params$nw_lag
  )
  L_midas_1 <- lag1$best_BIC$L
  p_ar_1    <- lag1$best_BIC$p_AR
  r_1       <- obj_est1$r
  
  cat("# Hyper #1:",
      "Lproxy =", Lproxy_1,
      "| L_midas =", L_midas_1,
      "| p_AR =", p_ar_1,
      "| r =", r_1, "\n")
  
  # --------------------------------------------
  # 2) Hyper selection #2 (post-COVID recalibration)
  #    la facciamo una volta sola: alla prima data >= covid_end
  # --------------------------------------------
  t_recalib <- which(dates >= params$covid_end)[1]
  do_recalib <- !is.na(t_recalib) && t_recalib <= t_end
  
  Lproxy_2 <- Lproxy_1
  L_midas_2 <- L_midas_1
  p_ar_2 <- p_ar_1
  r_2 <- r_1
  
  if (do_recalib) {
    cat("\n>>> HYPER #2 selection using data up to",
        as.character(dates[t_recalib]), "\n")
    
    obj_est2 <- build_XP_agg_up_to(
      X_full    = X_full,
      y_q       = y_q,
      dates     = dates,
      dates_q   = dates_q,
      Freq      = Freq,
      Unb       = Unb,
      current_t = t_recalib,
      agg_m     = agg_m,
      agg_q     = agg_q,
      Kmax      = params$Kmax
    )
    
    if (!is.null(obj_est2)) {
      pls2 <- select_L_autoproxy_3prf(obj_est2$X_lf[1:obj_est2$T_q, , drop=FALSE],
                                      obj_est2$y_q,
                                      Zmax = params$Zmax)
      Lproxy_2 <- pls2$L_opt
      
      lag2 <- choose_UMIDAS_lag(
        X_lf        = obj_est2$X_lf[1:obj_est2$T_q, , drop=FALSE],
        X_hf        = obj_est2$X_hf,
        y_q         = obj_est2$y_q,
        Lmax        = params$Lmax,
        Lproxy      = Lproxy_2,
        p_AR_max    = params$p_AR_max,
        Robust_F    = params$Robust_F,
        alpha       = params$alpha,
        robust_type = params$robust_type,
        nw_lag      = params$nw_lag
      )
      L_midas_2 <- lag2$best_BIC$L
      p_ar_2    <- lag2$best_BIC$p_AR
      r_2       <- obj_est2$r
      
      cat("# Hyper #2:",
          "Lproxy =", Lproxy_2,
          "| L_midas =", L_midas_2,
          "| p_AR =", p_ar_2,
          "| r =", r_2, "\n")
    } else {
      cat("(!) Hyper #2 skipped: not enough GDP published at recalib date.\n")
    }
  }
  
  # --------------------------------------------
  # 3) Real-time loop (expanding) with fixed hypers
  # --------------------------------------------
  now_M1 <- list()
  now_M2 <- list()
  now_M3 <- list()
  
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates[tt]
    
    # scegli regime hyper
    use_post <- do_recalib && tt >= t_recalib
    Lproxy_use <- if (use_post) Lproxy_2 else Lproxy_1
    L_midas_use <- if (use_post) L_midas_2 else L_midas_1
    p_ar_use <- if (use_post) p_ar_2 else p_ar_1
    
    cat("\n>>> REAL-TIME at", as.character(date_t),
        "| regime =", if (use_post) "post-COVID" else "pre-COVID",
        "| Lproxy =", Lproxy_use,
        "| L_midas =", L_midas_use,
        "| p_AR =", p_ar_use, "\n")
    
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
      Kmax      = params$Kmax
    )
    if (is.null(obj_rt) || length(obj_rt$y_q) < 2) next
    
    # taglia ai trimestri disponibili (coerente con GDP pubblicato)
    X_lf_cut <- obj_rt$X_lf[1:obj_rt$T_q, , drop = FALSE]
    y_q_cut  <- obj_rt$y_q
    X_hf_cut <- obj_rt$X_hf
    
    # MF-TPRF con hyper fissati
    out <- MF_TPRF(
      X_lf        = X_lf_cut,
      X_hf        = X_hf_cut,
      y_q         = y_q_cut,
      Lproxy      = Lproxy_use,
      L_midas     = L_midas_use,
      p_AR        = p_ar_use,
      Robust_F    = params$Robust_F,
      alpha       = params$alpha,
      robust_type = params$robust_type,
      nw_lag      = params$nw_lag
    )
    
    y_rt_last <- tail(out$y_nowcast, 1)
    
    m_tr <- compute_m_tr(date_t, dates_q)  # tua funzione: 1/2/3
    if (is.na(m_tr)) next
    
    key <- as.character(date_t)
    if (m_tr == 1) now_M1[[key]] <- y_rt_last
    if (m_tr == 2) now_M2[[key]] <- y_rt_last
    if (m_tr == 3) now_M3[[key]] <- y_rt_last
  }
  
  list(
    hyper_pre  = list(Lproxy = Lproxy_1, L_midas = L_midas_1, p_AR = p_ar_1, r = r_1),
    hyper_post = list(Lproxy = Lproxy_2, L_midas = L_midas_2, p_AR = p_ar_2, r = r_2,
                      t_recalib = if (do_recalib) dates[t_recalib] else NA),
    M1 = now_M1, M2 = now_M2, M3 = now_M3
  )
}