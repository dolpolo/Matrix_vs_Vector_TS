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
# EXPANDING NOWCAST
# ==============================================================================

# X_tens_full <- X
# Y_q_all <- y_q_all
# proxy_name = "EA"
# dates_m <- tensor$dates
# dates_q <- dates_EA
# tensor_info <- tensor_info_pred
# tt <- 304 # [M1]
# tt <- 305 # [M2]
# tt <- 306 # [M3]


pseudo_realtime_Tensor_MF_TPRF <- function(
    X_tens_full,
    Y_q_all,
    proxy_name,
    params,
    dates_m,
    dates_q,
    Unb,
    agg_m,
    agg_q,
    tensor_info
) {
  
  # -------------------------------------------------------------------
  # 0) META: dimensioni e controlli
  # -------------------------------------------------------------------
  stopifnot(length(dim(X_tens_full)) == 3)
  T_m <- dim(X_tens_full)[1]
  P   <- dim(X_tens_full)[2]
  K   <- dim(X_tens_full)[3]
  
  N_m <- tensor_info$n_M
  N_q <- tensor_info$n_Q
  
  if ((N_m + N_q) != K) stop("N_m + N_q != dim(X_tens_full)[3]. GDP deve essere escluso.")
  if (!all(dim(Unb) == c(P, K))) stop("Unb deve essere [P x (N_m+N_q)].")
  if (!all(dim(agg_m) == c(P, N_m))) stop("agg_m deve essere [P x N_m].")
  if (!all(dim(agg_q) == c(P, N_q))) stop("agg_q deve essere [P x N_q].")
  
  if (nrow(Y_q_all) != length(dates_q)) stop("Y_q_all e dates_q non allineati.")
  if (is.null(colnames(Y_q_all))) stop("Y_q_all deve avere colnames.")
  if (!(proxy_name %in% colnames(Y_q_all))) stop("proxy_name non in colnames(Y_q_all).")
  
  countries_eval <- setdiff(colnames(Y_q_all), proxy_name)
  if (length(countries_eval) < 1) stop("Nessun paese oltre alla proxy in Y_q_all.")
  
  # attributi (se ti servono)
  attr(X_tens_full, "N_m")       <- N_m
  attr(X_tens_full, "N_q")       <- N_q
  attr(X_tens_full, "countries") <- tensor_info$countries
  attr(X_tens_full, "var_names") <- tensor_info$vars
  
  agg_m_global <- as.vector(t(agg_m))
  agg_q_global <- as.vector(t(agg_q))
  
  # -------------------------------------------------------------------
  # 1) Evaluation window
  # -------------------------------------------------------------------
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  if (length(t_start) == 0 || length(t_end) == 0) stop("start_eval/end_eval non in dates_m.")
  if (t_start > t_end) stop("start_eval > end_eval.")
  
  # -------------------------------------------------------------------
  # 2) Estimation window end (per hyper)
  # -------------------------------------------------------------------
  t_est_end <- max(which(dates_m < params$start_eval))
  if (!is.finite(t_est_end) || t_est_end < 24) stop("Estimation sample troppo corto.")
  
  cat("\n>>> HYPER-PARAM SELECTION using data up to", as.character(dates_m[t_est_end]), "\n")
  
  # ==========================================================
  # 2.A) COSTRUISCI DATASET ESTIMATION a t_est_end
  # ==========================================================
  X_cut_est <- unbalancedness_tensor(X_full = X_tens_full, Unb = Unb, current_t = t_est_end)
  if (all(is.na(X_cut_est))) stop("Estimation cut all-NA.")
  W_cut_est <- ifelse(is.na(X_cut_est), 0L, 1L)
  
  std_est   <- standardize_mat_with_na(X_cut_est)
  X_std_est <- std_est$X_scaled
  
  flat_est <- flatten_tensor_to_matrix(
    X_tens = X_std_est,
    W_tens = W_cut_est,
    N_m    = N_m,
    N_q    = N_q,
    agg_q  = agg_q
  )
  X_obs_est <- flat_est$X_mat
  
  A_est <- A_list(
    X_na  = X_obs_est,
    N_q   = flat_est$N_q_global,
    agg_q = flat_est$agg_q_global
  )
  
  init_est <- init_CL_ER(X_std = X_std_est, W = W_cut_est, kmax = params$kmax, do_plot = FALSE)
  
  # >>> r scelto nell'estimation e fissato
  r_fix <- init_est$r                  # c(r1,r2)
  r_vec_fix <- prod(r_fix)
  
  X_init_est <- flatten_tensor_init(X_tens = init_est$X_init, N_m = N_m, N_q = N_q)
  
  EM_est <- EM_algorithm(
    X_init   = X_init_est,
    X_obs    = X_obs_est,
    A_list   = A_est,
    r        = r_vec_fix,
    max_iter = 100,
    tol      = 1e-4
  )
  X_em_est <- EM_est$X_completed
  for (nm in c("N_m","N_q","countries","var_names")) attr(X_em_est, nm) <- attr(X_obs_est, nm)
  
  # aggrego (LF)
  X_m_est <- X_em_est[, 1:flat_est$N_m_global, drop = FALSE]
  X_q_est <- X_em_est[, (flat_est$N_m_global + 1):flat_est$N_global, drop = FALSE]
  
  X_mq_est <- agg_mq(X_m_est, agg_m_global)
  X_qq_est <- agg_qq(X_q_est, agg_q_global)
  X_em_agg_est <- cbind(X_mq_est, X_qq_est)
  for (nm in c("N_m","N_q","countries","var_names")) attr(X_em_agg_est, nm) <- attr(X_em_est, nm)
  
  # unflatten
  X_em_est_tens     <- unflatten_matrix_to_tensor(X_em_est)       # [T_m_est x P x K]
  X_em_agg_est_tens <- unflatten_matrix_to_tensor(X_em_agg_est)   # [T_q_est x P x K]
  T_q_est <- dim(X_em_agg_est_tens)[1]
  T_m_est <- dim(X_em_est_tens)[1]
  
  # trimestri proxy disponibili fino a t_est_end
  idx_pub_est <- which(dates_q < dates_m[t_est_end])
  if (length(idx_pub_est) < 2) stop("Troppi pochi trimestri proxy nell'estimation sample.")
  T_q_use_est <- min(max(idx_pub_est), T_q_est, nrow(Y_q_all))
  
  y_proxy_est <- as.numeric(Y_q_all[1:T_q_use_est, proxy_name])
  
  X_lf_est <- X_em_agg_est_tens[1:T_q_use_est, , , drop = FALSE]
  X_hf_est <- X_em_est_tens[1:T_m_est,       , , drop = FALSE]
  
  # >>> NIENTE center_Y: X è già centrato/standardizzato in std_est
  
  # ==========================================================
  # 2.B) SELEZIONE hyper su estimation
  # ==========================================================
  X_lf_est_mat <- matrix(NA_real_, nrow = dim(X_lf_est)[1], ncol = dim(X_lf_est)[2] * dim(X_lf_est)[3])
  for (t in seq_len(nrow(X_lf_est_mat))) X_lf_est_mat[t, ] <- as.vector(X_lf_est[t, , ])
  
  pls_obj    <- select_L_autoproxy_3prf(X_lf = X_lf_est_mat, y_q = y_proxy_est, Zmax = params$Zmax)
  Lproxy_fix <- pls_obj$L_opt
  cat("# [Hyper] Lproxy_fix:", Lproxy_fix, "\n")
  
  lag_sel <- choose_UMIDAS_lag_tensor_MF(
    X_lf   = X_lf_est,
    X_hf   = X_hf_est,
    y_q    = y_proxy_est,
    r      = r_fix,
    Lmax   = params$Lmax,
    Lproxy = Lproxy_fix,
    p_AR   = params$p_AR
  )
  L_midas_fix <- lag_sel$lag_BIC
  cat("# [Hyper] L_midas_fix (BIC):", L_midas_fix, "\n")
  
  # -------------------------------------------------------------------
  # 3) Evaluation rolling con hyper FISSI
  # -------------------------------------------------------------------
  now_M1 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M2 <- setNames(vector("list", length(countries_eval)), countries_eval)
  now_M3 <- setNames(vector("list", length(countries_eval)), countries_eval)
  
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    cat("\n>>> REAL-TIME at", as.character(date_t),
        "| r_fix =", paste0("(", r_fix[1], ",", r_fix[2], ")"),
        "| Lproxy_fix =", Lproxy_fix,
        "| L_midas_fix =", L_midas_fix, "\n")
    
    X_cut <- unbalancedness_tensor(X_full = X_tens_full, Unb = Unb, current_t = tt)
    if (all(is.na(X_cut))) next
    W_cut <- ifelse(is.na(X_cut), 0L, 1L)
    
    std_out <- standardize_mat_with_na(X_cut)
    X_std <- std_out$X_scaled
    
    flat <- flatten_tensor_to_matrix(X_tens = X_std, W_tens = W_cut, N_m = N_m, N_q = N_q, agg_q = agg_q)
    X_obs <- flat$X_mat
    
    A <- A_list(X_na = X_obs, N_q = flat$N_q_global, agg_q = flat$agg_q_global)
    
    init <- init_CL_ER(X_std = X_std, W = W_cut, kmax = params$kmax, do_plot = FALSE)
    X_init <- flatten_tensor_init(X_tens = init$X_init, N_m = N_m, N_q = N_q)
    
    EM_out <- EM_algorithm(X_init = X_init, X_obs = X_obs, A_list = A,
                           r = r_vec_fix, max_iter = 100, tol = 1e-4)
    
    X_em <- EM_out$X_completed
    T_m_cur <- nrow(X_em)
    for (nm in c("N_m","N_q","countries","var_names")) attr(X_em, nm) <- attr(X_obs, nm)
    
    X_m_em <- X_em[, 1:flat$N_m_global, drop = FALSE]
    X_q_em <- X_em[, (flat$N_m_global + 1):flat$N_global, drop = FALSE]
    
    X_mq <- agg_mq(X_m_em, agg_m_global)
    X_qq <- agg_qq(X_q_em, agg_q_global)
    X_em_agg <- cbind(X_mq, X_qq)
    for (nm in c("N_m","N_q","countries","var_names")) attr(X_em_agg, nm) <- attr(X_em, nm)
    
    X_em_tens     <- unflatten_matrix_to_tensor(X_em)
    X_em_agg_tens <- unflatten_matrix_to_tensor(X_em_agg)
    T_q_em <- dim(X_em_agg_tens)[1]
    
    idx_pub <- which(dates_q < date_t)
    if (length(idx_pub) < 2) next
    T_q_use <- min(max(idx_pub), T_q_em, nrow(Y_q_all))
    if (T_q_use < 2) next
    
    Y_cut <- Y_q_all[1:T_q_use, , drop = FALSE]
    X_lf_cut <- X_em_agg_tens[1:T_q_use, , , drop = FALSE]
    X_hf_cut <- X_em_tens[1:T_m_cur, , , drop = FALSE]
    
    # >>> NIENTE center_Y anche qui
    
    out_rt <- Tensor_MF_TPRF(
      X_lf       = X_lf_cut,
      X_hf       = X_hf_cut,
      Y_q_all    = Y_cut,
      proxy_name = proxy_name,
      Lproxy     = Lproxy_fix,
      L_midas    = L_midas_fix,
      p_AR       = params$p_AR,
      r          = r_fix
    )
    
    m_tr <- compute_m_tr(date_t, dates_q)
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
    r_fix       = r_fix,
    Lproxy_fix  = Lproxy_fix,
    L_midas_fix = L_midas_fix,
    countries   = countries_eval,
    M1 = lapply(now_M1, function(x) unlist(x)),
    M2 = lapply(now_M2, function(x) unlist(x)),
    M3 = lapply(now_M3, function(x) unlist(x))
  )
}