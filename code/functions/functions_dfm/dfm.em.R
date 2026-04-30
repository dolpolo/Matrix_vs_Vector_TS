# ==============================================================================
# DFM_EM (coerente con InitialCond monolitica)
# Stato z_t = [ f_t, f_{t-1},..., f_{t-pC_f+1} ;  e_{1,t},...,e_{1,t-pC_e+1}, ..., e_{N,t},...,e_{N,t-pC_e+1} ]
# - fattori: VAR(p_factor) aggiornato (solo prima block-row di A_f e Q_f)
# - idio: AR(q_idio) aggiornato PER OGNI SERIE i=1..N (solo prima riga del blocco, shift fisso)
# - loadings: Lambda aggiornato (solo la parte fattoriale in C; la parte idio in C resta fissa)
# - restrizioni trimestrali (MM / stock_flow) entrano SOLO nelle righe di misura delle trimestrali
# ==============================================================================
# X <- X_std
DFM_EM <- function(
    X, Init,
    max_iter = 100, tol = 1e-4,
    kappa_fixed_R = TRUE,
    clip_phi = 0.99
) {
  X <- as.matrix(X)
  
  # --------- recupero meta coerente con InitialCond ----------
  meta <- Init$meta
  r       <- meta$r
  p       <- meta$p_factor
  q       <- meta$q_idio
  NM      <- meta$NM
  NQ      <- meta$NQ
  N       <- meta$N
  restr   <- meta$restr
  L       <- meta$L
  pC_f    <- meta$pC_f
  pC_e    <- meta$pC_e
  
  # sanity
  stopifnot(ncol(X) == N)
  if (!isTRUE(all.equal(restr, Init$meta$restr))) stop("Init$meta$restr mismatch.")
  
  A  <- Init$A
  C  <- Init$C
  Qm <- Init$Q
  Rm <- Init$R
  Z0 <- Init$Z0
  V0 <- Init$V0
  
  ll_old <- compute_loglik(X, A, C, Qm, Rm, Z0, V0)
  
  loglik_seq <- numeric(max_iter + 1)
  crit_seq   <- numeric(max_iter)
  loglik_seq[1] <- ll_old
  
  for (iter in 1:max_iter) {
    cat("EM iter:", iter, "\n")
    
    Step <- EMstep(
      X = X,
      A = A, C = C, Q = Qm, R = Rm, Z0 = Z0, V0 = V0,
      r = r, p = p, q = q,
      NM = NM, NQ = NQ, N = N,
      restr = restr, L = L, pC_f = pC_f, pC_e = pC_e,
      agg = if (!is.null(Init$meta$agg_used)) Init$meta$agg_used else NULL,
      kappa_fixed_R = kappa_fixed_R,
      clip_phi = clip_phi
    )
    
    A  <- Step$A
    C  <- Step$C
    Qm <- Step$Q
    Rm <- Step$R
    Z0 <- Step$Z0
    V0 <- Step$V0
    
    ll_new <- compute_loglik(X, A, C, Qm, Rm, Z0, V0)
    loglik_seq[iter + 1] <- ll_new
    
    denom <- (abs(ll_new) + abs(ll_old)) / 2
    c_j <- if (denom == 0) 0 else (ll_new - ll_old) / denom
    crit_seq[iter] <- c_j
    
    cat("  loglik =", ll_new, "  c_j =", c_j, "\n")
    
    if (abs(c_j) < tol) {
      cat("Converged at iter", iter, "with c_j =", c_j, "\n")
      loglik_seq <- loglik_seq[1:(iter + 1)]
      crit_seq   <- crit_seq[1:iter]
      break
    }
    
    ll_old <- ll_new
    if (iter == max_iter) cat("Max iter reached without satisfying convergence.\n")
  }
  
  list(A = A, C = C, Q = Qm, R = Rm, Z0 = Z0, V0 = V0,
       loglik = loglik_seq, crit = crit_seq)
}


# ==============================================================================
# EMstep coerente con InitialCond monolitica
# ==============================================================================
EMstep <- function(
    X, A, C, Q, R, Z0, V0,
    r, p, q,
    NM, NQ, N,
    restr = c("MM", "stock_flow"),
    L, pC_f, pC_e,
    agg = NULL,
    kappa_fixed_R = TRUE,
    clip_phi = 0.99
) {
  restr <- match.arg(restr)
  X <- as.matrix(X)
  Tt <- nrow(X)
  
  # ------------- helper: pesi trimestrali --------------
  if (restr == "MM") {
    stopifnot(L == 5L)
    get_wQ <- function(i) c(1,2,3,2,1)
  } else {
    stopifnot(L == 3L)
    if (is.null(agg)) stop("agg must be provided (in Init$meta$agg_used) when restr='stock_flow'.")
    if (length(agg) == NQ) {
      aggN <- c(rep(NA_integer_, NM), as.integer(agg))
    } else {
      stopifnot(length(agg) == N)
      aggN <- as.integer(agg)
    }
    get_wQ <- function(i) {
      ai <- aggN[i]
      if (is.na(ai)) stop("agg for a quarterly variable is NA.")
      if (ai == 2L) return(rep(1, L))     # flow: sum
      if (ai == 1L) return(rep(1/L, L))   # stock: mean
      stop("agg must be coded 1 (stock) or 2 (flow).")
    }
  }
  
  # ------------- dimensioni stato coerenti --------------
  k_fac <- r * pC_f
  k_e   <- N * pC_e
  k     <- k_fac + k_e
  stopifnot(nrow(A) == k, ncol(A) == k)
  
  # ------------- (opzionale) finestra: puoi partire da 1;
  # qui NON taglio più da L perché la gestione dei missing e del quarterly
  # la fai via O_t (osservati). Il modello di misura è già built per usare i lags.
  Xf   <- X
  Tobs <- nrow(Xf)
  
  # osservati per kalman.na
  Olist <- lapply(1:Tobs, function(tt) which(!is.na(Xf[tt, ])))
  
  # ------------- E-step: Kalman filter + smoother (tuo) --------------
  kf <- kalman.na(
    par   = list(t = A, l = C, q = Q, h = R),
    y     = t(Xf),              # N x T
    k     = k,
    start = list(f = Z0, P = V0),
    O     = Olist
  )
  
  Z_sm  <- kf$fs     # k x T
  P_sm  <- kf$Ps     # k x k x T
  P_lag <- kf$Cs     # k x k x T  (Cov(z_{t}, z_{t-1}) in your convention)
  
  # moments
  Ezz <- function(tt) {
    P_sm[,,tt] + Z_sm[,tt,drop=FALSE] %*% t(Z_sm[,tt,drop=FALSE])
  }
  Ezz1 <- function(tt) {
    # tt>=2: E[z_t z_{t-1}'] = Cov(z_t,z_{t-1}) + E[z_t]E[z_{t-1}]'
    P_lag_mat <- P_lag[,,tt-1]
    P_lag_mat + Z_sm[,tt,drop=FALSE] %*% t(Z_sm[,tt-1,drop=FALSE])
  }
  
  Tden <- Tobs - 1
  
  # =============================================================================
  # M-step (1): update VAR(p) sui fattori
  # =============================================================================
  idx_f_cur <- 1:r
  idx_reg   <- 1:(r*p)   # in z_{t-1}: [f_{t-1},...,f_{t-p}]
  
  S00 <- matrix(0, r*p, r*p)
  S10 <- matrix(0, r,   r*p)
  S11 <- matrix(0, r,   r)
  
  for (tt in 2:Tobs) {
    Et_t1  <- Ezz1(tt)
    Etm1m1 <- Ezz(tt - 1)
    Ett    <- Ezz(tt)
    
    S10 <- S10 + Et_t1[idx_f_cur, idx_reg, drop=FALSE]
    S00 <- S00 + Etm1m1[idx_reg, idx_reg, drop=FALSE]
    S11 <- S11 + Ett[idx_f_cur, idx_f_cur, drop=FALSE]
  }
  
  A_new <- A
  Q_new <- Q
  
  A_f <- S10 %*% solve(S00)  # r x (r*p)
  A_new[1:r, 1:(r*p)] <- A_f
  if (k_fac > r*p) A_new[1:r, (r*p + 1):k_fac] <- 0  # lags oltre p fino a pC_f
  
  Q_f <- (S11 - A_f %*% t(S10)) / Tden
  Q_f <- 0.5 * (Q_f + t(Q_f))
  Q_new[1:r, 1:r] <- Q_f
  
  # =============================================================================
  # M-step (2): update AR(q) idio PER OGNI SERIE i=1..N
  # Stato idio serie i in z_t: offset = k_fac + (i-1)*pC_e + 1
  # regressori in z_{t-1}: [e_{t-1},...,e_{t-q}] = prime q posizioni dello stesso blocco
  # =============================================================================
  if (q > 0) {
    for (i in 1:N) {
      off_i_t   <- k_fac + (i - 1) * pC_e + 1          # e_{i,t} in z_t
      off_i_tm1 <- k_fac + (i - 1) * pC_e + 1          # e_{i,t-1} in z_{t-1} è nello stesso offset
      
      idx_e_t   <- off_i_t
      idx_reg_e <- off_i_tm1:(off_i_tm1 + q - 1)       # [e_{t-1},...,e_{t-q}]
      
      S00e <- matrix(0, q, q)
      S10e <- matrix(0, 1, q)
      S11e <- 0
      
      for (tt in 2:Tobs) {
        Et_t1  <- Ezz1(tt)
        Etm1m1 <- Ezz(tt - 1)
        Ett    <- Ezz(tt)
        
        S10e <- S10e + Et_t1[idx_e_t, idx_reg_e, drop=FALSE]
        S00e <- S00e + Etm1m1[idx_reg_e, idx_reg_e, drop=FALSE]
        S11e <- S11e + Ett[idx_e_t, idx_e_t]
      }
      
      phi_i <- as.numeric(S10e %*% solve(S00e))   # length q
      # stabilizzazione "leggera"
      phi_i <- pmax(pmin(phi_i, clip_phi), -clip_phi)
      if (q > 1) phi_i <- 0.95 * phi_i
      
      sig2_i <- as.numeric((S11e - (S10e %*% matrix(phi_i, ncol = 1))) / Tden)
      sig2_i <- max(sig2_i, 1e-10)
      
      # aggiorna SOLO prima riga del blocco A_e per la serie i
      row_idx <- k_fac + (i - 1) * pC_e + 1
      col_idx <- k_fac + (i - 1) * pC_e + 1:q
      A_new[row_idx, col_idx] <- phi_i
      
      # shift sotto resta fisso (già in A)
      # Q: solo var innovazione sul primo stato idio
      Q_new[row_idx, row_idx] <- sig2_i
    }
  } else {
    # q=0: white noise idio -> A first row zero, Q updated by var
    for (i in 1:N) {
      idx <- k_fac + (i - 1) * pC_e + 1
      # stima var innovazione da E[e_t^2]
      s2 <- 0
      for (tt in 1:Tobs) s2 <- s2 + Ezz(tt)[idx, idx]
      s2 <- max(s2 / Tobs, 1e-10)
      Q_new[idx, idx] <- s2
      # A_new[idx, idx] deve essere 0 (già tale se costruito bene)
      A_new[idx, (k_fac + (i - 1) * pC_e + 1): (k_fac + (i - 1) * pC_e + pC_e)] <- 0
    }
  }
  
  # =============================================================================
  # M-step (3): update loadings Lambda (solo parte fattoriale in C)
  # - mensili: y = lambda_i' f_t + e_it (+ noise)
  # - trimestrali: y = sum_{l=0}^{L-1} w_l (lambda_i' f_{t-l} + e_{i,t-l})
  # idio in C resta fisso (1 sul current e per mensili; w su stack e per trimestrali)
  # =============================================================================
  C_new <- C
  
  # --- mensili ---
  for (i in 1:NM) {
    obs_t <- which(!is.na(Xf[, i]))
    if (length(obs_t) < 5) next
    
    idx_ei <- k_fac + (i - 1) * pC_e + 1
    
    Sff <- matrix(0, r, r)
    syf <- matrix(0, 1, r)
    
    for (tt in obs_t) {
      Ett <- Ezz(tt)
      
      Sff <- Sff + Ett[idx_f_cur, idx_f_cur, drop=FALSE]
      # y*f - E[e*f]
      Ef_row  <- matrix(Z_sm[idx_f_cur, tt], nrow = 1)                     # 1 x r
      Eef_row <- matrix(Ett[idx_ei, idx_f_cur], nrow = 1)                  # 1 x r
      syf <- syf + (as.numeric(Xf[tt, i]) * Ef_row - Eef_row)
    }
    
    lambda_i <- syf %*% solve(Sff)  # 1 x r
    C_new[i, 1:r] <- as.numeric(lambda_i)
    if (k_fac > r) C_new[i, (r + 1):k_fac] <- 0
  }
  
  # --- trimestrali ---
  idx_facLag <- 1:(r * L)
  
  for (j in 1:NQ) {
    i <- NM + j
    obs_t <- which(!is.na(Xf[, i]))
    if (length(obs_t) < 5) next
    
    w_i <- get_wQ(i)                 # length L
    stopifnot(length(w_i) == L)
    
    # g_t = sum w_l f_{t-l}  => implemento via S_mat
    S_mat <- kronecker(matrix(w_i, nrow = 1), diag(r))  # r x (rL)
    
    # idio stack col: e_t,...,e_{t-L+1} nello stato della serie i
    idx_e_stack <- k_fac + (i - 1) * pC_e + (1:L)
    c_idio      <- matrix(w_i, nrow = 1)               # 1 x L
    
    Sgg <- matrix(0, r, r)
    syg <- matrix(0, 1, r)
    
    for (tt in obs_t) {
      Ett <- Ezz(tt)
      
      Ezf_zf <- Ett[idx_facLag, idx_facLag, drop=FALSE]          # (rL)x(rL)
      Egg    <- S_mat %*% Ezf_zf %*% t(S_mat)                    # r x r
      Sgg    <- Sgg + Egg
      
      # E[idio * g]
      Ezq_zf <- Ett[idx_e_stack, idx_facLag, drop=FALSE]         # L x (rL)
      Ezq_g  <- Ezq_zf %*% t(S_mat)                              # L x r
      Eid_g  <- c_idio %*% Ezq_g                                 # 1 x r
      
      Eg_col <- S_mat %*% Z_sm[idx_facLag, tt, drop=FALSE]       # r x 1
      syg <- syg + (as.numeric(Xf[tt, i]) * t(Eg_col) - Eid_g)
    }
    
    lambda_i <- syg %*% solve(Sgg)  # 1 x r
    
    # riscrivi la parte fattoriale della riga trimestrale con i pesi w
    for (ell in 0:(L - 1)) {
      cols <- (ell * r + 1):((ell + 1) * r)
      C_new[i, cols] <- w_i[ell + 1] * as.numeric(lambda_i)
    }
    if (k_fac > r * L) C_new[i, (r * L + 1):k_fac] <- 0
  }
  
  # =============================================================================
  # R update (tipicamente fissata se idio è nello stato)
  # =============================================================================
  R_new <- R
  if (!kappa_fixed_R) {
    # se vuoi stimarla davvero: va derivata da residui attesi; per ora lasciamo invariata
    R_new <- R
  }
  
  list(
    A  = A_new,
    C  = C_new,
    Q  = Q_new,
    R  = R_new,
    Z0 = Z_sm[, 1],
    V0 = P_sm[,, 1]
  )
}


# ==============================================================================
# Log-likelihood via kalman.na
# ==============================================================================
compute_loglik <- function(X, A, C, Q, R, Z0, V0) {
  X <- as.matrix(X)
  Tobs <- nrow(X)
  k <- nrow(A)
  
  Olist <- lapply(1:Tobs, function(tt) which(!is.na(X[tt, ])))
  
  kf <- kalman.na(
    par   = list(t = A, l = C, q = Q, h = R),
    y     = t(X),
    k     = k,
    start = list(f = Z0, P = V0),
    O     = Olist
  )
  
  as.numeric(kf$loglik)
}
