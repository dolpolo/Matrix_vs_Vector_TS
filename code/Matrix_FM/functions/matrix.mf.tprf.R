# ==============================================================================
# PROXIES IN THE MF-3PRF (Matrix MF-TPRF) — SINGLE FUNCTION, NO NESTED FUNCTIONS
#
# INPUT:
#   X_lf   : array T x p1 x p2, already centered over time, COMPLETE (no NA)
#   y_q    : vector length T (GDP EA), can be non-centered
#   Lproxy : total proxies = 1 + number of residual proxies
#   r      : c(r1, r2)
#   p_AR   : fixed AR order in Step 3 (default 1)
#
# OUTPUT:
#   list(Z, R, C, F_hat, fit3)
#     Z     : T x Lproxy   (y_q, e1, ..., e_{Lproxy-1})
#     R, C  : last-iteration bilinear loadings
#     F_hat : last-iteration tensor factors
#     fit3  : last Step-3 lm object
# ==============================================================================
build_autoproxy_Tensor_3prf <- function(X_lf, y_q, Lproxy, r, p_AR = 1) {
  
  # =========================
  # STEP 0 — checks & setup
  # =========================
  T_obs <- dim(X_lf)[1]
  p1    <- dim(X_lf)[2]
  p2    <- dim(X_lf)[3]
  
  if (length(y_q) != T_obs) stop("Length(y_q) must match dim(X_lf)[1].")
  if (anyNA(X_lf)) stop("X_lf contains NA. This version expects complete X_lf.")
  if (length(r) != 2) stop("r must be c(r1, r2).")
  r1 <- r[1]; r2 <- r[2]
  if (r1 > p1 || r2 > p2) stop("r components must be <= (p1, p2).")
  
  if (length(Lproxy) != 1 || !is.finite(Lproxy) || Lproxy < 1) stop("Lproxy must be an integer >= 1.")
  Lproxy <- as.integer(Lproxy)
  
  if (length(p_AR) != 1 || !is.finite(p_AR) || p_AR < 1) stop("p_AR must be an integer >= 1.")
  p_AR <- as.integer(p_AR)
  if (T_obs <= p_AR) stop("Need T_obs > p_AR to run Step 3 with AR lags.")
  
  dn <- dimnames(X_lf)
  dates     <- if (!is.null(dn[[1]])) dn[[1]] else as.character(seq_len(T_obs))
  countries <- if (!is.null(dn[[2]])) dn[[2]] else paste0("c", seq_len(p1))
  vars      <- if (!is.null(dn[[3]])) dn[[3]] else paste0("v", seq_len(p2))
  
  col_names <- if (Lproxy == 1) "y_q" else c("y_q", paste0("e", 1:(Lproxy - 1)))
  Z <- matrix(NA_real_, nrow = T_obs, ncol = Lproxy,
              dimnames = list(dates, col_names))
  Z[, 1] <- y_q
  
  # Storage for last iteration
  R_last <- NULL; C_last <- NULL; F_last <- NULL; fit_last <- NULL
  
  # ======================================================
  # STEP 4 — Automatic-proxy loop (ell = 1..Lproxy-1)
  #   At iteration ell:
  #     use Z_curr = [y_q, e1, ..., e_{ell-1}]  (T x ell)
  #     run Step 1–2 to get (R,C,F)
  #     run Step 3 to get residual e_ell
  # ======================================================
  n_iter <- if (Lproxy == 1) 1 else (Lproxy - 1)
  
  for (ell in seq_len(n_iter)) {
    
    # Current proxy set
    Z_curr <- Z[, 1:ell, drop = FALSE]
    Z_targ <- sweep(Z_curr, 2, colMeans(Z_curr, na.rm = TRUE), FUN = "-")  # demean columns
    Lc <- ncol(Z_curr)
    
    # ======================================================
    # STEP 1 — Build M (single proxy) OR {G, Mlist} (multi proxy)
    # ======================================================
    if (Lc == 1) {
      # ---- single proxy: M = sum_t z_t X_t
      z <- as.numeric(Z_targ[, 1])
      M <- matrix(0, nrow = p1, ncol = p2)
      for (t in seq_len(T_obs)) {
        if (is.na(z[t])) next
        M <- M + z[t] * X_lf[t, , ]
      }
      
      # ======================================================
      # STEP 2 — Extract (R,C) and compute factors F_t
      #   single proxy: SVD(M) => R,C
      # ======================================================
      sv   <- svd(M)
      Rhat <- sv$u[, 1:r1, drop = FALSE]
      Chat <- sv$v[, 1:r2, drop = FALSE]
      
    } else {
      # ---- multi proxy: build G=Z'Z and moments M^(m)=sum_t z_{m,t} X_t
      G <- matrix(0, nrow = Lc, ncol = Lc)
      Mlist <- vector("list", Lc)
      for (m in seq_len(Lc)) Mlist[[m]] <- matrix(0, nrow = p1, ncol = p2)
      
      for (t in seq_len(T_obs)) {
        zt <- Z_targ[t, ]
        if (anyNA(zt)) next
        G <- G + tcrossprod(zt)
        Xt <- X_lf[t, , ]
        for (m in seq_len(Lc)) {
          Mlist[[m]] <- Mlist[[m]] + zt[m] * Xt
        }
      }
      
      if (qr(G)$rank < Lc) stop("G = Z'Z is singular (proxies collinear). Reduce Lproxy.")
      W <- solve(G)  # whitening
      
      # ======================================================
      # STEP 2 — Extract (R,C) via Srow/Scol and compute factors F_t
      #   Srow = sum_{m,ell} W_{ell,m} M^(m) (M^(ell))'
      #   Scol = sum_{m,ell} W_{ell,m} (M^(m))' M^(ell)
      # ======================================================
      Srow <- matrix(0, nrow = p1, ncol = p1)
      Scol <- matrix(0, nrow = p2, ncol = p2)
      
      for (m in seq_len(Lc)) {
        for (ell2 in seq_len(Lc)) {
          w_em <- W[ell2, m]
          if (w_em == 0) next
          Mm <- Mlist[[m]]
          Me <- Mlist[[ell2]]
          Srow <- Srow + w_em * (Mm %*% t(Me))
          Scol <- Scol + w_em * (t(Mm) %*% Me)
        }
      }
      
      er <- eigen(Srow, symmetric = TRUE)
      ec <- eigen(Scol, symmetric = TRUE)
      Rhat <- er$vectors[, 1:r1, drop = FALSE]
      Chat <- ec$vectors[, 1:r2, drop = FALSE]
    }
    
    rownames(Rhat) <- countries; colnames(Rhat) <- paste0("r", seq_len(r1))
    rownames(Chat) <- vars;      colnames(Chat) <- paste0("c", seq_len(r2))
    
    # Factors F_t = R' X_t C
    Fhat <- array(NA_real_, dim = c(T_obs, r1, r2),
                  dimnames = list(dates, paste0("f1_", seq_len(r1)), paste0("f2_", seq_len(r2))))
    for (t in seq_len(T_obs)) {
      Fhat[t, , ] <- t(Rhat) %*% X_lf[t, , ] %*% Chat
    }
    
    # Save last (even if Lproxy==1 we want these)
    R_last <- Rhat; C_last <- Chat; F_last <- Fhat
    
    # If Lproxy==1, we stop after first pass (no residual proxies needed)
    if (Lproxy == 1) break
    
    # ======================================================
    # STEP 3 — Forecasting regression with fixed AR(p_AR)
    #   y_t ~ [y_{t-1},...,y_{t-p}] + vec(F_{t-1})
    #   residual e_ell = y - y_hat
    # ======================================================
    F_vec <- matrix(NA_real_, nrow = T_obs, ncol = r1 * r2)
    for (t in seq_len(T_obs)) F_vec[t, ] <- as.vector(Fhat[t, , ])
    
    idx_t <- (p_AR + 1):T_obs
    y_dep <- y_q[idx_t]
    
    Y_lags <- sapply(1:p_AR, function(k) y_q[idx_t - k])
    colnames(Y_lags) <- paste0("y_lag", 1:p_AR)
    
    X_Flag <- F_vec[idx_t - 1, , drop = FALSE]
    
    df <- data.frame(y = y_dep, Y_lags, X_Flag)
    fit <- lm(y ~ ., data = df)
    fit_last <- fit
    
    yhat <- rep(NA_real_, T_obs)
    yhat[idx_t] <- as.numeric(predict(fit, newdata = df))
    
    e <- y_q - yhat
    e[1:p_AR] <- 0
    
    # add new proxy (residual) for next iteration
    Z[, ell + 1] <- e
  }
  
  list(Z = Z, R = R_last, C = C_last, F_hat = F_last, fit3 = fit_last)
}

# ==============================================================================
# NUMBER OF PROXIES IN THE MF-3PRF: 
#===============================================================================
# Richiede: library(plsdof)  # dove sta pls.ic()

select_L_autoproxy_3prf <- function(
    X_lf, y_q, Zmax,
    center_X = TRUE, scale_X = TRUE,
    center_y = TRUE, scale_y = TRUE,
    naive = FALSE,
    use.kernel = FALSE,
    compute.jacobian = FALSE,  # FALSE => Krylov
    verbose = TRUE,
    # gestione primo periodo (non prevedibile)
    residual_first = c("y1", "zero"),
    return_table = TRUE
) {
  residual_first <- match.arg(residual_first)
  stopifnot(is.matrix(X_lf), is.numeric(y_q), nrow(X_lf) == length(y_q))
  
  T_obs <- nrow(X_lf)
  
  # ---------- preprocessing coerente per DoF e per la tua yhat ----------
  Xp <- X_lf
  yp <- as.numeric(y_q)
  
  if (center_X || scale_X) Xp <- scale(Xp, center = center_X, scale = scale_X)
  if (center_y) yp <- yp - mean(yp)
  if (scale_y)  yp <- yp / sd(yp)
  
  # Zmax non può superare min(N, T-1)
  Zmax <- min(Zmax, ncol(Xp), nrow(Xp) - 1)
  if (Zmax < 1) stop("Zmax troppo piccolo o dati insufficienti.")
  
  # ---------- DoF(m) dal pacchetto (Krylov) ----------
  pls.object <- plsdof::pls.ic(
    X = Xp, y = yp, m = Zmax,
    criterion = "bic",
    naive = naive,
    use.kernel = use.kernel,
    compute.jacobian = compute.jacobian,
    verbose = verbose
  )
  
  DoF_raw <- as.numeric(pls.object$DoF)
  # robustezza: lunghezza può essere < Zmax in caso di crash/DoF negativi
  m_dof <- min(Zmax, length(DoF_raw))
  DoF_vec <- rep(NA_real_, Zmax)
  DoF_vec[1:m_dof] <- DoF_raw[1:m_dof]
  
  # ---------- loop sequenziale: costruisco proxy e calcolo RSS/BIC ----------
  # inizializzo Z con target-proxy
  Z <- matrix(yp, ncol = 1)
  colnames(Z) <- "y_proxy"
  
  T_eff <- T_obs - 1  # one-step ahead -> uso t=2..T
  
  RSS <- rep(NA_real_, Zmax)
  sigma2_hat <- rep(NA_real_, Zmax)
  BIC <- rep(NA_real_, Zmax)
  
  # helper: dato Z corrente (k colonne), calcola y_hat (one-step ahead) e residui
  compute_yhat_resid <- function(Z_current) {
    Z_tilde <- cbind(1, Z_current)
    
    A <- crossprod(Z_tilde)
    B <- crossprod(Z_tilde, Xp)
    B_full <- t(solve(A, B))
    Phi <- B_full[, -1, drop = FALSE]  # N x k
    
    Phi_cross <- crossprod(Phi)
    F_hat <- Xp %*% Phi %*% solve(Phi_cross)  # T x k
    
    y_t <- yp[-1]
    F_lag <- F_hat[-nrow(F_hat), , drop = FALSE]
    df_reg <- data.frame(y_t = y_t, F_lag)
    
    fit <- lm(y_t ~ ., data = df_reg)
    
    y_hat <- rep(NA_real_, T_obs)
    y_hat[-1] <- as.numeric(predict(fit, newdata = df_reg))
    
    e <- yp - y_hat
    if (residual_first == "y1")  e[1] <- yp[1]
    if (residual_first == "zero") e[1] <- 0
    
    list(y_hat = y_hat, e = e)
  }
  
  # m = 1 (solo target-proxy)
  tmp <- compute_yhat_resid(Z)
  e_use <- (yp - tmp$y_hat)[-1]
  RSS[1] <- sum(e_use^2)
  if (is.finite(DoF_vec[1]) && (T_eff - DoF_vec[1]) > 0) {
    sigma2_hat[1] <- RSS[1] / (T_eff - DoF_vec[1])
    BIC[1] <- (RSS[1] / T_eff) + log(T_eff) * sigma2_hat[1] * (DoF_vec[1] / T_eff)
  }
  
  # m = 2..Zmax: aggiungo proxy residuo e aggiorno
  if (Zmax >= 2) {
    for (m in 2:Zmax) {
      # aggiungo nuova proxy = residuo dall'ultima stima
      Z <- cbind(Z, tmp$e)
      colnames(Z)[ncol(Z)] <- paste0("e", m - 1)
      
      tmp <- compute_yhat_resid(Z)
      
      e_use <- (yp - tmp$y_hat)[-1]
      RSS[m] <- sum(e_use^2)
      
      dof_m <- DoF_vec[m]
      if (is.finite(dof_m) && (T_eff - dof_m) > 0) {
        sigma2_hat[m] <- RSS[m] / (T_eff - dof_m)
        BIC[m] <- (RSS[m] / T_eff) + log(T_eff) * sigma2_hat[m] * (dof_m / T_eff)
      }
    }
  }
  
  tab <- data.frame(
    m = 1:Zmax,
    DoF = DoF_vec,
    RSS = RSS,
    sigma2_hat = sigma2_hat,
    BIC = BIC
  )
  
  if (all(!is.finite(tab$BIC))) stop("Tutti i BIC sono NA/Inf: controlla DoF o Zmax.")
  L_opt <- tab$m[which.min(tab$BIC)]
  
  if (!return_table) return(L_opt)
  
  list(L_opt = L_opt, table = tab, pls_object = pls.object)
}



# ==============================================================================
# CHOOSE LAG FOR MF-UMIDAS USING MATRIX MF-TPRF FACTORS
#
# INPUT
#   X_lf : array (T_q x p1 x p2)  low-frequency panel, centered, complete
#   X_hf : array (T_m x p1 x p2)  high-frequency panel, centered, complete
#   y_q  : vector length T_q      quarterly target
#   r    : c(r1, r2)              bilinear ranks
#   Lproxy : number of proxies (y + residual proxies)
#   Lmax   : max UMIDAS factor lag (in quarters)
#   p_AR   : fixed AR order in y (quarterly)
#
# OUTPUT
#   list(results, lag_AIC, lag_BIC, Z, R, C)
# ==============================================================================
choose_UMIDAS_lag_tensor_MF <- function(X_lf, X_hf, y_q,
                                        r,
                                        Lmax   = 5,
                                        Lproxy = 1,
                                        p_AR   = 1) {
  
  # -------------------------
  # STEP 0 — checks
  # -------------------------
  if (anyNA(X_lf) || anyNA(X_hf)) stop("X_lf / X_hf contain NA. This version expects complete arrays.")
  if (length(dim(X_lf)) != 3 || length(dim(X_hf)) != 3) stop("X_lf and X_hf must be 3D arrays: T x p1 x p2.")
  if (length(r) != 2) stop("r must be c(r1, r2).")
  
  T_q <- dim(X_lf)[1]
  T_m <- dim(X_hf)[1]
  p1  <- dim(X_lf)[2]
  p2  <- dim(X_lf)[3]
  
  if (length(y_q) != T_q) stop("length(y_q) must match dim(X_lf)[1].")
  if (dim(X_hf)[2] != p1 || dim(X_hf)[3] != p2) stop("X_hf must have same (p1,p2) as X_lf.")
  if (Lmax < 1 || Lmax >= T_q) stop("Lmax must be >=1 and < T_q.")
  if (p_AR < 1) stop("p_AR must be >= 1 (fixed AR order).")
  if (T_m < 3 * T_q) stop("Need at least 3*T_q months in X_hf to cover all quarters.")
  
  r1 <- r[1]; r2 <- r[2]
  if (r1 > p1 || r2 > p2) stop("r components must be <= (p1, p2).")
  
  # -------------------------
  # STEP 1 — build Lproxy proxies Z (automatic-proxy MF-TPRF)
  #   Z: T_q x Lproxy  (y_q, e1, ..., e_{Lproxy-1})
  # -------------------------
  proxy_out <- build_autoproxy_Tensor_3prf(X_lf = X_lf, y_q = y_q, Lproxy = Lproxy, r = r, p_AR = p_AR)
  Z <- proxy_out$Z
  
  # Use ALL proxies for the final targeting step (common choice)
  Z_curr <- Z[, 1:Lproxy, drop = FALSE]
  
  # Optional but recommended: demean Z only for Step 2 targeting
  Z_targ <- sweep(Z_curr, 2, colMeans(Z_curr, na.rm = TRUE), FUN = "-")
  Lc <- ncol(Z_targ)
  
  # -------------------------
  # STEP 2 — MF-TPRF TARGETING: estimate R, C using X_lf and Z_targ
  #   Case Lc=1: SVD(M), M = sum_t z_t X_t
  #   Case Lc>1: whitening W=(Z'Z)^{-1}, build Srow/Scol, eigenvectors
  # -------------------------
  if (Lc == 1) {
    
    z <- as.numeric(Z_targ[, 1])
    M <- matrix(0, nrow = p1, ncol = p2)
    for (t in seq_len(T_q)) M <- M + z[t] * X_lf[t, , ]
    
    sv <- svd(M)
    Rhat <- sv$u[, 1:r1, drop = FALSE]
    Chat <- sv$v[, 1:r2, drop = FALSE]
    
  } else {
    
    G <- matrix(0, nrow = Lc, ncol = Lc)
    Mlist <- vector("list", Lc)
    for (m in seq_len(Lc)) Mlist[[m]] <- matrix(0, nrow = p1, ncol = p2)
    
    for (t in seq_len(T_q)) {
      zt <- Z_targ[t, ]
      if (anyNA(zt)) next
      G <- G + tcrossprod(zt)
      Xt <- X_lf[t, , ]
      for (m in seq_len(Lc)) Mlist[[m]] <- Mlist[[m]] + zt[m] * Xt
    }
    
    if (qr(G)$rank < Lc) stop("G = Z'Z is singular (proxies collinear). Reduce Lproxy.")
    W <- solve(G)
    
    Srow <- matrix(0, nrow = p1, ncol = p1)
    Scol <- matrix(0, nrow = p2, ncol = p2)
    for (m in seq_len(Lc)) {
      for (ell in seq_len(Lc)) {
        w_em <- W[ell, m]
        if (w_em == 0) next
        Mm <- Mlist[[m]]
        Me <- Mlist[[ell]]
        Srow <- Srow + w_em * (Mm %*% t(Me))
        Scol <- Scol + w_em * (t(Mm) %*% Me)
      }
    }
    
    er <- eigen(Srow, symmetric = TRUE)
    ec <- eigen(Scol, symmetric = TRUE)
    Rhat <- er$vectors[, 1:r1, drop = FALSE]
    Chat <- ec$vectors[, 1:r2, drop = FALSE]
  }
  
  # -------------------------
  # STEP 3 — Monthly factors from HF data: F_t = R' X_hf,t C
  #   F_hat_m: T_m x (r1*r2) after vectorization
  # -------------------------
  kF <- r1 * r2
  F_hat_m <- matrix(NA_real_, nrow = T_m, ncol = kF)
  for (t in seq_len(T_m)) {
    Ft <- t(Rhat) %*% X_hf[t, , ] %*% Chat   # r1 x r2
    F_hat_m[t, ] <- as.vector(Ft)
  }
  
  # -------------------------
  # STEP 3.1 — Map months to quarters (complete quarters only)
  #   F1,F2,F3: T_q x kF
  # -------------------------
  F1 <- F_hat_m[seq(1, 3 * T_q, by = 3), , drop = FALSE]
  F2 <- F_hat_m[seq(2, 3 * T_q, by = 3), , drop = FALSE]
  F3 <- F_hat_m[seq(3, 3 * T_q, by = 3), , drop = FALSE]
  
  # -------------------------
  # STEP 4 — Choose UMIDAS lag L by AIC/BIC
  #   y_tau ~ AR(p_AR) + [F1_tau,...,F3_tau, ..., lags up to L-1]
  # -------------------------
  results <- data.frame(L = integer(), AIC = numeric(), BIC = numeric())
  
  for (L in 1:Lmax) {
    
    start_tau <- max(L, p_AR) + 1
    y_dep <- y_q[start_tau:T_q]
    
    # AR block
    Y_lag <- NULL
    for (j in 1:p_AR) {
      Y_lag <- cbind(Y_lag, y_q[(start_tau - j):(T_q - j)])
    }
    colnames(Y_lag) <- paste0("y_lag", 1:p_AR)
    
    # UMIDAS factor block: for ell = 0..L-1 include (F1,F2,F3) at tau-ell
    Xreg_F <- NULL
    for (ell in 0:(L - 1)) {
      idx <- (start_tau - ell):(T_q - ell)
      Xreg_F <- cbind(Xreg_F,
                      F1[idx, , drop = FALSE],
                      F2[idx, , drop = FALSE],
                      F3[idx, , drop = FALSE])
    }
    
    Xreg <- cbind(Y_lag, Xreg_F)
    
    fit <- lm(y_dep ~ Xreg)
    results <- rbind(results,
                     data.frame(L = L, AIC = AIC(fit), BIC = BIC(fit)))
  }
  
  list(
    results = results,
    lag_AIC = results$L[which.min(results$AIC)],
    lag_BIC = results$L[which.min(results$BIC)],
    Z = Z,
    R = Rhat,
    C = Chat
  )
}


# ==============================================================================
# MF-3PRF: 
#===============================================================================
#' PROXY
#' LOADINGS 
#' FATTORI
#' TARGET PREDICTION
#' add the high friquency predictors
# ==============================================================================

# X_lf        = X_em_agg_tens
# X_hf        = X_em_tens
# Y_q_all     = y_q_all
# proxy_name = "EA"


Tensor_MF_TPRF <- function(X_lf, X_hf, Y_q_all,
                           proxy_name  = "EA",
                           Lproxy      = 1,
                           L_midas     = 1,
                           p_AR        = 1,
                           r           = c(1, 1)
                           ) {
  
  #--------------------------------------------------
  # STEP 0 — checks
  #--------------------------------------------------
  if (anyNA(X_lf) || anyNA(X_hf)) stop("X_lf / X_hf contain NA. This version expects complete arrays.")
  if (length(dim(X_lf)) != 3 || length(dim(X_hf)) != 3) stop("X_lf and X_hf must be 3D arrays: T x p1 x p2.")
  if (length(r) != 2) stop("r must be c(r1, r2).")
  
  y_proxy <- as.numeric(Y_q_all[, proxy_name])
  countries <- setdiff(colnames(Y_q_all), proxy_name)
  
  T_q <- dim(X_lf)[1]
  T_m <- dim(X_hf)[1]
  p1  <- dim(X_lf)[2]
  p2  <- dim(X_lf)[3]
  
  if (length(y_proxy) != T_q) stop("Proxy y length must match dim(X_lf)[1].")
  if (dim(X_hf)[2] != p1 || dim(X_hf)[3] != p2) stop("X_hf must have same (p1,p2) as X_lf.")
  if (L_midas < 1 || L_midas > T_q) stop("L_midas must be in 1..T_q.")
  if (p_AR < 1) stop("p_AR must be >= 1.")
  if (T_m < 3 * T_q) stop("Need at least 3*T_q months in X_hf to cover all quarters.")
  
  r1 <- r[1]; r2 <- r[2]
  if (r1 > p1 || r2 > p2) stop("r components must be <= (p1, p2).")
  K <- r1 * r2
  
  # dimnames helper
  dn_q <- dimnames(X_lf)[[1]]
  dn_m <- dimnames(X_hf)[[1]]
  if (is.null(dn_q)) dn_q <- as.character(seq_len(T_q))
  if (is.null(dn_m)) dn_m <- as.character(seq_len(T_m))
  
  #--------------------------------------------------
  # STEP 1 — Autoproxy construction at LF (Matrix MF-TPRF)
  #--------------------------------------------------
  auto_out <- build_autoproxy_Tensor_3prf(
    X_lf   = X_lf,
    y_q    = y_proxy,
    Lproxy = Lproxy,
    r      = r,
    p_AR   = p_AR
  )
  
  Z_q <- auto_out$Z   # T_q x Lproxy
  
  #--------------------------------------------------
  # STEP 2 — Recompute final R,C using ALL proxies (demeaned) for targeting
  #         (More stable / more "PLS-like" targeting step)
  #--------------------------------------------------
  Z_curr <- Z_q[, 1:Lproxy, drop = FALSE]
  Z_targ <- sweep(Z_curr, 2, colMeans(Z_curr, na.rm = TRUE), FUN = "-")
  Lc <- ncol(Z_targ)
  
  if (Lc == 1) {
    z <- as.numeric(Z_targ[, 1])
    M <- matrix(0, nrow = p1, ncol = p2)
    for (t in seq_len(T_q)) M <- M + z[t] * X_lf[t, , ]
    sv <- svd(M)
    Rhat <- sv$u[, 1:r1, drop = FALSE]
    Chat <- sv$v[, 1:r2, drop = FALSE]
  } else {
    G <- matrix(0, nrow = Lc, ncol = Lc)
    Mlist <- vector("list", Lc)
    for (m in seq_len(Lc)) Mlist[[m]] <- matrix(0, nrow = p1, ncol = p2)
    
    for (t in seq_len(T_q)) {
      zt <- Z_targ[t, ]
      if (anyNA(zt)) next
      G <- G + tcrossprod(zt)
      Xt <- X_lf[t, , ]
      for (m in seq_len(Lc)) Mlist[[m]] <- Mlist[[m]] + zt[m] * Xt
    }
    
    if (qr(G)$rank < Lc) stop("G = Z'Z is singular (proxies collinear). Reduce Lproxy.")
    W <- solve(G)
    
    Srow <- matrix(0, nrow = p1, ncol = p1)
    Scol <- matrix(0, nrow = p2, ncol = p2)
    for (m in seq_len(Lc)) {
      for (ell in seq_len(Lc)) {
        w_em <- W[ell, m]
        if (w_em == 0) next
        Mm <- Mlist[[m]]
        Me <- Mlist[[ell]]
        Srow <- Srow + w_em * (Mm %*% t(Me))
        Scol <- Scol + w_em * (t(Mm) %*% Me)
      }
    }
    
    er <- eigen(Srow, symmetric = TRUE)
    ec <- eigen(Scol, symmetric = TRUE)
    Rhat <- er$vectors[, 1:r1, drop = FALSE]
    Chat <- ec$vectors[, 1:r2, drop = FALSE]
  }
  
  #--------------------------------------------------
  # STEP 2bis — Monthly factors from HF tensor: F_t = R' X_t C
  #--------------------------------------------------
  F_hf <- array(NA_real_, dim = c(T_m, r1, r2),
                dimnames = list(dn_m, paste0("f1_", seq_len(r1)), paste0("f2_", seq_len(r2))))
  for (t in seq_len(T_m)) {
    F_hf[t, , ] <- t(Rhat) %*% X_hf[t, , ] %*% Chat
  }
  
  F_vec <- matrix(NA_real_, nrow = T_m, ncol = K)
  for (t in seq_len(T_m)) F_vec[t, ] <- as.vector(F_hf[t, , ])
  
  #--------------------------------------------------
  # STEP 3.1 — Quarterly blocks F1/F2/F3 (complete quarters + ragged edge)
  #--------------------------------------------------
  F1 <- F_vec[seq(1, 3 * T_q, by = 3), , drop = FALSE]
  F2 <- F_vec[seq(2, 3 * T_q, by = 3), , drop = FALSE]
  F3 <- F_vec[seq(3, 3 * T_q, by = 3), , drop = FALSE]
  
  rem <- T_m - 3 * T_q
  F_next1 <- if (rem >= 1) F_vec[3 * T_q + 1, , drop = FALSE] else NULL
  F_next2 <- if (rem >= 2) F_vec[3 * T_q + 2, , drop = FALSE] else NULL
  F_next3 <- if (rem >= 3) F_vec[3 * T_q + 3, , drop = FALSE] else NULL
  
  #--------------------------------------------------
  # STEP 3.2 — U-MIDAS per country + STEP 4 Nowcast mensile
  #--------------------------------------------------
  results_country <- list()
  
  for (cc in countries) {
    
    y_q <- as.numeric(Y_q_all[, cc])
    if (length(y_q) != T_q) stop(paste0("Country ", cc, ": y length mismatch."))
    
    start_tau <- max(L_midas, p_AR) + 1
    y_dep <- y_q[start_tau:T_q]
    
    # AR block
    Y_lag <- NULL
    for (j in 1:p_AR) {
      Y_lag <- cbind(Y_lag, y_q[(start_tau - j):(T_q - j)])
    }
    colnames(Y_lag) <- paste0("y_lag", 1:p_AR)
    
    # UMIDAS factor block (lags 0..L_midas-1)
    Xreg_F <- NULL
    for (ell_id in 1:L_midas) {
      ell <- ell_id - 1
      idx <- (start_tau - ell):(T_q - ell)
      Xreg_F <- cbind(Xreg_F,
                      F1[idx, , drop = FALSE],
                      F2[idx, , drop = FALSE],
                      F3[idx, , drop = FALSE])
    }
    
    Xreg <- cbind(Y_lag, Xreg_F)
    
    # OLS with intercept
    X_tilde <- cbind(1, Xreg)
    beta_hat <- solve(t(X_tilde) %*% X_tilde, t(X_tilde) %*% y_dep)
    
    beta0 <- as.numeric(beta_hat[1])
    rho_hat <- as.numeric(beta_hat[2:(1 + p_AR)])
    beta_vec <- as.numeric(beta_hat[(2 + p_AR):length(beta_hat)])
    
    beta_mat <- matrix(beta_vec, nrow = L_midas, ncol = 3 * K, byrow = TRUE)
    
    # ---- Nowcast monthly series ----
    y_nowcast <- rep(NA_real_, T_m)
    
    # 4.a complete quarters
    for (tau in start_tau:T_q) {
      
      month_idx <- ((tau - 1) * 3 + 1):(tau * 3)
      
      AR_part_tau <- sum(rho_hat * sapply(1:p_AR, function(j) y_q[tau - j]))
      
      for (m in 1:3) {
        contrib <- 0
        
        for (ell_id in 1:L_midas) {
          ell <- ell_id - 1
          lag_q <- tau - ell
          if (lag_q < 1) next
          
          for (mm in 1:3) {
            
            if (ell_id == 1 && mm > m) next  # ragged edge within quarter for ell=0
            
            F_qm <- switch(mm,
                           `1` = F1[lag_q, ],
                           `2` = F2[lag_q, ],
                           `3` = F3[lag_q, ])
            
            start_col <- (mm - 1) * K + 1
            end_col   <- mm * K
            beta_block <- beta_mat[ell_id, start_col:end_col]
            
            contrib <- contrib + sum(F_qm * beta_block)
          }
        }
        
        y_nowcast[month_idx[m]] <- beta0 + AR_part_tau + contrib
      }
    }
    
    # 4.b current quarter T_q+1 (if rem months)
    if (rem > 0) {
      
      tau_curr <- T_q + 1
      month_curr <- 3 * T_q + seq_len(rem)
      
      AR_part_curr <- sum(rho_hat * sapply(1:p_AR, function(j) y_q[tau_curr - j]))
      
      for (m in 1:rem) {
        contrib <- 0
        
        for (ell_id in 1:L_midas) {
          ell <- ell_id - 1
          lag_q <- tau_curr - ell
          
          if (ell_id >= 2 && (lag_q < 1 || lag_q > T_q)) next
          
          for (mm in 1:3) {
            
            if (ell_id == 1) {
              if (mm > m) next
              F_qm <- switch(mm, `1` = F_next1, `2` = F_next2, `3` = F_next3)
              if (is.null(F_qm)) next
            } else {
              F_qm <- switch(mm,
                             `1` = F1[lag_q, ],
                             `2` = F2[lag_q, ],
                             `3` = F3[lag_q, ])
            }
            
            start_col <- (mm - 1) * K + 1
            end_col   <- mm * K
            beta_block <- beta_mat[ell_id, start_col:end_col]
            
            contrib <- contrib + sum(F_qm * beta_block)
          }
        }
        
        y_nowcast[month_curr[m]] <- beta0 + AR_part_curr + contrib
      }
    }
    
    results_country[[cc]] <- list(
      beta0     = beta0,
      rho_hat   = rho_hat,
      beta_mat  = beta_mat,
      y_nowcast = y_nowcast,
      start_tau = start_tau
    )
  }
  
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  list(
    Z_q        = Z_q,
    R          = Rhat,
    C          = Chat,
    F_hf       = F_hf,
    F1         = F1,
    F2         = F2,
    F3         = F3,
    by_country = results_country
  )
}
