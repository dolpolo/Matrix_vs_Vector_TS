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
# GRID SEARCH (p_AR, L_midas) for Matrix MF-TPRF UMIDAS
#   Common effective sample across all (p_AR, L).
#
# INPUT
#   X_lf : array [T_q x p1 x p2] quarterly predictors (complete, centered)
#   X_hf : array [T_m x p1 x p2] monthly predictors (complete, centered)
#   y_q  : vector length T_q      quarterly target (country or EA)
#   r    : c(r1,r2)
#
#   Lproxy : fixed number of proxies (already chosen)
#   Lmax   : max MIDAS lag length (quarters)
#   p_AR_max: max AR order (allow 0)
#
#   use_autoproxy: if TRUE build Z via matrix_mf_tprf_autoproxy() using y_proxy
#   y_proxy      : vector length T_q used to build proxies (e.g. EA GDP). If NULL -> y_q.
#   Z_q          : optional precomputed proxy matrix/basis (T_q x Lproxy). Used if use_autoproxy=FALSE.
#
# OUTPUT
#   list(results, best_BIC, best_AIC, start_tau_global, Z_q, R, C)
# ==============================================================================

choose_UMIDAS_grid_tensor_MF <- function(
    X_lf, X_hf, y_q,
    r = c(1,1),
    Lproxy = 1,
    Lmax = 5,
    p_AR_max = 4,
    use_autoproxy = TRUE,
    y_proxy = NULL,
    Z_q = NULL,
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8
) {
  
  # ---- checks
  stopifnot(length(dim(X_lf)) == 3, length(dim(X_hf)) == 3)
  if (anyNA(X_lf) || anyNA(X_hf)) stop("X_lf / X_hf contain NA; expects complete arrays.")
  
  T_q <- dim(X_lf)[1]
  T_m <- dim(X_hf)[1]
  p1  <- dim(X_lf)[2]
  p2  <- dim(X_lf)[3]
  
  y_q <- as.numeric(y_q)
  if (length(y_q) != T_q) stop("length(y_q) must match dim(X_lf)[1].")
  if (dim(X_hf)[2] != p1 || dim(X_hf)[3] != p2) stop("X_hf must match (p1,p2) of X_lf.")
  if (T_m < 3 * T_q) stop("Need at least 3*T_q months in X_hf to cover all quarters.")
  
  if (length(r) != 2) stop("r must be c(r1,r2).")
  r1 <- as.integer(r[1]); r2 <- as.integer(r[2])
  if (r1 < 1 || r2 < 1 || r1 > p1 || r2 > p2) stop("Invalid r.")
  K <- r1 * r2
  
  if (!is.finite(Lproxy) || Lproxy < 1) stop("Lproxy must be integer >= 1.")
  Lproxy <- as.integer(Lproxy)
  
  if (!is.finite(Lmax) || Lmax < 1 || Lmax >= T_q) stop("Lmax must be in 1..(T_q-1).")
  Lmax <- as.integer(Lmax)
  
  if (!is.finite(p_AR_max) || p_AR_max < 0) stop("p_AR_max must be integer >= 0.")
  p_AR_max <- as.integer(p_AR_max)
  
  # --------------------------------------------------------
  # STEP A: proxies Z_q (fixed Lproxy)
  # --------------------------------------------------------
  if (use_autoproxy) {
    
    if (is.null(y_proxy)) y_proxy <- y_q
    y_proxy <- as.numeric(y_proxy)
    if (length(y_proxy) != T_q) stop("y_proxy must have length T_q.")
    
    auto <- matrix_mf_tprf_autoproxy(
      X_lf = X_lf, X_hf = X_hf, y_q = y_proxy,
      Lproxy = Lproxy, L_midas = 1,         # <-- temporary; only needed for residual loop,
      p_AR = 0,                              #     but if your proxy-builder needs L_midas/p_AR
      r = c(r1, r2),
      standardize_y = standardize_proxy,
      orthonormalize_each_iter = orthonormalize_each_iter,
      orthonormalize_final = orthonormalize_final_Z,
      ils_maxit = ils_maxit, ils_tol = ils_tol
    )
    
    Z_q_use <- auto$Z
    
  } else {
    
    if (is.null(Z_q)) stop("If use_autoproxy=FALSE you must supply Z_q.")
    Z_q_use <- as.matrix(Z_q)
    if (nrow(Z_q_use) != T_q || ncol(Z_q_use) != Lproxy)
      stop("Z_q must be T_q x Lproxy.")
  }
  
  # --------------------------------------------------------
  # STEP B: targeting once (R,C) using ALL proxies
  # --------------------------------------------------------
  pass1 <- mf_tprf_pass1_ils(
    X_lf = X_lf,
    Z_q  = Z_q_use,
    r = c(r1, r2),
    orthonormalize_Z = TRUE,
    center_Z = TRUE,
    maxit = ils_maxit,
    tol = ils_tol
  )
  
  # --------------------------------------------------------
  # STEP C: monthly factors once
  # --------------------------------------------------------
  pass2 <- mf_tprf_pass2_factors(
    X_hf = X_hf,
    R = pass1$R,
    C = pass1$C,
    T_q = T_q,
    scale_fac = 1 / (p1 * p2)
  )
  F1 <- pass2$F1; F2 <- pass2$F2; F3 <- pass2$F3
  
  # --------------------------------------------------------
  # STEP D: GRID SEARCH on (p_AR, L) with COMMON effective sample
  # global start: needs AR lags up to p_AR_max and MIDAS lag up to (Lmax-1)
  # --------------------------------------------------------
  start_tau_global <- 1L + max(p_AR_max, Lmax - 1L)
  if (start_tau_global > T_q) stop("Insufficient sample: start_tau_global > T_q. Reduce Lmax or p_AR_max.")
  n_eff <- T_q - start_tau_global + 1L
  
  results <- data.frame(
    p_AR = integer(),
    L    = integer(),
    n    = integer(),
    k    = integer(),
    sigma2 = numeric(),
    AIC  = numeric(),
    BIC  = numeric()
  )
  
  y_dep <- y_q[start_tau_global:T_q]
  
  for (p_AR in 0:p_AR_max) {
    for (L in 1:Lmax) {
      
      # ---- AR block aligned on common range
      X_AR <- NULL
      if (p_AR > 0) {
        X_AR <- sapply(1:p_AR, function(j) y_q[(start_tau_global - j):(T_q - j)])
        X_AR <- as.matrix(X_AR)
        colnames(X_AR) <- paste0("y_lag", 1:p_AR)
      }
      
      # ---- MIDAS block: ell=0..L-1, each adds (F1,F2,F3) at tau-ell
      X_MID <- NULL
      for (ell in 0:(L - 1)) {
        idx <- (start_tau_global - ell):(T_q - ell)  # length n_eff
        X_MID <- cbind(
          X_MID,
          F1[idx, , drop = FALSE],
          F2[idx, , drop = FALSE],
          F3[idx, , drop = FALSE]
        )
      }
      
      Xreg <- cbind(X_AR, X_MID)
      
      # OLS with intercept
      fit <- lm(y_dep ~ ., data = data.frame(y = y_dep, Xreg))
      
      # MLE-style residual variance (RSS / n)
      resid <- residuals(fit)
      sigma2_hat <- mean(resid^2)
      
      # number of parameters (intercept + AR + MIDAS blocks)
      k_par <- 1L + p_AR + 3L * K * L
      
      AIC_val <- log(sigma2_hat) + (2 * k_par) / n_eff
      BIC_val <- log(sigma2_hat) + (k_par * log(n_eff)) / n_eff
      
      results <- rbind(
        results,
        data.frame(
          p_AR = p_AR, L = L,
          n = n_eff, k = k_par,
          sigma2 = sigma2_hat,
          AIC = AIC_val, BIC = BIC_val
        )
      )
    }
  }
  
  best_AIC <- results[which.min(results$AIC), ]
  best_BIC <- results[which.min(results$BIC), ]
  
  list(
    results = results,
    start_tau_global = start_tau_global,
    n_eff = n_eff,
    best_AIC = best_AIC,
    best_BIC = best_BIC,
    Z_q = Z_q_use,
    R = pass1$R,
    C = pass1$C
  )
}


# ==============================================================================
# PASS 1 (Matrix MF-TPRF): ILS targeting given quarterly tensor X_lf and proxies Z_q
#   Inputs:
#     X_lf : array [Tq x p1 x p2]  (quarterly predictors)
#     Z_q  : matrix [Tq x L]       (proxies; will be centered+QR if requested)
#     r    : c(r1, r2)
#   Output:
#     list(R, C, Z_used, Mlist, obj_path)
# ==============================================================================

mf_tprf_pass1_ils <- function(
    X_lf, Z_q, r = c(1, 1),
    orthonormalize_Z = TRUE,
    center_Z = TRUE,
    maxit = 100,
    tol = 1e-8
) {
  
  # ---- checks
  stopifnot(length(dim(X_lf)) == 3)
  if (anyNA(X_lf)) stop("X_lf contains NA; Pass 1 expects complete X_lf.")
  Tq <- dim(X_lf)[1]
  p1 <- dim(X_lf)[2]
  p2 <- dim(X_lf)[3]
  
  Z_q <- as.matrix(Z_q)
  if (nrow(Z_q) != Tq) stop("nrow(Z_q) must match dim(X_lf)[1].")
  
  if (length(r) != 2) stop("r must be c(r1,r2).")
  r1 <- as.integer(r[1]); r2 <- as.integer(r[2])
  if (r1 < 1 || r2 < 1 || r1 > p1 || r2 > p2) {
    stop("Invalid r: must satisfy 1<=r1<=p1, 1<=r2<=p2.")
  }
  
  if (center_Z) {
    Zc <- scale(Z_q, center = TRUE, scale = FALSE)
  } else {
    Zc <- Z_q
  }
  
  if (orthonormalize_Z) {
    qrZ <- qr(Zc)
    if (qrZ$rank < ncol(Zc)) stop("Z_q columns are linearly dependent (after centering).")
    Z <- qr.Q(qrZ)
  } else {
    Z <- Zc
  }
  
  L <- ncol(Z)
  
  # ---- build moments M_l = (1/Tq) sum_t z_{t,l} X_t
  Mlist <- vector("list", L)
  for (l in 1:L) Mlist[[l]] <- matrix(0, nrow = p1, ncol = p2)
  
  for (t in 1:Tq) {
    Xt <- X_lf[t, , ]
    zt <- Z[t, ]
    for (l in 1:L) {
      Mlist[[l]] <- Mlist[[l]] + zt[l] * Xt
    }
  }
  for (l in 1:L) Mlist[[l]] <- Mlist[[l]] / Tq
  
  # ---- initialize C from eig( sum_l M_l' M_l )
  S_C0 <- matrix(0, nrow = p2, ncol = p2)
  for (l in 1:L) {
    Ml <- Mlist[[l]]
    S_C0 <- S_C0 + t(Ml) %*% Ml
  }
  ec0 <- eigen(S_C0, symmetric = TRUE)
  C_hat <- sqrt(p2) * ec0$vectors[, 1:r2, drop = FALSE]
  
  # ---- ILS loop
  obj_path <- numeric(0)
  obj_old <- -Inf
  
  for (it in 1:maxit) {
    
    # update R given C
    CCt <- C_hat %*% t(C_hat)
    S_R <- matrix(0, nrow = p1, ncol = p1)
    for (l in 1:L) {
      Ml <- Mlist[[l]]
      S_R <- S_R + Ml %*% CCt %*% t(Ml)
    }
    er <- eigen(S_R, symmetric = TRUE)
    R_hat <- sqrt(p1) * er$vectors[, 1:r1, drop = FALSE]
    
    # update C given R
    RRt <- R_hat %*% t(R_hat)
    S_C <- matrix(0, nrow = p2, ncol = p2)
    for (l in 1:L) {
      Ml <- Mlist[[l]]
      S_C <- S_C + t(Ml) %*% RRt %*% Ml
    }
    ec <- eigen(S_C, symmetric = TRUE)
    C_hat <- sqrt(p2) * ec$vectors[, 1:r2, drop = FALSE]
    
    # objective: proportional to reduced criterion under normalization
    obj <- 0
    for (l in 1:L) {
      A <- t(R_hat) %*% Mlist[[l]] %*% C_hat
      obj <- obj + sum(A * A)
    }
    obj_path <- c(obj_path, obj)
    
    if (abs(obj - obj_old) < tol * max(1, abs(obj_old))) break
    obj_old <- obj
  }
  
  # names
  dn <- dimnames(X_lf)
  countries <- if (!is.null(dn[[2]])) dn[[2]] else paste0("c", seq_len(p1))
  vars      <- if (!is.null(dn[[3]])) dn[[3]] else paste0("v", seq_len(p2))
  rownames(R_hat) <- countries; colnames(R_hat) <- paste0("r", 1:r1)
  rownames(C_hat) <- vars;      colnames(C_hat) <- paste0("c", 1:r2)
  
  list(
    R = R_hat,
    C = C_hat,
    Z_used = Z,
    Mlist = Mlist,
    obj_path = obj_path
  )
}

# ==============================================================================
# PASS 2 (Matrix MF-TPRF): monthly factors from HF tensor given R,C
#   Inputs:
#     X_hf : array [Tm x p1 x p2] (monthly predictors)
#     R    : [p1 x r1]
#     C    : [p2 x r2]
#     scale_fac : default 1/(p1*p2) as in main
#   Output:
#     list(F_hf, F_vec, F1, F2, F3, rem, F_next1, F_next2, F_next3)
# ==============================================================================

mf_tprf_pass2_factors <- function(
    X_hf, R, C,
    T_q = NULL,
    scale_fac = NULL
) {
  
  stopifnot(length(dim(X_hf)) == 3)
  if (anyNA(X_hf)) stop("X_hf contains NA; Pass 2 expects complete X_hf.")
  
  Tm <- dim(X_hf)[1]
  p1 <- dim(X_hf)[2]
  p2 <- dim(X_hf)[3]
  
  R <- as.matrix(R); C <- as.matrix(C)
  r1 <- ncol(R); r2 <- ncol(C)
  
  if (nrow(R) != p1) stop("nrow(R) must match dim(X_hf)[2].")
  if (nrow(C) != p2) stop("nrow(C) must match dim(X_hf)[3].")
  
  if (is.null(scale_fac)) scale_fac <- 1 / (p1 * p2)
  
  # factors tensor
  dnm <- dimnames(X_hf)[[1]]
  if (is.null(dnm)) dnm <- as.character(seq_len(Tm))
  
  F_hf <- array(NA_real_, dim = c(Tm, r1, r2),
                dimnames = list(dnm, paste0("f1_", 1:r1), paste0("f2_", 1:r2)))
  
  for (t in 1:Tm) {
    F_hf[t, , ] <- scale_fac * (t(R) %*% X_hf[t, , ] %*% C)
  }
  
  # vectorized factors
  K <- r1 * r2
  F_vec <- matrix(NA_real_, nrow = Tm, ncol = K)
  for (t in 1:Tm) F_vec[t, ] <- as.vector(F_hf[t, , ])
  
  # within-quarter blocks
  if (is.null(T_q)) {
    T_q <- floor(Tm / 3)  # best guess
  }
  if (Tm < 3 * T_q) stop("Pass 2: need Tm >= 3*T_q to form full-quarter blocks.")
  
  F1 <- F_vec[seq(1, 3 * T_q, by = 3), , drop = FALSE]
  F2 <- F_vec[seq(2, 3 * T_q, by = 3), , drop = FALSE]
  F3 <- F_vec[seq(3, 3 * T_q, by = 3), , drop = FALSE]
  
  rem <- Tm - 3 * T_q
  F_next1 <- if (rem >= 1) F_vec[3 * T_q + 1, , drop = FALSE] else NULL
  F_next2 <- if (rem >= 2) F_vec[3 * T_q + 2, , drop = FALSE] else NULL
  F_next3 <- if (rem >= 3) F_vec[3 * T_q + 3, , drop = FALSE] else NULL
  
  list(
    F_hf = F_hf,
    F_vec = F_vec,
    F1 = F1, F2 = F2, F3 = F3,
    rem = rem,
    F_next1 = F_next1, F_next2 = F_next2, F_next3 = F_next3
  )
}

# ==============================================================================
# PASS 3 (Matrix MF-TPRF): quarterly regression fit (in-sample fitted yhat_q)
#   Model:
#     y_tau = beta0 + sum_{j=1..p_AR} rho_j y_{tau-j}
#             + sum_{ell=0..L_midas-1} [ F1_{tau-ell}, F2_{tau-ell}, F3_{tau-ell} ] * beta_{ell}
#
#   Inputs:
#     y_q     : vector length T_q (target quarterly)
#     F1,F2,F3: matrices [T_q x K] (within-quarter blocks from monthly factors)
#     p_AR    : integer >= 0 (allow 0 => no AR terms)
#     L_midas : integer >= 1
#
#   Output:
#     list(fit, yhat_q, start_tau, beta0, rho_hat, beta_mat, df_used)
# ==============================================================================

mf_tprf_pass3_fit_quarterly <- function(
    y_q, F1, F2, F3,
    p_AR = 1,
    L_midas = 1
) {
  
  y_q <- as.numeric(y_q)
  T_q <- length(y_q)
  
  F1 <- as.matrix(F1); F2 <- as.matrix(F2); F3 <- as.matrix(F3)
  if (nrow(F1) != T_q || nrow(F2) != T_q || nrow(F3) != T_q)
    stop("F1/F2/F3 must have nrow = length(y_q).")
  
  if (!is.finite(p_AR) || p_AR < 0) stop("p_AR must be integer >= 0.")
  if (!is.finite(L_midas) || L_midas < 1) stop("L_midas must be integer >= 1.")
  p_AR <- as.integer(p_AR)
  L_midas <- as.integer(L_midas)
  
  K <- ncol(F1)
  if (ncol(F2) != K || ncol(F3) != K) stop("F1/F2/F3 must have same number of columns (K).")
  
  # common effective sample: need tau-ell >= 1 for ell=L_midas-1 and tau-j >= 1 for j=p_AR
  start_tau <- max(p_AR, L_midas - 1) + 1L
  if (start_tau > T_q) stop("Not enough quarters given p_AR and L_midas.")
  
  idx_tau <- start_tau:T_q
  y_dep <- y_q[idx_tau]
  
  # --- AR regressors (optional)
  X_AR <- NULL
  if (p_AR > 0) {
    X_AR <- sapply(1:p_AR, function(j) y_q[idx_tau - j])
    X_AR <- as.matrix(X_AR)
    colnames(X_AR) <- paste0("y_lag", 1:p_AR)
  }
  
  # --- MIDAS regressors (lags ell = 0..L_midas-1)
  X_MID <- NULL
  for (ell_id in 1:L_midas) {
    ell <- ell_id - 1L
    lag_idx <- idx_tau - ell
    
    block <- cbind(
      F1[lag_idx, , drop = FALSE],
      F2[lag_idx, , drop = FALSE],
      F3[lag_idx, , drop = FALSE]
    )
    X_MID <- cbind(X_MID, block)
  }
  
  Xreg <- cbind(X_AR, X_MID)
  df_used <- data.frame(y = y_dep, Xreg)
  
  # OLS with intercept
  fit <- lm(y ~ ., data = df_used)
  
  # in-sample fitted aligned to length T_q
  yhat_q <- rep(NA_real_, T_q)
  yhat_q[idx_tau] <- as.numeric(predict(fit, newdata = df_used))
  
  # extract parameters into your preferred format
  beta_hat <- coef(fit)
  beta0 <- as.numeric(beta_hat[1])
  
  if (p_AR > 0) {
    rho_hat <- as.numeric(beta_hat[2:(1 + p_AR)])
  } else {
    rho_hat <- numeric(0)
  }
  
  beta_vec <- as.numeric(beta_hat[(2 + p_AR):length(beta_hat)])
  if (length(beta_vec) != L_midas * (3 * K)) {
    stop("Coefficient length mismatch: expected L_midas*(3*K). Check regressor construction.")
  }
  beta_mat <- matrix(beta_vec, nrow = L_midas, ncol = 3 * K, byrow = TRUE)
  
  list(
    fit = fit,
    yhat_q = yhat_q,
    start_tau = start_tau,
    beta0 = beta0,
    rho_hat = rho_hat,
    beta_mat = beta_mat,
    df_used = df_used
  )
}

# ==============================================================================
# PASS 4 (Matrix MF-TPRF): monthly nowcast given estimated coefficients and factors
#
#   For each quarter tau, compute monthly nowcasts for m=1,2,3 using:
#     y_{tau|m} = beta0 + AR_part(tau) + sum_{ell=0..L_midas-1} sum_{mm<=m if ell=0}
#                  < beta_{ell,mm}, F_{mm, tau-ell} >
#
#   Inputs:
#     beta0, rho_hat, beta_mat : from mf_tprf_pass3_fit_quarterly()
#     y_q      : quarterly target (length T_q) needed for AR part
#     F1,F2,F3 : [T_q x K] blocks for complete quarters
#     F_next1/2/3 : [1 x K] blocks for current (incomplete) quarter if rem>0
#     p_AR, L_midas
#     T_q, T_m  : sizes (Tm can be > 3*Tq to include rem months)
#
#   Output:
#     y_nowcast : vector length T_m (monthly nowcasts)
# ==============================================================================

mf_tprf_nowcast_monthly <- function(
    beta0, rho_hat, beta_mat,
    y_q,
    F1, F2, F3,
    F_next1 = NULL, F_next2 = NULL, F_next3 = NULL,
    p_AR = 1,
    L_midas = 1,
    T_q = NULL,
    T_m = NULL
) {
  
  y_q <- as.numeric(y_q)
  if (is.null(T_q)) T_q <- length(y_q)
  if (length(y_q) != T_q) stop("y_q length mismatch with T_q.")
  
  F1 <- as.matrix(F1); F2 <- as.matrix(F2); F3 <- as.matrix(F3)
  K <- ncol(F1)
  if (nrow(F1) != T_q || nrow(F2) != T_q || nrow(F3) != T_q) stop("F1/F2/F3 must be T_q x K.")
  if (ncol(F2) != K || ncol(F3) != K) stop("F1/F2/F3 must have same K.")
  
  if (!is.finite(p_AR) || p_AR < 0) stop("p_AR must be integer >= 0.")
  if (!is.finite(L_midas) || L_midas < 1) stop("L_midas must be integer >= 1.")
  p_AR <- as.integer(p_AR)
  L_midas <- as.integer(L_midas)
  
  beta_mat <- as.matrix(beta_mat)
  if (nrow(beta_mat) != L_midas || ncol(beta_mat) != 3 * K)
    stop("beta_mat must be L_midas x (3*K).")
  
  # T_m default: full quarters plus possible rem months if you provided next blocks
  if (is.null(T_m)) {
    # if next blocks exist, assume up to 3 extra months, otherwise exact 3*T_q
    extra <- 0L
    if (!is.null(F_next1)) extra <- max(extra, 1L)
    if (!is.null(F_next2)) extra <- max(extra, 2L)
    if (!is.null(F_next3)) extra <- max(extra, 3L)
    T_m <- 3L * T_q + extra
  }
  T_m <- as.integer(T_m)
  if (T_m < 3L * T_q) stop("T_m must be >= 3*T_q.")
  
  rem <- T_m - 3L * T_q
  if (rem > 0L) {
    # basic consistency checks for next blocks
    if (rem >= 1L && is.null(F_next1)) stop("rem>=1 but F_next1 is NULL.")
    if (rem >= 2L && is.null(F_next2)) stop("rem>=2 but F_next2 is NULL.")
    if (rem >= 3L && is.null(F_next3)) stop("rem>=3 but F_next3 is NULL.")
  }
  
  # nowcast vector
  y_nowcast <- rep(NA_real_, T_m)
  
  # common effective quarter index for AR+MIDAS (must match Pass 3)
  start_tau <- max(p_AR, L_midas - 1) + 1L
  
  # helper: AR contribution at quarter tau
  AR_part <- function(tau) {
    if (p_AR == 0) return(0)
    sum(rho_hat * sapply(1:p_AR, function(j) y_q[tau - j]))
  }
  
  # ---- complete quarters: tau = start_tau..T_q
  for (tau in start_tau:T_q) {
    
    month_idx <- ((tau - 1L) * 3L + 1L):(tau * 3L)
    ar_tau <- AR_part(tau)
    
    for (m in 1:3) {
      
      contrib <- 0
      
      for (ell_id in 1:L_midas) {
        ell <- ell_id - 1L
        lag_q <- tau - ell
        if (lag_q < 1L) next
        
        for (mm in 1:3) {
          
          # ragged edge within-quarter only for ell=0
          if (ell_id == 1L && mm > m) next
          
          F_qm <- switch(mm,
                         `1` = F1[lag_q, ],
                         `2` = F2[lag_q, ],
                         `3` = F3[lag_q, ]
          )
          
          start_col <- (mm - 1L) * K + 1L
          end_col   <- mm * K
          beta_block <- beta_mat[ell_id, start_col:end_col]
          
          contrib <- contrib + sum(F_qm * beta_block)
        }
      }
      
      y_nowcast[month_idx[m]] <- beta0 + ar_tau + contrib
    }
  }
  
  # ---- current quarter (T_q + 1) if rem months are available
  if (rem > 0L) {
    
    tau_curr <- T_q + 1L
    month_curr <- 3L * T_q + seq_len(rem)
    
    # AR part: needs y_q up to T_q; for tau_curr it uses last observed quarters
    # only valid if tau_curr - j <= T_q, which holds
    ar_curr <- AR_part(tau_curr)
    
    for (m in 1:rem) {
      
      contrib <- 0
      
      for (ell_id in 1:L_midas) {
        ell <- ell_id - 1L
        lag_q <- tau_curr - ell
        
        for (mm in 1:3) {
          
          if (ell_id == 1L) {
            # ell=0: only months up to m are available
            if (mm > m) next
            
            F_qm <- switch(mm,
                           `1` = F_next1,
                           `2` = F_next2,
                           `3` = F_next3
            )
            if (is.null(F_qm)) next
            F_qm <- as.numeric(F_qm)
            
          } else {
            # ell>=1: use completed quarters only
            if (lag_q < 1L || lag_q > T_q) next
            
            F_qm <- switch(mm,
                           `1` = F1[lag_q, ],
                           `2` = F2[lag_q, ],
                           `3` = F3[lag_q, ]
            )
          }
          
          start_col <- (mm - 1L) * K + 1L
          end_col   <- mm * K
          beta_block <- beta_mat[ell_id, start_col:end_col]
          
          contrib <- contrib + sum(F_qm * beta_block)
        }
      }
      
      y_nowcast[month_curr[m]] <- beta0 + ar_curr + contrib
    }
  }
  
  y_nowcast
}

# ==============================================================================
# AUTO-PROXIES for Matrix MF-TPRF (NO nested functions)
# ==============================================================================

matrix_mf_tprf_autoproxy <- function(
    X_lf, X_hf, y_q,
    Lproxy = 1,
    L_midas = 1,
    p_AR = 1,              # allow 0
    r = c(1, 1),
    standardize_y = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8
) {
  
  # ---- checks
  stopifnot(length(dim(X_lf)) == 3, length(dim(X_hf)) == 3)
  if (anyNA(X_lf) || anyNA(X_hf)) stop("X_lf / X_hf contain NA. This version expects complete arrays.")
  
  Tq <- dim(X_lf)[1]
  Tm <- dim(X_hf)[1]
  p1 <- dim(X_lf)[2]
  p2 <- dim(X_lf)[3]
  
  if (length(y_q) != Tq) stop("length(y_q) must match dim(X_lf)[1].")
  if (dim(X_hf)[2] != p1 || dim(X_hf)[3] != p2) stop("X_hf must have same (p1,p2) as X_lf.")
  if (Tm < 3 * Tq) stop("Need at least 3*Tq months in X_hf to cover all quarters.")
  
  if (!is.finite(Lproxy) || Lproxy < 1) stop("Lproxy must be integer >= 1.")
  if (!is.finite(L_midas) || L_midas < 1) stop("L_midas must be integer >= 1.")
  Lproxy  <- as.integer(Lproxy)
  L_midas <- as.integer(L_midas)
  
  if (!is.finite(p_AR) || p_AR < 0) stop("p_AR must be integer >= 0.")
  p_AR <- as.integer(p_AR)
  
  if (length(r) != 2) stop("r must be c(r1,r2).")
  r1 <- as.integer(r[1]); r2 <- as.integer(r[2])
  if (r1 < 1 || r2 < 1 || r1 > p1 || r2 > p2) stop("Invalid r.")
  
  # dimnames
  dnq <- dimnames(X_lf)
  dates_q <- if (!is.null(dnq[[1]])) dnq[[1]] else as.character(seq_len(Tq))
  
  # ---- initialize y (EA proxy) and proxy matrix
  y1 <- as.numeric(y_q)
  if (standardize_y) y1 <- (y1 - mean(y1)) / sd(y1)
  
  coln <- if (Lproxy == 1) "y_ea" else c("y_ea", paste0("e", 1:(Lproxy - 1)))
  Z_raw <- matrix(NA_real_, nrow = Tq, ncol = Lproxy, dimnames = list(dates_q, coln))
  Z_raw[, 1] <- y1
  
  # storage
  yhat_list <- if (Lproxy > 1) vector("list", Lproxy - 1) else NULL
  e_list    <- if (Lproxy > 1) vector("list", Lproxy - 1) else NULL
  
  last_R <- NULL; last_C <- NULL; last_F_hf <- NULL; last_fit3 <- NULL
  
  # ---- if only one proxy, return immediately
  if (Lproxy == 1) {
    Z_out <- if (orthonormalize_final) qr.Q(qr(Z_raw)) else Z_raw
    colnames(Z_out) <- coln; rownames(Z_out) <- dates_q
    return(list(
      Z_raw = Z_raw, Z = Z_out,
      yhat_list = NULL, e_list = NULL,
      last_R = NULL, last_C = NULL, last_F_hf = NULL, last_fit3 = NULL
    ))
  }
  
  # ---- iterative residual proxies
  for (ell in 1:(Lproxy - 1)) {
    
    Z_curr <- Z_raw[, 1:ell, drop = FALSE]
    
    # PASS 1
    pass1 <- mf_tprf_pass1_ils(
      X_lf = X_lf,
      Z_q  = Z_curr,
      r    = c(r1, r2),
      orthonormalize_Z = orthonormalize_each_iter,
      center_Z = TRUE,
      maxit = ils_maxit,
      tol   = ils_tol
    )
    
    # PASS 2
    pass2 <- mf_tprf_pass2_factors(
      X_hf = X_hf,
      R    = pass1$R,
      C    = pass1$C,
      T_q  = Tq,
      scale_fac = 1 / (p1 * p2)
    )
    
    # PASS 3 (quarterly fitted for EA proxy)
    fit3 <- mf_tprf_pass3_fit_quarterly(
      y_q = y1,
      F1 = pass2$F1, F2 = pass2$F2, F3 = pass2$F3,
      p_AR = p_AR,
      L_midas = L_midas
    )
    
    yhat <- fit3$yhat_q
    e    <- y1 - yhat
    e[is.na(e)] <- y1[is.na(e)]  # same convention as your vector code
    
    yhat_list[[ell]] <- yhat
    e_list[[ell]]    <- e
    
    Z_raw[, ell + 1] <- e
    
    # save last
    last_R <- pass1$R
    last_C <- pass1$C
    last_F_hf <- pass2$F_hf
    last_fit3 <- fit3$fit
  }
  
  # ---- final orthonormal basis Z'Z=I (optional, as in main)
  Z_out <- if (orthonormalize_final) qr.Q(qr(scale(Z_raw, center = TRUE, scale = FALSE))) else Z_raw
  colnames(Z_out) <- coln; rownames(Z_out) <- dates_q
  
  list(
    Z_raw = Z_raw,
    Z     = Z_out,
    yhat_list = yhat_list,
    e_list    = e_list,
    last_R = last_R,
    last_C = last_C,
    last_F_hf = last_F_hf,
    last_fit3 = last_fit3,
    meta = list(
      Lproxy = Lproxy, L_midas = L_midas, p_AR = p_AR, r = c(r1, r2),
      standardize_y = standardize_y,
      orthonormalize_each_iter = orthonormalize_each_iter,
      orthonormalize_final = orthonormalize_final
    )
  )
}

# ==============================================================================
# Matrix MF-TPRF (clean): Autoproxy + final targeting + per-country UMIDAS + nowcast
# ==============================================================================

Tensor_MF_TPRF <- function(
    X_lf, X_hf, Y_q_all,
    proxy_name = "EA",
    Lproxy = 1,
    L_midas = 1,
    p_AR = 1,            # allow 0
    r = c(1, 1),
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8
) {
  
  # ---- checks
  stopifnot(length(dim(X_lf)) == 3, length(dim(X_hf)) == 3)
  if (anyNA(X_lf) || anyNA(X_hf)) stop("X_lf / X_hf contain NA. This version expects complete arrays.")
  
  T_q <- dim(X_lf)[1]
  T_m <- dim(X_hf)[1]
  p1  <- dim(X_lf)[2]
  p2  <- dim(X_lf)[3]
  
  if (dim(X_hf)[2] != p1 || dim(X_hf)[3] != p2) stop("X_hf must have same (p1,p2) as X_lf.")
  if (T_m < 3 * T_q) stop("Need at least 3*T_q months in X_hf to cover all quarters.")
  
  if (!proxy_name %in% colnames(Y_q_all)) stop("proxy_name not found in colnames(Y_q_all).")
  
  if (!is.finite(p_AR) || p_AR < 0) stop("p_AR must be integer >= 0.")
  p_AR <- as.integer(p_AR)
  
  r1 <- as.integer(r[1]); r2 <- as.integer(r[2])
  if (r1 < 1 || r2 < 1 || r1 > p1 || r2 > p2) stop("Invalid r.")
  
  # identify target countries (all except proxy column)
  countries <- setdiff(colnames(Y_q_all), proxy_name)
  
  y_proxy <- as.numeric(Y_q_all[, proxy_name])
  if (length(y_proxy) != T_q) stop("Proxy y length must match dim(X_lf)[1].")
  
  # dimnames
  dnq <- dimnames(X_lf)
  dnm <- dimnames(X_hf)
  dates_q <- if (!is.null(dnq[[1]])) dnq[[1]] else as.character(seq_len(T_q))
  dates_m <- if (!is.null(dnm[[1]])) dnm[[1]] else as.character(seq_len(T_m))
  
  # --------------------------------------------------
  # STEP 1: autoproxy Z_q (EA)
  # --------------------------------------------------
  auto <- matrix_mf_tprf_autoproxy(
    X_lf = X_lf, X_hf = X_hf, y_q = y_proxy,
    Lproxy = Lproxy, L_midas = L_midas, p_AR = p_AR,
    r = c(r1, r2),
    standardize_y = standardize_proxy,
    orthonormalize_each_iter = orthonormalize_each_iter,
    orthonormalize_final = orthonormalize_final_Z,
    ils_maxit = ils_maxit, ils_tol = ils_tol
  )
  Z_q <- auto$Z  # final proxy basis (possibly orthonormal)
  
  # --------------------------------------------------
  # STEP 2: final Pass 1 (targeting) using ALL proxies
  # --------------------------------------------------
  pass1 <- mf_tprf_pass1_ils(
    X_lf = X_lf,
    Z_q  = Z_q,
    r = c(r1, r2),
    orthonormalize_Z = TRUE,  # enforce clean basis in targeting
    center_Z = TRUE,
    maxit = ils_maxit,
    tol = ils_tol
  )
  Rhat <- pass1$R
  Chat <- pass1$C
  
  # --------------------------------------------------
  # STEP 3: factors from HF tensor + within-quarter blocks
  # --------------------------------------------------
  pass2 <- mf_tprf_pass2_factors(
    X_hf = X_hf,
    R = Rhat, C = Chat,
    T_q = T_q,
    scale_fac = 1 / (p1 * p2)
  )
  
  # --------------------------------------------------
  # STEP 4: per-country quarterly fit + monthly nowcast
  # --------------------------------------------------
  results_country <- vector("list", length(countries))
  names(results_country) <- countries
  
  for (cc in countries) {
    
    y_cc <- as.numeric(Y_q_all[, cc])
    if (length(y_cc) != T_q) stop(paste0("Country ", cc, ": y length mismatch."))
    
    fit3 <- mf_tprf_pass3_fit_quarterly(
      y_q = y_cc,
      F1 = pass2$F1, F2 = pass2$F2, F3 = pass2$F3,
      p_AR = p_AR,
      L_midas = L_midas
    )
    
    y_now <- mf_tprf_nowcast_monthly(
      beta0 = fit3$beta0,
      rho_hat = fit3$rho_hat,
      beta_mat = fit3$beta_mat,
      y_q = y_cc,
      F1 = pass2$F1, F2 = pass2$F2, F3 = pass2$F3,
      F_next1 = pass2$F_next1, F_next2 = pass2$F_next2, F_next3 = pass2$F_next3,
      p_AR = p_AR,
      L_midas = L_midas,
      T_q = T_q,
      T_m = T_m
    )
    
    results_country[[cc]] <- list(
      fit3 = fit3$fit,
      beta0 = fit3$beta0,
      rho_hat = fit3$rho_hat,
      beta_mat = fit3$beta_mat,
      start_tau = fit3$start_tau,
      yhat_q = fit3$yhat_q,
      y_nowcast = y_now
    )
  }
  
  # --------------------------------------------------
  # Output
  # --------------------------------------------------
  list(
    proxy_name = proxy_name,
    Z_q = Z_q,
    autoproxy = auto,
    R = Rhat,
    C = Chat,
    factors = list(
      F_hf = pass2$F_hf,
      F1 = pass2$F1, F2 = pass2$F2, F3 = pass2$F3,
      rem = pass2$rem,
      F_next1 = pass2$F_next1, F_next2 = pass2$F_next2, F_next3 = pass2$F_next3,
      dates_m = dates_m, dates_q = dates_q
    ),
    by_country = results_country,
    meta = list(
      Lproxy = Lproxy, L_midas = L_midas, p_AR = p_AR, r = c(r1, r2),
      standardize_proxy = standardize_proxy
    )
  )
}