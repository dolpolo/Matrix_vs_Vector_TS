# ==============================================================================
# InitialCond(): Banbura–Modugno style initialization (monolithic)
# - Input: X_std (T x N) standardized panel with NA
# - Mixed frequency: quarterly variables (last NQ cols) handled via:
#     restr = "MM"        -> L = 5, weights = (1,2,3,2,1)
#     restr = "stock_flow"-> L = 3, weights = sum (flow) or mean (stock) using agg
# - Factors: VAR(p_factor) on PCA factors from X_imp
# - Idiosyncratic: AR(q_idio) for EVERY series, same q for all
# - State padded so that factor/idiosyncratic stacks are at least length L
# - Output: A, C, Q, R, Z0, V0 + init objects
# ==============================================================================

# restr <- "stock_flow"
InitialCond <- function(
    X_std,                 # T x N with NA, already standardized
    r,                     # number of factors
    p_factor,              # VAR order for factors
    q_idio,                # AR order for idiosyncratic components (global)
    NM, NQ,                # N = NM + NQ, quarterly are last NQ columns
    restr = c("MM", "stock_flow"),
    agg = NULL,            # length N or length NQ: 1 stock (mean), 2 flow (sum), used if stock_flow
    kappa = 1e-4,          # measurement noise R = kappa I (small)
    force_stationary = TRUE,
    clip = 0.99,
    ridge = 1e-6           # ridge in XP factor WLS inversion
) {
  restr <- match.arg(restr)
  X_std <- as.matrix(X_std)
  Tt <- nrow(X_std)
  N  <- ncol(X_std)
  
  stopifnot(N == NM + NQ)
  stopifnot(r >= 1, p_factor >= 1, q_idio >= 0)
  stopifnot(NM >= 0, NQ >= 1)
  
  # -------------------------------
  # 0) Decide L and quarterly weights ω
  # -------------------------------
  if (restr == "MM") {
    L <- 5L
    w_MM <- c(1, 2, 3, 2, 1)
  } else { # stock_flow
    L <- 3L
  }
  
  # Normalize agg handling:
  # allow agg length N (all vars) or length NQ (only quarterly).
  if (restr == "stock_flow") {
    if (is.null(agg)) stop("agg must be provided when restr='stock_flow'.")
    if (length(agg) == NQ) {
      aggN <- c(rep(NA_integer_, NM), agg)
    } else {
      stopifnot(length(agg) == N)
      aggN <- agg
    }
    aggN <- as.integer(aggN)
  }
  
  # dopo aver creato aggN (lunghezza N) se restr == "stock_flow"
  if (restr == "stock_flow") {
    agg_used <- aggN
  } else {
    agg_used <- NULL
  }
  
  get_wQ <- function(i) {
    # i is global index (NM+1 ... N)
    if (restr == "MM") return(w_MM)
    
    ai <- aggN[i]
    if (is.na(ai)) stop("agg for a quarterly variable is NA: provide agg_q (length NQ) or full agg (length N).")
    if (ai == 2L) return(rep(1, L))        # flow: sum
    if (ai == 1L) return(rep(1 / L, L))    # stock: mean
    stop("agg must be coded 1 (stock) or 2 (flow).")
  }
  
  # -------------------------------
  # 1) XP-style initialization for missing data (Xiong–Pelger)
  #    -> returns: Lambda (N x r), F_hat (T x r), C_hat (T x N)
  # -------------------------------
  xp <- estimate_factors_XP(X_std, r = r, ridge = ridge)
  
  # Keep XP objects only for diagnostics / imputing:
  Lambda_XP <- xp$Lambda          # N x r, Lambda'Lambda/N = I
  F_XP      <- xp$F_hat           # T x r (may have NA if an entire row is missing)
  C_hat     <- xp$C_hat           # T x N common component (std scale)
  
  # -------------------------------
  # 2) Impute NA using XP common component
  # -------------------------------
  X_imp <- X_std
  miss  <- is.na(X_imp)
  if (any(miss)) X_imp[miss] <- C_hat[miss]
  
  # -------------------------------
  # 3) PCA init on X_imp (truncated SVD)
  #    then normalize to Lambda'Lambda/N = I (XP convention)
  # -------------------------------
  s <- svd(X_imp, nu = r, nv = r)
  
  U <- s$u[, 1:r, drop = FALSE]                 # T x r
  F0 <- sqrt(Tt) * U                            # so F0'F0/T = I
  Lambda0 <- crossprod(X_imp, F0) / Tt          # N x r
  
  # rescale to Lambda'Lambda/N = I (keep common component unchanged)
  D  <- crossprod(Lambda0) / N
  ee <- eigen(D, symmetric = TRUE)
  
  vals <- pmax(as.numeric(ee$values), 1e-12)
  r_eff <- length(vals)
  
  # diag "safe": NON usare diag(x) quando x è scalare
  D_inv_sqrt <- ee$vectors %*% diag(1 / sqrt(vals), nrow = r_eff, ncol = r_eff) %*% t(ee$vectors)
  D_sqrt     <- ee$vectors %*% diag(sqrt(vals),     nrow = r_eff, ncol = r_eff) %*% t(ee$vectors)
  
  Lambda_hat <- Lambda0 %*% D_inv_sqrt
  F_hat      <- F0      %*% D_sqrt
  
  # -------------------------------
  # 4) VAR(p_factor) on factors (OLS)
  # -------------------------------
  stopifnot(Tt > p_factor + 5)
  
  Y <- F_hat[(p_factor + 1):Tt, , drop = FALSE]  # (T-p) x r
  Z <- do.call(cbind, lapply(1:p_factor, function(l)
    F_hat[(p_factor + 1 - l):(Tt - l), , drop = FALSE]
  ))                                             # (T-p) x (r*p)
  
  # OLS multivariate
  B <- solve(crossprod(Z), crossprod(Z, Y))      # (r*p) x r
  Ures <- Y - Z %*% B
  Q_u  <- crossprod(Ures) / nrow(Ures)           # r x r (innovation covariance)
  
  # A_l matrices (r x r): f_t = A1 f_{t-1} + ... + Ap f_{t-p} + u_t
  A_list <- lapply(1:p_factor, function(l) {
    t(B[((l - 1) * r + 1):(l * r), , drop = FALSE])
  })
  
  # -------------------------------
  # 5) Residuals for idio init
  # -------------------------------
  Xhat_common <- F_hat %*% t(Lambda_hat)         # T x N
  Ehat <- X_imp - Xhat_common                    # T x N (complete)
  
  # -------------------------------
  # 6) Estimate AR(q_idio) for EACH series i=1..N (OLS), optional stationarity clip
  # -------------------------------
  if (q_idio == 0) {
    Phi  <- matrix(0, nrow = N, ncol = 0)
    Sig2 <- apply(Ehat, 2, function(e) mean(as.numeric(e)^2))
  } else {
    Phi  <- matrix(0, nrow = N, ncol = q_idio)
    Sig2 <- numeric(N)
    
    for (i in 1:N) {
      e <- as.numeric(Ehat[, i])
      if (length(e) <= q_idio + 5) {
        Phi[i, ] <- rep(0, q_idio)
        Sig2[i]  <- var(e, na.rm = TRUE)
        next
      }
      
      Yi <- e[(q_idio + 1):length(e)]
      Zi <- do.call(cbind, lapply(1:q_idio, function(l)
        e[(q_idio + 1 - l):(length(e) - l)]
      ))
      
      ph <- as.numeric(solve(crossprod(Zi), crossprod(Zi, Yi)))
      ui <- as.numeric(Yi - Zi %*% ph)
      s2 <- mean(ui^2)
      
      if (force_stationary) {
        # lightweight stabilization for init
        ph <- pmax(pmin(ph, clip), -clip)
        if (q_idio > 1) ph <- 0.95 * ph
      }
      
      Phi[i, ] <- ph
      Sig2[i]  <- s2
    }
  }
  
  # -------------------------------
  # 7) Build augmented state transition A and covariance Q
  #    - factors: companion length pC_f = max(p_factor, L)
  #    - idio: companion length pC_e = max(max(q_idio,1), L)
  # -------------------------------
  pC_f <- max(p_factor, L)
  pC_e <- max(max(q_idio, 1L), L)  # if q_idio=0 still need at least L lags for quarterly aggregation
  
  # 7a) Factor companion (r*pC_f)
  A_f <- matrix(0, r * pC_f, r * pC_f)
  for (l in 1:p_factor) {
    A_f[1:r, ((l - 1) * r + 1):(l * r)] <- A_list[[l]]
  }
  if (pC_f > 1) {
    A_f[(r + 1):(r * pC_f), 1:(r * (pC_f - 1))] <- diag(r * (pC_f - 1))
  }
  
  Q_f <- matrix(0, r * pC_f, r * pC_f)
  Q_f[1:r, 1:r] <- Q_u
  
  # 7b) Idio block: N independent AR(q_idio) padded to pC_e
  A_e <- matrix(0, N * pC_e, N * pC_e)
  Q_e <- matrix(0, N * pC_e, N * pC_e)
  
  for (i in 1:N) {
    idx <- ((i - 1) * pC_e + 1):(i * pC_e)
    
    Ablock <- matrix(0, pC_e, pC_e)
    if (q_idio > 0) {
      Ablock[1, 1:q_idio] <- Phi[i, ]
    }
    if (pC_e > 1) {
      Ablock[2:pC_e, 1:(pC_e - 1)] <- diag(pC_e - 1)
    }
    
    A_e[idx, idx] <- Ablock
    Q_e[idx[1], idx[1]] <- Sig2[i]   # innovation variance only on e_t
  }
  
  # 7c) Global A, Q (block diagonal)
  A <- as.matrix(Matrix::bdiag(A_f, A_e))
  Q <- as.matrix(Matrix::bdiag(Q_f, Q_e))
  
  # -------------------------------
  # 8) Measurement matrix C (N x k_state)
  #    Monthly:   x_{i,t} = lambda_i' f_t + e_{i,t}
  #    Quarterly: x_{i,t} = sum_{ell=0}^{L-1} ω_ell (lambda_i' f_{t-ell} + e_{i,t-ell})
  #              (to be used at quarter-end via missingness pattern in the Kalman filter)
  # -------------------------------
  k_fac   <- r * pC_f
  k_e     <- N * pC_e
  k_state <- k_fac + k_e
  
  C <- matrix(0, N, k_state)
  
  # 8a) Monthly variables (i = 1..NM)
  for (i in 1:NM) {
    C[i, 1:r] <- Lambda_hat[i, ]
    col_e_i <- k_fac + (i - 1) * pC_e + 1
    C[i, col_e_i] <- 1
  }
  
  # 8b) Quarterly variables (i = NM+1..N)
  for (j in 1:NQ) {
    i <- NM + j
    w <- get_wQ(i)                    # length L
    
    # factors: first L blocks of factor state
    for (ell in 0:(L - 1)) {
      cols_f <- (ell * r + 1):((ell + 1) * r)
      C[i, cols_f] <- C[i, cols_f] + w[ell + 1] * Lambda_hat[i, ]
    }
    
    # idio: first L elements of that series idio state
    idx_e_cols <- k_fac + (i - 1) * pC_e + (1:L)
    C[i, idx_e_cols] <- C[i, idx_e_cols] + w
  }
  
  # -------------------------------
  # 9) Measurement noise
  # -------------------------------
  R <- diag(kappa, N)
  
  # -------------------------------
  # 10) Initial state mean/cov
  # -------------------------------
  Z0 <- rep(0, k_state)
  V0 <- diag(1e-2, k_state)
  
  list(
    A = A, C = C, Q = Q, R = R,
    Z0 = Z0, V0 = V0,
    init = list(
      X_imp      = X_imp,
      F_hat      = F_hat,
      Lambda_hat = Lambda_hat,
      Ehat       = Ehat,
      XP         = list(Lambda = Lambda_XP, F_hat = F_XP, C_hat = C_hat),
      Phi        = Phi,
      Sig2       = Sig2
    ),
    meta = list(
      r = r, p_factor = p_factor, q_idio = q_idio,
      NM = NM, NQ = NQ, N = N,
      restr = restr, L = L,
      pC_f = pC_f, pC_e = pC_e,
      kappa = kappa,
      agg_used = agg_used 
    )
  )
}
