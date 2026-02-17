# ==============================================================================
# KALMAN (for arbitrary NA patterns)
# ==============================================================================

kalman.na <- function(par, y, k, start, O){
  
  # Size
  n <- ncol(y)   # numero di periodi (Tobs)
  d <- nrow(y)   # numero di serie (N)
  
  # Parameters
  tA <- par$t    # state transition (A)
  Z  <- par$l    # loadings (C)
  Q  <- par$q    # state variance
  H  <- par$h    # measurement variance
  
  # Initialize containers
  f.f <- f.u <- f.s <- matrix(0, k, n)
  v   <- matrix(0, d, n)
  P.f <- P.u <- P.s <- array(0, c(k, k, n))
  S.1 <- array(0, c(d, d, n))   # F_t^{-1}
  L   <- array(0, c(k, k, n))
  qK  <- array(0, c(k, d, n))   # quasi Kalman gain
  C.f <- C.s <- array(0, c(k, k, n))
  
  # Log-likelihood
  const  <- log(2*pi)
  loglik <- 0
  
  # --------------------------
  # Kalman filter
  # --------------------------
  f.f[,1]  <- f.u[,1]  <- start$f
  P.f[,,1] <- P.u[,,1] <- start$P
  
  # t = 1
  obs1 <- O[[1]]
  if (length(obs1) > 0) {
    Z1 <- Z[obs1, , drop = FALSE]      # m1 x k
    H1 <- H[obs1, obs1, drop = FALSE]  # m1 x m1
    P1 <- as.matrix(P.f[,,1])          # k x k
    
    v1 <- y[obs1, 1, drop = FALSE] - Z1 %*% f.f[,1, drop = FALSE]    # m1 x 1
    F1 <- Z1 %*% P1 %*% t(Z1) + H1                                   # m1 x m1
    S1 <- solve(F1)                                                  # m1 x m1
    
    v[obs1, 1]       <- v1
    S.1[obs1,obs1,1] <- S1
    qK[,obs1,1]      <- P1 %*% t(Z1) %*% S1                          # k x m1
    
    f.u[,1]  <- f.f[,1, drop = FALSE] + qK[,obs1,1] %*% v1
    P.u[,,1] <- P1 - qK[,obs1,1] %*% Z1 %*% P1
    L[,,1]   <- tA - tA %*% qK[,obs1,1] %*% Z1
    C.f[,,1] <- P1 %*% t(L[,,1])
    
    # loglik
    vi    <- as.matrix(v1)
    S1_i  <- as.matrix(S1)
    quad  <- as.numeric(t(vi) %*% S1_i %*% vi)
    detS1 <- determinant(S1_i, logarithm = TRUE)
    logdetF <- -as.numeric(detS1$modulus)
    loglik  <- loglik - 0.5 * (length(obs1)*const + logdetF + quad)
  } else {
    # nessuna osservazione al tempo 1
    f.u[,1]  <- f.f[,1, drop = FALSE]
    P.u[,,1] <- P.f[,,1]
    L[,,1]   <- tA
    C.f[,,1] <- P.f[,,1] %*% t(L[,,1])
  }
  
  # t = 2,...,n
  for (i in 2:n) {
    # Predizione
    P_prev    <- as.matrix(P.u[,,i-1])
    f.f[,i]   <- tA %*% f.u[,i-1, drop = FALSE]
    P.f[,,i]  <- tA %*% P_prev %*% t(tA) + Q
    
    obs <- O[[i]]
    if (length(obs) > 0) {
      Zi <- Z[obs, , drop = FALSE]      # m x k
      Hi <- H[obs, obs, drop = FALSE]   # m x m
      Pi <- as.matrix(P.f[,,i])         # k x k
      
      vi   <- y[obs, i, drop = FALSE] - Zi %*% f.f[,i, drop = FALSE]   # m x 1
      Fi   <- Zi %*% Pi %*% t(Zi) + Hi                                  # m x m
      S1_i <- solve(Fi)                                                 # m x m
      
      v[obs, i]       <- vi
      S.1[obs,obs,i]  <- S1_i
      qK[,obs,i]      <- Pi %*% t(Zi) %*% S1_i                          # k x m
      
      f.u[,i]  <- f.f[,i, drop = FALSE] + qK[,obs,i] %*% vi
      P.u[,,i] <- Pi - qK[,obs,i] %*% Zi %*% Pi
      
      L[,,i]   <- tA - tA %*% qK[,obs,i] %*% Zi
      C.f[,,i] <- Pi %*% t(L[,,i])
      
      # loglik
      vi_m   <- as.matrix(vi)                                           # m x 1
      S1_m   <- as.matrix(S1_i)                                         # m x m
      quad   <- as.numeric(t(vi_m) %*% S1_m %*% vi_m)
      detS1  <- determinant(S1_m, logarithm = TRUE)
      logdetF <- -as.numeric(detS1$modulus)
      loglik <- loglik - 0.5 * (length(obs)*const + logdetF + quad)
      
    } else {
      # nessuna osservazione al tempo i
      f.u[,i]  <- f.f[,i, drop = FALSE]
      P.u[,,i] <- P.f[,,i]
      L[,,i]   <- tA
      C.f[,,i] <- P.f[,,i] %*% t(L[,,i])
    }
  }
  
  # --------------------------
  # Kalman smoother
  # --------------------------
  C.s[,,n] <- t(C.f[,,n])
  r <- matrix(0, k, 1)
  N <- matrix(0, k, k)
  
  for (i in n:1) {
    if (i < n) {
      C.s[,,i] <- t(C.f[,,i] %*% (diag(k) - N %*% P.f[,,i+1]))
    }
    obs <- O[[i]]
    if (length(obs) > 0) {
      Zi  <- Z[obs, , drop = FALSE]                                # m x k
      S1i <- S.1[obs,obs,i]                                        # m x m (non usare drop=FALSE)
      S1i <- as.matrix(S1i)                                        # forzo 2D
      
      vi  <- v[obs, i, drop = FALSE]                               # m x 1
      
      r <- t(Zi) %*% S1i %*% vi + t(L[,,i]) %*% r                  # k x 1
      N <- t(Zi) %*% S1i %*% Zi + t(L[,,i]) %*% N %*% L[,,i]       # k x k
    } else {
      r <- t(L[,,i]) %*% r
      N <- t(L[,,i]) %*% N %*% L[,,i]
    }
    Pi       <- as.matrix(P.f[,,i])
    f.s[,i]  <- f.f[,i, drop = FALSE] + Pi %*% r
    P.s[,,i] <- Pi - Pi %*% N %*% Pi
  }
  
  return(list(
    ff     = f.f,
    Pf     = P.f,
    fs     = f.s,
    Ps     = P.s,
    Cs     = C.s,
    loglik = as.numeric(loglik)
  ))
}




## =========================================================
## 1. Kalman smoother + fitted values DFM (Last Run)
## =========================================================

dfm_kalman_fit <- function(X_in, res, n_lags_Q) {
  # X_in : T x N (standardizzata)
  # res  : output di DFM_EM (A, C, Q, R, Z0, V0)
  # n_lags_Q : numero di lag trimestrali (es. 5 per 1-2-3-2-1)
  
  X  <- as.matrix(X_in)
  Tt <- nrow(X)
  N  <- ncol(X)
  
  # stessa finestra usata nell'EM (a partire da n_lags_Q)
  Xf   <- X[n_lags_Q:Tt, , drop = FALSE]   # T_eff x N
  Tobs <- nrow(Xf)
  
  A  <- res$A
  C  <- res$C
  Q  <- res$Q
  R  <- res$R
  Z0 <- res$Z0
  V0 <- res$V0
  
  k <- nrow(A)   # dimensione dello stato
  
  # Kalman smoother con NA arbitrari
  kf <- kalman.na(
    par   = list(t = A, l = C, q = Q, h = R),
    y     = t(Xf),  # N x T_eff
    k     = k,
    start = list(f = Z0, P = V0),
    O     = lapply(1:Tobs, function(t) which(!is.na(Xf[t, ])))
  )
  
  Z_sm <- kf$fs  # k x T_eff
  
  # Fitted values: X_hat_t = C * Z_t   →  N x T_eff
  X_hat <- t(C %*% Z_sm)   # T_eff x N
  
  # Residui
  resid <- Xf - X_hat
  
  list(
    X_fit   = X_hat,                # T_eff x N
    resid   = resid,                # T_eff x N
    Z_sm    = Z_sm,                 # k x T_eff
    t_index = n_lags_Q:Tt           # indici temporali rispetto a X_in
  )
}

