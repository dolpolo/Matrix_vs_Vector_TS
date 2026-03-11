# ==============================================================================
# Autoproxy function
# ==============================================================================

Autoproxy_TPRF <- function(X, y, L) {
  stopifnot(is.matrix(X), length(y) == nrow(X), L >= 1)
  
  T <- nrow(X)
  N <- ncol(X)
  
  # Z: proxies (T x k)
  Z <- matrix(y, ncol = 1)
  colnames(Z) <- "y"
  
  if (L == 1) return(Z)
  
  for (ell in 1:(L - 1)) {
    
    # ------------------------------------------------------------
    # Step 1: for each i, OLS x_i,t on [1, Z_t]  -> loadings Phi (N x k)
    # ------------------------------------------------------------
    k <- ncol(Z)
    Phi <- matrix(NA_real_, nrow = N, ncol = k)
    
    for (i in 1:N) {
      df_i <- data.frame(xi = X[, i], Z)
      fit_i <- lm(xi ~ ., data = df_i)     # includes intercept
      Phi[i, ] <- coef(fit_i)[-1]          # drop intercept, keep slopes on Z
    }
    
    # ------------------------------------------------------------
    # Step 2: targeted factors F_t = X_t * Phi * (Phi'Phi)^{-1}
    # ------------------------------------------------------------
    F_hat <- X %*% Phi %*% solve(crossprod(Phi))  # (T x k)
    
    # ------------------------------------------------------------
    # Step 3: one-step-ahead regression y_{t} on [1, F_{t-1}]
    #        i.e. y[2:T] on F_hat[1:(T-1),]
    # ------------------------------------------------------------
    y_t   <- y[-1]
    F_lag <- F_hat[-T, , drop = FALSE]
    
    df_y <- data.frame(y = y_t, F_lag)
    fit_y <- lm(y ~ ., data = df_y)        # includes intercept
    
    y_hat <- rep(NA_real_, T)
    y_hat[-1] <- as.numeric(predict(fit_y, newdata = df_y))
    
    e <- y - y_hat
    e[1] <- 0                               # neutral handling of first NA
    
    Z <- cbind(Z, e)
    colnames(Z)[ncol(Z)] <- paste0("e", ell)
  }
  
  Z
}


# ==============================================================================
# TPRF function: 
#===============================================================================
#' PROXY
#' LOADINGS 
#' FATTORI
#' TARGET PREDICTION
# ==============================================================================

TPRF <- function(X, y, L) {
  stopifnot(is.matrix(X), length(y) == nrow(X), L >= 1)
  
  T <- nrow(X)
  N <- ncol(X)
  
  # Step 0: proxies
  Z <- Autoproxy_TPRF(X, y, L)
  k <- ncol(Z)
  
  if (is.null(k) || k < 1) {
    stop("Autoproxy_TPRF returned zero proxy columns.")
  }
  
  # Step 1: OLS x_i,t on [1, Z_t] -> Phi_hat (N x k)
  Phi_hat <- matrix(NA_real_, nrow = N, ncol = k)
  for (i in 1:N) {
    df_i <- data.frame(xi = X[, i], Z)
    fit_i <- lm(xi ~ ., data = df_i)
    Phi_hat[i, ] <- coef(fit_i)[-1]
  }
  
  # Step 2: targeted factors
  F_hat <- X %*% Phi_hat %*% solve(crossprod(Phi_hat))  # T x k
  F_hat <- as.matrix(F_hat)
  colnames(F_hat) <- paste0("F", 1:k)
  
  # Step 3: y[2:T] on F_hat[1:(T-1),]
  y_t   <- y[-1]
  F_lag <- F_hat[-T, , drop = FALSE]
  df_y  <- data.frame(y = y_t, F_lag)
  fit   <- lm(y ~ ., data = df_y)
  
  y_hat <- rep(NA_real_, T)
  y_hat[-1] <- as.numeric(predict(fit, newdata = df_y))
  
  # One-step-ahead forecast: y_{T+1} using F_T
  newdata_fc <- as.data.frame(F_hat[T, , drop = FALSE])
  y_forecast <- as.numeric(predict(fit, newdata = newdata_fc))
  
  list(
    Z          = Z,
    Phi_hat    = Phi_hat,
    F_hat      = F_hat,
    fit        = fit,
    y_hat      = y_hat,
    y_forecast = y_forecast,
    resid      = y - y_hat
  )
}
