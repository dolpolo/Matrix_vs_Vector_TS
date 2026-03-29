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

# ==============================================================================
# HELPERS
# ==============================================================================

safe_solve <- function(A, B = NULL) {
  if (is.null(B)) {
    tryCatch(
      solve(A),
      error = function(e) MASS::ginv(A)
    )
  } else {
    tryCatch(
      solve(A, B),
      error = function(e) MASS::ginv(A) %*% B
    )
  }
}

# ------------------------------------------------------------------------------
# Regress all columns of X on [1, Z] in one shot
# Returns only slopes on Z (intercept dropped), like your original code
# X: T x N
# Z: T x k
# output: N x k
# ------------------------------------------------------------------------------
ols_slopes_multivariate <- function(X, Z) {
  stopifnot(is.matrix(X), is.matrix(Z), nrow(X) == nrow(Z))
  
  W <- cbind(1, Z)                  # T x (k+1)
  XtX <- crossprod(W)               # (k+1) x (k+1)
  XtY <- crossprod(W, X)            # (k+1) x N
  
  B <- safe_solve(XtX, XtY)         # (k+1) x N
  Phi <- t(B[-1, , drop = FALSE])   # N x k
  
  Phi
}

# ------------------------------------------------------------------------------
# OLS of y on [1, X] using linear algebra
# X: T x k
# y: length T
# ------------------------------------------------------------------------------
ols_fit_with_intercept <- function(y, X) {
  stopifnot(is.numeric(y), is.matrix(X), length(y) == nrow(X))
  
  W <- cbind(1, X)                  # T x (k+1)
  XtX <- crossprod(W)
  Xty <- crossprod(W, y)
  
  beta <- as.numeric(safe_solve(XtX, Xty))
  fitted <- as.numeric(W %*% beta)
  resid  <- as.numeric(y - fitted)
  
  list(
    coefficients = beta,
    fitted.values = fitted,
    residuals = resid
  )
}

# ------------------------------------------------------------------------------
# Predict from coefficients of y ~ 1 + X
# beta includes intercept
# Xnew: m x k
# ------------------------------------------------------------------------------
ols_predict_with_intercept <- function(beta, Xnew) {
  Xnew <- as.matrix(Xnew)
  Wnew <- cbind(1, Xnew)
  as.numeric(Wnew %*% beta)
}


# ==============================================================================
# Autoproxy function with generic forecast horizon h (optimized)
# ==============================================================================

Autoproxy_TPRF_h <- function(X, y, L, h = 1) {
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(y))
  stopifnot(length(y) == nrow(X))
  stopifnot(L >= 1, h >= 1)
  
  T <- nrow(X)
  
  if (T <= h) {
    stop("Not enough observations for the requested forecast horizon h.")
  }
  
  y <- as.numeric(y)
  
  # Initial proxy set: current y
  Z <- matrix(y, ncol = 1)
  colnames(Z) <- "y"
  
  if (L == 1) return(Z)
  
  for (ell in 1:(L - 1)) {
    
    # ------------------------------------------------------------
    # Step 1: regress all predictors x_i,t on [1, Z_t] in one shot
    # ------------------------------------------------------------
    Phi <- ols_slopes_multivariate(X, Z)   # N x k
    
    # ------------------------------------------------------------
    # Step 2: targeted factors
    # ------------------------------------------------------------
    Ginv <- safe_solve(crossprod(Phi))
    F_hat <- X %*% Phi %*% Ginv
    F_hat <- as.matrix(F_hat)
    
    # ------------------------------------------------------------
    # Step 3: direct h-step regression y_t on F_{t-h}
    #         y[(h+1):T] on F_hat[1:(T-h), ]
    # ------------------------------------------------------------
    y_t   <- y[(h + 1):T]
    F_lag <- F_hat[1:(T - h), , drop = FALSE]
    
    fit_y <- ols_fit_with_intercept(y_t, F_lag)
    
    y_hat <- rep(NA_real_, T)
    y_hat[(h + 1):T] <- fit_y$fitted.values
    
    # Residual proxy
    e <- y - y_hat
    e[1:h] <- 0
    
    Z <- cbind(Z, e)
    colnames(Z)[ncol(Z)] <- paste0("e", ell)
  }
  
  Z
}


# ==============================================================================
# TPRF with generic forecast horizon h (optimized)
# ==============================================================================

TPRF_h <- function(X, y, L, h = 1) {
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(y))
  stopifnot(length(y) == nrow(X))
  stopifnot(L >= 1, h >= 1)
  
  T <- nrow(X)
  
  if (T <= h) {
    stop("Not enough observations for the requested forecast horizon h.")
  }
  
  y <- as.numeric(y)
  
  # Step 0: proxies
  Z <- Autoproxy_TPRF_h(X = X, y = y, L = L, h = h)
  k <- ncol(Z)
  
  if (is.null(k) || k < 1) {
    stop("Autoproxy_TPRF_h returned zero proxy columns.")
  }
  
  # Step 1: OLS x_i,t on [1, Z_t] -> Phi_hat, all at once
  Phi_hat <- ols_slopes_multivariate(X, Z)   # N x k
  
  # Step 2: targeted factors
  Ginv <- safe_solve(crossprod(Phi_hat))
  F_hat <- X %*% Phi_hat %*% Ginv
  F_hat <- as.matrix(F_hat)
  colnames(F_hat) <- paste0("F", 1:k)
  
  # Step 3: direct h-step regression y_t on F_{t-h}
  y_t   <- y[(h + 1):T]
  F_lag <- F_hat[1:(T - h), , drop = FALSE]
  fit   <- ols_fit_with_intercept(y_t, F_lag)
  
  y_hat <- rep(NA_real_, T)
  y_hat[(h + 1):T] <- fit$fitted.values
  
  # h-step-ahead forecast: y_{T+h} using F_T
  y_forecast <- ols_predict_with_intercept(
    beta = fit$coefficients,
    Xnew = F_hat[T, , drop = FALSE]
  )
  
  list(
    h          = h,
    Z          = Z,
    Phi_hat    = Phi_hat,
    F_hat      = F_hat,
    fit        = fit,
    y_hat      = y_hat,
    y_forecast = y_forecast,
    resid      = y - y_hat
  )
}