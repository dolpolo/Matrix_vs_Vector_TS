# ==============================================================================
# 1. All-purpose covariance estimator (Xiong & Pelger, 2020)
# ==============================================================================

all_purpose_covariance <- function(X_std) {
  # X_std: T x N standardized data (rows: time, cols: series)
  #        NA entries indicate missing observations.
  
  Tt <- nrow(X_std)
  N  <- ncol(X_std)
  
  # Missingness indicator: 1 = observed, 0 = missing
  M_obs <- ifelse(is.na(X_std), 0, 1)
  
  # Replace NA with 0 (used later for factor estimation)
  X_filled <- X_std
  X_filled[is.na(X_filled)] <- 0
  
  # Pairwise covariance using only common observed dates
  Sigma_hat <- matrix(0, nrow = N, ncol = N)
  
  for (i in 1:N) {
    for (j in i:N) {
      idx_ij <- which(M_obs[, i] == 1 & M_obs[, j] == 1)
      m_ij   <- length(idx_ij)
      
      if (m_ij > 0) {
        cov_ij <- sum(X_std[idx_ij, i] * X_std[idx_ij, j]) / m_ij
      } else {
        cov_ij <- 0
      }
      
      Sigma_hat[i, j] <- cov_ij
      Sigma_hat[j, i] <- cov_ij
    }
  }
  
  list(
    Sigma_tilde = Sigma_hat,  # N x N covariance estimator
    Missingness = M_obs,      # T x N, 0/1
    X_filled    = X_filled    # T x N, NA replaced by 0
  )
}




# ==============================================================================
# 3. Factor estimation à la Xiong & Pelger (given r)
#    Uses all-purpose covariance + missingness weights
# ==============================================================================

estimate_factors_XP <- function(X_std, r, ridge = 1e-6) {
  # X_std: T x N standardized data with NAs
  # r    : number of factors
  
  
  # result from all_purpose_covariance(X_std)
  cov_out <- all_purpose_covariance(X_std)
  
  Sigma    <- cov_out$Sigma_tilde
  W        <- cov_out$Missingness   # T x N
  X_filled <- cov_out$X_filled      # T x N
  
  Tt <- nrow(X_std)
  N  <- ncol(X_std)
  
  # PCA on (1/N) * Sigma
  Sigma_star <- Sigma / N
  eig        <- eigen(Sigma_star, symmetric = TRUE)
  V          <- eig$vectors[, 1:r, drop = FALSE]  # N x r
  
  # Loadings Λ = sqrt(N) * V_r
  Lambda <- sqrt(N) * V                            # N x r
  
  # Time-varying factors (XP formula)
  F_hat <- matrix(NA_real_, nrow = Tt, ncol = r)
  
  for (t in 1:Tt) {
    if (all(W[t, ] == 0)) {
      # No observed series at time t → no information
      next
    }
    
    W_t <- diag(W[t, ])          # N x N
    Y_t <- X_filled[t, ]         # N vector (no NAs)
    
    A_t <- t(Lambda) %*% W_t %*% Lambda         # r x r
    b_t <- t(Lambda) %*% W_t %*% Y_t            # r x 1
    
    # Small ridge term for numerical stability
    F_hat[t, ] <- solve(A_t + ridge * diag(r), b_t)
  }
  
  # Common component
  C_hat <- F_hat %*% t(Lambda)   # T x N
  
  list(
    Lambda = Lambda,
    F_hat  = F_hat,
    C_hat  = C_hat
  )
}


# ==============================================================================
# Initialization for EM (Hepenstrick & Marcellino)
# XP covariance + ER + XP factors → replace NAs with common component
# ==============================================================================

init_XP_ER <- function(X_std, Kmax = 15) {
  # 1. Covariance estimator
  cov_out <- all_purpose_covariance(X_std)
  
  # 2. Number of factors via ER
  ER_out <- select_num_factors_ER(cov_out$Sigma_tilde, Kmax = Kmax)
  r      <- ER_out$r
  
  # 3. Factors and common component (XP)
  fac_out <- estimate_factors_XP(X_std, r)
  C_hat   <- fac_out$C_hat
  
  # 4. Replace only missing entries with the common component
  X_init <- X_std
  na_idx <- is.na(X_std)
  X_init[na_idx] <- C_hat[na_idx]
  
  list(
    X_init = X_init,          # T x N, no NAs (for EM)
    r      = r,
    Sigma  = cov_out$Sigma_tilde,
    Lambda = fac_out$Lambda,
    F_hat  = fac_out$F_hat,
    C_hat  = C_hat,
    ER     = ER_out$ER
  )
}

