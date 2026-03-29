# ==============================================================================
# IMPUTATION OF NAN using Cen and Lam
# ==============================================================================

mfm.cl <- function(Y, W, r) {
  
  # Dimensions
  n  <- dim(Y)[1]
  d1 <- dim(Y)[2]
  d2 <- dim(Y)[3]
  r1 <- r[1]
  r2 <- r[2]
  
  # Initialize covariance matrices
  S_R <- matrix(0, d1, d1)
  S_C <- matrix(0, d2, d2)
  
  # Counters for valid entries in covariance sums
  count_R <- matrix(0, d1, d1)
  count_C <- matrix(0, d2, d2)
  
  # Compute mode-wise sample covariance matrices (handling missing values properly)
  for (t in 1:n) {
    W_t <- W[t,,]
    
    for (i in 1:d1) {
      for (j in i:d1) {
        valid_cols <- which(W_t[i,] * W_t[j,] == 1)  # Columns where both rows i, j are observed
        if (length(valid_cols) > 0) {
          S_R[i, j] <- S_R[i, j] + sum(Y[t, i, valid_cols] * Y[t, j, valid_cols])
          S_R[j, i] <- S_R[i, j]  # Symmetric
          count_R[i, j] <- count_R[i, j] + length(valid_cols)
          count_R[j, i] <- count_R[i, j]  # Symmetric
        }
      }
    }
    
    for (i in 1:d2) {
      for (j in i:d2) {
        valid_rows <- which(W_t[,i] * W_t[,j] == 1)  # Rows where both columns i, j are observed
        if (length(valid_rows) > 0) {
          S_C[i, j] <- S_C[i, j] + sum(Y[t, valid_rows, i] * Y[t, valid_rows, j])
          S_C[j, i] <- S_C[i, j]  # Symmetric
          count_C[i, j] <- count_C[i, j] + length(valid_rows)
          count_C[j, i] <- count_C[i, j]  # Symmetric
        }
      }
    }
  }
  
  # Normalize by valid observation counts to avoid division by zero
  S_R[count_R > 0] <- S_R[count_R > 0] / count_R[count_R > 0]
  S_C[count_C > 0] <- S_C[count_C > 0] / count_C[count_C > 0]
  
  # PCA: Extract first r1 and r2 eigenvectors as factor loadings
  eig_R <- eigen(S_R, symmetric = TRUE)
  R_hat <- eig_R$vectors[, 1:r1, drop = FALSE]
  
  eig_C <- eigen(S_C, symmetric = TRUE)
  C_hat <- eig_C$vectors[, 1:r2, drop = FALSE]
  
  # Estimate factor matrix F_t
  F_hat <- array(0, dim = c(n, r1, r2))
  
  for (t in 1:n) {
    W_t_vec <- as.vector(W[t,,])
    Y_t_vec <- as.vector(Y[t,,])
    
    if (sum(W_t_vec) > 0) {
      Q_tensor <- kronecker(C_hat, R_hat)
      Q_obs <- Q_tensor[W_t_vec == 1, , drop = FALSE]
      Y_obs <- Y_t_vec[W_t_vec == 1]
      
      if (nrow(Q_obs) >= max(r1 * r2, 1)) {
        F_hat[t,,] <- matrix(solve(t(Q_obs) %*% Q_obs) %*% (t(Q_obs) %*% Y_obs), nrow = r1, ncol = r2)
      }
    }
  }
  
  # Final imputation
  Y_hat <- Y_imputed <- array(0, dim = c(n, d1, d2))
  dimnames(Y_hat) <- dimnames(Y)
  dimnames(Y_imputed) <- dimnames(Y)
  for (t in 1:n) {
    Y_hat[t,,] <- R_hat %*% F_hat[t,,] %*% t(C_hat)
    Y_imputed[t,,] <- R_hat %*% F_hat[t,,] %*% t(C_hat)
    Y_imputed[t,,][W[t,,] == 1] <- Y[t,,][W[t,,] == 1]  # Keep observed values unchanged
  }
  
  return(list(R = R_hat, C = C_hat, f = F_hat, Y_hat = Y_hat, Y_imputed = Y_imputed))
}

# imputation <- mfm.cl(X_std, W_x, params$Kmax)

# ==============================================================================
# Initialization 
# CL covariance + iter_ER + CL factors → replace NAs with common component
# ==============================================================================

init_CL_Yu <- function(Y_std, W, kmax = c(1,10),
                       max_iter_dim = 20, do_plot = FALSE) {
  # 1) Dimension selection (Yu-style) using your mfm.f
  k_hat <- mfm.f(Y_std, W, kmax = kmax, max_iter = max_iter_dim, do_plot = do_plot)
  r1 <- k_hat[1]; r2 <- k_hat[2]
  
  # 2) Final Cen & Lam estimation + imputation given (r1,r2)
  out_cl <- mfm.cl(Y_std, W, r = c(r1, r2))
  
  # 3) Replace ONLY missing entries (recommended)
  Y_init <- Y_std
  Y_init[W == 0] <- out_cl$Y_hat[W == 0]
  
  list(
    r      = c(r1, r2),
    Y_init = Y_init,          # no missing (if W is 0/1 mask of obs)
    R_hat  = out_cl$R,
    C_hat  = out_cl$C,
    F_hat  = out_cl$f,
    Y_hat  = out_cl$Y_hat
  )
}
