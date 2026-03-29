# ==============================================================================
# 1. Eigenvalue Ratio estimator (Ahn & Horenstein, 2013)
#    Applied to a covariance matrix Sigma
# ==============================================================================

select_num_factors_ER <- function(Sigma, Kmax = 15) {
  # Sigma: N x N covariance matrix
  
  eigvals <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  eigvals <- sort(eigvals, decreasing = TRUE)
  
  max_k <- min(Kmax, length(eigvals) - 1)
  if (max_k < 1) {
    stop("Not enough eigenvalues to compute ER.")
  }
  
  ER_vec <- eigvals[1:max_k] / eigvals[2:(max_k + 1)]
  r      <- which.max(ER_vec)
  
  list(
    r     = r,
    ER    = ER_vec,
    eig   = eigvals
  )
}