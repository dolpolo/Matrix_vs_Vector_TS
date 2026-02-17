# ==============================================================================
# NUMBER OF FACTORS
# ==============================================================================

# Y_std <- X_std
# W     <- W_x
# kmax <- c(1,10)
mfm.f <- function(Y_std, W, kmax, max_iter = 20, do_plot = TRUE) {
  n  <- dim(Y_std)[1]
  p1 <- dim(Y_std)[2]
  p2 <- dim(Y_std)[3]
  
  k1_old <- kmax[1]
  k2_old <- kmax[2]
  
  c     <- 0.01
  delta <- max(1 / sqrt(n * p2), 1 / sqrt(n * p1), 1 / p1)
  
  for (iter in 1:max_iter) {
    k_old <- c(k1_old, k2_old)
    
    imputation <- mfm.cl(Y_std, W, k_old)
    X <- imputation$Y_imputed
    
    # Step 2: Column loadings
    M2.h <- matrix(0, p2, p2)
    for (j in 1:p1) {
      for (t in 1:n) {
        row_vec <- matrix(X[t, j, ], nrow = p2, ncol = 1)
        M2.h <- M2.h + row_vec %*% t(row_vec)
      }
    }
    M2.h <- M2.h / (n * p1 * p2)
    Q2.h <- eigen(M2.h)$vectors[, 1:k2_old]
    C.h  <- sqrt(p2) * Q2.h
    
    # Step 3: Project to Y.h (row-reduced)
    Y.h <- array(0, c(n, p1, k2_old))
    for (t in 1:n) {
      Y.h[t, , ] <- (1 / p2) * X[t, , ] %*% C.h
    }
    
    M1.t <- matrix(0, p1, p1)
    for (t in 1:n) {
      M1.t <- M1.t + Y.h[t, , ] %*% t(Y.h[t, , ])
    }
    M1.t <- M1.t / (n * p1)
    eigvals_M1 <- eigen(M1.t, symmetric = TRUE, only.values = TRUE)$values
    variance_explained_M1 <- eigvals_M1 / sum(eigvals_M1)
    
    ER1    <- eigvals_M1[-length(eigvals_M1)] / (eigvals_M1[-1] + c * delta)
    k1_new <- which.max(ER1)
    
    # Step 4: Row loadings
    M1.h <- matrix(0, p1, p1)
    for (j in 1:p2) {
      for (t in 1:n) {
        col_vec <- matrix(X[t, , j], nrow = p1, ncol = 1)
        M1.h <- M1.h + col_vec %*% t(col_vec)
      }
    }
    M1.h <- M1.h / (n * p1 * p2)
    Q1.h <- eigen(M1.h)$vectors[, 1:k1_new]
    R.h  <- sqrt(p1) * Q1.h
    
    # Step 5: Project to Z.h (column-reduced)
    Z.h <- array(0, c(n, p2, k1_new))
    for (t in 1:n) {
      Z.h[t, , ] <- (1 / p1) * t(X[t, , ]) %*% R.h
    }
    
    M2.t <- matrix(0, p2, p2)
    for (t in 1:n) {
      M2.t <- M2.t + Z.h[t, , ] %*% t(Z.h[t, , ])
    }
    M2.t <- M2.t / (n * p2)
    eigvals_M2 <- eigen(M2.t, symmetric = TRUE, only.values = TRUE)$values
    variance_explained_M2 <- eigvals_M2 / sum(eigvals_M2)
    
    ER2    <- eigvals_M2[-length(eigvals_M2)] / (eigvals_M2[-1] + c * delta)
    k2_new <- which.max(ER2)
    
    cat(sprintf("Iter %d: k1 = %d, k2 = %d\n", iter, k1_new, k2_new))
    
    if (k1_new == k1_old && k2_new == k2_old || iter == max_iter) {
      
      if (isTRUE(do_plot)) {
        par(mfrow = c(2, 2))
        
        # ROW FACTORS - Scree plot (bar + line)
        bar_x1 <- barplot(
          eigvals_M1,
          col  = "skyblue",
          main = "Scree Plot - Row Factors",
          xlab = "Factor",
          ylab = "Eigenvalue"
        )
        lines(bar_x1, eigvals_M1, type = "b", pch = 19, col = "blue")
        
        # ROW FACTORS - Cumulative variance (starts from 0)
        plot(
          0:length(variance_explained_M1),
          c(0, cumsum(variance_explained_M1)),
          type = "o",
          col  = "blue",
          main = "Cumulative Variance - Row Factors",
          xlab = "Number of Factors",
          ylab = "Cumulative Variance"
        )
        
        # COLUMN FACTORS - Scree plot (bar + line)
        bar_x2 <- barplot(
          eigvals_M2,
          col  = "lightcoral",
          main = "Scree Plot - Column Factors",
          xlab = "Factor",
          ylab = "Eigenvalue"
        )
        lines(bar_x2, eigvals_M2, type = "b", pch = 19, col = "red")
        
        # COLUMN FACTORS - Cumulative variance (starts from 0)
        plot(
          0:length(variance_explained_M2),
          c(0, cumsum(variance_explained_M2)),
          type = "o",
          col  = "red",
          main = "Cumulative Variance - Column Factors",
          xlab = "Number of Factors",
          ylab = "Cumulative Variance"
        )
        
        par(mfrow = c(1, 1))
      }
      
      break
    }
    
    k1_old <- k1_new
    k2_old <- k2_new
  }
  
  return(c(k1_new, k2_new))
}


# ==============================================================================
# IMPUTATION OF NAN
# ==============================================================================


# Y     <- X_std
# W     <- W_x
# r <- c(1,10)

# Can and Lam
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
# PROJECTED ESTIMATES
# ==============================================================================

# Projected estimator for matrix factor models developed in Yu et al. 
# (JoE, 2021)

mfm.pe <- function(X, k){
  
  n <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(X)[3]
  
  # Initial Estimators
  M1.h <- matrix(0, p1, p1)
  M2.h <- matrix(0, p2, p2)
  
  # covarianza tra righe (paesi) lungo tutte le variabili p2
  for (j in 1:p2){
    M1.h <- M1.h + t(X[,,j]) %*% X[,,j]
  }
  # covarianza tra colonne (variabili) lungo tutti i paesi in p1
  for (j in 1:p1){
    M2.h <- M2.h + t(X[,j,]) %*% X[,j,]
  }
  
  M1.h <- M1.h / (n * p1 * p2)
  M2.h <- M2.h / (n * p1 * p2)
  
  Q1.h <- as.matrix(eigen(M1.h)$vectors[, 1:k[1]])
  Q2.h <- as.matrix(eigen(M2.h)$vectors[, 1:k[2]])
  
  R.h <- sqrt(p1) * Q1.h
  C.h <- sqrt(p2) * Q2.h
  
  Y.h <- array(0, c(n, p1, k[2]))
  Z.h <- array(0, c(n, p2, k[1]))
  for (t in 1:n){
    Y.h[t,,] <- p2^(-1) * X[t,,] %*% C.h
    Z.h[t,,] <- p1^(-1) * t(X[t,,]) %*% R.h
  }
  
  # Recursive Estimation (Algorithm 1)
  M1.t <- matrix(0, p1, p1)
  M2.t <- matrix(0, p2, p2)
  for (t in 1:n){
    M1.t <- M1.t + matrix(Y.h[t,,], p1, k[2]) %*% t(matrix(Y.h[t,,], p1, k[2]))
    M2.t <- M2.t + matrix(Z.h[t,,], p2, k[1]) %*% t(matrix(Z.h[t,,], p2, k[1]))
  }
  M1.t <- M1.t / (n * p1)
  M2.t <- M2.t / (n * p2)
  
  Q1.t <- as.matrix(eigen(M1.t)$vectors[, 1:k[1]])
  Q2.t <- as.matrix(eigen(M2.t)$vectors[, 1:k[2]])
  
  R.t <- sqrt(p1) * Q1.t
  C.t <- sqrt(p2) * Q2.t
  
  # Aggiunta dei nomi se presenti in X
  if (!is.null(dimnames(X))) {
    dn <- dimnames(X)
    if (!is.null(dn[[2]])) {
      rownames(R.t) <- dn[[2]]
    }
    if (!is.null(dn[[3]])) {
      rownames(C.t) <- dn[[3]]
    }
  }
  colnames(R.t) <- paste0("F_row_", 1:k[1])
  colnames(C.t) <- paste0("F_col_", 1:k[2])
  
  # Fattori e stima ricostruita
  F.h <- array(0, c(n, k[1], k[2]))
  S.h <- array(0, c(n, p1, p2))
  for (t in 1:n){
    F.h[t,,] <- (p1 * p2)^(-1) * t(R.t) %*% X[t,,] %*% C.t
    S.h[t,,] <- (p1 * p2)^(-1) * R.t %*% t(R.t) %*% X[t,,] %*% C.t %*% t(C.t)
  }
  
  return(list(
    row.load = R.t,
    col.load = C.t,
    factor = F.h,
    fitted = S.h
  ))
}




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
# 2. Eigenvalue Ratio estimator (Ahn & Horenstein, 2013)
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

