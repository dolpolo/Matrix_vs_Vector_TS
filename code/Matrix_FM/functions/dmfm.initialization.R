# ==============================================================================
# NUMBER OF FACTORS
# ==============================================================================

mfm.f <- function(Y_std, W, kmax, max_iter = 20) {
  n <- dim(Y_std)[1]
  p1 <- dim(Y_std)[2]
  p2 <- dim(Y_std)[3]
  
  k1_old <- kmax[1]
  k2_old <- kmax[2]
  
  c <- 0.01
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
    C.h <- sqrt(p2) * Q2.h
    
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
    
    ER1 <- eigvals_M1[-length(eigvals_M1)] / (eigvals_M1[-1] + c * delta)
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
    R.h <- sqrt(p1) * Q1.h
    
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
    
    ER2 <- eigvals_M2[-length(eigvals_M2)] / (eigvals_M2[-1] + c * delta)
    k2_new <- which.max(ER2)
    
    cat(sprintf("Iter %d: k1 = %d, k2 = %d\n", iter, k1_new, k2_new))
    
    if (k1_new == k1_old && k2_new == k2_old || iter == max_iter) {
      par(mfrow = c(2, 2))
      
      # ROW FACTORS - Scree plot (bar + line)
      bar_x1 <- barplot(eigvals_M1, col = "skyblue", main = "Scree Plot - Row Factors",
                        xlab = "Factor", ylab = "Eigenvalue")
      lines(bar_x1, eigvals_M1, type = "b", pch = 19, col = "blue")
      
      # ROW FACTORS - Cumulative variance (starts from 0)
      plot(0:length(variance_explained_M1),
           c(0, cumsum(variance_explained_M1)), type = "o", col = "blue",
           main = "Cumulative Variance - Row Factors",
           xlab = "Number of Factors", ylab = "Cumulative Variance")
      
      # COLUMN FACTORS - Scree plot (bar + line)
      bar_x2 <- barplot(eigvals_M2, col = "lightcoral", main = "Scree Plot - Column Factors",
                        xlab = "Factor", ylab = "Eigenvalue")
      lines(bar_x2, eigvals_M2, type = "b", pch = 19, col = "red")
      
      # COLUMN FACTORS - Cumulative variance (starts from 0)
      plot(0:length(variance_explained_M2),
           c(0, cumsum(variance_explained_M2)), type = "o", col = "red",
           main = "Cumulative Variance - Column Factors",
           xlab = "Number of Factors", ylab = "Cumulative Variance")
      
      par(mfrow = c(1, 1))
      break
    }
    
    k1_old <- k1_new
    k2_old <- k2_new
  }
  
  return(c(k1_new, k2_new))
}


# ==============================================================================
# NUMBER OF FACTORS' LAGS
# ==============================================================================

mar_model_selection_auto <- function(Y_std, W, k_hat, max_lag = 10, verbose = TRUE) {
  imputation <- mfm.cl(Y_std, W, k_hat)
  pe <- mfm.pe(imputation$Y_imputed, k_hat)
  Factors <- pe$factor
  
  T <- dim(Factors)[1]
  k1 <- dim(Factors)[2]
  k2 <- dim(Factors)[3]
  k_total <- k1 * k2
  
  # Lista dei fattori vettorizzati
  F_vec_list <- lapply(1:T, function(t) as.vector(Factors[t, , , drop = FALSE]))
  
  bic_values <- numeric(max_lag)
  aic_values <- numeric(max_lag)
  loglik_values <- numeric(max_lag)
  
  for (p in 1:max_lag) {
    if ((T - p) < 2) next
    
    # Y = F_{t}
    Y_mat <- do.call(rbind, F_vec_list[(p+1):T])
    
    # X = [F_{t-1}, ..., F_{t-p}]
    X_mat <- do.call(rbind, lapply((p+1):T, function(t) {
      as.vector(do.call(cbind, F_vec_list[(t - 1):(t - p)]))
    }))
    
    # Verifica dimensioni compatibili
    if (nrow(Y_mat) != nrow(X_mat)) {
      warning(sprintf("Skipping lag p = %d due to dimension mismatch", p))
      next
    }
    
    # Regressione
    B_hat <- ginv(t(X_mat) %*% X_mat) %*% t(X_mat) %*% Y_mat
    Res <- Y_mat - X_mat %*% B_hat
    Sigma_hat <- cov(Res)
    
    n_obs <- T - p
    loglik <- -0.5 * n_obs * (k_total * log(2 * pi) + log(det(Sigma_hat)) + 1)
    loglik_values[p] <- loglik
    
    num_params <- p * k_total^2 + 0.5 * k_total * (k_total + 1)
    bic_values[p] <- -2 * loglik + log(n_obs) * num_params
    aic_values[p] <- -2 * loglik + 2 * num_params
  }
  
  best_p_bic <- which.min(bic_values)
  best_p_aic <- which.min(aic_values)
  
  if (verbose) {
    cat("Automatic MAR(p) - Model Selection\n")
    cat("BIC values:\n"); print(bic_values)
    cat("AIC values:\n"); print(aic_values)
    cat("Best lag (BIC):", best_p_bic, "\n")
    cat("Best lag (AIC):", best_p_aic, "\n")
    
    df <- data.frame(
      Lag = 1:max_lag,
      BIC = bic_values,
      AIC = aic_values
    )
    
    p1 <- ggplot(df, aes(x = Lag, y = BIC)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      geom_vline(xintercept = best_p_bic, linetype = "dashed", color = "blue") +
      scale_x_continuous(breaks = 1:max_lag) +
      ggtitle("BIC vs Lag") +
      theme_minimal()
    
    p2 <- ggplot(df, aes(x = Lag, y = AIC)) +
      geom_line(color = "red") +
      geom_point(color = "red") +
      geom_vline(xintercept = best_p_aic, linetype = "dashed", color = "red") +
      scale_x_continuous(breaks = 1:max_lag) +
      ggtitle("AIC vs Lag") +
      theme_minimal()
  
     print(p1)
   # print(p1 + p2)
  }
  
  return(list(
    BIC = bic_values,
    AIC = aic_values,
    logLik = loglik_values,
    best_lag_BIC = best_p_bic,
    best_lag_AIC = best_p_aic
  ))
}


# ==============================================================================
# IMPUTATION OF NAN
# ==============================================================================

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
# MAR PARAMETERS
# ==============================================================================

# Questa funzione implementa MLE per matrix autoregressive models in Chen et al.
# (JoE, 2021). Questa funzione è attualmente utilizzata per ottenere gli starting
# values per il modello autoregressivo del DMFM nella funzione dmfm.em.R. 
#
# Commento 31/07/23
# La stessa funzione (simile) è implementata nella cartella MAR e funziona. In 
# futuro sarebbe da uniformare se possibile e tenerne una. 

mar.fit <- function(Y){
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  ##### Qui ho dovuto mettere P e Q positive altrimenti non funziona quando k=c(1,1). 4/10/22 #######
  # Obtain projected estimates of A, B, P, Q
  . <- list()
  y <- matrix(apply(Y,1,vec),prod(p), n)
  .$BA    <- var.ols(y)$b
  B.A.t   <- rearrange(.$BA, p[1], p[2])
  B.A.t.u <- matrix(svd(B.A.t)$u[,1],p[1],p[1])
  B.A.t.v <- matrix(svd(B.A.t)$v[,1],p[2],p[2])
  B.A.t.d <- svd(B.A.t)$d[1]
  .$A  <- B.A.t.u*sign(B.A.t.u[1,1])
  .$B  <- B.A.t.d*B.A.t.v*sign(B.A.t.u[1,1])
  .$QP <- var.ols(y)$V
  Q.P.t <- rearrange(.$QP, p[1], p[2])
  Q.P.t.u <- matrix(svd(Q.P.t)$u[,1],p[1],p[1])
  Q.P.t.v <- matrix(svd(Q.P.t)$v[,1],p[2],p[2])
  Q.P.t.d <- svd(Q.P.t)$d[1]
  .$P <- Q.P.t.u*sign(Q.P.t.u[1,1])
  .$Q <- Q.P.t.d*Q.P.t.v*sign(Q.P.t.u[1,1])
  
  # Iterate
  criterion <- T
  eps <- 0.01
  while(criterion){
    
    # Save parameters previous iteration
    .. <- .
    
    # Update A
    a1.tmp <- a2.tmp <- 0
    for (i in 2:n){
      a1.tmp <- a1.tmp + matrix(Y[i,,],p[1],p[2])%*%solve(.$Q)%*%.$B%*%t(matrix(Y[i-1,,],p[1],p[2]))
      a2.tmp <- a2.tmp + matrix(Y[i-1,,],p[1],p[2])%*%t(.$B)%*%solve(.$Q)%*%.$B%*%t(matrix(Y[i-1,,],p[1],p[2]))
    }
    .$A <- a1.tmp%*%solve(a2.tmp)
    .$A <- .$A/norm(.$A,"F")
    
    # Update P
    p.tmp <- 0
    for (i in 2:n){
      p.tmp <- p.tmp + (matrix(Y[i,,],p[1],p[2])-.$A%*%matrix(Y[i-1,,],p[1],p[2])%*%t(.$B))%*%solve(.$Q)%*%t(matrix(Y[i,,],p[1],p[2])-.$A%*%matrix(Y[i-1,,],p[1],p[2])%*%t(.$B))
    }
    .$P <- p.tmp/(p[2]*(n-1))
    .$P <- .$P/norm(.$P, "F")
    
    # Update B
    b1.tmp <- b2.tmp <- 0
    for (i in 2:n){
      b1.tmp <- b1.tmp + t(matrix(Y[i,,],p[1],p[2]))%*%solve(.$P)%*%.$A%*%matrix(Y[i-1,,],p[1],p[2])
      b2.tmp <- b2.tmp + t(matrix(Y[i-1,,],p[1],p[2]))%*%t(.$A)%*%solve(.$P)%*%.$A%*%matrix(Y[i-1,,],p[1],p[2])
    }
    .$B <- b1.tmp%*%solve(b2.tmp)
    
    # Update Q
    q.tmp <- 0
    for (i in 2:n){
      q.tmp <- q.tmp + t(matrix(Y[i,,],p[1],p[2])-.$A%*%matrix(Y[i-1,,],p[1],p[2])%*%t(.$B))%*%solve(.$P)%*%(matrix(Y[i,,],p[1],p[2])-.$A%*%matrix(Y[i-1,,],p[1],p[2])%*%t(.$B))
    }
    .$Q <- q.tmp/(p[1]*(n-1))
    
    # Check convergence
    llk.old <- mvar.llk(.., Y)
    llk.new <- mvar.llk(., Y)
    
    delta <- 2*abs(llk.new-llk.old)/abs(llk.new+llk.old)
    if (delta < eps) {
      criterion <- F
    }
  }
  return(.)
}

rearrange <- function(x,p,q){
  y <- NULL
  for (j in 1:q){
    for(i in 1:q){
      y <- cbind(y, vec(x[(i-1)*p+(1:p),(j-1)*p+(1:p)]))
    }
  }
  return(y)
}

var.ols <- function(y){
  
  # Compute autoregressive parameter
  # p <- dim(y)[1]
  n <- dim(y)[2]
  # u <- d <- matrix(0, p, p)
  # for (i in 2:n){
  #   u <- u + y[,i]%*%t(y[,i-1])
  #   d <- d + y[,i-1]%*%t(y[,i-1])
  # }
  # b <- u%*%solve(d)
  
  b <- (y[,-1,drop=F]%*%t(y[,-n,drop=F]))%*%solve(y[,-n,drop=F]%*%t(y[,-n,drop=F]))
  
  # Compute residuals
  # e <- matrix(0, p, n)
  # for (i in 2:n){
  #   e[,i] <- y[,i] - b%*%y[,i-1]
  # }
  e <- y[,-1] - b%*%y[,-n]
  
  # Compute vcov
  # V <- matrix(0, p, p)
  # for (i in 1:n){
  #   V <- V + e[,i]%*%t(e[,i])/n
  # }
  V <- (e%*%t(e))/(n-1)
  
  out <- list(b=b, V=V)
  
  return(out)
}

mvar.llk <- function(., f){
  
  n <- dim(f)[1]
  k <- dim(f)[-1]
  
  llk <- rep(0, n)
  for (i in 2:n){
    llk[i] <- tr(solve(.$P)%*%(f[i,,]-.$A%*%f[i-1,,]%*%t(.$B))%*%solve(.$Q)%*%t(f[i,,]-.$A%*%f[i-1,,]%*%t(.$B))) 
  }
  
  out <- - (n-1)*k[1]*log(det(.$Q)) - (n-1)*k[2]*log(det(.$P)) - sum(llk)
  
  return(out)
  
}

# ==============================================================================
# INITIALIZE INPUTS
# ==============================================================================
# dmfm.na.sv----
dmfm.na.sv <- function(Y, X, k, W, t){

  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  # Obtain starting values R,C
  . <- list()
  par <- mfm.pe(X, k)
  R.pc <- par$row.load
  C.pc <- par$col.load
  Y.pc <- par$fitted
  F.pc <- par$factor
  if (t=="fv"){
    .$CR <- as.matrix(C.pc)%x%as.matrix(R.pc)
  }else{
    .$R <- as.matrix(R.pc)
    .$C <- as.matrix(C.pc)
  }
  
  # Obtain estimates of H and K
  k.tmp <- 0
  h.tmp <- 0
  for (i in 1:n){
    k.tmp <- k.tmp + t(Y[i,,]-Y.pc[i,,])%*%(Y[i,,]-Y.pc[i,,])
  }
  K.tmp <-  diag(diag(k.tmp/(n*p[1])))
  for (i in 1:n){
    h.tmp <- h.tmp + (Y[i,,]-Y.pc[i,,])%*%solve(K.tmp)%*%t(Y[i,,]-Y.pc[i,,])
  }
  if (t=="fv"){
    .$KH <- diag(diag(k.tmp/(n*p[1])))%x%diag(diag(h.tmp/(n*p[2])))
    # .$KH <- diag(apply(apply(Y-Y.pc,1,vec)^2, 1, mean))
  }else{
    .$K <- K.tmp
    .$H <- diag(diag(h.tmp/(n*p[2])))
  }
  
  # Obtain estimates of A, B, P, Q
  est <- mar.fit(F.pc)
  if (t=="fm"){
    .$A <- est$A
    .$B <- est$B
    .$P <- est$P
    .$Q <- est$Q
  } else {
    .$BA <- est$BA
    .$QP <- est$QP
  }
  
  # Kalman 
  O <- lapply(split(W,seq(nrow(W))), function(x){which(vec(x)==1)})
  if (t=="fm"){
    kout  <- kalman.na(list(t=.$B%x%.$A,
                            l=.$C%x%.$R,
                            q=.$Q%x%.$P,
                            h=.$K%x%.$H),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  } else if (t=="fv"){
    kout  <- kalman.na(list(t=.$BA,
                            l=.$CR,
                            q=.$QP,
                            h=.$KH),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }else{
    kout  <- kalman.na(list(t=.$BA,
                            l=.$C%x%.$R,
                            q=.$QP,
                            h=.$K%x%.$H),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }
  
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
  .$Y  <- array(0, c(n, p[1], p[2]))
  for (i in 1:n){
    .$ff[i,,] <- matrix(kout$ff[,i],k[1],k[2])
    .$fs[i,,] <- matrix(kout$fs[,i],k[1],k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    if (t=="fv"){
      .$Y[i,,] <- matrix(.$CR%*%kout$fs[,i], p[1], p[2])
    }else{
      .$Y[i,,] <- .$R%*%.$fs[i,,]%*%t(.$C)
    }
  }
  
  return(.)
  
}

# dmfm.na.2.sv----
# This function is proposed by Matteo and produce starting values obtained
# on a subset of the original matrix removing the NA's
dmfm.na.2.sv <- function(X, k, W, t){
  
  n <- dim(X)[1]
  p <- dim(X)[-1]
  
  # Remove NA's
  complete_time <- which(apply(X, 1, function(slice) !anyNA(slice)))
  Y <- X[complete_time,,]
  n <- dim(Y)[1]
  
  # Obtain starting values R,C
  . <- list()
  par <- mfm.pe(Y, k)
  R.pc <- par$row.load
  C.pc <- par$col.load
  Y.pc <- par$fitted
  F.pc <- par$factor
  if (t=="fv"){
    .$CR <- as.matrix(C.pc)%x%as.matrix(R.pc)
  }else{
    .$R <- as.matrix(R.pc)
    .$C <- as.matrix(C.pc)
  }
  
  # Obtain estimates of H and K
  k.tmp <- 0
  h.tmp <- 0
  for (i in 1:n){
    k.tmp <- k.tmp + t(Y[i,,]-Y.pc[i,,])%*%(Y[i,,]-Y.pc[i,,])
  }
  K.tmp <-  diag(diag(k.tmp/(n*p[1])))
  for (i in 1:n){
    h.tmp <- h.tmp + (Y[i,,]-Y.pc[i,,])%*%solve(K.tmp)%*%t(Y[i,,]-Y.pc[i,,])
  }
  if (t=="fv"){
    .$KH <- diag(diag(k.tmp/(n*p[1])))%x%diag(diag(h.tmp/(n*p[2])))
    # .$KH <- diag(apply(apply(Y-Y.pc,1,vec)^2, 1, mean))
  }else{
    .$K <- K.tmp
    .$H <- diag(diag(h.tmp/(n*p[2])))
  }
  
  # Obtain estimates of A, B, P, Q
  est <- mar.fit(F.pc)
  if (t=="fm"){
    .$A <- est$A
    .$B <- est$B
    .$P <- est$P
    .$Q <- est$Q
  } else {
    .$BA <- est$BA
    .$QP <- est$QP
  }
  
  # Kalman 
  O <- lapply(split(W,seq(nrow(W))), function(x){which(vec(x)==1)})
  if (t=="fm"){
    kout  <- kalman.na(list(t=.$B%x%.$A,
                            l=.$C%x%.$R,
                            q=.$Q%x%.$P,
                            h=.$K%x%.$H),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  } else if (t=="fv"){
    kout  <- kalman.na(list(t=.$BA,
                            l=.$CR,
                            q=.$QP,
                            h=.$KH),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }else{
    kout  <- kalman.na(list(t=.$BA,
                            l=.$C%x%.$R,
                            q=.$QP,
                            h=.$K%x%.$H),
                       apply(X,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }
  
  n <- dim(X)[1]
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
  .$Y  <- array(0, c(n, p[1], p[2]))
  for (i in 1:n){
    .$ff[i,,] <- matrix(kout$ff[,i],k[1],k[2])
    .$fs[i,,] <- matrix(kout$fs[,i],k[1],k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    if (t=="fv"){
      .$Y[i,,] <- matrix(.$CR%*%kout$fs[,i], p[1], p[2])
    }else{
      .$Y[i,,] <- .$R%*%.$fs[i,,]%*%t(.$C)
    }
  }
  
  return(.)
  
}


# X <- Y_std
# k <- k_hat
# t = "dmfm"

dmfm.na.2.sv.idio <- function(X, k, W, t,
                              kappa   = 1e-4,
                              max_als = 50,
                              tol_als = 1e-6,
                              clip_phi = 0.98,
                              eps     = 1e-8,
                              n_scale = 5){
  
  stopifnot(length(dim(X)) == 3)
  
  n <- dim(X)[1]
  p <- dim(X)[-1]    # p = c(p1,p2)
  p1 <- p[1]; p2 <- p[2]
  pp <- p1 * p2
  rr <- prod(k)
  
  # ------------------------------------------------------------
  # 1) Remove NA's (complete slices) for PC-based initialization
  # ------------------------------------------------------------
  complete_time <- which(apply(X, 1, function(slice) !anyNA(slice)))
  if (length(complete_time) < 3) {
    stop("Too few complete slices to initialize (need at least 3).")
  }
  Y  <- X[complete_time, , , drop = FALSE]
  nC <- dim(Y)[1]
  
  # ------------------------------------------------------------
  # 2) Starting values R,C and factors from MFM on complete slices
  # ------------------------------------------------------------
  . <- list()
  par  <- mfm.pe(Y, k)
  R.pc <- par$row.load
  C.pc <- par$col.load
  Y.pc <- par$fitted
  F.pc <- par$factor
  
  if (t == "fv"){
    .$CR <- as.matrix(C.pc) %x% as.matrix(R.pc)     # (pp) x (rr)
  } else {
    .$R  <- as.matrix(R.pc)
    .$C  <- as.matrix(C.pc)
  }
  
  # ------------------------------------------------------------
  # 3) Measurement nugget: Xi_t ~ N(0, kappa I_{pp})
  #    (ONLY in measurement covariance, not in the state)
  # ------------------------------------------------------------
  .$kappa  <- kappa
  .$h_meas <- diag(kappa, pp)
  
  # ------------------------------------------------------------
  # 4) Starting values for persistent idiosyncratic block
  #    Proxy: E_hat ≈ \tilde E on complete times (kappa small)
  # ------------------------------------------------------------
  E.hat <- Y - Y.pc
  .$Ehat <- E.hat
  
  # Entrywise sufficient statistics (proxy)
  S10 <- matrix(0, p1, p2)
  S00 <- matrix(0, p1, p2)
  for (tt in 2:nC){
    Et   <- E.hat[tt,,]
    Etm1 <- E.hat[tt-1,,]
    S10 <- S10 + (Et * Etm1)
    S00 <- S00 + (Etm1 * Etm1)
  }
  S00 <- pmax(S00, eps)
  
  # ALS init
  a  <- rep(0, p1)
  b  <- rep(1, p2)
  hE <- rep(1, p1)
  kE <- rep(1, p2)
  
  # ------------------------------------------------------------
  # 5) ALS for (a,b): phi_ij = a_i b_j
  # ------------------------------------------------------------
  for (iter in 1:max_als){
    a_old <- a; b_old <- b
    
    for (i in 1:p1){
      num <- sum((b / pmax(kE, eps)) * S10[i,])
      den <- sum(((b^2) / pmax(kE, eps)) * S00[i,])
      a[i] <- num / pmax(den, eps)
    }
    for (j in 1:p2){
      num <- sum((a / pmax(hE, eps)) * S10[,j])
      den <- sum(((a^2) / pmax(hE, eps)) * S00[,j])
      b[j] <- num / pmax(den, eps)
    }
    
    # Stability clamp: max_ij |a_i b_j| <= clip_phi
    phi_max <- max(abs(outer(a, b, "*")))
    if (is.finite(phi_max) && phi_max > clip_phi && phi_max > 0){
      s <- clip_phi / phi_max
      a <- a * sqrt(s)
      b <- b * sqrt(s)
    }
    
    if (max(abs(a-a_old), abs(b-b_old)) < tol_als) break
  }
  
  # ------------------------------------------------------------
  # 6) Innovation variances (hE,kE) from RSS, alternating scaling
  # ------------------------------------------------------------
  phi <- outer(a, b, "*")
  .$phi <- phi
  
  RSS <- matrix(0, p1, p2)
  for (tt in 2:nC){
    Et   <- E.hat[tt,,]
    Etm1 <- E.hat[tt-1,,]
    res  <- Et - (phi * Etm1)
    RSS  <- RSS + (res * res)
  }
  RSS <- pmax(RSS, eps)
  
  for (iter in 1:n_scale){
    for (i in 1:p1){
      hE[i] <- (1 / ((nC-1) * p2)) * sum(RSS[i,] / pmax(kE, eps))
    }
    hE <- pmax(hE, eps)
    for (j in 1:p2){
      kE[j] <- (1 / ((nC-1) * p1)) * sum(RSS[,j] / pmax(hE, eps))
    }
    kE <- pmax(kE, eps)
  }
  
  # Normalize mean(hE)=1 (scale identification)
  mh <- mean(hE)
  if (is.finite(mh) && mh > 0){
    hE <- hE / mh
    kE <- kE * mh
  }
  
  .$A_E <- diag(a,  p1, p1)
  .$B_E <- diag(b,  p2, p2)
  .$H_E <- diag(hE, p1, p1)
  .$K_E <- diag(kE, p2, p2)
  
  # ------------------------------------------------------------
  # 7) Starting values for factor dynamics (A,B,P,Q) from F.pc
  # ------------------------------------------------------------
  est <- mar.fit(F.pc)
  if (t == "fm"){
    .$A <- est$A; .$B <- est$B; .$P <- est$P; .$Q <- est$Q
  } else {
    .$BA <- est$BA; .$QP <- est$QP
  }
  
  # ============================================================
  # 8) Kalman on AUGMENTED state: s_t = [ f_t ; e_t ]
  #    x_t = [Lf I] s_t + Xi_t,   Xi_t ~ N(0, kappa I)
  # ============================================================
  m_state <- rr + pp
  
  # Missing selector
  O <- lapply(split(W, seq(nrow(W))), function(x){ which(vec(x) == 1) })
  
  # Factor block
  if (t == "fm") {
    Tf <- .$B %x% .$A
    Lf <- .$C %x% .$R
    Qf <- .$Q %x% .$P
  } else if (t == "fv") {
    Tf <- .$BA
    Lf <- .$CR
    Qf <- .$QP
  } else {
    Tf <- .$BA
    Lf <- .$C %x% .$R
    Qf <- .$QP
  }
  
  # Persistent idiosyncratic block
  Te <- .$B_E %x% .$A_E
  Qe <- (diag(diag(.$K_E)) %x% diag(diag(.$H_E)))  # dg(K_E) ⊗ dg(H_E)
  
  # Augmented transition and innovations
  Ts <- rbind(
    cbind(Tf, matrix(0, rr, pp)),
    cbind(matrix(0, pp, rr), Te)
  )
  
  Qs <- rbind(
    cbind(Qf, matrix(0, rr, pp)),
    cbind(matrix(0, pp, rr), Qe)
  )
  
  # Augmented measurement matrix
  Ls <- cbind(Lf, diag(pp))
  
  # Measurement covariance (nugget only)
  Hx <- .$h_meas
  
  kout <- kalman.na(
    list(t = Ts, l = Ls, q = Qs, h = Hx),
    apply(X, 1, vec),
    m_state,
    list(f = rep(0, m_state),
         P = diag(m_state)),
    O
  )
  
  # ============================================================
  # 9) Save outputs: factors + persistent idiosyncratic state
  # ============================================================
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, m_state, m_state))
  
  .$Y      <- array(0, c(n, p1, p2))   # common component
  .$Etilde <- array(0, c(n, p1, p2))   # persistent idiosyncratic state
  
  for (i in 1:n){
    
    s_ff <- kout$ff[, i]
    s_fs <- kout$fs[, i]
    
    f_ff <- s_ff[1:rr]
    f_fs <- s_fs[1:rr]
    e_fs <- s_fs[(rr+1):(rr+pp)]
    
    .$ff[i,,] <- matrix(f_ff, k[1], k[2])
    .$fs[i,,] <- matrix(f_fs, k[1], k[2])
    
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    
    .$Etilde[i,,] <- matrix(e_fs, p1, p2)
    
    if (t == "fv"){
      .$Y[i,,] <- matrix(Lf %*% f_fs, p1, p2)
    } else {
      .$Y[i,,] <- .$R %*% .$fs[i,,] %*% t(.$C)
    }
  }
  
  return(.)
}
# ----

dmfm.na.2.sv.vec <- function(X, k, W, t = "dmfm") {
  n <- dim(X)[1]
  p <- dim(X)[-1]
  
  # Remove NA's
  complete_time <- which(apply(X, 1, function(slice) !anyNA(slice)))
  Y <- X[complete_time, , , drop = FALSE]
  n_clean <- dim(Y)[1]
  
  . <- list()
  
  # === 1. Obtain Projected Estimates ===
  par <- mfm.pe.vec(Y, k)
  R.pc <- par$row.load  # Should be 1 x 1 identity
  C.pc <- par$col.load
  Y.pc <- par$fitted
  F.pc <- par$factor
  
  if (t == "fv") {
    .$CR <- as.matrix(C.pc)  # R is identity
  } else {
    .$R <- as.matrix(R.pc)
    .$C <- as.matrix(C.pc)
  }
  
  # === 2. Estimate H and K ===
  k.tmp <- 0
  h.tmp <- 0
  for (i in 1:n_clean) {
    res <- Y[i, , ] - Y.pc[i, , ]
    k.tmp <- k.tmp + t(res) %*% res
  }
  K.tmp <- matrix(k.tmp / (n_clean * p[1]), nrow = 1, ncol = 1)  # scalare se p1 = 1
  for (i in 1:n_clean) {
    res <- Y[i, , ] - Y.pc[i, , ]  # 1 × p2
    h.tmp <- h.tmp + (1 / K.tmp[1,1]) * res %*% t(res)
  }
  if (t == "fv") {
    .$KH <- diag(diag(k.tmp / (n_clean * p[1])))
  } else {
    .$K <- K.tmp
    .$H <- diag(diag(h.tmp / (n_clean * p[2])))
  }
  
  # === 3. Estimate Dynamics ===
  est <- mar.fit(F.pc)
  if (t == "fm") {
    .$A <- est$A
    .$B <- est$B
    .$P <- est$P
    .$Q <- est$Q
  } else {
    .$BA <- est$BA
    .$QP <- est$QP
  }
  
  # === 4. Kalman Filter & Smoother ===
  O <- lapply(split(W, seq(nrow(W))), function(x) which(vec(x) == 1))
  y_vec <- apply(X, 1, vec)
  init <- list(f = rep(0, prod(k)), P = diag(prod(k)))
  
  if (t == "fm") {
    kout <- kalman.na(list(
      t = .$B %x% .$A,
      l = .$C %x% .$R,
      q = .$Q %x% .$P,
      h = .$K %x% .$H
    ), y_vec, prod(k), init, O)
  } else if (t == "fv") {
    kout <- kalman.na(list(
      t = .$BA,
      l = .$CR,
      q = .$QP,
      h = .$KH
    ), y_vec, prod(k), init, O)
  } else {
    kout <- kalman.na(list(
      t = .$BA,
      l = .$C %x% .$R,
      q = .$QP,
      h = .$K %x% .$H
    ), y_vec, prod(k), init, O)
  }
  
  # === 5. Populate Output ===
  n <- dim(X)[1]
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
  .$Y <- array(0, c(n, p[1], p[2]))
  
  for (i in 1:n) {
    .$ff[i,,] <- matrix(kout$ff[,i], k[1], k[2])
    .$fs[i,,] <- matrix(kout$fs[,i], k[1], k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    
    if (t == "fv") {
      .$Y[i,,] <- matrix(.$CR %*% kout$fs[,i], p[1], p[2])
    } else {
      .$Y[i,,] <- .$R %*% .$fs[i,,] %*% t(.$C)
    }
  }
  
  return(.)
}

