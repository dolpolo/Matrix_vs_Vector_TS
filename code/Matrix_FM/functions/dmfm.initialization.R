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

# ----

# Modified version of mfm.f to select k2 (column factors) only.
# - Assumes k1 = 1 (row factor) is fixed.
# - Performs eigen decomposition only on column covariances.
# - Estimates column loadings and projects onto them to compute reduced data.
# - Computes eigenvalue ratios to choose optimal k2.

dfm.f.vec <- function(Y_std, W, kmax, max_iter = 100) {
  # Assumiamo che Y_std abbia dimensione T x 1 x p2
  n <- dim(Y_std)[1]
  p2 <- dim(Y_std)[3]
  k2_old <- kmax[2]
  
  c <- 0.01
  delta <- max(1 / sqrt(n * p2), 1 / p2)
  
  for (iter in 1:max_iter) {
    imputation <- mfm.cl.vec(Y_std, W, c(1, k2_old))
    X <- imputation$Y_imputed
    
    # Step: Covarianza sulle colonne
    M2.h <- matrix(0, p2, p2)
    for (t in 1:n) {
      row_vec <- as.matrix(X[t, 1, ])
      M2.h <- M2.h + row_vec %*% t(row_vec)
    }
    M2.h <- M2.h / (n * p2)
    
    eig <- eigen(M2.h, symmetric = TRUE)
    eigvals_M2 <- eig$values
    variance_explained_M2 <- eigvals_M2 / sum(eigvals_M2)
    
    ER2 <- eigvals_M2[-length(eigvals_M2)] / (eigvals_M2[-1] + c * delta)
    if (length(ER2) == 0) {
      warning("ER2 vuoto: k2 impostato a 1")
      k2_new <- 1
    } else {
      k2_new <- which.max(ER2)
    }
    
    cat(sprintf("Iter %d: k2 = %d\n", iter, k2_new))
    
    if (k2_new == k2_old || iter == max_iter) {
      par(mfrow = c(1, 2))
      
      # Scree plot
      bar_x2 <- barplot(eigvals_M2, col = "lightcoral", main = "Scree Plot - Column Factors",
                        xlab = "Factor", ylab = "Eigenvalue")
      lines(bar_x2, eigvals_M2, type = "b", pch = 19, col = "red")
      
      # Cumulative variance
      plot(0:length(variance_explained_M2),
           c(0, cumsum(variance_explained_M2)), type = "o", col = "red",
           main = "Cumulative Variance - Column Factors",
           xlab = "Number of Factors", ylab = "Cumulative Variance")
      
      par(mfrow = c(1, 1))
      break
    }
    
    k2_old <- k2_new
  }
  
  return(c(1, k2_new))
}



# ==============================================================================
# NUMBER OF FACTORS' LAGS
# ==============================================================================

#  somma i fattori ai lag precedenti e usa una sola matrice B per rappresentare 
# l’effetto combinato dei p ritardi.

factorize_bic_matrix_loglik_mar <- function(Y_std, W, k_hat, max_lag = 10, verbose = TRUE) {
  imputation <- mfm.cl(Y_std, W, k_hat)
  pe <- mfm.pe(imputation$Y_imputed, k_hat)
  Factors <- pe$factor
  
  T <- dim(Factors)[1]
  k1 <- dim(Factors)[2]
  k2 <- dim(Factors)[3]
  F_list <- lapply(1:T, function(t) Factors[t,,])
  
  bic_values <- numeric(max_lag)
  aic_values <- numeric(max_lag)
  loglik_values <- numeric(max_lag)
  
  for (p in 1:max_lag) {
    n_obs <- T - p
    A_list <- list()
    B_list <- list()
    P_sum <- matrix(0, k1, k1)
    Q_sum <- matrix(0, k2, k2)
    
    # MAR(1) stimato tra F_t e F_{t-p} (solo un lag alla volta)
    for (i in (p+1):T) {
      Y <- F_list[[i]]
      X <- F_list[[i - p]]
      
      svd_tmp <- svd(Y %*% X)
      A_hat <- svd_tmp$u
      B_hat <- svd_tmp$v
      D_hat <- svd_tmp$d
      
      A <- A_hat %*% diag(sign(A_hat[1,]))  # orientamento
      B <- diag(D_hat) %*% t(B_hat)
      
      E <- Y - A %*% X %*% t(B)
      P_sum <- P_sum + E %*% solve(Q_sum + diag(k2)) %*% t(E)  # iniziale Q unit
      Q_sum <- Q_sum + t(E) %*% solve(P_sum + diag(k1)) %*% E
      
      A_list[[i - p]] <- A
      B_list[[i - p]] <- B
    }
    
    # Stima finale P, Q
    P_hat <- P_sum / n_obs
    Q_hat <- Q_sum / n_obs
    
    # Calcolo log-likelihood matriciale
    llk <- 0
    for (i in (p+1):T) {
      Y <- F_list[[i]]
      X <- F_list[[i - p]]
      A <- A_list[[i - p]]
      B <- B_list[[i - p]]
      E <- Y - A %*% X %*% t(B)
      term <- tr(solve(P_hat) %*% E %*% solve(Q_hat) %*% t(E))
      llk <- llk + term
    }
    
    loglik <- - (n_obs * k1 * log(det(Q_hat)) + n_obs * k2 * log(det(P_hat)) + llk + n_obs * k1 * k2 * log(2*pi)) / 2
    loglik_values[p] <- loglik
    
    # Numero parametri (A, B, P, Q)
    num_params <- 2 * k1 * k2 + 0.5 * k1 * (k1 + 1) + 0.5 * k2 * (k2 + 1)
    bic_values[p] <- -2 * loglik + log(n_obs) * num_params
    aic_values[p] <- -2 * loglik + 2 * num_params
  }
  
  best_p_bic <- which.min(bic_values)
  best_p_aic <- which.min(aic_values)
  
  if (verbose) {
    cat("BIC (matricial MAR):\n")
    print(bic_values)
    cat("Best lag (BIC):", best_p_bic, "\n")
    cat("Best lag (AIC):", best_p_aic, "\n")
    
    df <- data.frame(
      Lag = 1:max_lag,
      BIC = bic_values,
      AIC = aic_values
    )
    
    # Plot
    p1 <- ggplot(df, aes(x = Lag, y = BIC)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      geom_vline(xintercept = best_p_bic, linetype = "dashed", color = "blue") +
      ggtitle("BIC (MAR)") +
      theme_minimal()
    
    p2 <- ggplot(df, aes(x = Lag, y = AIC)) +
      geom_line(color = "red") +
      geom_point(color = "red") +
      geom_vline(xintercept = best_p_aic, linetype = "dashed", color = "red") +
      ggtitle("AIC (MAR)") +
      theme_minimal()
    
    print(p1 + p2)
  }
  
  return(list(
    Factors = Factors,
    BIC = bic_values,
    AIC = aic_values,
    logLik = loglik_values,
    best_lag_BIC = best_p_bic,
    best_lag_AIC = best_p_aic
  ))
}

# Qui ogni ritardo ha una matrice di regressione distinta. Il vettore dei fattori 
# è spiegato da una concatenazione dei p vettori ritardati

mar_model_selection <- function(Y_std, W, k_hat, max_lag = 10, verbose = TRUE) {
  
  imputation <- mfm.cl(Y_std, W, k_hat)
  pe <- mfm.pe(imputation$Y_imputed, k_hat)
  Factors <- pe$factor
  
  T <- dim(Factors)[1]
  k1 <- dim(Factors)[2]
  k2 <- dim(Factors)[3]
  
  F_list <- lapply(1:T, function(t) Factors[t,,])
  
  bic_values <- numeric(max_lag)
  aic_values <- numeric(max_lag)
  loglik_values <- numeric(max_lag)
  
  for (p in 1:max_lag) {
    Y_list <- F_list[(p+1):T]
    
    # Costruisci la lista delle matrici X_t = [F_{t-1}, ..., F_{t-p}]
    X_list <- lapply((p+1):T, function(t) {
      do.call(cbind, lapply(1:p, function(lag) as.vector(F_list[[t - lag]])))
    })
    
    Y_mat <- do.call(rbind, lapply(Y_list, function(m) as.vector(m)))
    X_mat <- do.call(rbind, X_list)
    
    # Regressione multivariata: vec(F_t) = B * [vec(F_{t-1}); ... ; vec(F_{t-p})] + e
    B_hat <- solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% as.matrix(Y_mat)
    Res <- Y_mat - X_mat %*% B_hat
    Sigma_hat <- cov(Res)
    
    n_obs <- T - p
    loglik <- -0.5 * n_obs * (k1 * k2 * log(2 * pi) + log(det(Sigma_hat)) + 1)
    loglik_values[p] <- loglik
    
    # Numero di parametri: p*(k1*k2)^2 + covarianza residuo (simm. pos def)
    num_params <- p * (k1 * k2)^2 + 0.5 * k1 * k2 * (k1 * k2 + 1)
    bic_values[p] <- -2 * loglik + log(n_obs) * num_params
    aic_values[p] <- -2 * loglik + 2 * num_params
  }
  
  best_p_bic <- which.min(bic_values)
  best_p_aic <- which.min(aic_values)
  
  if (verbose) {
    cat("MAR(p) - Model Selection\n")
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
    
    print(p1 + p2)
  }
  
  return(list(
    BIC = bic_values,
    AIC = aic_values,
    logLik = loglik_values,
    best_lag_BIC = best_p_bic,
    best_lag_AIC = best_p_aic
  ))
}

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

# Versione adattata per modelli vettoriali (k1 = 1)

# Modified version of factorize_bic_matrix_loglik_mar for vector DFM.
# - Uses mfm.cl.vec() and mfm.pe.vec() to obtain factors.
# - Each F_t is treated as a 1 x k2 matrix.
# - Matrix autoregression A_t, B_t estimated via SVD between F_t and F_{t-p}.
# - BIC/AIC computed for varying MAR(p) lags.

# === Funzione 1: MAR(1) matrix factor model - stima singolo lag ===
factorize_bic_mar_vec <- function(Y_std, W, k_hat, max_lag = 10, verbose = TRUE) {
  imputation <- mfm.cl.vec(Y_std, W, k_hat)
  pe <- mfm.pe.vec(imputation$Y_imputed, k_hat)
  Factors <- pe$factor
  
  T <- dim(Factors)[1]
  k2 <- dim(Factors)[3]  # solo colonna
  F_list <- lapply(1:T, function(t) Factors[t, 1, ])
  
  bic_values <- numeric(max_lag)
  aic_values <- numeric(max_lag)
  loglik_values <- numeric(max_lag)
  
  for (p in 1:max_lag) {
    n_obs <- T - p
    B_list <- list()
    Q_sum <- matrix(0, k2, k2)
    SSR <- 0
    
    for (i in (p+1):T) {
      Y <- matrix(F_list[[i]], ncol = 1)
      X <- matrix(F_list[[i - p]], ncol = 1)
      
      B_hat <- Y %*% solve(t(X) %*% X) %*% t(X)
      E <- Y - B_hat %*% X
      SSR <- SSR + t(E) %*% E
      
      B_list[[i - p]] <- B_hat
    }
    
    Q_hat <- SSR / n_obs
    
    llk <- -0.5 * n_obs * (k2 * log(2 * pi) + log(det(Q_hat)) + 1)
    loglik_values[p] <- llk
    
    num_params <- k2^2 + 0.5 * k2 * (k2 + 1)
    bic_values[p] <- -2 * llk + log(n_obs) * num_params
    aic_values[p] <- -2 * llk + 2 * num_params
  }
  
  best_p_bic <- which.min(bic_values)
  best_p_aic <- which.min(aic_values)
  
  if (verbose) {
    df <- data.frame(Lag = 1:max_lag, BIC = bic_values, AIC = aic_values)
    p1 <- ggplot(df, aes(x = Lag, y = BIC)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      geom_vline(xintercept = best_p_bic, linetype = "dashed", color = "blue") +
      ggtitle("BIC (DFM)") +
      theme_minimal()
    
    p2 <- ggplot(df, aes(x = Lag, y = AIC)) +
      geom_line(color = "red") +
      geom_point(color = "red") +
      geom_vline(xintercept = best_p_aic, linetype = "dashed", color = "red") +
      ggtitle("AIC (DFM)") +
      theme_minimal()
    
    print(p1 + p2)
  }
  
  return(list(
    Factors = Factors,
    BIC = bic_values,
    AIC = aic_values,
    logLik = loglik_values,
    best_lag_BIC = best_p_bic,
    best_lag_AIC = best_p_aic
  ))
}

# === Funzione 2: modello MAR(p) vettoriale ===

# Modified version of mar_model_selection for vector DFM.
# - Regresses vec(F_t) on stacked vec(F_{t-1}, ..., F_{t-p}) as in a MAR(p) model.
# - Handles factor arrays of shape T x 1 x k2.
# - Computes BIC and AIC to select optimal lag p.

mar_model_selection_vec <- function(Y_std, W, k_hat, max_lag = 10, verbose = TRUE) {
  imputation <- mfm.cl.vec(Y_std, W, k_hat)
  pe <- mfm.pe.vec(imputation$Y_imputed, k_hat)
  Factors <- pe$factor
  
  T <- dim(Factors)[1]
  k2 <- dim(Factors)[3]
  F_list <- lapply(1:T, function(t) Factors[t, 1, ])
  
  bic_values <- numeric(max_lag)
  aic_values <- numeric(max_lag)
  loglik_values <- numeric(max_lag)
  
  for (p in 1:max_lag) {
    Y_list <- F_list[(p+1):T]
    X_list <- lapply((p+1):T, function(t) {
      do.call(cbind, lapply(1:p, function(lag) as.vector(F_list[[t - lag]])))
    })
    
    Y_mat <- do.call(rbind, lapply(Y_list, matrix, nrow = 1))
    X_mat <- do.call(rbind, lapply(X_list, matrix, nrow = 1))
    
    B_hat <- solve(t(X_mat) %*% X_mat) %*% t(X_mat) %*% Y_mat
    Res <- Y_mat - X_mat %*% B_hat
    Sigma_hat <- cov(Res)
    
    n_obs <- T - p
    loglik <- -0.5 * n_obs * (k2 * log(2 * pi) + log(det(Sigma_hat)) + 1)
    loglik_values[p] <- loglik
    
    num_params <- p * k2^2 + 0.5 * k2 * (k2 + 1)
    bic_values[p] <- -2 * loglik + log(n_obs) * num_params
    aic_values[p] <- -2 * loglik + 2 * num_params
  }
  
  best_p_bic <- which.min(bic_values)
  best_p_aic <- which.min(aic_values)
  
  if (verbose) {
    df <- data.frame(Lag = 1:max_lag, BIC = bic_values, AIC = aic_values)
    p1 <- ggplot(df, aes(x = Lag, y = BIC)) +
      geom_line(color = "blue") +
      geom_point(color = "blue") +
      geom_vline(xintercept = best_p_bic, linetype = "dashed", color = "blue") +
      ggtitle("BIC vs Lag (DFM)") +
      theme_minimal()
    
    p2 <- ggplot(df, aes(x = Lag, y = AIC)) +
      geom_line(color = "red") +
      geom_point(color = "red") +
      geom_vline(xintercept = best_p_aic, linetype = "dashed", color = "red") +
      ggtitle("AIC vs Lag (DFM)") +
      theme_minimal()
    
    print(p1 + p2)
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

# ----

# Modified version of mfm.cl for T x 1 x p2 data (vector DFM).
# - Computes only column covariances (S_C), ignoring rows (S_R).
# - Avoids slicing issues when p1 = 1.
# - Estimates column loadings (C_hat) and factor matrix (F_hat).
# - Performs imputation using: Y_hat[t,,] = R_hat %*% F_hat[t,,] %*% t(C_hat)
#   with R_hat set as a 1-row identity matrix.

mfm.cl.vec <- function(Y, W, r) {
  # Dimensioni
  n <- dim(Y)[1]   # T
  d1 <- 1          # un solo paese
  d2 <- dim(Y)[3]  # numero di variabili
  r1 <- r[1]
  r2 <- r[2]
  
  # Covarianze (S_R sarà scalare)
  S_R <- matrix(0, 1, 1)
  S_C <- matrix(0, d2, d2)
  count_R <- matrix(0, 1, 1)
  count_C <- matrix(0, d2, d2)
  
  for (t in 1:n) {
    Y_t <- Y[t, 1, ]
    W_t <- W[t, 1, ]
    
    # Mode-1: solo diagonale (unica riga)
    S_R[1, 1] <- S_R[1, 1] + sum((Y_t * W_t)^2)
    count_R[1, 1] <- count_R[1, 1] + sum(W_t)
    
    # Mode-2: colonne (variabili)
    for (i in 1:d2) {
      for (j in i:d2) {
        if (W_t[i] == 1 && W_t[j] == 1) {
          S_C[i, j] <- S_C[i, j] + Y_t[i] * Y_t[j]
          count_C[i, j] <- count_C[i, j] + 1
        }
      }
    }
  }
  
  # Simmetrizza
  S_C <- (S_C + t(S_C)) - diag(diag(S_C))
  count_C <- (count_C + t(count_C)) - diag(diag(count_C))
  
  # Normalizza
  if (count_R[1, 1] > 0) {
    S_R[1, 1] <- S_R[1, 1] / count_R[1, 1]
  }
  S_C[count_C > 0] <- S_C[count_C > 0] / count_C[count_C > 0]
  
  # Autovettori
  R_hat <- matrix(1, 1, r1)  # unico valore per r1 = 1
  C_hat <- eigen(S_C, symmetric = TRUE)$vectors[, 1:r2, drop = FALSE]
  
  # Stima dei fattori
  F_hat <- array(0, dim = c(n, r1, r2))
  for (t in 1:n) {
    Y_t <- as.vector(Y[t, 1, ])
    W_t <- as.vector(W[t, 1, ])
    Y_vec <- Y_t
    W_vec <- W_t
    
    Q <- kronecker(C_hat, R_hat)
    Q_obs <- Q[W_vec == 1, , drop = FALSE]
    Y_obs <- Y_vec[W_vec == 1]
    
    if (length(Y_obs) >= r1 * r2) {
      f_vec <- solve(t(Q_obs) %*% Q_obs) %*% t(Q_obs) %*% Y_obs
      F_hat[t,,] <- matrix(f_vec, r1, r2)
    }
  }
  
  # Imputazione
  Y_hat <- array(0, dim = c(n, 1, d2))
  Y_imputed <- Y
  for (t in 1:n) {
    Y_hat[t, 1, ] <- R_hat %*% F_hat[t,,] %*% t(C_hat)
    missing <- which(W[t, 1, ] == 0)
    if (length(missing) > 0) {
      Y_imputed[t, 1, missing] <- Y_hat[t, 1, missing]
    }
  }
  
  return(list(
    R = R_hat,
    C = C_hat,
    f = F_hat,
    Y_hat = Y_hat,
    Y_imputed = Y_imputed
  ))
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

# ----

# Projected estimates for vector DFM.
# - Assumes p1 = 1 (one row).
# - Estimates column loadings (C.t) only.
# - Factor matrix F.h[t,,] is computed as: F_t = X_t %*% C
# - Reconstructed signal: S.h[t,,] = R %*% F_t %*% t(C), with R = 1-row identity.

mfm.pe.vec <- function(X, k) {
  n <- dim(X)[1]
  p1 <- dim(X)[2]  # dovrebbe essere 1
  p2 <- dim(X)[3]
  
  # Step 1: Stima della covarianza tra variabili
  M2.h <- matrix(0, p2, p2)
  for (t in 1:n) {
    x_t <- matrix(X[t, 1, ], nrow = p2)
    M2.h <- M2.h + x_t %*% t(x_t)
  }
  M2.h <- M2.h / n
  
  Q2.h <- eigen(M2.h, symmetric = TRUE)$vectors[, 1:k[2], drop = FALSE]
  C.h <- sqrt(p2) * Q2.h
  
  # Fattori
  F.h <- array(0, dim = c(n, 1, k[2]))
  S.h <- array(0, dim = c(n, 1, p2))
  for (t in 1:n) {
    x_t <- matrix(X[t, 1, ], nrow = 1)  # 1 × p2
    f_t <- x_t %*% C.h / p2
    F.h[t, 1, ] <- f_t
    S.h[t, 1, ] <- f_t %*% t(C.h)
  }
  
  return(list(
    row.load = matrix(1, nrow = 1, ncol = 1),  # costante nel DFM vettoriale
    col.load = C.h,
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

# ---------------------------- Dolp Method -------------------------------------

initialize_dmfm_em_inputs <- function(Y, k) {
  # Dimensioni
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  # Step 0: Costruisci matrice W (1 = osservato, 0 = missing)
  W <- array(1, dim = dim(Y))
  W[is.na(Y)] <- 0
  
  # Step 1: Imputazione iniziale (Chen et al.)
  imputazione <- mfm.cl(Y, W, k)
  Y_imp <- imputazione$Y_imputed
  
  # Step 2: Projected Estimator (Yu et al.)
  pe <- mfm.pe(Y_imp, k)
  R0 <- pe$row.load
  C0 <- pe$col.load
  F0 <- pe$factor
  
  # Step 3: Calcolo E_t^(0)
  E0 <- array(0, dim = c(T, p1, p2))
  for (t in 1:T) {
    E0[t,,] <- Y_imp[t,,] - R0 %*% F0[t,,] %*% t(C0)
  }
  
  # Step 4: Calcolo K^(0)
  K0 <- matrix(0, p2, p2)
  for (t in 1:T) {
    K0 <- K0 + t(E0[t,,]) %*% E0[t,,]
  }
  K0 <- diag(diag(K0)) / (T * p1)  # Solo diagonale
  
  # Step 5: Calcolo H^(0)
  K0inv <- diag(1 / diag(K0))  # inverso della diagonale di K0
  H0 <- matrix(0, p1, p1)
  for (t in 1:T) {
    H0 <- H0 + E0[t,,] %*% K0inv %*% t(E0[t,,])
  }
  H0 <- diag(diag(H0)) / (T * p1)  # Solo diagonale
  
  # Step 6: Costruzione ftilde
  ftilde <- t(sapply(1:T, function(t) {
    yt <- as.vector(Y_imp[t,,])
    (1 / (p1 * p2)) * t(kronecker(C0, R0)) %*% yt
  }))
  
  # Step 7: Calcolo BA^(0)
  BA_num <- matrix(0, prod(k), prod(k))
  BA_den <- matrix(0, prod(k), prod(k))
  for (t in 2:T) {
    BA_num <- BA_num + tcrossprod(ftilde[t,], ftilde[t-1,])
    BA_den <- BA_den + tcrossprod(ftilde[t-1,], ftilde[t-1,])
  }
  BA0 <- BA_num %*% solve(BA_den)
  
  # Step 8: Calcolo QP^(0)
  QP0 <- matrix(0, prod(k), prod(k))
  for (t in 2:T) {
    diff <- ftilde[t,] - BA0 %*% ftilde[t-1,]
    QP0 <- QP0 + tcrossprod(diff)
  }
  QP0 <- QP0 / (T-1)
  
  # Step 9: Inizializzazione Ps0, Cs0, Pf0, ff0
  Ps0 <- array(0, c(T, prod(k), prod(k)))
  Cs0 <- array(0, c(T, prod(k), prod(k)))
  Pf0 <- array(0, c(T, prod(k), prod(k)))
  ff0 <- array(0, c(T, k[1], k[2]))
  for (t in 1:T) {
    Ps0[t,,] <- diag(prod(k))  # Identity matrix
    Cs0[t,,] <- matrix(0, prod(k), prod(k))  # Zero matrix
    Pf0[t,,] <- diag(prod(k))  # Identity for prediction error
    ff0[t,,] <- matrix(0, k[1], k[2])  # Factors prediction mean at 0
  }
  
  # Output lista finale
  init <- list(
    R  = R0,
    C  = C0,
    H  = H0,
    K  = K0,
    BA = BA0,
    QP = QP0,
    fs = F0,   # T x k1 x k2
    ff = ff0,  # T x k1 x k2 (predictions)
    Ps = Ps0,  # T x prod(k) x prod(k) smoothing error
    Pf = Pf0,  # T x prod(k) x prod(k) prediction error
    Cs = Cs0   # T x prod(k) x prod(k) cross covariances
  )
  
  return(list(
    init = init,
    W = W,
    Y = Y,
    k = k
  ))
}


# ---------------------------- Trapin Method -----------------------------------


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

