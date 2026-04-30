# ==============================================================================
# SELECT THE RIGHT NUMBER OF FACTORS OF X (BAI-NG 2013)
# ==============================================================================

# compute the number of factors just on monthly variables
# Xm <- X[,1:N_m]

# cols_complete <- which(colSums(is.na(Xm)) == 0)
# X_m_completed <- Xm[, cols_complete]

IC_PC_criteria <- function(Xm, k_max = 10) {
  
  # 1. Standardizzazione opzionale (consigliata)
  X <- scale(Xm, center = FALSE, scale = FALSE)
  
  N <- ncol(X)
  T <- nrow(X)
  
  # 2. PCA
  pca <- prcomp(X, center = FALSE, scale. = FALSE)
  
  # Otteniamo loadings e fattori stimati
  F_hat <- pca$x               # T × N
  Lambda_hat <- pca$rotation   # N × N
  
  # 3. Pre-allocazione risultati
  V <- numeric(k_max)
  
  IC1 <- numeric(k_max)
  IC2 <- numeric(k_max)
  IC3 <- numeric(k_max)
  
  PC1 <- numeric(k_max)
  PC2 <- numeric(k_max)
  PC3 <- numeric(k_max)
  
  # Penalizzazioni corrette Bai-Ng 2003
  g1 <- log(N*T)/(N*T)
  g2 <- log(min(N,T)) / min(N,T)
  g3 <- log(min(N,T)) / (N*T)
  
  # 4. Loop sui possibili numeri di fattori
  for(k in 1:k_max){
    
    # ricostruzione con k fattori
    F_k <- F_hat[, 1:k]                
    L_k <- Lambda_hat[, 1:k]           
    X_hat <- F_k %*% t(L_k)            
    
    # V(k)
    resid <- X - X_hat
    V[k] <- sum(resid^2) / (N * T)
    
    # Criteri IC
    IC1[k] <- log(V[k]) + k * g1
    IC2[k] <- log(V[k]) + k * g2
    IC3[k] <- log(V[k]) + k * g3
    
    # Criteri PC
    PC1[k] <- V[k] + k * g1
    PC2[k] <- V[k] + k * g2
    PC3[k] <- V[k] + k * g3
  }
  
  # 5. Tabella finale
  results <- data.frame(
    k = 1:k_max,
    V = V,
    IC1 = IC1,
    IC2 = IC2,
    IC3 = IC3,
    PC1 = PC1,
    PC2 = PC2,
    PC3 = PC3
  )
  
  return(results)
}


# K_sel <- IC_PC_criteria(X_m_completed, 10)

# 6. Trovo i minimi dei criteri
# which.min(K_sel$IC1)
# which.min(K_sel$IC2)
# which.min(K_sel$IC3)
# which.min(K_sel$PC1)
# which.min(K_sel$PC2)
# which.min(K_sel$PC3)




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
  
  eps <- 1e-12
  ER_vec <- eigvals[1:max_k] / pmax(eigvals[2:(max_k + 1)], eps)
  
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

# ==============================================================================
# PLOT FACTOR SELECTED
# ==============================================================================

plot_factors_XP <- function(F_hat, dates = NULL, center = TRUE,
                            ncol = 1, scales = "free_y",
                            title = "XP Estimated Factors") {
  stopifnot(is.matrix(F_hat) || is.data.frame(F_hat))
  F_hat <- as.matrix(F_hat)
  
  # center factors (aesthetic only)
  if (center) F_hat <- scale(F_hat, center = TRUE, scale = FALSE)
  
  df <- as.data.frame(F_hat)
  colnames(df) <- paste0("Factor_", seq_len(ncol(df)))
  
  df$Time <- if (!is.null(dates)) as.Date(dates) else seq_len(nrow(df))
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = dplyr::starts_with("Factor_"),
    names_to = "Factor",
    values_to = "Value"
  )
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = Time, y = Value)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ Factor, ncol = ncol, scales = scales) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = title, x = NULL, y = NULL)
}

plot_factor_selection_diag <- function(ER_out, Kmax = NULL, r_sel = NULL,
                                       clip_neg = TRUE, eps = 1e-12) {
  eig <- sort(as.numeric(ER_out$eig), decreasing = TRUE)
  if (clip_neg) eig <- pmax(eig, 0)
  
  if (is.null(Kmax)) Kmax <- min(length(eig) - 1, length(ER_out$ER))
  Kmax <- min(Kmax, length(eig) - 1)
  if (is.null(r_sel)) r_sel <- ER_out$r
  
  ER_vec <- eig[1:Kmax] / pmax(eig[2:(Kmax + 1)], eps)
  df_er  <- data.frame(k = 1:Kmax, ER = ER_vec)
  
  df_var <- data.frame(
    k = 1:length(eig),
    cum_var = cumsum(eig) / max(sum(eig), eps)
  )
  
  p1 <- ggplot2::ggplot(df_er, ggplot2::aes(x = k, y = ER)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_vline(xintercept = r_sel, linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = 1:Kmax) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = "Eigenvalue Ratio", x = "k", y = expression(lambda[k]/lambda[k+1]))
  
  p2 <- ggplot2::ggplot(df_var, ggplot2::aes(x = k, y = cum_var)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_vline(xintercept = r_sel, linetype = "dashed") +
    ggplot2::scale_x_continuous(breaks = seq(1, max(df_var$k), by = 5)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::labs(title = "Cumulative Explained Variance",
                  x = "Number of components", y = "Cumulative variance")
  
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p1 + p2 + patchwork::plot_layout(widths = c(1, 1))
  } else {
    list(p1 = p1, p2 = p2)
  }
}



# ==============================================================================
# LAG of FACTOR SELECTION
# ==============================================================================

select_p_var_ic <- function(F_hat, pmax = 9, include_const = TRUE) {
  F_hat <- as.matrix(F_hat)
  if (anyNA(F_hat)) stop("F_hat contains NA: remove or impute before VAR order selection.")
  
  Tt <- nrow(F_hat)
  r  <- ncol(F_hat)
  
  out <- data.frame(p = 1:pmax, AIC = NA_real_, BIC = NA_real_, HQ = NA_real_)
  
  for (p in 1:pmax) {
    Y <- F_hat[(p+1):Tt, , drop = FALSE]            # (T-p) x r
    Z <- NULL
    
    if (include_const) {
      Z <- matrix(1, nrow = Tt - p, ncol = 1)
    } else {
      Z <- matrix(nrow = Tt - p, ncol = 0)
    }
    
    for (lag in 1:p) {
      Z <- cbind(Z, F_hat[(p+1-lag):(Tt-lag), , drop = FALSE])
    }
    
    Tp <- nrow(Y)
    k  <- r^2 * p + (include_const * r)
    
    Bhat <- solve(crossprod(Z), crossprod(Z, Y))    # OLS multivariato
    E    <- Y - Z %*% Bhat
    Sigma <- crossprod(E) / Tp
    
    ld <- determinant(Sigma, logarithm = TRUE)$modulus
    ld <- as.numeric(ld)
    
    out$AIC[out$p == p] <- ld + (2 / Tp) * k
    out$BIC[out$p == p] <- ld + (log(Tp) / Tp) * k
    out$HQ[out$p == p]  <- ld + (2 * log(log(Tp)) / Tp) * k
  }
  
  list(
    table = out,
    p_AIC = out$p[which.min(out$AIC)],
    p_BIC = out$p[which.min(out$BIC)],
    p_HQ  = out$p[which.min(out$HQ)]
  )
}


# ==============================================================================
# GLOBAL AR(p) ORDER SELECTION FOR IDIOSYNCRATIC RESIDUALS
# - Ehat: T x N matrix (can contain NA)
# - Selects ONE common p for all series by minimizing summed IC across series
# - IC style matches select_p_var_ic (log sigma^2 + penalty/T)
# ==============================================================================

select_p_ar_ic_global <- function(Ehat, pmax = 9, include_const = TRUE, min_T = 30) {
  Ehat <- as.matrix(Ehat)
  Tt <- nrow(Ehat)
  N  <- ncol(Ehat)
  
  out <- data.frame(p = 0:pmax, AIC = NA_real_, BIC = NA_real_, HQ = NA_real_,
                    n_used = NA_integer_)
  
  # pre-count usable obs per series
  n_obs <- colSums(!is.na(Ehat))
  keep  <- which(n_obs >= min_T)
  if (length(keep) == 0) stop("No series has enough observations (min_T) to select AR order.")
  
  for (p in 0:pmax) {
    aic_sum <- 0
    bic_sum <- 0
    hq_sum  <- 0
    used    <- 0L
    
    for (j in keep) {
      e <- Ehat[, j]
      e <- e[!is.na(e)]
      Tj <- length(e)
      
      # need enough data for this p
      if (Tj <= p + 5) next
      
      if (p == 0) {
        # white noise, optional intercept
        mu   <- if (include_const) mean(e) else 0
        u    <- e - mu
        sig2 <- mean(u^2)
        Tp   <- Tj
        k    <- 1L + as.integer(include_const)   # variance + intercept(optional)
      } else {
        Y <- e[(p + 1):Tj]
        Z <- do.call(cbind, lapply(1:p, function(l) e[(p + 1 - l):(Tj - l)]))
        if (include_const) Z <- cbind(1, Z)
        
        # OLS
        B <- solve(crossprod(Z), crossprod(Z, Y))
        u <- as.numeric(Y - Z %*% B)
        sig2 <- mean(u^2)
        Tp <- length(u)
        
        # parameters: AR(p) + variance + intercept(optional)
        k <- p + 1L + as.integer(include_const)
      }
      
      # guard
      if (!is.finite(sig2) || sig2 <= 0) next
      
      # IC per observation (same structure as your VAR function)
      ld <- log(sig2)
      aic_sum <- aic_sum + ld + (2 / Tp) * k
      bic_sum <- bic_sum + ld + (log(Tp) / Tp) * k
      hq_sum  <- hq_sum  + ld + (2 * log(log(Tp)) / Tp) * k
      
      used <- used + 1L
    }
    
    out$AIC[out$p == p] <- aic_sum
    out$BIC[out$p == p] <- bic_sum
    out$HQ[out$p == p]  <- hq_sum
    out$n_used[out$p == p] <- used
  }
  
  # If some p used zero series, they will be 0 sums; set to Inf so they can't win
  out$AIC[out$n_used == 0] <- Inf
  out$BIC[out$n_used == 0] <- Inf
  out$HQ[out$n_used  == 0] <- Inf
  
  list(
    table = out,
    p_AIC = out$p[which.min(out$AIC)],
    p_BIC = out$p[which.min(out$BIC)],
    p_HQ  = out$p[which.min(out$HQ)]
  )
}
