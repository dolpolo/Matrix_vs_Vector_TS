# ==============================================================================
# STANDARDIZATION with NA
# ==============================================================================

# Standardization of variables in the tensor
standardize_mat_with_na <- function(X) {
  T <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(X)[3]
  
  X_std <- array(NA, dim = dim(X), dimnames = dimnames(X))
  mean_X <- matrix(0, p1, p2)
  sd_X <- matrix(1, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      X_ij <- X[, i, j]
      mu <- mean(X_ij, na.rm = TRUE)
      sigma <- sd(X_ij, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      X_std[, i, j] <- (X_ij - mu) / sigma
      mean_X[i, j] <- mu
      sd_X[i, j] <- sigma
    }
  }
  
  return(list(X_scaled = X_std, mean = mean_X, sd = sd_X))
}


center_Y <- function(Y) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  Y_centered <- array(NA, dim = dim(Y))
  dimnames(Y_centered) <- dimnames(Y) 
  mean_Y <- matrix(0, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu <- mean(y_ij, na.rm = TRUE)
      Y_centered[, i, j] <- y_ij - mu
      mean_Y[i, j] <- mu
    }
  }
  
  return(list(Y_centered = Y_centered, mean = mean_Y))
}

decenter_Y <- function(Y_centered, mean_Y) {
  T <- dim(Y_centered)[1]
  p1 <- dim(Y_centered)[2]
  p2 <- dim(Y_centered)[3]
  
  Y_original <- array(NA, dim = dim(Y_centered))
  dimnames(Y_original) <- dimnames(Y_centered)
  
  for (t in 1:T) {
    Y_original[t,,] <- Y_centered[t,,] + mean_Y
  }
  
  return(Y_original)
}


# ==============================================================================
# UTILS
# ==============================================================================


make_global_col_names <- function(countries, var_names, N_m, N_q) {
  P1 <- length(countries)
  V  <- length(var_names)
  stopifnot(V == N_m + N_q)
  
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  
  col_names <- character(N_global)
  col_country <- character(N_global)
  col_var     <- character(N_global)
  col_freq    <- character(N_global)   # "M" / "Q"
  
  idx <- 1
  # blocco mensili
  for (p in seq_len(P1)) {
    for (k in seq_len(N_m)) {
      v <- k
      col_names[idx]   <- paste0(countries[p], "_", var_names[v])
      col_country[idx] <- countries[p]
      col_var[idx]     <- var_names[v]
      col_freq[idx]    <- "M"
      idx <- idx + 1
    }
  }
  # blocco trimestrali
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      col_names[idx]   <- paste0(countries[p], "_", var_names[v])
      col_country[idx] <- countries[p]
      col_var[idx]     <- var_names[v]
      col_freq[idx]    <- "Q"
      idx <- idx + 1
    }
  }
  
  list(
    names    = col_names,
    country  = col_country,
    var_base = col_var,
    freq     = col_freq
  )
}



# ==============================================================================
# Aggregation rule
# ==============================================================================

aggregate_tensor_to_quarterly <- function(X_tens, agg, N_m, N_q) {
  # X_tens: [T_m x P1 x (N_m+N_q)] già imputato (niente NA idealmente)
  # agg   : [P1 x (N_m+N_q)] con codici 1=stock(media), 2=flow(somma)
  
  T_m <- dim(X_tens)[1]
  P1  <- dim(X_tens)[2]
  V   <- dim(X_tens)[3]
  stopifnot(V == (N_m + N_q))
  stopifnot(all(dim(agg) == c(P1, V)))
  
  T_q <- floor(T_m / 3)
  idx_m3 <- seq(3, by = 3, length.out = T_q)
  
  dn <- dimnames(X_tens)
  time_q <- if (!is.null(dn[[1]])) dn[[1]][idx_m3] else as.character(idx_m3)
  
  X_q <- array(
    NA_real_,
    dim = c(T_q, P1, V),
    dimnames = list(time_q, dn[[2]], dn[[3]])
  )
  
  for (tau in seq_len(T_q)) {
    m3 <- 3 * tau
    m2 <- m3 - 1
    m1 <- m3 - 2
    
    block <- X_tens[c(m1, m2, m3), , , drop = FALSE]  # [3 x P1 x V]
    
    # per ogni (p,v) applica stock/flow sul blocco di 3 mesi
    for (p in seq_len(P1)) {
      for (v in seq_len(V)) {
        if (agg[p, v] == 2) {
          X_q[tau, p, v] <- sum(block[, p, v], na.rm = TRUE)
        } else if (agg[p, v] == 1) {
          X_q[tau, p, v] <- mean(block[, p, v], na.rm = TRUE)
        } else {
          stop("agg must be 1 (stock) or 2 (flow).")
        }
      }
    }
  }
  
  X_q
}

# ==============================================================================
# VECTORIZATION  (Tensor -> vector)
# ==============================================================================
tensor_to_vector <- function(X_tens, N_m, N_q) {
  stopifnot(length(dim(X_tens)) == 3)
  T_gl <- dim(X_tens)[1]
  P1   <- dim(X_tens)[2]
  V    <- dim(X_tens)[3]
  stopifnot(V == (N_m + N_q))
  
  dn <- dimnames(X_tens)
  time_names <- dn[[1]]
  countries  <- dn[[2]]
  var_names  <- dn[[3]]
  
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  
  X_mat <- matrix(NA_real_, nrow = T_gl, ncol = N_global)
  
  # rownames
  if (!is.null(time_names)) rownames(X_mat) <- time_names
  
  # colnames coerenti (se hai make_global_col_names usalo)
  col_meta <- make_global_col_names(
    countries = countries,
    var_names = var_names,
    N_m       = N_m,
    N_q       = N_q
  )
  colnames(X_mat) <- col_meta$names
  
  # fill: M block then Q block, country-major
  col_idx <- 1
  
  # Monthly block
  for (p in seq_len(P1)) {
    for (k in seq_len(N_m)) {
      v <- k
      X_mat[, col_idx] <- X_tens[, p, v]
      col_idx <- col_idx + 1
    }
  }
  
  # Quarterly block
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      X_mat[, col_idx] <- X_tens[, p, v]
      col_idx <- col_idx + 1
    }
  }
  
  # store metadata for perfect inversion
  attr(X_mat, "N_m")        <- N_m
  attr(X_mat, "N_q")        <- N_q
  attr(X_mat, "countries")  <- countries
  attr(X_mat, "var_names")  <- var_names
  
  X_mat
}


# ==============================================================================
# TENSORIZATION  (vector -> Tensor)
# ==============================================================================
vector_to_tensor <- function(X_mat,
                             N_m       = attr(X_mat, "N_m"),
                             N_q       = attr(X_mat, "N_q"),
                             countries = attr(X_mat, "countries"),
                             var_names = attr(X_mat, "var_names")) {
  stopifnot(is.matrix(X_mat))
  stopifnot(!is.null(N_m), !is.null(N_q), !is.null(countries), !is.null(var_names))
  
  T_gl <- nrow(X_mat)
  P1   <- length(countries)
  V    <- N_m + N_q
  
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  stopifnot(ncol(X_mat) == N_global)
  
  dates <- rownames(X_mat)
  if (is.null(dates)) dates <- as.character(seq_len(T_gl))
  
  X_tens <- array(
    NA_real_,
    dim = c(T_gl, P1, V),
    dimnames = list(dates, countries, var_names)
  )
  
  col_idx <- 1
  
  # Monthly block
  for (p in seq_len(P1)) {
    for (k in seq_len(N_m)) {
      v <- k
      X_tens[, p, v] <- X_mat[, col_idx]
      col_idx <- col_idx + 1
    }
  }
  
  # Quarterly block
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      X_tens[, p, v] <- X_mat[, col_idx]
      col_idx <- col_idx + 1
    }
  }
  
  X_tens
}
