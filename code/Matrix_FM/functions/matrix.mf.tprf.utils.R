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


# X_tens <- X_std
# W_tens <- W_x
# var_type <- Freq
# agg_q
flatten_tensor_to_matrix <- function(X_tens, W_tens,
                                     N_m, N_q, agg_q) {
  
  # --------------------------------------------------
  # 1. METADATA DAL TENSORE
  # --------------------------------------------------
  countries <- dimnames(X_tens)[[2]]
  var_names <- dimnames(X_tens)[[3]]
  time_names <- dimnames(X_tens)[[1]]
  
  T_gl <- dim(X_tens)[1]
  P1   <- dim(X_tens)[2]
  V    <- dim(X_tens)[3]
  
  stopifnot(V == N_m + N_q)
  stopifnot(length(countries) == P1)
  
  # --------------------------------------------------
  # 2. DIMENSIONI GLOBALI
  # --------------------------------------------------
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  
  # --------------------------------------------------
  # 3. ALLOCAZIONE OUTPUT
  # --------------------------------------------------
  X_mat <- matrix(NA_real_, nrow = T_gl, ncol = N_global)
  W_mat <- matrix(0L,       nrow = T_gl, ncol = N_global)
  
  if (!is.null(time_names)) {
    stopifnot(length(time_names) == T_gl)
    rownames(X_mat) <- time_names
    rownames(W_mat) <- time_names
  }
  
  # --------------------------------------------------
  # 4. NOMI DELLE COLONNE (COERENTI CON L’ORDINE)
  # --------------------------------------------------
  col_meta <- make_global_col_names(
    countries = countries,
    var_names = var_names,
    N_m       = N_m,
    N_q       = N_q
  )
  
  colnames(X_mat) <- col_meta$names
  colnames(W_mat) <- col_meta$names
  
  # --------------------------------------------------
  # 5. RIEMPIMENTO MATRICE (STESSO ORDINE DEI NOMI)
  # --------------------------------------------------
  col_idx <- 1
  
  ## Mensili
  for (p in seq_len(P1)) {
    for (v in seq_len(N_m)) {
      X_mat[, col_idx] <- X_tens[, p, v]
      W_mat[, col_idx] <- W_tens[, p, v]
      col_idx <- col_idx + 1
    }
  }
  
  ## Trimestrali
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      X_mat[, col_idx] <- X_tens[, p, v]
      W_mat[, col_idx] <- W_tens[, p, v]
      col_idx <- col_idx + 1
    }
  }
  
  # --------------------------------------------------
  # 6. AGGREGATION CODE PER Q (STESSO ORDINE)
  # --------------------------------------------------
  agg_q_global <- character(N_q_global)
  idx <- 1
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      agg_q_global[idx] <- as.character(agg_q[p, k])
      idx <- idx + 1
    }
  }
  
  # --------------------------------------------------
  # 7. ATTRIBUTI UTILI A VALLE
  # --------------------------------------------------
  attr(X_mat, "N_m")       <- N_m
  attr(X_mat, "N_q")       <- N_q
  attr(X_mat, "countries")<- countries
  attr(X_mat, "var_names")<- var_names
  
  # --------------------------------------------------
  # 8. OUTPUT
  # --------------------------------------------------
  list(
    X_mat        = X_mat,
    W_mat        = W_mat,
    N_m_global   = N_m_global,
    N_q_global   = N_q_global,
    N_global     = N_global,
    agg_q_global = agg_q_global
  )
}

################################################################################
##################################### X Vettorizzata ###########################
################################################################################
flatten_tensor_init <- function(
    X_tens,
    N_m,
    N_q,
    countries = dimnames(X_tens)[[2]],
    var_names = dimnames(X_tens)[[3]]
) {
  T_gl <- dim(X_tens)[1]
  P1   <- dim(X_tens)[2]
  V    <- dim(X_tens)[3]
  stopifnot(V == N_m + N_q)
  
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  
  X_mat <- matrix(NA_real_, nrow = T_gl, ncol = N_global)
  
  # nomi del tempo
  time_names <- dimnames(X_tens)[[1]]
  if (!is.null(time_names)) {
    stopifnot(length(time_names) == T_gl)
    rownames(X_mat) <- time_names
  }
  
  # stessi colnames di flatten_tensor_to_matrix
  col_meta  <- make_global_col_names(countries, var_names, N_m, N_q)
  col_names <- col_meta$names
  
  # riempi i dati nello stesso ordine
  col_idx <- 1
  # blocco M
  for (p in seq_len(P1)) {
    for (k in seq_len(N_m)) {
      v <- k
      X_mat[ , col_idx] <- X_tens[ , p, v]
      col_idx <- col_idx + 1
    }
  }
  # blocco Q
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      X_mat[ , col_idx] <- X_tens[ , p, v]
      col_idx <- col_idx + 1
    }
  }
  
  colnames(X_mat) <- col_names
  
  X_mat
}




################################################################################
##################################### Tensorize ###########################
################################################################################
unflatten_matrix_to_tensor <- function(X_mat,
                                       N_m      = attr(X_mat, "N_m"),
                                       N_q      = attr(X_mat, "N_q"),
                                       countries= attr(X_mat, "countries"),
                                       var_names= attr(X_mat, "var_names"),
                                       dates    = rownames(X_mat)) {
  T_gl <- nrow(X_mat)
  P1   <- length(countries)
  V    <- N_m + N_q
  
  N_m_global <- P1 * N_m
  N_q_global <- P1 * N_q
  N_global   <- N_m_global + N_q_global
  stopifnot(ncol(X_mat) == N_global)
  
  # inizializza tensore
  X_tens <- array(
    NA_real_,
    dim = c(T_gl, P1, V),
    dimnames = list(
      dates,
      countries,
      var_names
    )
  )
  
  # blocco mensile
  col_idx <- 1
  for (p in seq_len(P1)) {
    for (k in seq_len(N_m)) {
      v <- k
      X_tens[ , p, v] <- X_mat[ , col_idx]
      col_idx <- col_idx + 1
    }
  }
  
  # blocco trimestrale
  for (p in seq_len(P1)) {
    for (k in seq_len(N_q)) {
      v <- N_m + k
      X_tens[ , p, v] <- X_mat[ , col_idx]
      col_idx <- col_idx + 1
    }
  }
  
  X_tens
}

