# ==============================================================================
# STANDARDIZATION HANDLING NA
# ==============================================================================

standardize_with_na <- function(X) {
  # X : matrix T x N with NA
  
  # Compute column means and sds (ignoring NAs)
  col_means <- apply(X, 2, mean, na.rm = TRUE)
  col_sds   <- apply(X, 2, sd,   na.rm = TRUE)
  
  # Replace zero or NA std devs with 1 (prevents division by zero)
  col_sds[col_sds == 0 | is.na(col_sds)] <- 1
  
  # Standardize
  X_centered <- sweep(X, 2, col_means, "-")
  X_std      <- sweep(X_centered, 2, col_sds, "/")
  
  return(list(
    X_std   = X_std,
    mean    = col_means,
    sd      = col_sds
  ))
}



destandardize <- function(X_std, mean, sd) {
  # X_std: matrix T x N in standardized scale
  sweep(sweep(X_std, 2, sd, "*"), 2, mean, "+")
}

