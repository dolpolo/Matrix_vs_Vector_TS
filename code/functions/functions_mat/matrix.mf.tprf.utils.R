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


# ==============================================================================
# UTILITIES FOR RESULTS
# ==============================================================================

get_size_tag <- function(N_m, N_q) {
  total <- N_m + N_q
  
  if (total <= 25) {
    "small"
  } else if (total <= 50) {
    "medium"
  } else {
    "large"
  }
}

build_result_filename <- function(path_out,
                                  model,
                                  stage,
                                  Size,
                                  sel,
                                  countries = NULL,
                                  N_m = NA,
                                  N_q = NA,
                                  Lproxy = NA,
                                  L_midas = NA,
                                  p_ar = NA,
                                  r1 = NA,
                                  r2 = NA,
                                  robust_f = NA,
                                  covid_m = NA,
                                  covid_q = NA,
                                  ext = "rds",
                                  timestamp = TRUE,
                                  include_details = TRUE) {
  
  dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
  
  base_name <- paste0(
    stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel
  )
  
  if (isTRUE(include_details)) {
    cc_lab <- if (!is.null(countries) && length(countries) > 0) {
      paste(countries, collapse = "-")
    } else {
      "NA"
    }
    
    base_name <- paste0(
      base_name,
      "_cc-", cc_lab,
      "_Nm-", N_m,
      "_Nq-", N_q,
      "_Lproxy-", Lproxy,
      "_Lmidas-", L_midas,
      "_pAR-", p_ar,
      "_r1-", r1,
      "_r2-", r2,
      "_RobustF-", robust_f,
      "_CovidM-", covid_m,
      "_CovidQ-", covid_q
    )
  }
  
  if (isTRUE(timestamp)) {
    base_name <- paste0(
      base_name,
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  }
  
  file.path(path_out, paste0(base_name, ".", ext))
}


list_to_df_country <- function(lst_by_country, tag) {
  if (is.null(lst_by_country) || length(lst_by_country) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      country          = character(),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  
  out <- lapply(names(lst_by_country), function(cc) {
    x <- lst_by_country[[cc]]
    if (is.null(x) || length(x) == 0L) return(NULL)
    
    data.frame(
      date             = as.Date(names(x)),
      country          = cc,
      nowcast          = as.numeric(unlist(x)),
      month_in_quarter = tag
    )
  })
  
  dplyr::bind_rows(out)
}

find_result_file <- function(path, model, stage, Size, sel, ext = "rds") {
  pattern <- paste0(
    "^", stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel,
    ".*\\.", ext, "$"
  )
  
  files <- list.files(path, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0L) {
    stop(
      "No file found for model = ", model,
      ", stage = ", stage,
      ", Size = ", Size,
      ", sel = ", sel
    )
  }
  
  if (length(files) > 1L) {
    message("Multiple files found. Using the most recent one.")
    files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  }
  
  files[1]
}


rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

list_to_latex_table <- function(df, caption, label, digits = 4) {
  stopifnot(all(c("country", "period", "M1", "M2", "M3") %in% names(df)))
  
  rows <- sprintf(
    "%s & %s & %.*f & %.*f & %.*f \\\\",
    df$country, df$period,
    digits, df$M1,
    digits, df$M2,
    digits, df$M3
  )
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{llccc}\n",
    "\\toprule\n",
    "Country & Period & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(rows, collapse = "\n"),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}


# ==============================================================================
# 2. HELPERS
# ==============================================================================

normalize_summary_table <- function(df, model_name, default_period = "Full sample") {
  
  if (!is.data.frame(df)) {
    stop("Input object must be a data.frame.")
  }
  
  required_basic <- c("country", "M1", "M2", "M3")
  if (!all(required_basic %in% names(df))) {
    stop(
      "Summary table must contain at least these columns: ",
      paste(required_basic, collapse = ", ")
    )
  }
  
  if (!"period" %in% names(df)) {
    df$period <- default_period
  }
  
  df %>%
    mutate(
      model   = model_name,
      country = as.character(country),
      period  = as.character(period)
    ) %>%
    select(country, period, M1, M2, M3, model)
}

build_comparison_table <- function(df_matrix, df_vector, df_vectensor) {
  bind_rows(df_matrix, df_vector, df_vectensor) %>%
    pivot_longer(
      cols = c("M1", "M2", "M3"),
      names_to = "horizon",
      values_to = "RMSFE"
    ) %>%
    mutate(col_id = paste0(model, "_", horizon)) %>%
    select(country, period, col_id, RMSFE) %>%
    pivot_wider(names_from = col_id, values_from = RMSFE) %>%
    arrange(
      country,
      factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID"))
    )
}

comparison_to_latex <- function(df_comp, caption, label) {
  
  needed <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  for (nm in needed) {
    if (!nm %in% names(df_comp)) df_comp[[nm]] <- NA_real_
  }
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{llccc ccc ccc}\n",
    "\\toprule\n",
    " & & \\multicolumn{3}{c}{Matrix MF-TPRF} & \\multicolumn{3}{c}{Vector MF-TPRF} & \\multicolumn{3}{c}{Vec-from-Tensor MF-TPRF} \\\\\n",
    "Country & Period & M1 & M2 & M3 & M1 & M2 & M3 & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf(
        "%s & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
        df_comp$country,
        df_comp$period,
        df_comp$Matrix_M1, df_comp$Matrix_M2, df_comp$Matrix_M3,
        df_comp$Vector_M1, df_comp$Vector_M2, df_comp$Vector_M3,
        df_comp$VecTensor_M1, df_comp$VecTensor_M2, df_comp$VecTensor_M3
      ),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

build_relative_table <- function(df_comp) {
  
  needed <- c(
    "country", "period",
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  missing_cols <- setdiff(needed, names(df_comp))
  if (length(missing_cols) > 0) {
    stop("Missing columns in comparison table: ", paste(missing_cols, collapse = ", "))
  }
  
  df_comp %>%
    mutate(
      rel_MV_M1 = Matrix_M1 / Vector_M1,
      rel_MV_M2 = Matrix_M2 / Vector_M2,
      rel_MV_M3 = Matrix_M3 / Vector_M3,
      rel_MT_M1 = Matrix_M1 / VecTensor_M1,
      rel_MT_M2 = Matrix_M2 / VecTensor_M2,
      rel_MT_M3 = Matrix_M3 / VecTensor_M3,
      rel_VT_M1 = Vector_M1 / VecTensor_M1,
      rel_VT_M2 = Vector_M2 / VecTensor_M2,
      rel_VT_M3 = Vector_M3 / VecTensor_M3
    ) %>%
    select(
      country, period,
      rel_MV_M1, rel_MV_M2, rel_MV_M3,
      rel_MT_M1, rel_MT_M2, rel_MT_M3,
      rel_VT_M1, rel_VT_M2, rel_VT_M3
    ) %>%
    arrange(
      country,
      factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID"))
    )
}

relative_to_latex <- function(df_rel, caption, label) {
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{llccc ccc ccc}\n",
    "\\toprule\n",
    " & & \\multicolumn{3}{c}{Matrix / Vector} & \\multicolumn{3}{c}{Matrix / VecTensor} & \\multicolumn{3}{c}{Vector / VecTensor} \\\\\n",
    "Country & Period & M1 & M2 & M3 & M1 & M2 & M3 & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf(
        "%s & %s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\\\",
        df_rel$country,
        df_rel$period,
        df_rel$rel_MV_M1, df_rel$rel_MV_M2, df_rel$rel_MV_M3,
        df_rel$rel_MT_M1, df_rel$rel_MT_M2, df_rel$rel_MT_M3,
        df_rel$rel_VT_M1, df_rel$rel_VT_M2, df_rel$rel_VT_M3
      ),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

build_saved_filename <- function(stage, model, Size, sel, ext = "rds", timestamp = TRUE) {
  
  base_name <- paste0(
    stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel
  )
  
  if (timestamp) {
    base_name <- paste0(base_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  
  paste0(base_name, ".", ext)
}


# ==============================================================================
# 3. LOCAL HELPERS
# ==============================================================================

build_run_tag <- function(params) {
  paste0(
    "eval-", format(params$start_eval, "%Y%m"), "_", format(params$end_eval, "%Y%m"),
    "_sel-", params$sel_method,
    "_Nm-", params$n_m,
    "_Nq-", params$n_q,
    "_KmaxVec-", params$Kmax_vec,
    "_Zmax-", params$Zmax,
    "_Lmax-", params$Lmax,
    "_pARmax-", params$p_AR_max,
    "_CovidM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-", as.integer(isTRUE(params$covid_mask_q))
  )
}

list_to_df_nowcast <- function(lst, tag) {
  if (length(lst) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  
  data.frame(
    date             = as.Date(names(lst)),
    nowcast          = as.numeric(unlist(lst)),
    month_in_quarter = tag,
    row.names        = NULL
  )
}
