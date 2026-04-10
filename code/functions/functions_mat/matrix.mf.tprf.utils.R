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

format_month_year <- function(x) {
  format(x, "%b %Y")
}

build_period_labels <- function(params) {
  pre_end    <- seq(params$covid_start, by = "-1 month", length.out = 2)[2]
  post_start <- seq(params$covid_end,   by = "+1 month", length.out = 2)[2]
  
  list(
    "Pre-COVID" = c(
      title = "Pre-COVID",
      date  = paste0(format_month_year(params$start_eval), " -- ", format_month_year(pre_end))
    ),
    "COVID period" = c(
      title = "COVID",
      date  = paste0(format_month_year(params$covid_start), " -- ", format_month_year(params$covid_end))
    ),
    "Post-COVID" = c(
      title = "Post-COVID",
      date  = paste0(format_month_year(post_start), " -- ", format_month_year(params$end_eval))
    )
  )
}

build_final_large_style_latex <- function(
    summary_matrix,
    summary_vector,
    summary_vectensor,
    params,
    dm_vecp_wide = NULL,
    dm_vecc_wide = NULL,
    Size = "large",
    sel = "LASSO"
) {
  
  country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  period_order  <- c("Pre-COVID", "COVID period", "Post-COVID")
  period_labels <- build_period_labels(params)
  
  # ---------------------------------------------------------------------------
  # 1. Extract rolling RMSFE tables
  # ---------------------------------------------------------------------------
  mat_rt <- summary_matrix$tab_rt_all %>%
    dplyr::filter(country %in% country_order, period %in% period_order) %>%
    dplyr::mutate(
      country = factor(country, levels = country_order),
      period  = factor(period, levels = period_order)
    ) %>%
    dplyr::arrange(period, country)
  
  vec_rt <- summary_vector$tab_rt_all %>%
    dplyr::filter(country %in% country_order, period %in% period_order) %>%
    dplyr::mutate(
      country = factor(country, levels = country_order),
      period  = factor(period, levels = period_order)
    ) %>%
    dplyr::arrange(period, country)
  
  vtp_rt <- summary_vectensor$tab_rt_all %>%
    dplyr::filter(country %in% country_order, period %in% period_order) %>%
    dplyr::mutate(
      country = factor(country, levels = country_order),
      period  = factor(period, levels = period_order)
    ) %>%
    dplyr::arrange(period, country)
  
  # ---------------------------------------------------------------------------
  # 2. Merge and compute relative RMSFE
  # ---------------------------------------------------------------------------
  df <- mat_rt %>%
    dplyr::rename(
      Matrix_M1 = M1,
      Matrix_M2 = M2,
      Matrix_M3 = M3
    ) %>%
    dplyr::left_join(
      vtp_rt %>%
        dplyr::rename(
          VecP_M1 = M1,
          VecP_M2 = M2,
          VecP_M3 = M3
        ),
      by = c("country", "period")
    ) %>%
    dplyr::left_join(
      vec_rt %>%
        dplyr::rename(
          VecC_M1 = M1,
          VecC_M2 = M2,
          VecC_M3 = M3
        ),
      by = c("country", "period")
    ) %>%
    dplyr::mutate(
      Rel_VecP_M1 = Matrix_M1 / VecP_M1,
      Rel_VecP_M2 = Matrix_M2 / VecP_M2,
      Rel_VecP_M3 = Matrix_M3 / VecP_M3,
      Rel_VecC_M1 = Matrix_M1 / VecC_M1,
      Rel_VecC_M2 = Matrix_M2 / VecC_M2,
      Rel_VecC_M3 = Matrix_M3 / VecC_M3
    )
  
  if (!is.null(dm_vecp_wide)) {
    df <- df %>%
      dplyr::left_join(dm_vecp_wide, by = c("country", "period"))
  }
  
  if (!is.null(dm_vecc_wide)) {
    df <- df %>%
      dplyr::left_join(dm_vecc_wide, by = c("country", "period"))
  }
  
  df <- df %>%
    dplyr::mutate(
      star_VecP_M1 = ifelse(Rel_VecP_M1 < 1, p_to_stars(p_VecP_M1), ""),
      star_VecP_M2 = ifelse(Rel_VecP_M2 < 1, p_to_stars(p_VecP_M2), ""),
      star_VecP_M3 = ifelse(Rel_VecP_M3 < 1, p_to_stars(p_VecP_M3), ""),
      star_VecC_M1 = ifelse(Rel_VecC_M1 < 1, p_to_stars(p_VecC_M1), ""),
      star_VecC_M2 = ifelse(Rel_VecC_M2 < 1, p_to_stars(p_VecC_M2), ""),
      star_VecC_M3 = ifelse(Rel_VecC_M3 < 1, p_to_stars(p_VecC_M3), "")
    )
  
  # ---------------------------------------------------------------------------
  # 3. Hyperparameters
  # ---------------------------------------------------------------------------
  hp_mat_pre  <- summary_matrix$hyper$pre
  hp_mat_post <- summary_matrix$hyper$post
  
  hp_vec_pre <- if (!is.null(summary_vector$hyper_rt_pre_all)) {
    summary_vector$hyper_rt_pre_all %>%
      dplyr::filter(country %in% country_order) %>%
      dplyr::mutate(country = factor(country, levels = country_order)) %>%
      dplyr::arrange(country)
  } else {
    NULL
  }
  
  hp_vec_post <- if (!is.null(summary_vector$hyper_rt_post_all)) {
    summary_vector$hyper_rt_post_all %>%
      dplyr::filter(country %in% country_order) %>%
      dplyr::mutate(country = factor(country, levels = country_order)) %>%
      dplyr::arrange(country)
  } else {
    NULL
  }
  
  hp_vtp_pre <- if (!is.null(summary_vectensor$hyper_rt_pre_all)) {
    summary_vectensor$hyper_rt_pre_all %>%
      dplyr::filter(country %in% country_order) %>%
      dplyr::mutate(country = factor(country, levels = country_order)) %>%
      dplyr::arrange(country)
  } else {
    NULL
  }
  
  hp_vtp_post <- if (!is.null(summary_vectensor$hyper_rt_post_all)) {
    summary_vectensor$hyper_rt_post_all %>%
      dplyr::filter(country %in% country_order) %>%
      dplyr::mutate(country = factor(country, levels = country_order)) %>%
      dplyr::arrange(country)
  } else {
    NULL
  }
  
  matrix_hp_pre_txt <- paste0(
    "$L=", hp_mat_pre$Lproxy,
    ",\\; P_f=", hp_mat_pre$L_midas - 1,
    ",\\; P_\\rho=", hp_mat_pre$p_AR,
    ",\\; (\\hat r_1,\\hat r_2)=(",
    hp_mat_pre$r_targeted[1], ",", hp_mat_pre$r_targeted[2], ")$"
  )
  
  matrix_hp_post_txt <- paste0(
    "$L=", hp_mat_post$Lproxy,
    ",\\; P_f=", hp_mat_post$L_midas - 1,
    ",\\; P_\\rho=", hp_mat_post$p_AR,
    ",\\; (\\hat r_1,\\hat r_2)=(",
    hp_mat_post$r_targeted[1], ",", hp_mat_post$r_targeted[2], ")$"
  )
  
  # ---------------------------------------------------------------------------
  # 4. Helpers
  # ---------------------------------------------------------------------------
  fmt <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
  
  fmt_star <- function(x, s = "") {
    if (length(x) == 0 || is.na(x)) return("")
    s <- ifelse(is.na(s) || s == "", "", paste0("\\sym{", s, "}"))
    paste0(sprintf("%.3f", x), s)
  }
  
  hp_cells <- function(hp_df, cc) {
    if (is.null(hp_df)) return(c("", "", "", ""))
    row <- hp_df %>% dplyr::filter(country == cc)
    if (nrow(row) == 0) return(c("", "", "", ""))
    c(
      as.character(row$Lproxy[1]),
      as.character(row$L_midas[1] - 1),
      as.character(row$p_AR[1]),
      as.character(row$r[1])
    )
  }
  
  block_rows <- function(period_name, period_label = NULL,
                         hp_matrix_txt = NULL, hp_vp = NULL, hp_vc = NULL) {
    out <- c()
    
    if (is.null(period_label)) {
      line1 <- "{\\fontsize{6.0}{6.4}\\selectfont\\textbf{Pre-COVID}}"
      line2 <- ""
    } else {
      line1 <- paste0(
        "{\\fontsize{6.0}{6.4}\\selectfont\\textbf{",
        period_label["title"],
        "}}"
      )
      line2 <- paste0(
        "{\\fontsize{4.6}{4.8}\\selectfont ",
        period_label["date"],
        "}"
      )
    }
    
    out <- c(
      out,
      "\\specialrule{0.08em}{0.15em}{0.08em}",
      paste0("\\multicolumn{18}{l}{", line1, "}\\\\[-0.60em]"),
      paste0("\\multicolumn{18}{l}{", line2, "}\\\\[-0.32em]"),
      "\\cmidrule{1-18}"
    )
    
    if (!is.null(hp_matrix_txt)) {
      out <- c(
        out,
        paste0(
          "\\rowcolor{hypergray}\n",
          "\\textit{\\fontsize{5.1}{6}\\selectfont Hyperparameters}",
          " & \\multicolumn{3}{@{}>{\\columncolor{hypergray}}c@{\\hspace{0.22cm}}}{\\fontsize{5.0}{5.8}\\selectfont ",
          hp_matrix_txt,
          "}",
          " & \\multicolumn{3}{c}{}",
          " & \\cellcolor{hypergray}\\hphdr $L$",
          " & \\cellcolor{hypergray}\\hphdr $P_f$",
          " & \\cellcolor{hypergray}\\hphdr $P_\\rho$",
          " & \\cellcolor{hypergray}\\hphdr $r$",
          " & \\multicolumn{3}{c}{}",
          " & \\cellcolor{hypergray}\\hphdr $L$",
          " & \\cellcolor{hypergray}\\hphdr $P_f$",
          " & \\cellcolor{hypergray}\\hphdr $P_\\rho$",
          " & \\cellcolor{hypergray}\\hphdr $r$",
          "\\\\[0.25em]"
        )
      )
    }
    
    df_p <- df %>% dplyr::filter(period == period_name)
    
    for (cc in country_order) {
      row <- df_p %>% dplyr::filter(country == cc)
      if (nrow(row) == 0) next
      
      hpvp <- if (!is.null(hp_vp)) hp_cells(hp_vp, cc) else c("", "", "", "")
      hpvc <- if (!is.null(hp_vc)) hp_cells(hp_vc, cc) else c("", "", "", "")
      
      out <- c(
        out,
        paste0(
          cc, " & ",
          "\\matcell{", fmt(row$Matrix_M1), "} & ",
          "\\matcell{", fmt(row$Matrix_M2), "} & ",
          "\\matcell{", fmt(row$Matrix_M3), "} & ",
          fmt_star(row$Rel_VecP_M1, row$star_VecP_M1), " & ",
          fmt_star(row$Rel_VecP_M2, row$star_VecP_M2), " & ",
          fmt_star(row$Rel_VecP_M3, row$star_VecP_M3), " & ",
          "\\hpstyle ", hpvp[1], " & ", "\\hpstyle ", hpvp[2], " & ", "\\hpstyle ", hpvp[3], " & ", "\\hpstyle ", hpvp[4], " & ",
          fmt_star(row$Rel_VecC_M1, row$star_VecC_M1), " & ",
          fmt_star(row$Rel_VecC_M2, row$star_VecC_M2), " & ",
          fmt_star(row$Rel_VecC_M3, row$star_VecC_M3), " & ",
          "\\hpstyle ", hpvc[1], " & ", "\\hpstyle ", hpvc[2], " & ", "\\hpstyle ", hpvc[3], " & ", "\\hpstyle ", hpvc[4], " \\\\"
        )
      )
    }
    
    out
  }
  
  size_caption <- paste0(tolower(Size), " info set")
  sel_caption  <- ifelse(
    toupper(sel) == "LASSO",
    "LASSO preselection",
    "correlation preselection"
  )
  
  table_note <- paste0(
    "The first shaded block reports the Matrix MF--TPRF RMSFE for M1, M2, and M3 nowcasts. ",
    "Relative RMSFE with respect to the pooled vector benchmark (VEC-P) and the country-specific vector benchmark (VEC-C) ",
    "are reported in the two blocks on the right; values below one favour Matrix MF--TPRF. ",
    "Asterisks denote one-sided Diebold--Mariano significance in favour of Matrix MF--TPRF, ",
    "based on squared forecast errors and a Newey--West HAC variance estimator over the common evaluation sample ",
    "($^{*}\\, p\\!<\\!0.10$, $^{**}\\, p\\!<\\!0.05$, $^{***}\\, p\\!<\\!0.01$)."
  )
  
  latex_lines <- c(
    "%=========================================================",
    paste0("% ", toupper(Size), " INFORMATION SET -- ", toupper(sel)),
    "%=========================================================",
    "\\begin{table}[htbp]",
    "\\centering",
    "\\scriptsize",
    "\\renewcommand{\\arraystretch}{1.08}",
    "\\setlength{\\tabcolsep}{1.7pt}",
    "",
    "\\newcommand{\\hpstyle}{\\fontsize{4.6}{5.2}\\selectfont}",
    "\\newcommand{\\hphdr}{\\fontsize{4.7}{5.4}\\selectfont}",
    "\\newcommand{\\matcell}[1]{\\cellcolor{matrixgray}#1}",
    "\\newcommand{\\sym}[1]{\\rlap{\\hspace{0.05em}\\textsuperscript{\\fontsize{4.0}{4.0}\\selectfont#1}}}",
    "",
    "\\definecolor{topgray}{gray}{0.86}",
    "\\definecolor{hypergray}{gray}{0.91}",
    "\\definecolor{matrixgray}{gray}{0.965}",
    "",
    paste0(
      "\\caption{Matrix MF--TPRF and relative RMSFE by country: ",
      size_caption, " \\& ", sel_caption, "}"
    ),
    paste0("\\label{tab:EA_nowcast_", tolower(Size), "_", tolower(sel), "}"),
    "",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{",
    ">{\\centering\\arraybackslash}p{1.15cm}",
    "@{\\hspace{0.22cm}}",
    ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.97cm}",
    ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.97cm}",
    ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.97cm}",
    "@{\\hspace{0.22cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{0.86cm}}",
    "*{4}{>{\\centering\\arraybackslash}p{0.14cm}}",
    "@{\\hspace{0.22cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{0.86cm}}",
    "*{4}{>{\\centering\\arraybackslash}p{0.14cm}}",
    "}",
    "\\toprule",
    paste0(
      "\\multicolumn{1}{>{\\cellcolor{topgray}}c@{\\hspace{0.22cm}}}{\\textbf{", toupper(Size), "}}",
      " & \\multicolumn{3}{@{}>{\\columncolor{matrixgray}}c@{\\hspace{0.22cm}}}{\\textbf{Matrix MF--TPRF}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-P MF--TPRF}}$}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-C MF--TPRF}}$}}"
    ),
    "\\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-11}\\cmidrule(lr){12-18}",
    "\\textbf{Country}",
    "& \\cellcolor{matrixgray}M1 & \\cellcolor{matrixgray}M2 & \\cellcolor{matrixgray}M3",
    "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
    "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
    "\\\\",
    block_rows("Pre-COVID",
               period_label = period_labels[["Pre-COVID"]],
               hp_matrix_txt = matrix_hp_pre_txt,
               hp_vp = hp_vtp_pre,
               hp_vc = hp_vec_pre),
    block_rows("COVID period",
               period_label = period_labels[["COVID period"]],
               hp_matrix_txt = NULL,
               hp_vp = hp_vtp_pre,
               hp_vc = hp_vec_pre),
    block_rows("Post-COVID",
               period_label = period_labels[["Post-COVID"]],
               hp_matrix_txt = matrix_hp_post_txt,
               hp_vp = hp_vtp_post,
               hp_vc = hp_vec_post),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.15cm}",
    paste0(
      "\\parbox{0.96\\textwidth}{\\footnotesize ",
      table_note,
      "}"
    ),
    "\\end{table}"
  )
  
  paste(latex_lines, collapse = "\n")
}

save_selection_wide_to_latex <- function(df_selection_wide,
                                         file,
                                         caption = "Variables selected by country",
                                         label = "tab:variables_selected_by_country",
                                         country_cols = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
                                         check_symbol = "\\checkmark") {
  
  df <- df_selection_wide
  
  # Controllo colonne minime
  required_cols <- c("base_name", "frequency", "model_size", country_cols)
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing columns in df_selection_wide: ",
         paste(missing_cols, collapse = ", "))
  }
  
  # Ordina
  df <- df[order(df$frequency, df$base_name), ]
  
  # Trasforma 1/0 in checkmark / vuoto
  for (cc in country_cols) {
    df[[cc]] <- ifelse(is.na(df[[cc]]) | df[[cc]] == 0, "", check_symbol)
  }
  
  # Header tabella
  align <- paste(c("l", "c", "c", rep("c", length(country_cols))), collapse = " ")
  
  header_1 <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\small\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{", align, "}\n",
    "\\toprule\n"
  )
  
  header_2 <- paste(
    c("\\textbf{Variable}", "\\textbf{Fr.}", "\\textbf{Size}",
      paste0("\\textbf{", country_cols, "}")),
    collapse = " & "
  )
  header_2 <- paste0(header_2, " \\\\ \n\\midrule\n")
  
  # Corpo tabella
  body <- apply(df[, c("base_name", "frequency", "model_size", country_cols)], 1, function(row) {
    paste0(paste(row, collapse = " & "), " \\\\")
  })
  
  footer <- paste0(
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  latex_code <- paste0(
    header_1,
    header_2,
    paste(body, collapse = "\n"),
    footer
  )
  
  writeLines(latex_code, con = file)
  
  invisible(latex_code)
}
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

normalize_selection_ids <- function(df) {
  df %>%
    dplyr::mutate(
      base_name = dplyr::recode(
        base_name,
        "TASS.LBD" = "TASS.LDB",
        "TLB.LBD"  = "TLB.LDB"
      )
    )
}

escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([&_#%$])", "\\\\\\1", x, perl = TRUE)
  x
}

combine_size_codes <- function(s, m, l) {
  out <- character(length(s))
  for (i in seq_along(s)) {
    parts <- c(
      if (isTRUE(s[i] == 1)) "S" else NULL,
      if (isTRUE(m[i] == 1)) "M" else NULL,
      if (isTRUE(l[i] == 1)) "L" else NULL
    )
    out[i] <- if (length(parts) == 0) "" else paste(parts, collapse = "/")
  }
  out
}

build_selection_size_table <- function(sel_small, sel_medium, sel_large,
                                       country_cols = c("DE","FR","IT","ES","NL","BE","AT","PT")) {
  
  key_cols <- c("base_name", "frequency")
  
  sel_small  <- normalize_selection_ids(sel_small)
  sel_medium <- normalize_selection_ids(sel_medium)
  sel_large  <- normalize_selection_ids(sel_large)
  
  s_tbl <- sel_small %>%
    dplyr::select(dplyr::all_of(key_cols), dplyr::all_of(country_cols)) %>%
    dplyr::rename_with(~ paste0(., "_S"), dplyr::all_of(country_cols))
  
  m_tbl <- sel_medium %>%
    dplyr::select(dplyr::all_of(key_cols), dplyr::all_of(country_cols)) %>%
    dplyr::rename_with(~ paste0(., "_M"), dplyr::all_of(country_cols))
  
  l_tbl <- sel_large %>%
    dplyr::select(dplyr::all_of(key_cols), dplyr::all_of(country_cols)) %>%
    dplyr::rename_with(~ paste0(., "_L"), dplyr::all_of(country_cols))
  
  out <- dplyr::full_join(s_tbl, m_tbl, by = key_cols) %>%
    dplyr::full_join(l_tbl, by = key_cols)
  
  for (cc in country_cols) {
    out[[cc]] <- combine_size_codes(
      out[[paste0(cc, "_S")]],
      out[[paste0(cc, "_M")]],
      out[[paste0(cc, "_L")]]
    )
  }
  
  out %>%
    dplyr::select(dplyr::all_of(key_cols), dplyr::all_of(country_cols)) %>%
    dplyr::arrange(frequency, base_name)
}

# ------------------------------------------------------------
# METADATA COMPLETO
# ------------------------------------------------------------
build_metadata_table <- function() {
  tribble(
    ~ID,          ~Series,                                                                 ~Cl, ~Cat, ~Tr, ~Fr, ~Del, ~Group,
    
    # (1) National Accounts / Real Economy
    "GDP",        "Real Gross Domestic Product",                                           "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "EXPGS",      "Real Export Goods and services",                                        "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "IMPGS",      "Real Import Goods and services",                                        "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFCE",       "Real Government Final consumption expenditure",                         "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "HFCE",       "Real Households consumption expenditure",                               "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSD",      "Real Households consumption expenditure: Durable Goods",                "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSSD",     "Real Households consumption expenditure: Semi-Durable Goods",           "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSSV",     "Real Households consumption expenditure: Services",                     "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSND",     "Real Households consumption expenditure: Non-Durable Goods",            "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GCF",        "Real Gross capital formation",                                          "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFCF",       "Real Gross fixed capital formation",                                    "R", "H", 1, "Q", 45, "National Accounts / Real Economy",  # alias usato nei tuoi selection
    "GFACON",     "Real Gross Fixed Capital Formation: Construction",                      "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFAMG",      "Real Gross Fixed Capital Formation: Machinery and Equipment",           "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "DFGDP",      "Real Gross Domestic Product Deflator",                                  "N", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "HPRC",       "Residential Property Prices (BIS)",                                     "N", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GNFCPS",     "Gross Profit Share of Non-Financial Corporations",                      "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GNFCIR",     "Gross Investment Share of Non-Financial Corporations",                  "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GHSR",       "Gross Households Savings Rate",                                         "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    
    # (2) Labor Market
    "TEMP",       "Total Employment (domestic concept)",                                   "R", "H", 1, "Q", 45, "Labor Market",
    "EMP",        "Employees (domestic concept)",                                          "R", "H", 1, "Q", 45, "Labor Market",
    "SEMP",       "Self Employment (domestic concept)",                                    "R", "H", 1, "Q", 45, "Labor Market",
    "THOURS",     "Hours Worked: Total",                                                   "R", "H", 1, "Q", 45, "Labor Market",
    "EMPAG",      "Quarterly Employment: Agriculture, Forestry, Fishing",                  "R", "H", 1, "Q", 45, "Labor Market",
    "EMPIN",      "Quarterly Employment: Industry",                                        "R", "H", 1, "Q", 45, "Labor Market",
    "EMPMN",      "Quarterly Employment: Manufacturing",                                   "R", "H", 1, "Q", 45, "Labor Market",
    "EMPCON",     "Quarterly Employment: Construction",                                    "R", "H", 1, "Q", 45, "Labor Market",
    "EMPRT",      "Quarterly Employment: Wholesale/Retail trade, transport, food",         "R", "H", 1, "Q", 45, "Labor Market",
    "EMPIT",      "Quarterly Employment: Information and Communication",                   "R", "H", 1, "Q", 45, "Labor Market",
    "EMPFC",      "Quarterly Employment: Financial and Insurance activities",              "R", "H", 1, "Q", 45, "Labor Market",
    "EMPRE",      "Quarterly Employment: Real Estate",                                     "R", "H", 1, "Q", 45, "Labor Market",
    "EMPPR",      "Quarterly Employment: Professional, Scientific, Technical activities",  "R", "H", 1, "Q", 45, "Labor Market",
    "EMPPA",      "Quarterly Employment: PA, education, health and social services",       "R", "H", 1, "Q", 45, "Labor Market",
    "EMPENT",     "Quarterly Employment: Arts and recreational activities",                "R", "H", 1, "Q", 45, "Labor Market",
    "UNETOT",     "Unemployment: Total",                                                   "R", "H", 0, "M", 45, "Labor Market",
    "UNEO25",     "Unemployment: Over 25 years",                                           "R", "H", 0, "M", 45, "Labor Market",
    "UNEU25",     "Unemployment: Under 25 years",                                          "R", "H", 0, "M", 45, "Labor Market",
    "RPRP",       "Real Labour Productivity (person)",                                     "R", "H", 1, "Q", 45, "Labor Market",
    "WS",         "Wages and salaries",                                                    "N", "H", 1, "Q", 45, "Labor Market",
    "ESC",        "Employers' Social Contributions",                                       "N", "H", 1, "Q", 45, "Labor Market",
    
    # (3) Credit Aggregates
    "TASS.SDB",   "Total Economy - Assets: Short-Term Debt Securities",                    "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.LDB",   "Total Economy - Assets: Long-Term Debt Securities",                     "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.SLN",   "Total Economy - Assets: Short-Term Loans",                              "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.LLN",   "Total Economy - Assets: Long-Term Loans",                               "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.SDB",    "Total Economy - Liabilities: Short-Term Debt Securities",               "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.LDB",    "Total Economy - Liabilities: Long-Term Debt Securities",                "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.SLN",    "Total Economy - Liabilities: Short-Term Loans",                         "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.LLN",    "Total Economy - Liabilities: Long-Term Loans",                          "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB.SLN",  "Non-Financial Corporations - Liabilities - Short-Term Loans",           "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB.LLN",  "Non-Financial Corporations - Liabilities - Long-Term Loans",            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS",      "General Government: Total Financial Assets",                            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.SLN",  "General Government - Assets: Short-Term Loans",                         "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.LLN",  "General Government - Assets: Long-Term Loans",                          "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB",       "General Government: Total Financial Liabilities",                       "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.SLN",   "General Government - Liabilities: Short-Term Loans",                    "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.LLN",   "General Government - Liabilities: Long-Term Loans",                     "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB",       "Households: Total Financial Liabilities",                               "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.SLN",   "Households - Liabilities: Short-Term Loans",                            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.LLN",   "Households - Liabilities: Long-Term Loans",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    
    # (4) Labor Costs
    "ULCCON",     "Nominal Unit Labor Costs: Construction",                                "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCIN",      "Nominal Unit Labor Costs: Industry",                                    "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCMN",      "Nominal Unit Labor Costs: Manufacturing",                               "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCFC",      "Nominal Unit Labor Costs: Financial Activities",                        "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRE",      "Nominal Unit Labor Costs: Real Estate",                                 "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCPR",      "Nominal Unit Labor Costs: Professional, Scientific, Technical activities","N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRT",      "Nominal Unit Labor Costs: Wholesale/Retail Trade, Transport, Food, IT", "N", "H", 1, "Q", 45, "Labor Costs",
    
    # (5) Financial Markets
    "REER42",     "Real Exchange Rate (42 main industrial countries)",                     "F", "H", 1, "M", 35, "Financial Markets",
    "SHIX",       "Stock Price Index",                                                     "F", "S", 1, "M", 1,  "Financial Markets",
    
    # (6) Interest Rates
    "LTIRT",      "Long-Term Interest Rates (EMU Criterion)",                              "F", "H", 2, "M", 35, "Interest Rates",
    
    # (7) Industrial Production and Turnover
    "IPMN",       "Industrial Production Index: Manufacturing",                            "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPCAG",      "Industrial Production Index: Capital Goods",                            "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPCOG",      "Industrial Production Index: Consumer Goods",                           "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPDCOG",     "Industrial Production Index: Durable Consumer Goods",                   "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPNDCOG",    "Industrial Production Index: Non Durable Consumer Goods",               "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPING",      "Industrial Production Index: Intermediate Goods",                       "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPNRG",      "Industrial Production Index: Energy",                                   "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNMN",      "Turnover Index: Manufacturing",                                         "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNCAG",     "Turnover Index: Capital Goods",                                         "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNCOG",     "Turnover Index: Consumer Goods",                                        "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNDCOG",    "Turnover Index: Durable Consumer Goods",                                "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNNDCOG",   "Turnover Index: Non Durable Consumer Goods",                            "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNING",     "Turnover Index: Intermediate Goods",                                    "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNNRG",     "Turnover Index: Energy",                                                "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    
    # (8) Prices
    "PPICAG",     "Producer Price Index: Capital Goods",                                   "N", "H", 1, "M", 40, "Prices",
    "PPICOG",     "Producer Price Index: Consumer Goods",                                  "N", "H", 1, "M", 40, "Prices",
    "PPIDCOG",    "Producer Price Index: Durable Consumer Goods",                          "N", "H", 1, "M", 40, "Prices",
    "PPINDCOG",   "Producer Price Index: Non Durable Consumer Goods",                      "N", "H", 1, "M", 40, "Prices",
    "PPIING",     "Producer Price Index: Intermediate Goods",                              "N", "H", 1, "M", 40, "Prices",
    "PPINRG",     "Producer Price Index: Energy",                                          "N", "H", 1, "M", 40, "Prices",
    "HICPOV",     "Harmonized Index of Consumer Prices: Overall Index",                    "N", "H", 1, "M", 40, "Prices",
    "HICPNEF",    "Harmonized Index of Consumer Prices: All Items, no Energy/Food",        "N", "H", 1, "M", 40, "Prices",
    "HICPG",      "Harmonized Index of Consumer Prices: Goods",                            "N", "H", 1, "M", 40, "Prices",
    "HICPIN",     "Harmonized Index of Consumer Prices: Industrial Goods",                 "N", "H", 1, "M", 40, "Prices",
    "HICPSV",     "Harmonized Index of Consumer Prices: Services",                         "N", "H", 1, "M", 40, "Prices",
    "HICPNG",     "Harmonized Index of Consumer Prices: Energy",                           "N", "H", 1, "M", 40, "Prices",
    
    # (9) Confidence Indicators
    "ICONFIX",    "Industrial Confidence Indicator",                                       "C", "S", 0, "M", 5,  "Confidence Indicators",
    "CCONFIX",    "Consumer Confidence Indicator",                                         "C", "S", 0, "M", 5,  "Confidence Indicators",
    "KCONFIX",    "Construction Confidence Indicator",                                     "C", "S", 0, "M", 5,  "Confidence Indicators",
    "SCONFIX",    "Services Confidence Indicator",                                         "C", "S", 0, "M", 5,  "Confidence Indicators",
    "ESENTIX",    "Economic Sentiment Indicator",                                          "C", "S", 0, "M", 5,  "Confidence Indicators",
    "RTCONFIX",   "Retail Confidence Indicator",                                           "C", "S", 0, "M", 5,  "Confidence Indicators",
    "BCI",        "Business Confidence Index",                                             "C", "S", 1, "M", 5,  "Confidence Indicators",
    "CCI",        "Consumer Confidence Index",                                             "C", "S", 1, "M", 5,  "Confidence Indicators"
  )
}

to_setcell <- function(x) {
  ifelse(is.na(x) | x == "", "", paste0("\\setcell{", x, "}"))
}

build_latex_selection_table <- function(sel_small, sel_medium, sel_large,
                                        metadata = build_metadata_table(),
                                        country_cols = c("DE","FR","IT","ES","NL","BE","AT","PT"),
                                        caption = "Macroeconomic predictors: country-specific selection across information-set sizes",
                                        label = "tab:variables_selected_country_size") {
  
  sel_all <- build_selection_size_table(
    sel_small  = sel_small,
    sel_medium = sel_medium,
    sel_large  = sel_large,
    country_cols = country_cols
  ) %>%
    rename(ID = base_name, Fr = frequency) %>%
    filter(ID != "GDP")
  
  # controllo diagnostico: se qualcosa manca nel metadata, fermati subito
  missing_in_meta <- setdiff(sel_all$ID, metadata$ID)
  if (length(missing_in_meta) > 0) {
    stop(
      "Missing metadata for these IDs: ",
      paste(sort(missing_in_meta), collapse = ", ")
    )
  }
  
  tab <- metadata %>%
    filter(ID != "GDP") %>%
    inner_join(sel_all, by = c("ID", "Fr")) %>%
    mutate(across(all_of(country_cols), to_setcell)) %>%
    mutate(
      ID     = escape_latex(ID),
      Series = escape_latex(Series),
      Cl     = escape_latex(Cl),
      Cat    = escape_latex(Cat),
      Fr     = escape_latex(Fr)
    )
  
  groups_order <- c(
    "National Accounts / Real Economy",
    "Labor Market",
    "Credit Aggregates",
    "Labor Costs",
    "Financial Markets",
    "Interest Rates",
    "Industrial Production and Turnover",
    "Prices",
    "Confidence Indicators"
  )
  
  tab$Group <- factor(tab$Group, levels = groups_order)
  
  tab <- tab %>%
    arrange(Group, Fr, ID) %>%
    mutate(N = row_number())
  
  header <- paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\tiny\n",
    "\\renewcommand{\\arraystretch}{0.82}\n",
    "\\setlength{\\tabcolsep}{1.6pt}\n",
    "\\newcommand{\\setcell}[1]{{\\fontsize{4.1}{4.4}\\selectfont #1}}\n\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "{\\fontsize{4.0}{4.4}\\selectfont\n",
    "\\begin{tabular}{c p{1.25cm} p{5.6cm} c c c c c c c c c c c c c}\n",
    "\\toprule\n",
    "\\textbf{N} & \\textbf{ID} & \\textbf{Series} & \\textbf{Cl.} & \\textbf{Cat.} & \\textbf{Tr.} & \\textbf{Fr.} & \\textbf{Del.} & ",
    paste0("\\textbf{", country_cols, "}", collapse = " & "),
    " \\\\\n",
    "\\midrule\n",
    "\\midrule\n"
  )
  
  body_lines <- c()
  
  for (grp in groups_order) {
    subtab <- tab %>% filter(Group == grp)
    if (nrow(subtab) == 0) next
    
    ncols_total <- 8 + length(country_cols)
    grp_line <- paste0("\\multicolumn{", ncols_total, "}{c}{\\textbf{", grp, "}} \\\\")
    
    rows <- subtab %>%
      select(N, ID, Series, Cl, Cat, Tr, Fr, Del, all_of(country_cols))
    
    row_lines <- apply(rows, 1, function(r) paste0(paste(r, collapse = " & "), " \\\\"))
    
    body_lines <- c(
      body_lines,
      grp_line,
      "\\midrule",
      row_lines,
      "\\midrule"
    )
  }
  
  if (length(body_lines) > 0) {
    body_lines <- body_lines[-length(body_lines)]
  }
  
  footer <- paste0(
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "}\n",
    "}\n",
    "\\parbox{0.98\\textwidth}{\\tiny\n",
    "\\textit{Notes:} Country entries report the information sets in which the predictor is selected: ",
    "\\textbf{S} = small, \\textbf{M} = medium, \\textbf{L} = large. ",
    "Combined labels such as \\textbf{S/M}, \\textbf{S/L}, \\textbf{M/L}, and \\textbf{S/M/L} indicate that the predictor is selected in more than one information set for the corresponding country. ",
    "\\textbf{Cl.} denotes the variable class, where $R$ indicates real activity variables, $N$ nominal variables, $C$ confidence indicators, and $F$ financial variables. ",
    "\\textbf{Cat.} distinguishes hard ($H$) from soft ($S$) indicators. ",
    "\\textbf{Tr.} reports the transformation code, \\textbf{Fr.} the sampling frequency, and \\textbf{Del.} the approximate release delay in days. ",
    "The target variable GDP is excluded from the table.}\n",
    "\\end{table}\n"
  )
  
  paste0(header, paste(body_lines, collapse = "\n"), footer)
}

# ==============================================================================
# DM TEST 
# ==============================================================================

make_quarter_id <- function(x) {
  paste0(year(x), "Q", quarter(x))
}

build_eval_df <- function(df_rt, df_y) {
  df_rt %>%
    mutate(
      quarter_id = make_quarter_id(date)
    ) %>%
    left_join(
      df_y %>%
        select(country, quarter_id, GDP, period),
      by = c("country", "quarter_id")
    ) %>%
    mutate(
      error = GDP - nowcast,
      se    = error^2
    )
}

p_to_stars <- function(p) {
  case_when(
    is.na(p)   ~ "",
    p < 0.01   ~ "***",
    p < 0.05   ~ "**",
    p < 0.10   ~ "*",
    TRUE       ~ ""
  )
}

fmt_star <- function(x, s = "") {
  if (length(x) == 0 || is.na(x)) return("")
  s <- ifelse(is.na(s), "", s)
  paste0(sprintf("%.3f", x), s)
}

dm_test_hac <- function(d, lag = NULL, alternative = c("less", "two.sided", "greater")) {
  alternative <- match.arg(alternative)
  d <- d[is.finite(d)]
  n <- length(d)
  
  if (n < 4) {
    return(data.frame(
      n = n,
      mean_d = NA_real_,
      stat = NA_real_,
      p_value = NA_real_
    ))
  }
  
  if (is.null(lag)) {
    lag <- floor(n^(1/3))
  }
  
  fit <- lm(d ~ 1)
  vc  <- sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = TRUE)
  
  mean_d  <- coef(fit)[1]
  se_mean <- sqrt(vc[1, 1])
  stat    <- mean_d / se_mean
  
  p_value <- switch(
    alternative,
    "less"      = pnorm(stat),
    "greater"   = 1 - pnorm(stat),
    "two.sided" = 2 * pnorm(-abs(stat))
  )
  
  data.frame(
    n = n,
    mean_d = mean_d,
    stat = stat,
    p_value = p_value
  )
}