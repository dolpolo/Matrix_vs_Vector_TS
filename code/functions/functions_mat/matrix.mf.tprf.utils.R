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
# utils_results_01_core.R
# Core utilities for final results
# ==============================================================================

# ==============================================================================
# 0. REQUIRED PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(lubridate)
library(sandwich)

# ==============================================================================
# 1. GLOBAL ORDERS AND DEFAULT LABELS
# ==============================================================================

country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")

period_order <- c(
  "Pre-COVID",
  "COVID period",
  "Post-COVID"
)

period_order_full <- c(
  "Full sample",
  "Pre-COVID",
  "COVID period",
  "Post-COVID"
)

month_order <- c("M1", "M2", "M3")

country_labels <- c(
  "DE" = "Germany",
  "FR" = "France",
  "IT" = "Italy",
  "ES" = "Spain",
  "NL" = "Netherlands",
  "BE" = "Belgium",
  "AT" = "Austria",
  "PT" = "Portugal"
)

model_labels <- c(
  matrix    = "Matrix MF-TPRF",
  vector    = "VEC-C",
  vectensor = "VEC-P",
  dfm       = "DFM"
)

# ==============================================================================
# 2. SIZE AND RUN IDENTIFIERS
# ==============================================================================

get_size_tag <- function(N_m, N_q) {
  total <- N_m + N_q
  
  if (total <= 25) {
    return("small")
  }
  
  if (total <= 50) {
    return("medium")
  }
  
  "large"
}

build_run_tag <- function(params) {
  required_names <- c(
    "start_eval", "end_eval", "sel_method",
    "n_m", "n_q", "Kmax_vec", "Zmax", "Lmax",
    "p_AR_max", "covid_mask_m", "covid_mask_q"
  )
  
  missing_names <- setdiff(required_names, names(params))
  if (length(missing_names) > 0) {
    stop(
      "Missing entries in params: ",
      paste(missing_names, collapse = ", ")
    )
  }
  
  paste0(
    "eval-", format(params$start_eval, "%Y%m"),
    "_", format(params$end_eval, "%Y%m"),
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

# ==============================================================================
# 3. FILE NAME BUILDERS AND FILE SEARCH
# ==============================================================================

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
    cc_lab <- if (!is.null(countries) && length(countries) > 0L) {
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

build_saved_filename <- function(stage,
                                 model,
                                 Size,
                                 sel,
                                 ext = "rds",
                                 timestamp = TRUE) {
  
  base_name <- paste0(
    stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel
  )
  
  if (isTRUE(timestamp)) {
    base_name <- paste0(
      base_name,
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  }
  
  paste0(base_name, ".", ext)
}

find_result_file <- function(path,
                             model,
                             stage,
                             Size,
                             sel,
                             ext = "rds") {
  
  pattern <- paste0(
    "^", stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel,
    ".*\\.", ext, "$"
  )
  
  files <- list.files(
    path       = path,
    pattern    = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0L) {
    stop(
      "No file found for model = ", model,
      ", stage = ", stage,
      ", Size = ", Size,
      ", sel = ", sel,
      ", path = ", path
    )
  }
  
  if (length(files) > 1L) {
    message("Multiple files found. Using the most recent one.")
    files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  }
  
  files[1]
}

find_dfm_summary_file <- function(path,
                                  Size,
                                  sel,
                                  ext = "rds") {
  
  pattern <- paste0(
    "^VECTOR_DFM_summary_Size-", Size,
    "_sel-", sel,
    ".*\\.", ext, "$"
  )
  
  files <- list.files(
    path       = path,
    pattern    = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0L) {
    stop(
      "No DFM summary file found for Size = ",
      Size,
      ", sel = ",
      sel,
      ", path = ",
      path
    )
  }
  
  if (length(files) > 1L) {
    message("Multiple DFM files found. Using the most recent one.")
    files <- files[order(file.info(files)$mtime, decreasing = TRUE)]
  }
  
  files[1]
}

# ==============================================================================
# 4. DATE AND PERIOD HELPERS
# ==============================================================================

format_month_year <- function(x) {
  format(as.Date(x), "%b %Y")
}

build_period_labels <- function(params) {
  required_names <- c("start_eval", "end_eval", "covid_start", "covid_end")
  
  missing_names <- setdiff(required_names, names(params))
  if (length(missing_names) > 0) {
    stop(
      "Missing entries in params: ",
      paste(missing_names, collapse = ", ")
    )
  }
  
  pre_end <- seq(
    from       = params$covid_start,
    by         = "-1 month",
    length.out = 2
  )[2]
  
  post_start <- seq(
    from       = params$covid_end,
    by         = "+1 month",
    length.out = 2
  )[2]
  
  list(
    "Pre-COVID" = c(
      title = "Pre-COVID",
      date  = paste0(
        format_month_year(params$start_eval),
        " -- ",
        format_month_year(pre_end)
      )
    ),
    "COVID period" = c(
      title = "COVID",
      date  = paste0(
        format_month_year(params$covid_start),
        " -- ",
        format_month_year(params$covid_end)
      )
    ),
    "Post-COVID" = c(
      title = "Post-COVID",
      date  = paste0(
        format_month_year(post_start),
        " -- ",
        format_month_year(params$end_eval)
      )
    )
  )
}

make_quarter_id <- function(x) {
  paste0(lubridate::year(x), "Q", lubridate::quarter(x))
}

# ==============================================================================
# 5. FORMATTERS
# ==============================================================================

escape_latex <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

fmt_num <- function(x, digits = 3, missing = "") {
  ifelse(
    is.na(x),
    missing,
    sprintf(paste0("%.", digits, "f"), x)
  )
}

fmt_int <- function(x, missing = "--") {
  ifelse(
    is.na(x),
    missing,
    as.character(as.integer(round(x)))
  )
}

fmt_pct <- function(x, digits = 1, missing = "--") {
  ifelse(
    is.na(x),
    missing,
    format(round(x, digits), nsmall = digits, trim = TRUE)
  )
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.01  ~ "***",
    p < 0.05  ~ "**",
    p < 0.10  ~ "*",
    TRUE      ~ ""
  )
}

fmt_star <- function(x, s = "", digits = 3, latex_sym = TRUE) {
  if (length(x) == 0L || is.na(x)) {
    return("")
  }
  
  s <- ifelse(is.na(s), "", s)
  
  if (isTRUE(latex_sym) && nzchar(s)) {
    s <- paste0("\\sym{", s, "}")
  }
  
  paste0(sprintf(paste0("%.", digits, "f"), x), s)
}

bold_if_good_share <- function(x, digits = 1, threshold = 50) {
  val <- fmt_pct(x, digits = digits)
  
  ifelse(
    !is.na(x) & x > threshold,
    paste0("\\textbf{", val, "}"),
    val
  )
}

fmt_small_share <- function(x, digits = 1, threshold = 50) {
  val <- fmt_pct(x, digits = digits)
  
  ifelse(
    !is.na(x) & x > threshold,
    paste0("\\textbf{", val, "}"),
    val
  )
}

# ==============================================================================
# 6. LIST AND NOWCAST CONVERSION HELPERS
# ==============================================================================

list_to_df_country <- function(lst_by_country, tag) {
  if (is.null(lst_by_country) || length(lst_by_country) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      country          = character(),
      nowcast          = numeric(),
      month_in_quarter = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  out <- lapply(names(lst_by_country), function(cc) {
    x <- lst_by_country[[cc]]
    
    if (is.null(x) || length(x) == 0L) {
      return(NULL)
    }
    
    data.frame(
      date             = as.Date(names(x)),
      country          = cc,
      nowcast          = as.numeric(unlist(x)),
      month_in_quarter = tag,
      row.names        = NULL,
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out)
}

list_to_df_nowcast <- function(lst, tag) {
  if (is.null(lst) || length(lst) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      nowcast          = numeric(),
      month_in_quarter = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    date             = as.Date(names(lst)),
    nowcast          = as.numeric(unlist(lst)),
    month_in_quarter = tag,
    row.names        = NULL,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 7. RMSFE HELPERS
# ==============================================================================

rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q, na.rm = TRUE)) {
    return(NA_real_)
  }
  
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

# ==============================================================================
# 8. SUMMARY TABLE NORMALIZATION
# ==============================================================================

normalize_summary_table <- function(df,
                                    model_name,
                                    default_period = "Full sample") {
  
  if (!is.data.frame(df)) {
    stop("Input object must be a data.frame.")
  }
  
  required_basic <- c("country", "M1", "M2", "M3")
  missing_basic  <- setdiff(required_basic, names(df))
  
  if (length(missing_basic) > 0L) {
    stop(
      "Summary table is missing columns: ",
      paste(missing_basic, collapse = ", ")
    )
  }
  
  if (!"period" %in% names(df)) {
    df$period <- default_period
  }
  
  df %>%
    dplyr::mutate(
      model   = model_name,
      country = as.character(country),
      period  = as.character(period)
    ) %>%
    dplyr::select(country, period, M1, M2, M3, model)
}

# ==============================================================================
# 9. COMPARISON TABLES
# ==============================================================================

build_comparison_table <- function(df_matrix,
                                   df_vector,
                                   df_vectensor,
                                   df_dfm = NULL) {
  
  df_all <- dplyr::bind_rows(
    df_matrix,
    df_vector,
    df_vectensor
  )
  
  if (!is.null(df_dfm)) {
    df_all <- dplyr::bind_rows(df_all, df_dfm)
  }
  
  df_all %>%
    tidyr::pivot_longer(
      cols      = c("M1", "M2", "M3"),
      names_to  = "horizon",
      values_to = "RMSFE"
    ) %>%
    dplyr::mutate(
      col_id = paste0(model, "_", horizon)
    ) %>%
    dplyr::select(country, period, col_id, RMSFE) %>%
    tidyr::pivot_wider(
      names_from  = col_id,
      values_from = RMSFE
    ) %>%
    dplyr::arrange(
      country,
      factor(period, levels = period_order_full)
    )
}

comparison_to_latex <- function(df_comp,
                                caption,
                                label,
                                include_dfm = TRUE,
                                digits = 4) {
  
  needed <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  if (isTRUE(include_dfm)) {
    needed <- c(needed, "DFM_M1", "DFM_M2", "DFM_M3")
  }
  
  for (nm in needed) {
    if (!nm %in% names(df_comp)) {
      df_comp[[nm]] <- NA_real_
    }
  }
  
  fmt <- function(x) {
    ifelse(
      is.na(x),
      "",
      sprintf(paste0("%.", digits, "f"), x)
    )
  }
  
  if (isTRUE(include_dfm)) {
    tabular <- "llccc ccc ccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{Matrix MF--TPRF}",
      " & \\multicolumn{3}{c}{VEC-C}",
      " & \\multicolumn{3}{c}{VEC-P}",
      " & \\multicolumn{3}{c}{DFM} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df_comp$country,
        df_comp$period,
        fmt(df_comp$Matrix_M1), fmt(df_comp$Matrix_M2), fmt(df_comp$Matrix_M3),
        fmt(df_comp$Vector_M1), fmt(df_comp$Vector_M2), fmt(df_comp$Vector_M3),
        fmt(df_comp$VecTensor_M1), fmt(df_comp$VecTensor_M2), fmt(df_comp$VecTensor_M3),
        fmt(df_comp$DFM_M1), fmt(df_comp$DFM_M2), fmt(df_comp$DFM_M3)
      ),
      collapse = "\n"
    )
  } else {
    tabular <- "llccc ccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{Matrix MF--TPRF}",
      " & \\multicolumn{3}{c}{VEC-C}",
      " & \\multicolumn{3}{c}{VEC-P} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df_comp$country,
        df_comp$period,
        fmt(df_comp$Matrix_M1), fmt(df_comp$Matrix_M2), fmt(df_comp$Matrix_M3),
        fmt(df_comp$Vector_M1), fmt(df_comp$Vector_M2), fmt(df_comp$Vector_M3),
        fmt(df_comp$VecTensor_M1), fmt(df_comp$VecTensor_M2), fmt(df_comp$VecTensor_M3)
      ),
      collapse = "\n"
    )
  }
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header_1,
    header_2,
    "\\midrule\n",
    rows,
    "\n\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# 10. RELATIVE RMSFE TABLES
# ==============================================================================

build_relative_table <- function(df_comp,
                                 include_dfm = TRUE) {
  
  needed <- c(
    "country", "period",
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  if (isTRUE(include_dfm)) {
    needed <- c(needed, "DFM_M1", "DFM_M2", "DFM_M3")
  }
  
  missing_cols <- setdiff(needed, names(df_comp))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in comparison table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df_comp %>%
    dplyr::mutate(
      rel_MV_M1 = Matrix_M1 / Vector_M1,
      rel_MV_M2 = Matrix_M2 / Vector_M2,
      rel_MV_M3 = Matrix_M3 / Vector_M3,
      
      rel_MT_M1 = Matrix_M1 / VecTensor_M1,
      rel_MT_M2 = Matrix_M2 / VecTensor_M2,
      rel_MT_M3 = Matrix_M3 / VecTensor_M3,
      
      rel_VT_M1 = Vector_M1 / VecTensor_M1,
      rel_VT_M2 = Vector_M2 / VecTensor_M2,
      rel_VT_M3 = Vector_M3 / VecTensor_M3
    )
  
  if (isTRUE(include_dfm)) {
    out <- out %>%
      dplyr::mutate(
        rel_MD_M1 = Matrix_M1 / DFM_M1,
        rel_MD_M2 = Matrix_M2 / DFM_M2,
        rel_MD_M3 = Matrix_M3 / DFM_M3
      ) %>%
      dplyr::select(
        country, period,
        rel_MV_M1, rel_MV_M2, rel_MV_M3,
        rel_MT_M1, rel_MT_M2, rel_MT_M3,
        rel_VT_M1, rel_VT_M2, rel_VT_M3,
        rel_MD_M1, rel_MD_M2, rel_MD_M3
      )
  } else {
    out <- out %>%
      dplyr::select(
        country, period,
        rel_MV_M1, rel_MV_M2, rel_MV_M3,
        rel_MT_M1, rel_MT_M2, rel_MT_M3,
        rel_VT_M1, rel_VT_M2, rel_VT_M3
      )
  }
  
  out %>%
    dplyr::arrange(
      country,
      factor(period, levels = period_order_full)
    )
}

relative_to_latex <- function(df_rel,
                              caption,
                              label,
                              include_dfm = TRUE,
                              digits = 3) {
  
  fmt <- function(x) {
    ifelse(
      is.na(x),
      "",
      sprintf(paste0("%.", digits, "f"), x)
    )
  }
  
  if (isTRUE(include_dfm)) {
    tabular <- "llccc ccc ccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{Matrix / VEC-C}",
      " & \\multicolumn{3}{c}{Matrix / VEC-P}",
      " & \\multicolumn{3}{c}{VEC-C / VEC-P}",
      " & \\multicolumn{3}{c}{Matrix / DFM} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df_rel$country,
        df_rel$period,
        fmt(df_rel$rel_MV_M1), fmt(df_rel$rel_MV_M2), fmt(df_rel$rel_MV_M3),
        fmt(df_rel$rel_MT_M1), fmt(df_rel$rel_MT_M2), fmt(df_rel$rel_MT_M3),
        fmt(df_rel$rel_VT_M1), fmt(df_rel$rel_VT_M2), fmt(df_rel$rel_VT_M3),
        fmt(df_rel$rel_MD_M1), fmt(df_rel$rel_MD_M2), fmt(df_rel$rel_MD_M3)
      ),
      collapse = "\n"
    )
  } else {
    tabular <- "llccc ccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{Matrix / VEC-C}",
      " & \\multicolumn{3}{c}{Matrix / VEC-P}",
      " & \\multicolumn{3}{c}{VEC-C / VEC-P} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df_rel$country,
        df_rel$period,
        fmt(df_rel$rel_MV_M1), fmt(df_rel$rel_MV_M2), fmt(df_rel$rel_MV_M3),
        fmt(df_rel$rel_MT_M1), fmt(df_rel$rel_MT_M2), fmt(df_rel$rel_MT_M3),
        fmt(df_rel$rel_VT_M1), fmt(df_rel$rel_VT_M2), fmt(df_rel$rel_VT_M3)
      ),
      collapse = "\n"
    )
  }
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header_1,
    header_2,
    "\\midrule\n",
    rows,
    "\n\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# 10.5 VEC-C BENCHMARK RMSFE TABLE
# ==============================================================================

build_vecc_benchmark_table <- function(df_comp,
                                       include_dfm = TRUE) {
  
  needed <- c(
    "country", "period",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  if (isTRUE(include_dfm)) {
    needed <- c(needed, "DFM_M1", "DFM_M2", "DFM_M3")
  }
  
  missing_cols <- setdiff(needed, names(df_comp))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in comparison table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df_comp %>%
    dplyr::mutate(
      vecc_M1 = Vector_M1,
      vecc_M2 = Vector_M2,
      vecc_M3 = Vector_M3,
      
      rel_VecC_VecP_M1 = Vector_M1 / VecTensor_M1,
      rel_VecC_VecP_M2 = Vector_M2 / VecTensor_M2,
      rel_VecC_VecP_M3 = Vector_M3 / VecTensor_M3
    )
  
  if (isTRUE(include_dfm)) {
    out <- out %>%
      dplyr::mutate(
        rel_VecC_DFM_M1 = Vector_M1 / DFM_M1,
        rel_VecC_DFM_M2 = Vector_M2 / DFM_M2,
        rel_VecC_DFM_M3 = Vector_M3 / DFM_M3
      ) %>%
      dplyr::select(
        country, period,
        vecc_M1, vecc_M2, vecc_M3,
        rel_VecC_VecP_M1, rel_VecC_VecP_M2, rel_VecC_VecP_M3,
        rel_VecC_DFM_M1, rel_VecC_DFM_M2, rel_VecC_DFM_M3
      )
  } else {
    out <- out %>%
      dplyr::select(
        country, period,
        vecc_M1, vecc_M2, vecc_M3,
        rel_VecC_VecP_M1, rel_VecC_VecP_M2, rel_VecC_VecP_M3
      )
  }
  
  out %>%
    dplyr::arrange(
      country,
      factor(period, levels = period_order_full)
    )
}

vecc_benchmark_to_latex <- function(df,
                                    caption,
                                    label,
                                    include_dfm = TRUE,
                                    digits_rmsfe = 4,
                                    digits_rel = 3) {
  
  fmt_rmsfe <- function(x) {
    ifelse(
      is.na(x),
      "",
      sprintf(paste0("%.", digits_rmsfe, "f"), x)
    )
  }
  
  fmt_rel <- function(x) {
    ifelse(
      is.na(x),
      "",
      sprintf(paste0("%.", digits_rel, "f"), x)
    )
  }
  
  if (isTRUE(include_dfm)) {
    
    tabular <- "llccc ccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{VEC-C RMSFE}",
      " & \\multicolumn{3}{c}{VEC-C / VEC-P}",
      " & \\multicolumn{3}{c}{VEC-C / DFM} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df$country,
        df$period,
        fmt_rmsfe(df$vecc_M1), fmt_rmsfe(df$vecc_M2), fmt_rmsfe(df$vecc_M3),
        fmt_rel(df$rel_VecC_VecP_M1), fmt_rel(df$rel_VecC_VecP_M2), fmt_rel(df$rel_VecC_VecP_M3),
        fmt_rel(df$rel_VecC_DFM_M1), fmt_rel(df$rel_VecC_DFM_M2), fmt_rel(df$rel_VecC_DFM_M3)
      ),
      collapse = "\n"
    )
    
  } else {
    
    tabular <- "llccc ccc"
    
    header_1 <- paste0(
      " & & \\multicolumn{3}{c}{VEC-C RMSFE}",
      " & \\multicolumn{3}{c}{VEC-C / VEC-P} \\\\\n"
    )
    
    header_2 <- paste0(
      "Country & Period",
      " & M1 & M2 & M3",
      " & M1 & M2 & M3 \\\\\n"
    )
    
    rows <- paste(
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s & %s \\\\",
        df$country,
        df$period,
        fmt_rmsfe(df$vecc_M1), fmt_rmsfe(df$vecc_M2), fmt_rmsfe(df$vecc_M3),
        fmt_rel(df$rel_VecC_VecP_M1), fmt_rel(df$rel_VecC_VecP_M2), fmt_rel(df$rel_VecC_VecP_M3)
      ),
      collapse = "\n"
    )
  }
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header_1,
    header_2,
    "\\midrule\n",
    rows,
    "\n\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# 11. SIMPLE LATEX TABLE EXPORT
# ==============================================================================

list_to_latex_table <- function(df,
                                caption,
                                label,
                                digits = 4) {
  
  required_cols <- c("country", "period", "M1", "M2", "M3")
  missing_cols  <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  rows <- sprintf(
    "%s & %s & %.*f & %.*f & %.*f \\\\",
    df$country,
    df$period,
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
# utils_results_02_selection_dm_factor.R
# Selection tables, metadata, DM tests, factor helpers
# ==============================================================================

# ==============================================================================
# 1. SELECTION TABLE HELPERS
# ==============================================================================

normalize_selection_ids <- function(df) {
  if (!is.data.frame(df)) {
    stop("Input object must be a data.frame.")
  }
  
  if (!"base_name" %in% names(df)) {
    stop("Column 'base_name' is missing.")
  }
  
  df %>%
    dplyr::mutate(
      base_name = dplyr::recode(
        base_name,
        "TASS.LBD" = "TASS.LDB",
        "TLB.LBD"  = "TLB.LDB"
      )
    )
}

combine_size_codes <- function(s, m, l) {
  out <- character(length(s))
  
  for (i in seq_along(s)) {
    parts <- c(
      if (isTRUE(s[i] == 1)) "S" else NULL,
      if (isTRUE(m[i] == 1)) "M" else NULL,
      if (isTRUE(l[i] == 1)) "L" else NULL
    )
    
    out[i] <- if (length(parts) == 0L) {
      ""
    } else {
      paste(parts, collapse = "/")
    }
  }
  
  out
}

build_selection_size_table <- function(sel_small,
                                       sel_medium,
                                       sel_large,
                                       country_cols = country_order) {
  
  key_cols <- c("base_name", "frequency")
  
  for (obj_name in c("sel_small", "sel_medium", "sel_large")) {
    obj <- get(obj_name)
    missing_cols <- setdiff(c(key_cols, country_cols), names(obj))
    
    if (length(missing_cols) > 0L) {
      stop(
        "Missing columns in ", obj_name, ": ",
        paste(missing_cols, collapse = ", ")
      )
    }
  }
  
  sel_small  <- normalize_selection_ids(sel_small)
  sel_medium <- normalize_selection_ids(sel_medium)
  sel_large  <- normalize_selection_ids(sel_large)
  
  s_tbl <- sel_small %>%
    dplyr::select(
      dplyr::all_of(key_cols),
      dplyr::all_of(country_cols)
    ) %>%
    dplyr::rename_with(
      .fn   = ~ paste0(.x, "_S"),
      .cols = dplyr::all_of(country_cols)
    )
  
  m_tbl <- sel_medium %>%
    dplyr::select(
      dplyr::all_of(key_cols),
      dplyr::all_of(country_cols)
    ) %>%
    dplyr::rename_with(
      .fn   = ~ paste0(.x, "_M"),
      .cols = dplyr::all_of(country_cols)
    )
  
  l_tbl <- sel_large %>%
    dplyr::select(
      dplyr::all_of(key_cols),
      dplyr::all_of(country_cols)
    ) %>%
    dplyr::rename_with(
      .fn   = ~ paste0(.x, "_L"),
      .cols = dplyr::all_of(country_cols)
    )
  
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
    dplyr::select(
      dplyr::all_of(key_cols),
      dplyr::all_of(country_cols)
    ) %>%
    dplyr::arrange(frequency, base_name)
}

save_selection_wide_to_latex <- function(df_selection_wide,
                                         file,
                                         caption = "Variables selected by country",
                                         label = "tab:variables_selected_by_country",
                                         country_cols = NULL,
                                         check_symbol = "\\checkmark",
                                         exclude_cols = "EA") {
  
  if (!is.data.frame(df_selection_wide)) {
    stop("df_selection_wide must be a data.frame.")
  }
  
  df <- df_selection_wide
  
  meta_cols <- c("base_name", "frequency", "model_size")
  
  missing_meta <- setdiff(meta_cols, names(df))
  if (length(missing_meta) > 0L) {
    stop(
      "Missing metadata columns in df_selection_wide: ",
      paste(missing_meta, collapse = ", ")
    )
  }
  
  if (is.null(country_cols)) {
    
    preferred_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
    
    country_cols <- intersect(preferred_order, names(df))
    
    if (length(country_cols) == 0L) {
      country_cols <- setdiff(names(df), c(meta_cols, exclude_cols))
    }
    
  } else {
    
    country_cols <- setdiff(country_cols, exclude_cols)
    country_cols <- intersect(country_cols, names(df))
  }
  
  if (length(country_cols) == 0L) {
    stop("No country-specific columns found in df_selection_wide.")
  }
  
  df <- df[, c(meta_cols, country_cols), drop = FALSE]
  
  df <- df[order(df$frequency, df$base_name), ]
  
  for (cc in country_cols) {
    df[[cc]] <- ifelse(
      is.na(df[[cc]]) | df[[cc]] == 0,
      "",
      check_symbol
    )
  }
  
  df <- df %>%
    dplyr::mutate(
      base_name  = escape_latex(base_name),
      frequency  = escape_latex(frequency),
      model_size = escape_latex(model_size)
    )
  
  align <- paste(
    c("l", "c", "c", rep("c", length(country_cols))),
    collapse = " "
  )
  
  header <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\small\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{", align, "}\n",
    "\\toprule\n",
    paste(
      c(
        "\\textbf{Variable}",
        "\\textbf{Fr.}",
        "\\textbf{Size}",
        paste0("\\textbf{", country_cols, "}")
      ),
      collapse = " & "
    ),
    " \\\\\n",
    "\\midrule\n"
  )
  
  body <- apply(
    df[, c("base_name", "frequency", "model_size", country_cols), drop = FALSE],
    1,
    function(row) paste0(paste(row, collapse = " & "), " \\\\")
  )
  
  footer <- paste0(
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  latex_code <- paste0(
    header,
    paste(body, collapse = "\n"),
    footer
  )
  
  writeLines(latex_code, con = file)
  
  invisible(latex_code)
}

# ==============================================================================
# DATED PREDICTOR-SELECTION REPORT
# ==============================================================================

load_fit_by_size <- function(
    Size,
    path_results,
    model_name,
    sel_method
) {
  
  file_fit <- find_result_file(
    path  = path_results,
    model = model_name,
    stage = "fit",
    Size  = Size,
    sel   = sel_method
  )
  
  message("\nLoaded ", Size, " fit from:\n", file_fit)
  
  fit_i <- readRDS(file_fit)
  
  if (is.null(fit_i$params) || is.null(fit_i$countries)) {
    stop(
      "Fit object for Size = ", Size,
      " must contain both `params` and `countries`."
    )
  }
  
  fit_method <- tolower(
    as.character(fit_i$params$sel_method)[1L]
  )
  
  if (!identical(fit_method, tolower(sel_method))) {
    stop(
      "Selection method mismatch for Size = ", Size, ". ",
      "Expected `", sel_method,
      "`, found `", fit_i$params$sel_method, "`."
    )
  }
  
  fit_i
}


check_common_fit_settings <- function(fit_results) {
  
  fields_to_check <- c(
    "start_eval",
    "covid_end",
    "target",
    "target_cc",
    "sel_method"
  )
  
  for (field_i in fields_to_check) {
    
    values_i <- vapply(
      fit_results,
      function(x) {
        
        if (!field_i %in% names(x$params)) {
          return(NA_character_)
        }
        
        paste(
          as.character(x$params[[field_i]]),
          collapse = "|"
        )
      },
      character(1)
    )
    
    if (anyNA(values_i)) {
      stop(
        "Missing parameter `", field_i,
        "` in at least one fitted object."
      )
    }
    
    if (length(unique(values_i)) != 1L) {
      stop(
        "Fitted objects differ in parameter `", field_i, "`: ",
        paste(
          names(values_i),
          values_i,
          sep = " = ",
          collapse = "; "
        )
      )
    }
  }
  
  reference_countries <- as.character(
    fit_results[[1L]]$countries
  )
  
  same_countries <- vapply(
    fit_results,
    function(x) {
      identical(
        as.character(x$countries),
        reference_countries
      )
    },
    logical(1)
  )
  
  if (!all(same_countries)) {
    stop("Fitted objects do not use the same country set.")
  }
  
  invisible(TRUE)
}


rebuild_selection_regime <- function(
    res_fit,
    Size,
    selection_end,
    active_from,
    regime_name,
    path_data_raw,
    path_data_adj
) {
  
  if (!"selection_end" %in% names(formals(prepare_all_countries))) {
    stop(
      "`prepare_all_countries()` must include a `selection_end` argument. ",
      "Update this function in matrix.mf.tprf.prep.R before running ",
      "dated selection reporting."
    )
  }
  
  params_i    <- res_fit$params
  countries_i <- res_fit$countries
  
  all_countries_i <- prepare_all_countries(
    countries     = countries_i,
    params        = params_i,
    path_raw      = path_data_raw,
    path_adj      = path_data_adj,
    covid_mask_m  = params_i$covid_mask_m,
    covid_mask_q  = params_i$covid_mask_q,
    selection_end = as.Date(selection_end)
  )
  
  selections_i <- selection_to_df(
    all_countries_i,
    params_i
  )
  
  wide_i <- selection_to_wide(
    selections_i
  )
  
  predictors_i <- selections_i |>
    dplyr::filter(
      tolower(as.character(base_name)) !=
        tolower(as.character(params_i$target))
    )
  
  n_monthly_predictors <- predictors_i |>
    dplyr::filter(frequency == "M") |>
    dplyr::distinct(base_name) |>
    nrow()
  
  n_quarterly_predictors <- predictors_i |>
    dplyr::filter(frequency == "Q") |>
    dplyr::distinct(base_name) |>
    nrow()
  
  dimensions_i <- tibble::tibble(
    regime               = regime_name,
    size                 = Size,
    selection_end        = as.Date(selection_end),
    active_from          = as.Date(active_from),
    monthly_predictors   = n_monthly_predictors,
    quarterly_predictors = n_quarterly_predictors
  )
  
  list(
    selections    = selections_i,
    wide          = wide_i,
    dimensions    = dimensions_i,
    all_countries = all_countries_i
  )
}


rebuild_regime_all_sizes <- function(
    fit_results,
    selection_end,
    active_from,
    regime_name,
    path_data_raw,
    path_data_adj
) {
  
  out <- lapply(
    names(fit_results),
    function(Size) {
      
      rebuild_selection_regime(
        res_fit       = fit_results[[Size]],
        Size          = Size,
        selection_end = selection_end,
        active_from   = active_from,
        regime_name   = regime_name,
        path_data_raw = path_data_raw,
        path_data_adj = path_data_adj
      )
    }
  )
  
  names(out) <- names(fit_results)
  
  out
}


compare_two_selections <- function(
    sel_initial,
    sel_updated,
    Size,
    target_name
) {
  
  initial_tbl <- sel_initial |>
    dplyr::filter(
      tolower(as.character(base_name)) !=
        tolower(as.character(target_name))
    ) |>
    dplyr::distinct(
      country,
      base_name,
      frequency
    ) |>
    dplyr::mutate(
      selected_initial = 1L
    )
  
  updated_tbl <- sel_updated |>
    dplyr::filter(
      tolower(as.character(base_name)) !=
        tolower(as.character(target_name))
    ) |>
    dplyr::distinct(
      country,
      base_name,
      frequency
    ) |>
    dplyr::mutate(
      selected_updated = 1L
    )
  
  dplyr::full_join(
    initial_tbl,
    updated_tbl,
    by = c("country", "base_name", "frequency")
  ) |>
    dplyr::mutate(
      selected_initial = dplyr::coalesce(
        selected_initial,
        0L
      ),
      
      selected_updated = dplyr::coalesce(
        selected_updated,
        0L
      ),
      
      status = dplyr::case_when(
        selected_initial == 1L &
          selected_updated == 1L ~ "Retained",
        
        selected_initial == 0L &
          selected_updated == 1L ~ "Added",
        
        selected_initial == 1L &
          selected_updated == 0L ~ "Dropped"
      ),
      
      size = Size
    )
}


build_dated_selection_report <- function(config) {
  
  required_config <- c(
    "model",
    "sel_method",
    "sizes",
    "countries",
    "path_results",
    "path_data_raw",
    "path_data_adj"
  )
  
  missing_config <- setdiff(
    required_config,
    names(config)
  )
  
  if (length(missing_config) > 0L) {
    stop(
      "Missing config entries: ",
      paste(missing_config, collapse = ", ")
    )
  }
  
  size_order <- c("small", "medium", "large")
  
  sizes <- tolower(
    as.character(config$sizes)
  )
  
  if (!setequal(sizes, size_order)) {
    stop(
      "`config$sizes` must contain exactly: ",
      paste(size_order, collapse = ", ")
    )
  }
  
  fit_results <- lapply(
    size_order,
    function(Size) {
      
      load_fit_by_size(
        Size         = Size,
        path_results = config$path_results,
        model_name   = config$model,
        sel_method   = config$sel_method
      )
    }
  )
  
  names(fit_results) <- size_order
  
  check_common_fit_settings(fit_results)
  
  params_ref <- fit_results$small$params
  
  model_countries <- setdiff(
    as.character(fit_results$small$countries),
    as.character(params_ref$target_cc)
  )
  
  if (!setequal(
    model_countries,
    as.character(config$countries)
  )) {
    stop(
      "`config$countries` does not match the non-target countries ",
      "stored in the fitted object."
    )
  }
  
  selection_end_initial <- as.Date(
    params_ref$start_eval %m-% lubridate::period(1, "month")
  )
  
  selection_end_updated <- as.Date(
    params_ref$covid_end
  )
  
  initial_regime_start <- as.Date(
    params_ref$start_eval
  )
  
  updated_regime_start <- as.Date(
    selection_end_updated %m+% lubridate::period(1, "month")
  )
  
  
  message(
    "\nInitial screening sample ends: ",
    format(selection_end_initial, "%Y-%m-%d"),
    "\nInitial regime begins: ",
    format(initial_regime_start, "%Y-%m-%d"),
    "\n\nPost-COVID screening sample ends: ",
    format(selection_end_updated, "%Y-%m-%d"),
    "\nUpdated regime begins: ",
    format(updated_regime_start, "%Y-%m-%d"),
    "\n"
  )
  
  initial <- rebuild_regime_all_sizes(
    fit_results   = fit_results,
    selection_end = selection_end_initial,
    active_from   = initial_regime_start,
    regime_name   = "Pre-evaluation",
    path_data_raw = config$path_data_raw,
    path_data_adj = config$path_data_adj
  )
  
  updated <- rebuild_regime_all_sizes(
    fit_results   = fit_results,
    selection_end = selection_end_updated,
    active_from   = updated_regime_start,
    regime_name   = "Post-COVID",
    path_data_raw = config$path_data_raw,
    path_data_adj = config$path_data_adj
  )
  
  # ---------------------------------------------------------------------------
  # Predictor-set dimensions.
  # ---------------------------------------------------------------------------
  
  dimension_list <- unname(
    c(
      lapply(initial, function(x) x$dimensions),
      lapply(updated, function(x) x$dimensions)
    )
  )
  
  if (any(vapply(dimension_list, is.null, logical(1)))) {
    stop(
      "At least one rebuilt selection regime does not contain a `dimensions` table."
    )
  }
  
  predictor_dimensions <- dplyr::bind_rows(dimension_list)
  
  required_dimension_cols <- c(
    "regime",
    "size",
    "selection_end",
    "active_from",
    "monthly_predictors",
    "quarterly_predictors"
  )
  
  missing_dimension_cols <- setdiff(
    required_dimension_cols,
    names(predictor_dimensions)
  )
  
  if (length(missing_dimension_cols) > 0L) {
    stop(
      "The predictor-dimension table is missing: ",
      paste(missing_dimension_cols, collapse = ", ")
    )
  }
  
  predictor_dimensions <- predictor_dimensions |>
    dplyr::mutate(
      regime = factor(
        as.character(regime),
        levels = c("Pre-evaluation", "Post-COVID")
      ),
      size = factor(
        as.character(size),
        levels = size_order
      )
    ) |>
    dplyr::arrange(
      regime,
      size
    )
  
  selection_transition <- dplyr::bind_rows(
    lapply(
      size_order,
      function(Size) {
        
        compare_two_selections(
          sel_initial = initial[[Size]]$selections,
          sel_updated = updated[[Size]]$selections,
          Size        = Size,
          target_name = params_ref$target
        )
      }
    )
  ) |>
    dplyr::mutate(
      size = factor(
        size,
        levels = size_order
      ),
      
      frequency = factor(
        frequency,
        levels = c("M", "Q")
      )
    ) |>
    dplyr::arrange(
      size,
      frequency,
      country,
      base_name
    )
  
  selection_turnover <- selection_transition |>
    dplyr::group_by(
      size,
      frequency
    ) |>
    dplyr::summarise(
      pre_selected  = sum(selected_initial),
      retained      = sum(status == "Retained"),
      added         = sum(status == "Added"),
      dropped       = sum(status == "Dropped"),
      post_selected = sum(selected_updated),
      
      jaccard = retained / (
        retained + added + dropped
      ),
      
      .groups = "drop"
    ) |>
    dplyr::arrange(
      size,
      frequency
    )
  
  list(
    config = config,
    
    dates = list(
      selection_end_initial = selection_end_initial,
      selection_end_updated = selection_end_updated,
      initial_regime_start  = initial_regime_start,
      updated_regime_start  = updated_regime_start
    ),
    
    fit_results = fit_results,
    
    initial = initial,
    updated = updated,
    
    predictor_dimensions = predictor_dimensions,
    selection_transition = selection_transition,
    selection_turnover   = selection_turnover
  )
}

# ==============================================================================
# 1. COMMON HELPERS FOR SELECTION REPORTING
# ==============================================================================

country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")


find_result_file <- function(
    path,
    model,
    stage,
    Size,
    sel,
    ext = "rds"
) {
  
  pattern <- paste0(
    "^", stage,
    "_", model,
    "_Size-", Size,
    "_sel-", sel,
    ".*\\.", ext, "$"
  )
  
  files <- list.files(
    path       = path,
    pattern    = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0L) {
    stop(
      "No file found for model = ", model,
      ", stage = ", stage,
      ", Size = ", Size,
      ", sel = ", sel,
      ", path = ", path
    )
  }
  
  if (length(files) > 1L) {
    message("Multiple matching files found. Using the most recent one.")
    
    files <- files[
      order(
        file.info(files)$mtime,
        decreasing = TRUE
      )
    ]
  }
  
  files[1L]
}


format_month_year <- function(x) {
  format(as.Date(x), "%b %Y")
}


escape_latex <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "\\\\",
    "\\\\textbackslash{}",
    x
  )
  
  gsub(
    "([#$%&_{}])",
    "\\\\\\1",
    x,
    perl = TRUE
  )
}


normalize_selection_ids <- function(df) {
  
  if (!is.data.frame(df)) {
    stop("Input object must be a data.frame.")
  }
  
  if (!"base_name" %in% names(df)) {
    stop("Column `base_name` is missing.")
  }
  
  df |>
    dplyr::mutate(
      base_name = dplyr::recode(
        as.character(base_name),
        "TASS.LBD" = "TASS.LDB",
        "TLB.LBD"  = "TLB.LDB",
        .default   = as.character(base_name)
      )
    )
}


combine_size_codes <- function(s, m, l) {
  
  vapply(
    seq_along(s),
    function(i) {
      
      labels_i <- c(
        if (isTRUE(s[i] == 1)) "S",
        if (isTRUE(m[i] == 1)) "M",
        if (isTRUE(l[i] == 1)) "L"
      )
      
      if (length(labels_i) == 0L) {
        ""
      } else {
        paste(labels_i, collapse = "/")
      }
    },
    character(1)
  )
}


build_selection_size_table <- function(
    sel_small,
    sel_medium,
    sel_large,
    country_cols = country_order
) {
  
  key_cols <- c("base_name", "frequency")
  
  selections <- list(
    small  = sel_small,
    medium = sel_medium,
    large  = sel_large
  )
  
  for (name_i in names(selections)) {
    
    missing_cols <- setdiff(
      c(key_cols, country_cols),
      names(selections[[name_i]])
    )
    
    if (length(missing_cols) > 0L) {
      stop(
        "Missing columns in `sel_", name_i, "`: ",
        paste(missing_cols, collapse = ", ")
      )
    }
  }
  
  sel_small  <- normalize_selection_ids(sel_small)
  sel_medium <- normalize_selection_ids(sel_medium)
  sel_large  <- normalize_selection_ids(sel_large)
  
  make_size_table <- function(df, suffix) {
    
    df |>
      dplyr::select(
        dplyr::all_of(key_cols),
        dplyr::all_of(country_cols)
      ) |>
      dplyr::rename_with(
        ~ paste0(.x, suffix),
        dplyr::all_of(country_cols)
      )
  }
  
  out <- make_size_table(sel_small, "_S") |>
    dplyr::full_join(
      make_size_table(sel_medium, "_M"),
      by = key_cols
    ) |>
    dplyr::full_join(
      make_size_table(sel_large, "_L"),
      by = key_cols
    )
  
  for (cc in country_cols) {
    
    out[[cc]] <- combine_size_codes(
      out[[paste0(cc, "_S")]],
      out[[paste0(cc, "_M")]],
      out[[paste0(cc, "_L")]]
    )
  }
  
  out |>
    dplyr::select(
      dplyr::all_of(key_cols),
      dplyr::all_of(country_cols)
    ) |>
    dplyr::arrange(
      frequency,
      base_name
    )
}

# ==============================================================================
# 2. METADATA TABLE
# ==============================================================================

build_metadata_table <- function() {
  tibble::tribble(
    ~ID,          ~Series,                                                                  ~Cl, ~Cat, ~Tr, ~Fr, ~Del, ~Group,
    
    # National Accounts / Real Economy
    "GDP",        "Real Gross Domestic Product",                                            "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "EXPGS",      "Real Export Goods and services",                                         "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "IMPGS",      "Real Import Goods and services",                                         "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFCE",       "Real Government Final consumption expenditure",                          "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "HFCE",       "Real Households consumption expenditure",                                "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSD",      "Real Households consumption expenditure: Durable Goods",                 "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSSD",     "Real Households consumption expenditure: Semi-Durable Goods",            "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSSV",     "Real Households consumption expenditure: Services",                      "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "CONSND",     "Real Households consumption expenditure: Non-Durable Goods",             "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GCF",        "Real Gross capital formation",                                           "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFCF",       "Real Gross fixed capital formation",                                     "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFACON",     "Real Gross Fixed Capital Formation: Construction",                       "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GFAMG",      "Real Gross Fixed Capital Formation: Machinery and Equipment",            "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "DFGDP",      "Real Gross Domestic Product Deflator",                                   "N", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "HPRC",       "Residential Property Prices (BIS)",                                      "N", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GNFCPS",     "Gross Profit Share of Non-Financial Corporations",                       "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GNFCIR",     "Gross Investment Share of Non-Financial Corporations",                   "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    "GHIR",       "Gross Investment Rate of Households",                                    "R", "H", 5, "Q", 45, "National Accounts / Real Economy",
    "GHSR",       "Gross Households Savings Rate",                                          "R", "H", 1, "Q", 45, "National Accounts / Real Economy",
    
    # Labor Market
    "TEMP",       "Total Employment (domestic concept)",                                    "R", "H", 1, "Q", 45, "Labor Market",
    "EMP",        "Employees (domestic concept)",                                           "R", "H", 1, "Q", 45, "Labor Market",
    "SEMP",       "Self Employment (domestic concept)",                                     "R", "H", 1, "Q", 45, "Labor Market",
    "THOURS",     "Hours Worked: Total",                                                    "R", "H", 1, "Q", 45, "Labor Market",
    "EMPAG",      "Quarterly Employment: Agriculture, Forestry, Fishing",                   "R", "H", 1, "Q", 45, "Labor Market",
    "EMPIN",      "Quarterly Employment: Industry",                                         "R", "H", 1, "Q", 45, "Labor Market",
    "EMPMN",      "Quarterly Employment: Manufacturing",                                    "R", "H", 1, "Q", 45, "Labor Market",
    "EMPCON",     "Quarterly Employment: Construction",                                     "R", "H", 1, "Q", 45, "Labor Market",
    "EMPRT",      "Quarterly Employment: Wholesale/Retail trade, transport, food",          "R", "H", 1, "Q", 45, "Labor Market",
    "EMPIT",      "Quarterly Employment: Information and Communication",                    "R", "H", 1, "Q", 45, "Labor Market",
    "EMPFC",      "Quarterly Employment: Financial and Insurance activities",               "R", "H", 1, "Q", 45, "Labor Market",
    "EMPRE",      "Quarterly Employment: Real Estate",                                      "R", "H", 1, "Q", 45, "Labor Market",
    "EMPPR",      "Quarterly Employment: Professional, Scientific, Technical activities",   "R", "H", 1, "Q", 45, "Labor Market",
    "EMPPA",      "Quarterly Employment: PA, education, health and social services",        "R", "H", 1, "Q", 45, "Labor Market",
    "EMPENT",     "Quarterly Employment: Arts and recreational activities",                 "R", "H", 1, "Q", 45, "Labor Market",
    "UNETOT",     "Unemployment: Total",                                                    "R", "H", 0, "M", 45, "Labor Market",
    "UNEO25",     "Unemployment: Over 25 years",                                            "R", "H", 0, "M", 45, "Labor Market",
    "UNEU25",     "Unemployment: Under 25 years",                                           "R", "H", 0, "M", 45, "Labor Market",
    "RPRP",       "Real Labour Productivity (person)",                                      "R", "H", 1, "Q", 45, "Labor Market",
    "WS",         "Wages and salaries",                                                     "N", "H", 1, "Q", 45, "Labor Market",
    "ESC",        "Employers' Social Contributions",                                        "N", "H", 1, "Q", 45, "Labor Market",
    
    # Credit Aggregates
    "TASS.SDB",   "Total Economy - Assets: Short-Term Debt Securities",                     "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.LDB",   "Total Economy - Assets: Long-Term Debt Securities",                      "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.SLN",   "Total Economy - Assets: Short-Term Loans",                               "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TASS.LLN",   "Total Economy - Assets: Long-Term Loans",                                "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.SDB",    "Total Economy - Liabilities: Short-Term Debt Securities",                "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.LDB",    "Total Economy - Liabilities: Long-Term Debt Securities",                 "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.SLN",    "Total Economy - Liabilities: Short-Term Loans",                          "F", "H", 1, "Q", 45, "Credit Aggregates",
    "TLB.LLN",    "Total Economy - Liabilities: Long-Term Loans",                           "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCASS",     "Non-Financial Corporations: Total Financial Assets",              "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB.SLN",  "Non-Financial Corporations - Liabilities - Short-Term Loans",            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB.LLN",  "Non-Financial Corporations - Liabilities - Long-Term Loans",             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCASS.SLN", "Non-Financial Corporations - Assets: Short-Term Loans",           "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCASS.LLN", "Non-Financial Corporations - Assets: Long-Term Loans",            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB",      "Non-Financial Corporations: Total Financial Liabilities",         "F", "H", 1, "Q", 45, "Credit Aggregates",
    
    "GGASS",      "General Government: Total Financial Assets",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.SLN",  "General Government - Assets: Short-Term Loans",                          "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.LLN",  "General Government - Assets: Long-Term Loans",                           "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB",       "General Government: Total Financial Liabilities",                        "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.SLN",   "General Government - Liabilities: Short-Term Loans",                     "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.LLN",   "General Government - Liabilities: Long-Term Loans",                      "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHASS",      "Households: Total Financial Assets",                              "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHASS.SLN",  "Households - Assets: Short-Term Loans",                            "F", "H", 1, "Q", 45, "Credit Aggregates","HHASS.LLN",  "Households - Assets: Long-Term Loans",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB",       "Households: Total Financial Liabilities",                                "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.SLN",   "Households - Liabilities: Short-Term Loans",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.LLN",   "Households - Liabilities: Long-Term Loans",                              "F", "H", 1, "Q", 45, "Credit Aggregates",

    # Labor Costs
    "ULCIN",      "Nominal Unit Labor Costs: Industry",                                     "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCMQ",      "Nominal Unit Labor Costs: Mining and Quarrying",                  "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCMN",      "Nominal Unit Labor Costs: Manufacturing",                                "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCCON",     "Nominal Unit Labor Costs: Construction",                                 "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRT",      "Nominal Unit Labor Costs: Wholesale/Retail Trade, Transport, Food, IT",  "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCFC",      "Nominal Unit Labor Costs: Financial Activities",                         "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRE",      "Nominal Unit Labor Costs: Real Estate",                                  "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCPR",      "Nominal Unit Labor Costs: Professional, Scientific, Technical activities","N", "H", 1, "Q", 45, "Labor Costs",
    
    # Financial Markets
    "REER42",     "Real Exchange Rate (42 main industrial countries)",                      "F", "H", 1, "M", 35, "Financial Markets",
    "SHIX",       "Stock Price Index",                                                      "F", "S", 1, "M", 1,  "Financial Markets",
    
    # Interest Rates
    "LTIRT",      "Long-Term Interest Rates (EMU Criterion)",                               "F", "H", 2, "M", 35, "Interest Rates",
    
    # Industrial Production and Turnover
    "IPMN",       "Industrial Production Index: Manufacturing",                             "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPCAG",      "Industrial Production Index: Capital Goods",                             "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPCOG",      "Industrial Production Index: Consumer Goods",                            "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPDCOG",     "Industrial Production Index: Durable Consumer Goods",                    "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPNDCOG",    "Industrial Production Index: Non Durable Consumer Goods",                "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPING",      "Industrial Production Index: Intermediate Goods",                        "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "IPNRG",      "Industrial Production Index: Energy",                                    "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNMN",      "Turnover Index: Manufacturing",                                          "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNCAG",     "Turnover Index: Capital Goods",                                          "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNCOG",     "Turnover Index: Consumer Goods",                                         "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNDCOG",    "Turnover Index: Durable Consumer Goods",                                 "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNNDCOG",   "Turnover Index: Non Durable Consumer Goods",                             "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNING",     "Turnover Index: Intermediate Goods",                                     "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    "TRNNRG",     "Turnover Index: Energy",                                                 "R", "H", 1, "M", 45, "Industrial Production and Turnover",
    
    # Prices
    "PPICAG",     "Producer Price Index: Capital Goods",                                    "N", "H", 1, "M", 40, "Prices",
    "PPICOG",     "Producer Price Index: Consumer Goods",                                   "N", "H", 1, "M", 40, "Prices",
    "PPIDCOG",    "Producer Price Index: Durable Consumer Goods",                           "N", "H", 1, "M", 40, "Prices",
    "PPINDCOG",   "Producer Price Index: Non Durable Consumer Goods",                       "N", "H", 1, "M", 40, "Prices",
    "PPIING",     "Producer Price Index: Intermediate Goods",                               "N", "H", 1, "M", 40, "Prices",
    "PPINRG",     "Producer Price Index: Energy",                                           "N", "H", 1, "M", 40, "Prices",
    "HICPOV",     "Harmonized Index of Consumer Prices: Overall Index",                     "N", "H", 1, "M", 40, "Prices",
    "HICPNEF",    "Harmonized Index of Consumer Prices: All Items, no Energy/Food",         "N", "H", 1, "M", 40, "Prices",
    "HICPG",      "Harmonized Index of Consumer Prices: Goods",                             "N", "H", 1, "M", 40, "Prices",
    "HICPIN",     "Harmonized Index of Consumer Prices: Industrial Goods",                  "N", "H", 1, "M", 40, "Prices",
    "HICPSV",     "Harmonized Index of Consumer Prices: Services",                          "N", "H", 1, "M", 40, "Prices",
    "HICPNG",     "Harmonized Index of Consumer Prices: Energy",                            "N", "H", 1, "M", 40, "Prices",
    
    # Confidence Indicators
    "ICONFIX",    "Industrial Confidence Indicator",                                        "C", "S", 0, "M", 5,  "Confidence Indicators",
    "CCONFIX",    "Consumer Confidence Indicator",                                          "C", "S", 0, "M", 5,  "Confidence Indicators",
    "KCONFIX",    "Construction Confidence Indicator",                                      "C", "S", 0, "M", 5,  "Confidence Indicators",
    "SCONFIX",    "Services Confidence Indicator",                                          "C", "S", 0, "M", 5,  "Confidence Indicators",
    "ESENTIX",    "Economic Sentiment Indicator",                                           "C", "S", 0, "M", 5,  "Confidence Indicators",
    "RTCONFIX",   "Retail Confidence Indicator",                                            "C", "S", 0, "M", 5,  "Confidence Indicators",
    "BCI",        "Business Confidence Index",                                              "C", "S", 1, "M", 5,  "Confidence Indicators",
    "CCI",        "Consumer Confidence Index",                                              "C", "S", 1, "M", 5,  "Confidence Indicators"
  )
}

to_setcell <- function(x) {
  ifelse(
    is.na(x) | x == "",
    "",
    paste0("\\setcell{", x, "}")
  )
}

# ==============================================================================
# LATEX TABLE: COUNTRY-SPECIFIC PREDICTOR SELECTION
# ==============================================================================

build_latex_selection_table <- function(
    sel_small,
    sel_medium,
    sel_large,
    metadata = build_metadata_table(),
    country_cols = country_order,
    caption = "Country-specific predictor selection",
    label = "tab:variables_selected_country_size",
    note = NULL
) {
  
  if (is.null(note)) {
    note <- paste0(
      "\\textit{Notes:} Country entries report the information sets in which ",
      "the predictor is selected: \\textbf{S} = small, \\textbf{M} = medium, ",
      "and \\textbf{L} = large. Combined labels indicate selection in more ",
      "than one information set for the corresponding country. ",
      "\\textbf{Cl.} denotes variable class: $R$ = real activity, ",
      "$N$ = nominal, $C$ = confidence, and $F$ = financial. ",
      "\\textbf{Cat.} distinguishes hard ($H$) from soft ($S$) indicators. ",
      "\\textbf{Tr.} reports the transformation code, \\textbf{Fr.} the ",
      "sampling frequency, and \\textbf{Del.} the approximate release delay ",
      "in days. GDP is excluded from the table."
    )
  }
  
  sel_all <- build_selection_size_table(
    sel_small    = sel_small,
    sel_medium   = sel_medium,
    sel_large    = sel_large,
    country_cols = country_cols
  ) |>
    dplyr::rename(
      ID = base_name,
      Fr = frequency
    ) |>
    dplyr::filter(ID != "GDP")
  
  missing_in_meta <- setdiff(sel_all$ID, metadata$ID)
  
  if (length(missing_in_meta) > 0L) {
    stop(
      "Missing metadata for these IDs: ",
      paste(sort(missing_in_meta), collapse = ", ")
    )
  }
  
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
  
  tab <- metadata |>
    dplyr::filter(ID != "GDP") |>
    dplyr::inner_join(
      sel_all,
      by = c("ID", "Fr")
    ) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(country_cols),
        to_setcell
      ),
      ID     = escape_latex(ID),
      Series = escape_latex(Series),
      Cl     = escape_latex(Cl),
      Cat    = escape_latex(Cat),
      Fr     = escape_latex(Fr),
      Group  = factor(Group, levels = groups_order)
    ) |>
    dplyr::arrange(Group, Fr, ID) |>
    dplyr::mutate(
      N = dplyr::row_number()
    )
  
  header <- paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\tiny\n",
    "\\renewcommand{\\arraystretch}{0.82}\n",
    "\\setlength{\\tabcolsep}{1.6pt}\n",
    "\\providecommand{\\setcell}[1]{{\\fontsize{4.1}{4.4}\\selectfont #1}}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "{\\fontsize{4.0}{4.4}\\selectfont\n",
    "\\begin{tabular}{c p{1.25cm} p{5.6cm} c c c c c c c c c c c c c}\n",
    "\\toprule\n",
    "\\textbf{N} & \\textbf{ID} & \\textbf{Series} & ",
    "\\textbf{Cl.} & \\textbf{Cat.} & \\textbf{Tr.} & ",
    "\\textbf{Fr.} & \\textbf{Del.} & ",
    paste0("\\textbf{", country_cols, "}", collapse = " & "),
    " \\\\\n",
    "\\midrule\n"
  )
  
  body_lines <- character(0)
  ncols_total <- 8L + length(country_cols)
  
  for (grp in groups_order) {
    
    subtab <- tab |>
      dplyr::filter(Group == grp)
    
    if (nrow(subtab) == 0L) {
      next
    }
    
    grp_line <- paste0(
      "\\multicolumn{",
      ncols_total,
      "}{c}{\\textbf{",
      escape_latex(grp),
      "}} \\\\"
    )
    
    rows <- subtab |>
      dplyr::select(
        N,
        ID,
        Series,
        Cl,
        Cat,
        Tr,
        Fr,
        Del,
        dplyr::all_of(country_cols)
      )
    
    row_lines <- apply(
      rows,
      1,
      function(r) {
        paste0(
          paste(r, collapse = " & "),
          " \\\\"
        )
      }
    )
    
    body_lines <- c(
      body_lines,
      grp_line,
      "\\midrule",
      row_lines,
      "\\midrule"
    )
  }
  
  if (length(body_lines) > 0L) {
    body_lines <- body_lines[-length(body_lines)]
  }
  
  footer <- paste0(
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "}\n",
    "}\n",
    "\\parbox{0.98\\textwidth}{\\tiny\n",
    note,
    "}\n",
    "\\end{table}\n"
  )
  
  paste0(
    header,
    paste(body_lines, collapse = "\n"),
    footer
  )
}


# ==============================================================================
# LATEX TABLE: PREDICTOR-SET COMPOSITION ACROSS DATED REGIMES
# ==============================================================================
build_latex_predictor_dimensions_table <- function(
    tbl,
    caption,
    label,
    evaluation_end = NULL,
    note = paste0(
      "\\textit{Notes:} Entries are counts of non-GDP predictors in the ",
      "country-union information set. Within each information set, ",
      "Total $N=N_m+N_q$ is the predictor-column dimension."
    )
) {
  
  required_cols <- c(
    "regime",
    "size",
    "active_from",
    "monthly_predictors",
    "quarterly_predictors"
  )
  
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in predictor-dimension table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  size_order <- c("small", "medium", "large")
  
  size_columns <- c(
    "small_M",  "small_Q",  "small_N",
    "medium_M", "medium_Q", "medium_N",
    "large_M",  "large_Q",  "large_N"
  )
  
  tbl_clean <- tbl |>
    dplyr::mutate(
      regime      = as.character(regime),
      size        = tolower(as.character(size)),
      active_from = as.Date(active_from)
    ) |>
    dplyr::distinct(
      regime,
      size,
      active_from,
      monthly_predictors,
      quarterly_predictors
    )
  
  tab <- tbl_clean |>
    dplyr::select(
      regime,
      size,
      monthly_predictors,
      quarterly_predictors
    ) |>
    tidyr::pivot_longer(
      cols      = c(monthly_predictors, quarterly_predictors),
      names_to  = "frequency",
      values_to = "n_predictors"
    ) |>
    dplyr::mutate(
      frequency = dplyr::recode(
        frequency,
        monthly_predictors   = "M",
        quarterly_predictors = "Q"
      )
    ) |>
    tidyr::pivot_wider(
      names_from  = c(size, frequency),
      values_from = n_predictors,
      names_glue  = "{size}_{frequency}"
    )
  
  for (size_i in size_order) {
    
    col_M <- paste0(size_i, "_M")
    col_Q <- paste0(size_i, "_Q")
    col_N <- paste0(size_i, "_N")
    
    if (!col_M %in% names(tab)) {
      tab[[col_M]] <- NA_integer_
    }
    
    if (!col_Q %in% names(tab)) {
      tab[[col_Q]] <- NA_integer_
    }
    
    tab[[col_N]] <- as.integer(
      tab[[col_M]] + tab[[col_Q]]
    )
  }
  
  latex_month_year <- function(x) {
    
    month_names <- c(
      "Jan.", "Feb.", "Mar.", "Apr.",
      "May", "Jun.", "Jul.", "Aug.",
      "Sep.", "Oct.", "Nov.", "Dec."
    )
    
    x <- as.Date(x)
    
    paste0(
      month_names[as.integer(format(x, "%m"))],
      "\\ ",
      format(x, "%Y")
    )
  }
  
  regime_labels <- c(
    "Pre-evaluation" = "Initial regime",
    "Post-COVID"     = "Post-COVID update"
  )
  
  if (!is.null(evaluation_end)) {
    
    regime_dates <- tbl_clean |>
      dplyr::group_by(regime) |>
      dplyr::summarise(
        active_from = dplyr::first(active_from),
        .groups = "drop"
      )
    
    initial_start <- regime_dates |>
      dplyr::filter(regime == "Pre-evaluation") |>
      dplyr::pull(active_from)
    
    update_start <- regime_dates |>
      dplyr::filter(regime == "Post-COVID") |>
      dplyr::pull(active_from)
    
    if (length(initial_start) == 1L && length(update_start) == 1L) {
      
      initial_end <- seq.Date(
        from       = update_start,
        by         = "-1 month",
        length.out = 2L
      )[2L]
      
      regime_labels["Pre-evaluation"] <- paste0(
        latex_month_year(initial_start),
        "--",
        latex_month_year(initial_end)
      )
      
      regime_labels["Post-COVID"] <- paste0(
        latex_month_year(update_start),
        "--",
        latex_month_year(as.Date(evaluation_end))
      )
    }
  }
  
  tab <- tab |>
    dplyr::mutate(
      regime_order = match(
        regime,
        c("Pre-evaluation", "Post-COVID")
      ),
      Regime = unname(regime_labels[regime])
    ) |>
    dplyr::arrange(regime_order) |>
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(size_columns),
        ~ ifelse(
          is.na(.x),
          "--",
          as.character(as.integer(.x))
        )
      )
    )
  
  body <- vapply(
    seq_len(nrow(tab)),
    function(i) {
      
      row_prefix <- if (tab$regime[i] == "Post-COVID") {
        "\\rowcolor{matrixgray}\n"
      } else {
        ""
      }
      
      values_i <- c(
        tab$Regime[i],
        unlist(
          tab[i, size_columns],
          use.names = FALSE
        )
      )
      
      paste0(
        row_prefix,
        paste(values_i, collapse = " & "),
        " \\\\"
      )
    },
    character(1)
  )
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{1.3}\n",
    "\\setlength{\\tabcolsep}{7.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n",
    "\\definecolor{hypergray}{gray}{0.91}\n",
    "\\definecolor{matrixgray}{gray}{0.965}\n",
    
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    
    "\\begin{adjustbox}{max width=\\textwidth,center}\n",
    "\\begin{tabular}{@{}c ccc ccc ccc@{}}\n",
    "\\toprule\n",
    
    "\\rowcolor{hypergray}\n",
    " & \\multicolumn{3}{c}{\\textbf{Small}} & ",
    "\\multicolumn{3}{c}{\\textbf{Medium}} & ",
    "\\multicolumn{3}{c}{\\textbf{Large}} \\\\\n",
    
    "\\cmidrule(lr){2-4}",
    "\\cmidrule(lr){5-7}",
    "\\cmidrule(lr){8-10}\n",
    
    "\\rowcolor{topgray}\n",
    "\\textbf{Regime} & ",
    "\\textbf{Monthly} & \\textbf{Quarterly} & \\textbf{Total $N$} & ",
    "\\textbf{Monthly} & \\textbf{Quarterly} & \\textbf{Total $N$} & ",
    "\\textbf{Monthly} & \\textbf{Quarterly} & \\textbf{Total $N$} \\\\\n",
    
    "\\midrule\n",
    paste(body, collapse = "\n"),
    
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{adjustbox}\n",
    
    "\\vspace{-0.3cm}\n",
    "\\begin{center}\n",
    "\\parbox{0.88\\textwidth}{\\centering\\scriptsize\n",
    note,
    "}\n",
    "\\end{center}\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# LATEX TABLE: COUNTRY-SPECIFIC SELECTION TURNOVER
# ==============================================================================

build_latex_selection_turnover_table <- function(
    tbl,
    caption,
    label,
    note = paste0(
      "\\textit{Notes:} Entries are country--predictor pairs. Retained pairs ",
      "are selected in both regimes; added and removed pairs are selected ",
      "only after and only before the update, respectively. Overlap is ",
      "computed relative to the union of pairs selected across regimes."
    )
) {
  
  required_cols <- c(
    "size",
    "frequency",
    "pre_selected",
    "retained",
    "added",
    "dropped",
    "post_selected",
    "jaccard"
  )
  
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in selection-turnover table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  tab <- tbl |>
    dplyr::mutate(
      size      = tolower(as.character(size)),
      frequency = toupper(as.character(frequency))
    )
  
  format_int <- function(x) {
    
    if (length(x) == 0L || is.na(x)) {
      return("--")
    }
    
    as.character(as.integer(round(x)))
  }
  
  format_pct <- function(x) {
    
    if (length(x) == 0L || is.na(x)) {
      return("--")
    }
    
    sprintf("%.1f", 100 * x)
  }
  
  build_row <- function(size_i, frequency_i, label_i) {
    
    row_i <- tab |>
      dplyr::filter(
        size == size_i,
        frequency == frequency_i
      )
    
    if (nrow(row_i) == 0L) {
      return(
        paste(
          c(label_i, rep("--", 6L)),
          collapse = " & "
        )
      )
    }
    
    if (nrow(row_i) > 1L) {
      stop(
        "More than one turnover row for size = ",
        size_i,
        " and frequency = ",
        frequency_i,
        "."
      )
    }
    
    paste(
      c(
        label_i,
        format_int(row_i$pre_selected),
        format_int(row_i$retained),
        format_int(row_i$added),
        format_int(row_i$dropped),
        format_int(row_i$post_selected),
        format_pct(row_i$jaccard)
      ),
      collapse = " & "
    )
  }
  
  size_labels <- c(
    small  = "Small set",
    medium = "Medium set",
    large  = "Large set"
  )
  
  blocks <- unlist(
    lapply(
      c("small", "medium", "large"),
      function(size_i) {
        
        prefix_i <- if (size_i == "small") {
          character(0)
        } else {
          "\\addlinespace[0.12em]"
        }
        
        c(
          prefix_i,
          "\\rowcolor{hypergray}",
          paste0(
            "\\multicolumn{7}{c}{\\textbf{",
            size_labels[size_i],
            "}} \\\\"
          ),
          "\\cmidrule(lr){1-7}",
          paste0(build_row(size_i, "M", "Monthly"), " \\\\"),
          paste0(build_row(size_i, "Q", "Quarterly"), " \\\\")
        )
      }
    )
  )
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{0.90}\n",
    "\\setlength{\\tabcolsep}{2.0pt}\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{@{}c c c c c c c@{}}\n",
    "\\toprule\n",
    "\\rowcolor{topgray}\n",
    "\\textbf{Frequency} & ",
    "\\textbf{Initial} & ",
    "\\textbf{Retained} & ",
    "\\textbf{Added} & ",
    "\\textbf{Removed} & ",
    "\\textbf{Updated} & ",
    "\\textbf{Overlap (\\%)} \\\\\n",
    "\\midrule\n\n",
    paste(blocks, collapse = "\n"),
    "\n\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\vspace{-0.3cm}\n",
    "\\begin{center}\n",
    "\\parbox{0.88\\textwidth}{\\centering\\scriptsize\n",
    note,
    "}\n",
    "\\end{center}\n",
    "\\end{table}\n"
  )
}
# ==============================================================================
# LATEX COMPACT TABLE: COUNTRY-SPECIFIC SELECTION CHANGES
# =============================================================================

build_latex_selection_changes_by_status_table <- function(
    selection_transition,
    caption,
    label,
    country_order = NULL
) {
  
  required_cols <- c(
    "country",
    "base_name",
    "frequency",
    "size",
    "status"
  )
  
  missing_cols <- setdiff(required_cols, names(selection_transition))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in selection-transition table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  size_levels <- c("small", "medium", "large")
  
  change_columns <- c(
    "small_added",
    "small_removed",
    "medium_added",
    "medium_removed",
    "large_added",
    "large_removed"
  )
  
  order_countries <- function(x) {
    
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    
    if (length(x) == 0L) {
      return(character(0))
    }
    
    if (is.null(country_order)) {
      return(sort(x))
    }
    
    country_rank <- match(x, country_order)
    
    x[
      order(
        is.na(country_rank),
        ifelse(is.na(country_rank), Inf, country_rank),
        x
      )
    ]
  }
  
  collapse_countries <- function(x) {
    
    x <- order_countries(x)
    
    if (length(x) == 0L) {
      return("--")
    }
    
    paste(x, collapse = ", ")
  }
  
  changes <- selection_transition |>
    dplyr::transmute(
      country = as.character(country),
      
      base_name = dplyr::recode(
        as.character(base_name),
        "TASS.LBD" = "TASS.LDB",
        "TLB.LBD"  = "TLB.LDB",
        .default   = as.character(base_name)
      ),
      
      frequency = toupper(as.character(frequency)),
      size      = tolower(as.character(size)),
      status_raw = tolower(as.character(status))
    ) |>
    dplyr::mutate(
      status = dplyr::case_when(
        status_raw == "added" ~ "Added",
        status_raw %in% c("dropped", "removed") ~ "Removed",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(
      !is.na(status),
      size %in% size_levels,
      frequency %in% c("M", "Q")
    ) |>
    dplyr::select(-status_raw)
  
  if (nrow(changes) == 0L) {
    return(
      paste0(
        "\\begin{table}[htbp]\n",
        "\\centering\n",
        "\\small\n",
        "\\caption{", escape_latex(caption), "}\n",
        "\\label{", label, "}\n",
        "\\begin{tabular}{l}\n",
        "\\toprule\n",
        "No country-specific changes in predictor selection. \\\\\n",
        "\\bottomrule\n",
        "\\end{tabular}\n",
        "\\end{table}\n"
      )
    )
  }
  
  change_cells <- changes |>
    dplyr::group_by(
      base_name,
      frequency,
      size,
      status
    ) |>
    dplyr::summarise(
      countries = collapse_countries(country),
      .groups   = "drop"
    ) |>
    dplyr::mutate(
      direction = tolower(status)
    ) |>
    dplyr::select(-status) |>
    tidyr::pivot_wider(
      names_from  = c(size, direction),
      values_from = countries,
      names_glue  = "{size}_{direction}",
      values_fill = list(countries = "--")
    )
  
  for (column_i in change_columns) {
    
    if (!column_i %in% names(change_cells)) {
      change_cells[[column_i]] <- "--"
    }
    
    change_cells[[column_i]] <- as.character(change_cells[[column_i]])
    
    change_cells[[column_i]][
      is.na(change_cells[[column_i]]) |
        !nzchar(change_cells[[column_i]])
    ] <- "--"
  }
  
  tab <- change_cells |>
    dplyr::mutate(
      frequency_order = match(frequency, c("M", "Q")),
      
      Predictor = paste0(
        "\\texttt{",
        escape_latex(base_name),
        "}"
      )
    ) |>
    dplyr::arrange(
      frequency_order,
      base_name
    ) |>
    dplyr::select(
      frequency,
      Predictor,
      dplyr::all_of(change_columns)
    )
  
  build_frequency_panel <- function(
    frequency_code,
    panel_title
  ) {
    
    panel <- tab |>
      dplyr::filter(frequency == frequency_code) |>
      dplyr::select(-frequency)
    
    if (nrow(panel) == 0L) {
      return(character(0))
    }
    
    rows <- vapply(
      seq_len(nrow(panel)),
      function(i) {
        
        row_colour <- if (i %% 2L == 1L) {
          "\\rowcolor{selectionrowa}\n"
        } else {
          "\\rowcolor{selectionrowb}\n"
        }
        
        row_i <- unlist(
          panel[i, ],
          use.names = FALSE
        )
        
        paste0(
          row_colour,
          paste(row_i, collapse = " & "),
          " \\\\"
        )
      },
      character(1)
    )
    
    separator <- c(
      "\\addlinespace[0.15em]",
      "\\specialrule{0.65pt}{0.20em}{0pt}",
      "\\specialrule{0.25pt}{0.05em}{0.18em}"
    )
    
    c(
      separator,
      "\\rowcolor{matrixgray}",
      paste0(
        "\\multicolumn{7}{c}{\\textbf{",
        panel_title,
        "}} \\\\[-0.15em]"
      ),
      "\\cmidrule(lr){1-7}",
      rows
    )
  }
  
  column_header <- paste0(
    "\\rowcolor{topgray}\n",
    "\\multicolumn{1}{c}{\\textbf{Predictor}} & ",
    "\\multicolumn{2}{c}{\\textbf{Small}} & ",
    "\\multicolumn{2}{c}{\\textbf{Medium}} & ",
    "\\multicolumn{2}{c}{\\textbf{Large}} \\\\\n",
    "\\cmidrule(lr){2-3}",
    "\\cmidrule(lr){4-5}",
    "\\cmidrule(lr){6-7}\n",
    "\\rowcolor{hypergray}\n",
    " & \\textbf{Added} & \\textbf{Removed} & ",
    "\\textbf{Added} & \\textbf{Removed} & ",
    "\\textbf{Added} & \\textbf{Removed} \\\\\n"
  )
  
  body <- paste(
    c(
      build_frequency_panel("M", "Monthly predictors"),
      build_frequency_panel("Q", "Quarterly predictors")
    ),
    collapse = "\n"
  )
  
  header <- paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n",
    "\\vspace{0.08cm}\n",
    "\\begingroup\n",
    "\\fontsize{6.5}{7.1}\\selectfont\n",
    "\\renewcommand{\\arraystretch}{0.86}\n",
    "\\setlength{\\tabcolsep}{0.45pt}\n",
    "\\begin{adjustbox}{",
    "max width=\\textwidth,",
    "max totalheight=0.85\\textheight,",
    "center",
    "}\n",
    "\\begin{tabular}{@{}",
    ">{\\centering\\arraybackslash}p{1.45cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    ">{\\centering\\arraybackslash}p{2.55cm}",
    "@{}}\n",
    "\\toprule\n",
    column_header,
    "\\midrule\n"
  )
  
  footer <- paste0(
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{adjustbox}\n",
    "\\par\\vspace{0.04cm}\\noindent",
    "\\parbox{0.96\\textwidth}{\\scriptsize ",
    "\\textit{Notes:} Added and removed countries are reported by information ",
    "set; \\texttt{--} denotes no change.}\n",
    "\\endgroup\n",
    "\\end{table}\n"
  )
  
  paste0(
    header,
    body,
    footer
  )
}



# ==============================================================================
# 3. DIEBOLD-MARIANO TEST HELPERS
# ==============================================================================

build_eval_df <- function(df_rt, df_y) {
  required_rt <- c("date", "country", "nowcast", "type")
  required_y  <- c("country", "quarter_id", "GDP", "period")
  
  missing_rt <- setdiff(required_rt, names(df_rt))
  missing_y  <- setdiff(required_y, names(df_y))
  
  if (length(missing_rt) > 0L) {
    stop("df_rt is missing columns: ", paste(missing_rt, collapse = ", "))
  }
  
  if (length(missing_y) > 0L) {
    stop("df_y is missing columns: ", paste(missing_y, collapse = ", "))
  }
  
  df_rt %>%
    dplyr::mutate(
      date       = as.Date(date),
      quarter_id = make_quarter_id(date)
    ) %>%
    dplyr::left_join(
      df_y %>%
        dplyr::select(country, quarter_id, GDP, period),
      by = c("country", "quarter_id")
    ) %>%
    dplyr::mutate(
      error = GDP - nowcast,
      se    = error^2
    )
}

dm_test_hac <- function(d,
                        lag = NULL,
                        alternative = c("less", "two.sided", "greater")) {
  
  alternative <- match.arg(alternative)
  d <- d[is.finite(d)]
  n <- length(d)
  
  if (n < 4L) {
    return(data.frame(
      n       = n,
      mean_d  = NA_real_,
      stat    = NA_real_,
      p_value = NA_real_
    ))
  }
  
  if (is.null(lag)) {
    lag <- floor(n^(1 / 3))
  }
  
  fit <- stats::lm(d ~ 1)
  
  vc <- sandwich::NeweyWest(
    fit,
    lag      = lag,
    prewhite = FALSE,
    adjust   = TRUE
  )
  
  mean_d  <- stats::coef(fit)[1]
  se_mean <- sqrt(vc[1, 1])
  stat    <- mean_d / se_mean
  
  p_value <- switch(
    alternative,
    "less"      = stats::pnorm(stat),
    "greater"   = 1 - stats::pnorm(stat),
    "two.sided" = 2 * stats::pnorm(-abs(stat))
  )
  
  data.frame(
    n       = n,
    mean_d  = unname(mean_d),
    stat    = unname(stat),
    p_value = unname(p_value)
  )
}

build_dm_wide <- function(eval_matrix,
                          eval_benchmark,
                          benchmark_tag,
                          alternative = "less") {
  
  p_prefix <- paste0("p_", benchmark_tag, "_")
  
  dm_base <- eval_matrix %>%
    dplyr::select(country, quarter_id, type, period, se_matrix = se) %>%
    dplyr::inner_join(
      eval_benchmark %>%
        dplyr::select(country, quarter_id, type, se_benchmark = se),
      by = c("country", "quarter_id", "type")
    ) %>%
    dplyr::filter(!is.na(se_matrix), !is.na(se_benchmark)) %>%
    dplyr::mutate(d = se_matrix - se_benchmark)
  
  dm_long <- dm_base %>%
    dplyr::group_by(country, period, type) %>%
    dplyr::group_modify(
      ~ dm_test_hac(.x$d, alternative = alternative)
    ) %>%
    dplyr::ungroup()
  
  dm_wide <- dm_long %>%
    dplyr::select(country, period, type, p_value) %>%
    tidyr::pivot_wider(
      names_from   = type,
      values_from  = p_value,
      names_prefix = p_prefix
    )
  
  list(
    base = dm_base,
    long = dm_long,
    wide = dm_wide
  )
}

# ==============================================================================
# 4. PLOT THEMES AND FACTOR HELPERS
# ==============================================================================

theme_factor_compare <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25")
    )
}

theme_country_compare <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25")
    )
}

standardize_by_series <- function(df, value_col = "value") {
  if (!"series" %in% names(df)) {
    stop("Column 'series' is required.")
  }
  
  if (!value_col %in% names(df)) {
    stop("Column '", value_col, "' is missing.")
  }
  
  df %>%
    dplyr::group_by(series) %>%
    dplyr::mutate(
      value_std = as.numeric(scale(.data[[value_col]]))
    ) %>%
    dplyr::ungroup()
}

is_factor_case <- function(Size, sel) {
  identical(Size, "small") && sel %in% c("corr", "LASSO")
}

build_fit_stats <- function(df,
                            y_col,
                            x_cols,
                            spec_name,
                            model_name) {
  
  missing_cols <- setdiff(c(y_col, x_cols), names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in regression data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  form <- stats::as.formula(
    paste(y_col, "~", paste(x_cols, collapse = " + "))
  )
  
  fit <- stats::lm(form, data = df)
  s   <- summary(fit)
  
  data.frame(
    model          = model_name,
    specification  = spec_name,
    R2             = unname(s$r.squared),
    Adj_R2         = unname(s$adj.r.squared),
    stringsAsFactors = FALSE
  )
}

make_quarter_design <- function(f_q, dates_q, y_q) {
  data.frame(
    date   = as.Date(dates_q),
    GDP    = as.numeric(y_q),
    Factor = as.numeric(f_q),
    stringsAsFactors = FALSE
  )
}

make_monthly_design <- function(f_m, dates_q, y_q) {
  T_q <- length(y_q)
  
  if (length(f_m) < 3 * T_q) {
    stop("f_m must have length at least 3 * length(y_q).")
  }
  
  f_use <- as.numeric(f_m[seq_len(3 * T_q)])
  
  data.frame(
    date = as.Date(dates_q),
    GDP  = as.numeric(y_q),
    M1   = f_use[seq(1, 3 * T_q, by = 3)],
    M2   = f_use[seq(2, 3 * T_q, by = 3)],
    M3   = f_use[seq(3, 3 * T_q, by = 3)],
    stringsAsFactors = FALSE
  )
}

make_factor_fit_latex <- function(df_long,
                                  caption,
                                  label) {
  
  required_cols <- c("model", "specification", "R2", "Adj_R2")
  missing_cols  <- setdiff(required_cols, names(df_long))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df_long: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df_r2 <- df_long %>%
    dplyr::mutate(R2 = round(100 * R2, 1)) %>%
    dplyr::select(model, specification, R2) %>%
    tidyr::pivot_wider(names_from = model, values_from = R2)
  
  df_adj <- df_long %>%
    dplyr::mutate(Adj_R2 = round(100 * Adj_R2, 1)) %>%
    dplyr::select(model, specification, Adj_R2) %>%
    tidyr::pivot_wider(names_from = model, values_from = Adj_R2)
  
  model_cols <- setdiff(colnames(df_r2), "specification")
  header_models <- paste(model_cols, collapse = " & ")
  
  rows_r2 <- paste(
    apply(
      df_r2,
      1,
      function(x) paste(c(x["specification"], x[model_cols]), collapse = " & ")
    ),
    collapse = " \\\\\n"
  )
  
  rows_adj <- paste(
    apply(
      df_adj,
      1,
      function(x) paste(c(x["specification"], x[model_cols]), collapse = " & ")
    ),
    collapse = " \\\\\n"
  )
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{l", paste(rep("c", length(model_cols)), collapse = ""), "}\n",
    "\\toprule\n",
    " & ", header_models, " \\\\\n",
    "\\midrule\n",
    "\\multicolumn{", length(model_cols) + 1, "}{l}{\\textit{$R^2$}} \\\\\n",
    rows_r2, " \\\\\n",
    "\\midrule\n",
    "\\multicolumn{", length(model_cols) + 1, "}{l}{\\textit{Adjusted $R^2$}} \\\\\n",
    rows_adj, " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# 5. COUNTRY EXTRACTION HELPERS
# ==============================================================================

extract_country_rt <- function(df_rt,
                               country_code,
                               model_label) {
  
  df_rt %>%
    dplyr::filter(country == country_code) %>%
    dplyr::mutate(
      model = model_label,
      type  = factor(type, levels = month_order)
    ) %>%
    dplyr::select(date, country, nowcast, type, model)
}

extract_country_gdp <- function(df_yq,
                                country_code) {
  
  df_yq %>%
    dplyr::filter(country == country_code) %>%
    dplyr::select(date, country, GDP)
}

extract_country_rt_multi <- function(df_rt,
                                     country_codes,
                                     model_label) {
  
  df_rt %>%
    dplyr::filter(country %in% country_codes) %>%
    dplyr::mutate(
      model   = model_label,
      type    = factor(type, levels = month_order),
      country = factor(country, levels = country_codes)
    ) %>%
    dplyr::select(date, country, nowcast, type, model)
}

extract_country_gdp_multi <- function(df_yq,
                                      country_codes) {
  
  df_yq %>%
    dplyr::filter(country %in% country_codes) %>%
    dplyr::mutate(
      country = factor(country, levels = country_codes)
    ) %>%
    dplyr::select(date, country, GDP)
}

# ==============================================================================
# utils_results_03_wins_latex.R
# Synthetic summary tables, forecast-level win analysis, LaTeX generators
# ==============================================================================

# ==============================================================================
# 1. RMSFE-LEVEL SYNTHETIC TABLES
# ==============================================================================

reshape_comp_to_long <- function(df_comp,
                                 include_dfm = TRUE) {
  
  base_cols <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  dfm_cols <- c("DFM_M1", "DFM_M2", "DFM_M3")
  
  cols <- base_cols
  
  if (isTRUE(include_dfm) && all(dfm_cols %in% names(df_comp))) {
    cols <- c(cols, dfm_cols)
  }
  
  missing_cols <- setdiff(cols, names(df_comp))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df_comp: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df_comp %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(cols),
      names_to  = c("model", "month"),
      names_sep = "_",
      values_to = "rmsfe"
    ) %>%
    tidyr::pivot_wider(
      names_from  = model,
      values_from = rmsfe
    ) %>%
    dplyr::mutate(
      month   = factor(month, levels = month_order),
      period  = factor(
        period,
        levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")
      ),
      country = factor(country, levels = country_order)
    )
}

compute_summary_metrics <- function(df_long) {
  
  required_cols <- c("Matrix", "Vector", "VecTensor")
  missing_cols  <- setdiff(required_cols, names(df_long))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df_long: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df_long %>%
    dplyr::mutate(
      win_vecp = as.integer(Matrix < VecTensor),
      win_vecc = as.integer(Matrix < Vector)
    )
  
  if ("DFM" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(
        win_dfm = as.integer(Matrix < DFM)
      )
  }
  
  out
}

build_period_month_table <- function(df_eval,
                                     include_dfm = TRUE) {
  
  required_cols <- c("period", "month", "win_vecp", "win_vecc")
  missing_cols  <- setdiff(required_cols, names(df_eval))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df_eval: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df_eval %>%
    dplyr::filter(period %in% period_order) %>%
    dplyr::group_by(period, month) %>%
    dplyr::summarise(
      share_vecp = 100 * mean(win_vecp, na.rm = TRUE),
      share_vecc = 100 * mean(win_vecc, na.rm = TRUE),
      share_dfm  = if ("win_dfm" %in% names(.)) {
        100 * mean(win_dfm, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      period = factor(period, levels = period_order),
      month  = factor(month, levels = month_order)
    ) %>%
    dplyr::arrange(period, month) %>%
    tidyr::pivot_wider(
      names_from  = month,
      values_from = c(share_vecp, share_vecc, share_dfm),
      names_glue  = "{.value}_{month}"
    ) %>%
    dplyr::rename(
      Period      = period,
      `VEC-P M1`  = share_vecp_M1,
      `VEC-P M2`  = share_vecp_M2,
      `VEC-P M3`  = share_vecp_M3,
      `VEC-C M1`  = share_vecc_M1,
      `VEC-C M2`  = share_vecc_M2,
      `VEC-C M3`  = share_vecc_M3,
      `DFM M1`    = share_dfm_M1,
      `DFM M2`    = share_dfm_M2,
      `DFM M3`    = share_dfm_M3
    ) %>%
    dplyr::mutate(
      dplyr::across(-Period, ~ round(.x, 1))
    )
  
  if (!isTRUE(include_dfm)) {
    out <- out %>%
      dplyr::select(-dplyr::starts_with("DFM "))
  }
  
  out
}

# ==============================================================================
# 2. FORECAST-LEVEL WIN TABLES
# ==============================================================================

build_forecast_win_table <- function(df,
                                     group_vars = NULL,
                                     include_dfm = TRUE) {
  
  has_dfm <- isTRUE(include_dfm) && "win_dfm" %in% names(df)
  
  if (is.null(group_vars)) {
    
    out <- df %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Group = "Overall", .before = 1)
    
    if (has_dfm) {
      out$`DFM wins`      <- sum(df$win_dfm, na.rm = TRUE)
      out$`DFM win share` <- 100 * mean(df$win_dfm, na.rm = TRUE)
    }
    
  } else {
    
    out <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (has_dfm) {
      dfm_part <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
          `DFM wins`      = sum(win_dfm, na.rm = TRUE),
          `DFM win share` = 100 * mean(win_dfm, na.rm = TRUE),
          .groups = "drop"
        )
      
      out <- out %>%
        dplyr::left_join(dfm_part, by = group_vars)
    }
  }
  
  out %>%
    dplyr::mutate(
      Total             = as.integer(Total),
      `VEC-P wins`      = as.integer(`VEC-P wins`),
      `VEC-C wins`      = as.integer(`VEC-C wins`),
      `VEC-P win share` = round(`VEC-P win share`, 1),
      `VEC-C win share` = round(`VEC-C win share`, 1),
      dplyr::across(dplyr::any_of("DFM wins"), as.integer),
      dplyr::across(dplyr::any_of("DFM win share"), ~ round(.x, 1))
    )
}

build_country_period_month_table <- function(df,
                                             include_dfm = TRUE) {
  
  required_cols <- c("period", "country", "type", "win_vecp", "win_vecc")
  missing_cols  <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- df %>%
    dplyr::filter(period %in% period_order) %>%
    dplyr::group_by(period, country, type) %>%
    dplyr::summarise(
      share_vecp = 100 * mean(win_vecp, na.rm = TRUE),
      share_vecc = 100 * mean(win_vecc, na.rm = TRUE),
      share_dfm  = if ("win_dfm" %in% names(.)) {
        100 * mean(win_dfm, na.rm = TRUE)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      period  = factor(period, levels = period_order),
      country = factor(country, levels = country_order),
      type    = factor(type, levels = month_order)
    ) %>%
    dplyr::arrange(period, country, type) %>%
    tidyr::pivot_wider(
      names_from  = type,
      values_from = c(share_vecp, share_vecc, share_dfm),
      names_glue  = "{.value}_{type}"
    ) %>%
    dplyr::rename(
      Period      = period,
      Country     = country,
      `VEC-P M1`  = share_vecp_M1,
      `VEC-P M2`  = share_vecp_M2,
      `VEC-P M3`  = share_vecp_M3,
      `VEC-C M1`  = share_vecc_M1,
      `VEC-C M2`  = share_vecc_M2,
      `VEC-C M3`  = share_vecc_M3,
      `DFM M1`    = share_dfm_M1,
      `DFM M2`    = share_dfm_M2,
      `DFM M3`    = share_dfm_M3
    ) %>%
    dplyr::mutate(
      dplyr::across(-c(Period, Country), ~ round(.x, 1))
    )
  
  if (!isTRUE(include_dfm)) {
    out <- out %>%
      dplyr::select(-dplyr::starts_with("DFM "))
  }
  
  out
}

# ==============================================================================
# 3. FORMAT HELPERS FOR LATEX WIN TABLES
# ==============================================================================

fmt_int <- function(x) {
  ifelse(
    is.na(x),
    "--",
    as.character(as.integer(round(x)))
  )
}

fmt_pct <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "--",
    format(round(x, digits), nsmall = digits, trim = TRUE)
  )
}

bold_if_good_share <- function(x, digits = 1) {
  val <- fmt_pct(x, digits)
  
  ifelse(
    !is.na(x) & x > 50,
    paste0("\\textbf{", val, "}"),
    val
  )
}

fmt_small_share <- function(x, digits = 1) {
  val <- fmt_pct(x, digits)
  
  ifelse(
    !is.na(x) & x > 50,
    paste0("\\textbf{", val, "}"),
    val
  )
}

format_win_table_values <- function(df_tex, original_df = df_tex) {
  
  if ("Total" %in% names(df_tex)) {
    df_tex$Total <- fmt_int(original_df$Total)
  }
  
  win_cols <- grep(" wins$", names(df_tex), value = TRUE)
  share_cols <- grep(" win share$", names(df_tex), value = TRUE)
  
  for (nm in win_cols) {
    df_tex[[nm]] <- fmt_int(original_df[[nm]])
  }
  
  for (nm in share_cols) {
    df_tex[[nm]] <- bold_if_good_share(original_df[[nm]])
  }
  
  df_tex
}

# ==============================================================================
# 4. LATEX: PERIOD x MONTH WIN SHARES
# ==============================================================================

make_latex_period_month_wins <- function(
    df,
    caption = "Win shares of Matrix MF--TPRF across periods and nowcast vintages",
    label   = "tab:summary_wins_period_month",
    note    = "Share of cases in which Matrix MF--TPRF outperforms each benchmark."
) {
  
  df_tex <- df
  df_tex$Period <- escape_latex(df_tex$Period)
  
  for (nm in names(df_tex)[-1]) {
    df_tex[[nm]] <- bold_if_good_share(df[[nm]])
  }
  
  has_dfm <- all(c("DFM M1", "DFM M2", "DFM M3") %in% names(df_tex))
  
  if (has_dfm) {
    ordered_cols <- c(
      "Period",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3",
      "DFM M1", "DFM M2", "DFM M3"
    )
    
    tabular <- "cccccccccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P win share (\\%)}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C win share (\\%)}}",
      " & \\multicolumn{3}{c}{\\textbf{DFM win share (\\%)}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Period}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  } else {
    ordered_cols <- c(
      "Period",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3"
    )
    
    tabular <- "ccccccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P win share (\\%)}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C win share (\\%)}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Period}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  }
  
  missing_cols <- setdiff(ordered_cols, names(df_tex))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  rows <- apply(
    df_tex[, ordered_cols, drop = FALSE],
    1,
    function(r) paste(r, collapse = " & ")
  )
  
  body <- paste(rows, collapse = " \\\\\n")
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\renewcommand{\\arraystretch}{1.04}\n",
    "\\setlength{\\tabcolsep}{3.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{\\small ", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    body, " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n\n",
    "\\vspace{0.03cm}\n",
    "\\parbox{0.82\\linewidth}{\\centering\\small ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}

# ==============================================================================
# 5. LATEX: FORECAST-LEVEL SUMMARY
# ==============================================================================

make_latex_forecast_wins_summary <- function(
    overall_df,
    period_df,
    vintage_df,
    caption = "Forecast-level wins of Matrix MF--TPRF against the benchmarks",
    label   = "tab:forecast_wins_summary",
    note    = "Forecast-level wins and win shares (\\%) of Matrix MF--TPRF against each benchmark."
) {
  
  overall_tex <- overall_df
  period_tex  <- period_df
  vintage_tex <- vintage_df
  
  names(overall_tex)[1] <- "Group"
  names(period_tex)[1]  <- "Group"
  names(vintage_tex)[1] <- "Group"
  
  has_dfm <- all(c("DFM wins", "DFM win share") %in% names(overall_tex))
  
  format_block <- function(df_block) {
    df_block[[1]] <- escape_latex(as.character(df_block[[1]]))
    df_block <- format_win_table_values(df_block, df_block)
    
    if (has_dfm) {
      use_cols <- c(
        "Group", "Total",
        "VEC-P wins", "VEC-P win share",
        "VEC-C wins", "VEC-C win share",
        "DFM wins", "DFM win share"
      )
    } else {
      use_cols <- c(
        "Group", "Total",
        "VEC-P wins", "VEC-P win share",
        "VEC-C wins", "VEC-C win share"
      )
    }
    
    apply(
      df_block[, use_cols, drop = FALSE],
      1,
      function(r) paste(r, collapse = " & ")
    )
  }
  
  rows_overall <- format_block(overall_tex)
  rows_period  <- format_block(period_tex)
  rows_vintage <- format_block(vintage_tex)
  
  if (has_dfm) {
    tabular <- "c@{\\hspace{1.0em}}ccccccc"
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{2}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Group} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share} \\\\\n"
    )
    ncols <- 8
  } else {
    tabular <- "c@{\\hspace{1.0em}}ccccc"
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Group} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share} \\\\\n"
    )
    ncols <- 6
  }
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\renewcommand{\\arraystretch}{1.08}\n",
    "\\setlength{\\tabcolsep}{3.8pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{\\small ", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n\n",
    rows_overall[1], " \\\\\n",
    "\\addlinespace[0.20em]\n\n",
    "\\multicolumn{", ncols, "}{@{}l}{\\fontsize{8.4}{9.2}\\selectfont\\textbf{By evaluation period}}\\\\[-0.28em]\n",
    "\\cmidrule(lr){1-", ncols, "}\n",
    paste(rows_period, collapse = " \\\\\n"), " \\\\\n",
    "\\addlinespace[0.20em]\n\n",
    "\\multicolumn{", ncols, "}{@{}l}{\\fontsize{8.4}{9.2}\\selectfont\\textbf{By nowcast vintage}}\\\\[-0.28em]\n",
    "\\cmidrule(lr){1-", ncols, "}\n",
    paste(rows_vintage, collapse = " \\\\\n"), " \\\\\n\n",
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.03cm}\n",
    "\\parbox{0.82\\linewidth}{\\centering\\small ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}

# ==============================================================================
# 6. LATEX: GENERIC FORECAST WIN TABLE
# ==============================================================================

make_latex_forecast_win_table <- function(
    df,
    caption,
    label,
    first_col_name = NULL,
    note = "Entries report forecast-level wins of Matrix MF--TPRF against each benchmark. Win shares denote the percentage of forecasts for which Matrix MF--TPRF yields a lower squared forecast error than the corresponding benchmark."
) {
  
  df_tex <- df
  
  if (!is.null(first_col_name)) {
    names(df_tex)[1] <- first_col_name
  }
  
  first_col <- names(df_tex)[1]
  
  df_tex[[1]] <- escape_latex(as.character(df_tex[[1]]))
  df_tex <- format_win_table_values(df_tex, df)
  
  has_dfm <- all(c("DFM wins", "DFM win share") %in% names(df_tex))
  
  if (has_dfm) {
    use_cols <- c(
      first_col, "Total",
      "VEC-P wins", "VEC-P win share",
      "VEC-C wins", "VEC-C win share",
      "DFM wins", "DFM win share"
    )
    
    tabular <- "lccccccc"
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{2}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{", escape_latex(first_col), "} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share} \\\\\n"
    )
  } else {
    use_cols <- c(
      first_col, "Total",
      "VEC-P wins", "VEC-P win share",
      "VEC-C wins", "VEC-C win share"
    )
    
    tabular <- "lccccc"
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{", escape_latex(first_col), "} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Win share}",
      " & \\textbf{Wins} & \\textbf{Win share} \\\\\n"
    )
  }
  
  missing_cols <- setdiff(use_cols, names(df_tex))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  rows <- apply(
    df_tex[, use_cols, drop = FALSE],
    1,
    function(r) paste(r, collapse = " & ")
  )
  
  body <- paste(rows, collapse = " \\\\\n")
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\footnotesize\n",
    "\\renewcommand{\\arraystretch}{1.04}\n",
    "\\setlength{\\tabcolsep}{3.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{\\small ", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    body, " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.03cm}\n",
    "\\parbox{0.72\\linewidth}{\\centering\\fontsize{6.0}{6.8}\\selectfont ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}

make_latex_country_period_month_wins <- function(
    df,
    caption = "Forecast-level win shares of Matrix MF--TPRF by country, period, and nowcast vintage",
    label   = "tab:forecast_wins_country_period_month",
    note    = "Entries report the share of forecasts for which Matrix MF--TPRF has a lower squared forecast error than the benchmark.",
    period_date_map = NULL
) {
  
  df_tex <- df
  
  if (is.null(period_date_map)) {
    period_date_map <- stats::setNames(
      rep("", length(period_order)),
      period_order
    )
  }
  
  df_tex$Period  <- as.character(df_tex$Period)
  df_tex$Country <- escape_latex(as.character(df_tex$Country))
  
  share_cols <- setdiff(names(df_tex), c("Period", "Country"))
  
  for (nm in share_cols) {
    df_tex[[nm]] <- fmt_small_share(df[[nm]])
  }
  
  has_dfm <- all(c("DFM M1", "DFM M2", "DFM M3") %in% names(df_tex))
  
  if (has_dfm) {
    ordered_cols <- c(
      "Country",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3",
      "DFM M1", "DFM M2", "DFM M3"
    )
    
    tabular <- "c@{\\hspace{0.45cm}}ccc@{\\hspace{0.45cm}}ccc@{\\hspace{0.45cm}}ccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{3}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Country}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
    
    ncols <- 10
  } else {
    ordered_cols <- c(
      "Country",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3"
    )
    
    tabular <- "c@{\\hspace{0.55cm}}ccc@{\\hspace{0.55cm}}ccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Country}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
    
    ncols <- 7
  }
  
  make_block <- function(period_name) {
    
    block_df <- df_tex %>%
      dplyr::filter(Period == period_name)
    
    if (nrow(block_df) == 0L) {
      return("")
    }
    
    date_line <- period_date_map[[period_name]]
    if (is.null(date_line) || is.na(date_line)) {
      date_line <- ""
    }
    
    rows <- apply(
      block_df[, ordered_cols, drop = FALSE],
      1,
      function(r) paste(r, collapse = " & ")
    )
    
    paste0(
      "\\addlinespace[0.25em]\n",
      "\\specialrule{0.05em}{0.08em}{0.27em}\n",
      "\\multicolumn{", ncols, "}{@{}l}{\\scriptsize\\textbf{",
      escape_latex(period_name),
      "}} \\\\\n",
      "\\multicolumn{", ncols, "}{@{}l}{\\fontsize{6.2}{6.8}\\selectfont ",
      escape_latex(date_line),
      "} \\\\[-0.27em]\n",
      "\\cmidrule{1-", ncols, "}\n",
      paste(rows, collapse = " \\\\\n"),
      " \\\\\n"
    )
  }
  
  body <- paste(
    make_block("Pre-COVID"),
    make_block("COVID period"),
    make_block("Post-COVID"),
    sep = "\n"
  )
  
  paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{0.8}\n",
    "\\setlength{\\tabcolsep}{3.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\resizebox{0.98\\textwidth}{!}{%\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    body,
    "\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n\n",
    "\\vspace{0.04cm}\n",
    "\\parbox{0.86\\linewidth}{\\footnotesize\\centering ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}

# ==============================================================================
# utils_results_04_final_large_table.R
# Final large LaTeX table: Matrix MF-TPRF vs VEC-P vs VEC-C vs DFM
# ==============================================================================

get_rt_table_for_final <- function(summary_obj,
                                   country_order,
                                   period_order) {
  
  if (is.null(summary_obj$tab_rt_all)) {
    stop("summary_obj does not contain tab_rt_all.")
  }
  
  summary_obj$tab_rt_all %>%
    dplyr::filter(country %in% country_order, period %in% period_order) %>%
    dplyr::mutate(
      country = factor(country, levels = country_order),
      period  = factor(period, levels = period_order)
    ) %>%
    dplyr::arrange(period, country)
}

get_hyper_country_table <- function(summary_obj,
                                    which = c("pre", "post"),
                                    country_order) {
  
  which <- match.arg(which)
  
  obj_name <- if (which == "pre") {
    "hyper_rt_pre_all"
  } else {
    "hyper_rt_post_all"
  }
  
  if (is.null(summary_obj[[obj_name]])) {
    return(NULL)
  }
  
  summary_obj[[obj_name]] %>%
    dplyr::filter(country %in% country_order) %>%
    dplyr::mutate(country = factor(country, levels = country_order)) %>%
    dplyr::arrange(country)
}

hp_cells_tprf <- function(hp_df, cc) {
  
  if (is.null(hp_df)) {
    return(c("", "", "", ""))
  }
  
  row <- hp_df %>%
    dplyr::filter(country == cc)
  
  if (nrow(row) == 0L) {
    return(c("", "", "", ""))
  }
  
  c(
    as.character(row$Lproxy[1]),
    as.character(row$L_midas[1] - 1),
    as.character(row$p_AR[1]),
    as.character(row$r[1])
  )
}

hp_cells_dfm <- function(hp_df, cc, include_pq = TRUE) {
  
  if (is.null(hp_df)) {
    return(if (isTRUE(include_pq)) c("", "", "") else c("", ""))
  }
  
  row <- hp_df %>%
    dplyr::filter(country == cc)
  
  if (nrow(row) == 0L) {
    return(if (isTRUE(include_pq)) c("", "", "") else c("", ""))
  }
  
  r_val <- if ("r" %in% names(row)) {
    row$r[1]
  } else if ("K" %in% names(row)) {
    row$K[1]
  } else {
    NA
  }
  
  pf_val <- if ("p" %in% names(row)) {
    row$p[1]
  } else if ("p_hat" %in% names(row)) {
    row$p_hat[1]
  } else if ("P_f" %in% names(row)) {
    row$P_f[1]
  } else {
    NA
  }
  
  pq_val <- if ("q" %in% names(row)) {
    row$q[1]
  } else if ("q_hat" %in% names(row)) {
    row$q_hat[1]
  } else if ("P_q" %in% names(row)) {
    row$P_q[1]
  } else {
    NA
  }
  
  if (isTRUE(include_pq)) {
    c(
      ifelse(is.na(pf_val), "", as.character(pf_val)),
      ifelse(is.na(pq_val), "", as.character(pq_val)),
      ifelse(is.na(r_val),  "", as.character(r_val))
    )
  } else {
    c(
      ifelse(is.na(pf_val), "", as.character(pf_val)),
      ifelse(is.na(r_val),  "", as.character(r_val))
    )
  }
}

matrix_hp_text <- function(hp) {
  
  if (is.null(hp)) {
    return("")
  }
  
  r1 <- if (!is.null(hp$r_targeted)) hp$r_targeted[1] else NA
  r2 <- if (!is.null(hp$r_targeted)) hp$r_targeted[2] else NA
  
  paste0(
    "$L=", hp$Lproxy,
    ",\\; P_f=", hp$L_midas - 1,
    ",\\; P_\\rho=", hp$p_AR,
    ",\\; (\\hat r_1,\\hat r_2)=(",
    r1, ",", r2, ")$"
  )
}

build_final_large_style_latex <- function(
    summary_matrix,
    summary_vector,
    summary_vectensor,
    summary_dfm = NULL,
    params,
    dm_vecp_wide = NULL,
    dm_vecc_wide = NULL,
    dm_dfm_wide  = NULL,
    Size = "large",
    sel = "LASSO",
    include_dfm = !is.null(summary_dfm)
) {
  
  country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  period_order  <- c("Pre-COVID", "COVID period", "Post-COVID")
  period_labels <- build_period_labels(params)
  
  dfm_is_q0 <- isTRUE(include_dfm) && (
    (!is.null(summary_dfm$idio_spec) && identical(summary_dfm$idio_spec, "q0")) ||
      (!is.null(summary_dfm$params$idio_spec) && identical(summary_dfm$params$idio_spec, "q0"))
  )
  
  dfm_include_pq <- isTRUE(include_dfm) && !isTRUE(dfm_is_q0)
  
  mat_rt <- get_rt_table_for_final(summary_matrix, country_order, period_order)
  vec_rt <- get_rt_table_for_final(summary_vector, country_order, period_order)
  vtp_rt <- get_rt_table_for_final(summary_vectensor, country_order, period_order)
  
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
  
  if (isTRUE(include_dfm)) {
    dfm_rt <- get_rt_table_for_final(summary_dfm, country_order, period_order)
    
    df <- df %>%
      dplyr::left_join(
        dfm_rt %>%
          dplyr::rename(
            DFM_M1 = M1,
            DFM_M2 = M2,
            DFM_M3 = M3
          ),
        by = c("country", "period")
      ) %>%
      dplyr::mutate(
        Rel_DFM_M1 = Matrix_M1 / DFM_M1,
        Rel_DFM_M2 = Matrix_M2 / DFM_M2,
        Rel_DFM_M3 = Matrix_M3 / DFM_M3
      )
  }
  
  if (!is.null(dm_vecp_wide)) {
    df <- df %>%
      dplyr::left_join(dm_vecp_wide, by = c("country", "period"))
  }
  
  if (!is.null(dm_vecc_wide)) {
    df <- df %>%
      dplyr::left_join(dm_vecc_wide, by = c("country", "period"))
  }
  
  if (isTRUE(include_dfm) && !is.null(dm_dfm_wide)) {
    df <- df %>%
      dplyr::left_join(dm_dfm_wide, by = c("country", "period"))
  }
  
  add_missing_pcols <- function(data, prefix) {
    for (m in c("M1", "M2", "M3")) {
      nm <- paste0("p_", prefix, "_", m)
      if (!nm %in% names(data)) {
        data[[nm]] <- NA_real_
      }
    }
    data
  }
  
  df <- df %>%
    add_missing_pcols("VecP") %>%
    add_missing_pcols("VecC")
  
  if (isTRUE(include_dfm)) {
    df <- df %>%
      add_missing_pcols("DFM")
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
  
  if (isTRUE(include_dfm)) {
    df <- df %>%
      dplyr::mutate(
        star_DFM_M1 = ifelse(Rel_DFM_M1 < 1, p_to_stars(p_DFM_M1), ""),
        star_DFM_M2 = ifelse(Rel_DFM_M2 < 1, p_to_stars(p_DFM_M2), ""),
        star_DFM_M3 = ifelse(Rel_DFM_M3 < 1, p_to_stars(p_DFM_M3), "")
      )
  }
  
  hp_mat_pre  <- summary_matrix$hyper$pre
  hp_mat_post <- summary_matrix$hyper$post
  
  hp_vec_pre  <- get_hyper_country_table(summary_vector, "pre",  country_order)
  hp_vec_post <- get_hyper_country_table(summary_vector, "post", country_order)
  
  hp_vtp_pre  <- get_hyper_country_table(summary_vectensor, "pre",  country_order)
  hp_vtp_post <- get_hyper_country_table(summary_vectensor, "post", country_order)
  
  if (isTRUE(include_dfm)) {
    hp_dfm_pre  <- get_hyper_country_table(summary_dfm, "pre",  country_order)
    hp_dfm_post <- get_hyper_country_table(summary_dfm, "post", country_order)
  } else {
    hp_dfm_pre  <- NULL
    hp_dfm_post <- NULL
  }
  
  matrix_hp_pre_txt  <- matrix_hp_text(hp_mat_pre)
  matrix_hp_post_txt <- matrix_hp_text(hp_mat_post)
  
  fmt <- function(x) {
    ifelse(is.na(x), "", sprintf("%.3f", x))
  }
  
  fmt_star <- function(x, s = "") {
    if (length(x) == 0L || is.na(x)) {
      return("")
    }
    
    s <- ifelse(
      is.na(s) || s == "",
      "",
      paste0("\\sym{", s, "}")
    )
    
    paste0(sprintf("%.3f", x), s)
  }
  
  block_rows <- function(period_name,
                         period_label,
                         hp_matrix_txt = NULL,
                         hp_vp = NULL,
                         hp_vc = NULL,
                         hp_dfm = NULL) {
    
    out <- c()
    
    line1 <- paste0(
      "{\\fontsize{5.8}{6.3}\\selectfont\\textbf{",
      period_label["title"],
      "}}"
    )
    
    line2 <- paste0(
      "{\\fontsize{4.5}{4.8}\\selectfont ",
      period_label["date"],
      "}"
    )
    
    ncols <- if (isTRUE(include_dfm)) {
      if (isTRUE(dfm_include_pq)) 24 else 23
    } else {
      18
    }
    
    out <- c(
      out,
      "\\specialrule{0.08em}{0.15em}{0.08em}",
      paste0("\\multicolumn{", ncols, "}{l}{", line1, "}\\\\[-0.60em]"),
      paste0("\\multicolumn{", ncols, "}{l}{", line2, "}\\\\[-0.32em]"),
      paste0("\\cmidrule{1-", ncols, "}")
    )
    
    if (!is.null(hp_matrix_txt)) {
      
      if (isTRUE(include_dfm)) {
        
        dfm_hp_header <- if (isTRUE(dfm_include_pq)) {
          paste0(
            " & \\cellcolor{hypergray}\\hphdr $P_f$",
            " & \\cellcolor{hypergray}\\hphdr $r$",
            " & \\cellcolor{hypergray}\\hphdr $P_q$"
          )
        } else {
          paste0(
            " & \\cellcolor{hypergray}\\hphdr $P_f$",
            " & \\cellcolor{hypergray}\\hphdr $r$"
          )
        }
        
        hp_line <- paste0(
          "\\rowcolor{hypergray}\n",
          "\\textit{\\fontsize{5.0}{5.8}\\selectfont Hyper-params}",
          " & \\multicolumn{3}{@{}>{\\columncolor{hypergray}}c@{\\hspace{0.18cm}}}{\\fontsize{4.8}{5.5}\\selectfont ",
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
          " & \\multicolumn{3}{c}{}",
          dfm_hp_header,
          "\\\\[0.25em]"
        )
        
      } else {
        
        hp_line <- paste0(
          "\\rowcolor{hypergray}\n",
          "\\textit{\\fontsize{5.0}{5.8}\\selectfont Hyper-params}",
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
      }
      
      out <- c(out, hp_line)
    }
    
    df_p <- df %>%
      dplyr::filter(period == period_name)
    
    for (cc in country_order) {
      
      row <- df_p %>%
        dplyr::filter(country == cc)
      
      if (nrow(row) == 0L) {
        next
      }
      
      hpvp  <- hp_cells_tprf(hp_vp, cc)
      hpvc  <- hp_cells_tprf(hp_vc, cc)
      hpdfm <- hp_cells_dfm(hp_dfm, cc, include_pq = dfm_include_pq)
      
      row_line <- paste0(
        cc, " & ",
        "\\matcell{", fmt(row$Matrix_M1), "} & ",
        "\\matcell{", fmt(row$Matrix_M2), "} & ",
        "\\matcell{", fmt(row$Matrix_M3), "} & ",
        
        fmt_star(row$Rel_VecP_M1, row$star_VecP_M1), " & ",
        fmt_star(row$Rel_VecP_M2, row$star_VecP_M2), " & ",
        fmt_star(row$Rel_VecP_M3, row$star_VecP_M3), " & ",
        "\\hpstyle ", hpvp[1], " & ",
        "\\hpstyle ", hpvp[2], " & ",
        "\\hpstyle ", hpvp[3], " & ",
        "\\hpstyle ", hpvp[4], " & ",
        
        fmt_star(row$Rel_VecC_M1, row$star_VecC_M1), " & ",
        fmt_star(row$Rel_VecC_M2, row$star_VecC_M2), " & ",
        fmt_star(row$Rel_VecC_M3, row$star_VecC_M3), " & ",
        "\\hpstyle ", hpvc[1], " & ",
        "\\hpstyle ", hpvc[2], " & ",
        "\\hpstyle ", hpvc[3], " & ",
        "\\hpstyle ", hpvc[4]
      )
      
      if (isTRUE(include_dfm)) {
        if (isTRUE(dfm_include_pq)) {
          row_line <- paste0(
            row_line, " & ",
            fmt_star(row$Rel_DFM_M1, row$star_DFM_M1), " & ",
            fmt_star(row$Rel_DFM_M2, row$star_DFM_M2), " & ",
            fmt_star(row$Rel_DFM_M3, row$star_DFM_M3), " & ",
            "\\hpstyle ", hpdfm[1], " & ",
            "\\hpstyle ", hpdfm[2], " & ",
            "\\hpstyle ", hpdfm[3]
          )
        } else {
          row_line <- paste0(
            row_line, " & ",
            fmt_star(row$Rel_DFM_M1, row$star_DFM_M1), " & ",
            fmt_star(row$Rel_DFM_M2, row$star_DFM_M2), " & ",
            fmt_star(row$Rel_DFM_M3, row$star_DFM_M3), " & ",
            "\\hpstyle ", hpdfm[1], " & ",
            "\\hpstyle ", hpdfm[2]
          )
        }
      }
      
      row_line <- paste0(row_line, " \\\\")
      
      out <- c(out, row_line)
    }
    
    out
  }
  
  size_caption <- paste0(tolower(Size), " information set")
  
  sel_caption <- ifelse(
    toupper(sel) == "LASSO",
    "LASSO preselection",
    "correlation preselection"
  )
  
  if (isTRUE(include_dfm)) {
    
    if (isTRUE(dfm_include_pq)) {
      dfm_tabular_spec <- c(
        "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
        "*{3}{>{\\centering\\arraybackslash}p{0.13cm}}"
      )
      dfm_top_header <- " & \\multicolumn{6}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{DFM}}$}}"
      dfm_second_header <- "& M1 & M2 & M3 & \\multicolumn{3}{c}{HP}"
      dfm_cmidrule <- "\\cmidrule(lr){19-24}"
    } else {
      dfm_tabular_spec <- c(
        "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
        "*{2}{>{\\centering\\arraybackslash}p{0.13cm}}"
      )
      dfm_top_header <- " & \\multicolumn{5}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{DFM}}$}}"
      dfm_second_header <- "& M1 & M2 & M3 & \\multicolumn{2}{c}{HP}"
      dfm_cmidrule <- "\\cmidrule(lr){19-23}"
    }
    
    table_note <- paste0(
      "The shaded block reports the Matrix MF--TPRF RMSFE for M1, M2, and M3 nowcasts. ",
      "The remaining blocks report relative RMSFE with respect to VEC-P, VEC-C, and DFM; ",
      "values below one favour Matrix MF--TPRF. ",
      "Asterisks denote one-sided Diebold--Mariano significance in favour of Matrix MF--TPRF ",
      "based on squared forecast errors and a Newey--West HAC variance estimator ",
      "($^{*}\\, p\\!<\\!0.10$, $^{**}\\, p\\!<\\!0.05$, $^{***}\\, p\\!<\\!0.01$)."
    )
    
    tabular_spec <- c(
      ">{\\centering\\arraybackslash}p{1.05cm}",
      "@{\\hspace{0.16cm}}",
      ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.82cm}",
      ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.82cm}",
      ">{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.82cm}",
      "@{\\hspace{0.16cm}}",
      "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
      "*{4}{>{\\centering\\arraybackslash}p{0.13cm}}",
      "@{\\hspace{0.16cm}}",
      "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
      "*{4}{>{\\centering\\arraybackslash}p{0.13cm}}",
      "@{\\hspace{0.16cm}}",
      dfm_tabular_spec
    )
    
    top_header <- paste0(
      "\\multicolumn{1}{>{\\cellcolor{topgray}}c@{\\hspace{0.16cm}}}{\\textbf{", toupper(Size), "}}",
      " & \\multicolumn{3}{@{}>{\\columncolor{matrixgray}}c@{\\hspace{0.16cm}}}{\\textbf{Matrix MF--TPRF}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-P MF--TPRF}}$}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-C MF--TPRF}}$}}",
      dfm_top_header
    )
    
    cmidrule_line <- paste0(
      "\\cmidrule(lr){2-4}",
      "\\cmidrule(lr){5-11}",
      "\\cmidrule(lr){12-18}",
      dfm_cmidrule
    )
    
    second_header <- paste0(
      "\\textbf{Country}",
      "& \\cellcolor{matrixgray}M1 & \\cellcolor{matrixgray}M2 & \\cellcolor{matrixgray}M3",
      "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
      "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
      dfm_second_header
    )
    
  } else {
    
    table_note <- paste0(
      "The shaded block reports the Matrix MF--TPRF RMSFE for M1, M2, and M3 nowcasts. ",
      "The remaining blocks report relative RMSFE with respect to VEC-P and VEC-C; ",
      "values below one favour Matrix MF--TPRF. ",
      "Asterisks denote one-sided Diebold--Mariano significance in favour of Matrix MF--TPRF ",
      "based on squared forecast errors and a Newey--West HAC variance estimator ",
      "($^{*}\\, p\\!<\\!0.10$, $^{**}\\, p\\!<\\!0.05$, $^{***}\\, p\\!<\\!0.01$)."
    )
    
    tabular_spec <- c(
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
      "*{4}{>{\\centering\\arraybackslash}p{0.14cm}}"
    )
    
    top_header <- paste0(
      "\\multicolumn{1}{>{\\cellcolor{topgray}}c@{\\hspace{0.22cm}}}{\\textbf{", toupper(Size), "}}",
      " & \\multicolumn{3}{@{}>{\\columncolor{matrixgray}}c@{\\hspace{0.22cm}}}{\\textbf{Matrix MF--TPRF}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-P MF--TPRF}}$}}",
      " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{Matrix MF--TPRF}}{\\text{VEC-C MF--TPRF}}$}}"
    )
    
    cmidrule_line <- "\\cmidrule(lr){2-4}\\cmidrule(lr){5-11}\\cmidrule(lr){12-18}"
    
    second_header <- paste0(
      "\\textbf{Country}",
      "& \\cellcolor{matrixgray}M1 & \\cellcolor{matrixgray}M2 & \\cellcolor{matrixgray}M3",
      "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
      "& M1 & M2 & M3 & \\multicolumn{4}{c}{HP}"
    )
  }
  
  latex_lines <- c(
    "%=========================================================",
    paste0("% ", toupper(Size), " INFORMATION SET -- ", toupper(sel)),
    "%=========================================================",
    "\\begin{table}[!ht]",
    "\\centering",
    "\\scriptsize",
    "\\renewcommand{\\arraystretch}{1.06}",
    "\\setlength{\\tabcolsep}{1.25pt}",
    "",
    "\\newcommand{\\hpstyle}{\\fontsize{4.2}{4.8}\\selectfont}",
    "\\newcommand{\\hphdr}{\\fontsize{4.3}{5.0}\\selectfont}",
    "\\newcommand{\\matcell}[1]{\\cellcolor{matrixgray}#1}",
    "\\newcommand{\\sym}[1]{\\textsuperscript{\\fontsize{3.8}{3.8}\\selectfont#1}}",
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
    paste(tabular_spec, collapse = ""),
    "}",
    "\\toprule",
    top_header,
    "\\\\",
    cmidrule_line,
    second_header,
    "\\\\",
    block_rows(
      "Pre-COVID",
      period_label  = period_labels[["Pre-COVID"]],
      hp_matrix_txt = matrix_hp_pre_txt,
      hp_vp         = hp_vtp_pre,
      hp_vc         = hp_vec_pre,
      hp_dfm        = hp_dfm_pre
    ),
    block_rows(
      "COVID period",
      period_label  = period_labels[["COVID period"]],
      hp_matrix_txt = NULL,
      hp_vp         = hp_vtp_pre,
      hp_vc         = hp_vec_pre,
      hp_dfm        = hp_dfm_pre
    ),
    block_rows(
      "Post-COVID",
      period_label  = period_labels[["Post-COVID"]],
      hp_matrix_txt = matrix_hp_post_txt,
      hp_vp         = hp_vtp_post,
      hp_vc         = hp_vec_post,
      hp_dfm        = hp_dfm_post
    ),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.15cm}",
    paste0("\\parbox{0.96\\textwidth}{\\footnotesize ", table_note, "}"),
    "\\end{table}"
  )
  
  paste(latex_lines, collapse = "\n")
}


build_final_vecc_style_latex <- function(
    summary_vector,
    summary_vectensor,
    summary_dfm,
    params,
    dm_vecp_wide = NULL,
    dm_dfm_wide  = NULL,
    Size = "large",
    sel = "LASSO"
) {
  
  country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  period_order  <- c("Pre-COVID", "COVID period", "Post-COVID")
  period_labels <- build_period_labels(params)
  
  dfm_is_q0 <- (
    (!is.null(summary_dfm$idio_spec) && identical(summary_dfm$idio_spec, "q0")) ||
      (!is.null(summary_dfm$params$idio_spec) && identical(summary_dfm$params$idio_spec, "q0"))
  )
  
  dfm_include_pq <- !isTRUE(dfm_is_q0)
  
  vecc_rt <- get_rt_table_for_final(summary_vector,    country_order, period_order)
  vecp_rt <- get_rt_table_for_final(summary_vectensor, country_order, period_order)
  dfm_rt  <- get_rt_table_for_final(summary_dfm,       country_order, period_order)
  
  df <- vecc_rt %>%
    dplyr::rename(
      VecC_M1 = M1,
      VecC_M2 = M2,
      VecC_M3 = M3
    ) %>%
    dplyr::left_join(
      vecp_rt %>%
        dplyr::rename(
          VecP_M1 = M1,
          VecP_M2 = M2,
          VecP_M3 = M3
        ),
      by = c("country", "period")
    ) %>%
    dplyr::left_join(
      dfm_rt %>%
        dplyr::rename(
          DFM_M1 = M1,
          DFM_M2 = M2,
          DFM_M3 = M3
        ),
      by = c("country", "period")
    ) %>%
    dplyr::mutate(
      Rel_VecP_M1 = VecC_M1 / VecP_M1,
      Rel_VecP_M2 = VecC_M2 / VecP_M2,
      Rel_VecP_M3 = VecC_M3 / VecP_M3,
      
      Rel_DFM_M1 = VecC_M1 / DFM_M1,
      Rel_DFM_M2 = VecC_M2 / DFM_M2,
      Rel_DFM_M3 = VecC_M3 / DFM_M3
    )
  
  if (!is.null(dm_vecp_wide)) {
    df <- df %>%
      dplyr::left_join(dm_vecp_wide, by = c("country", "period"))
  }
  
  if (!is.null(dm_dfm_wide)) {
    df <- df %>%
      dplyr::left_join(dm_dfm_wide, by = c("country", "period"))
  }
  
  add_missing_pcols <- function(data, prefix) {
    for (m in c("M1", "M2", "M3")) {
      nm <- paste0("p_", prefix, "_", m)
      if (!nm %in% names(data)) {
        data[[nm]] <- NA_real_
      }
    }
    data
  }
  
  df <- df %>%
    add_missing_pcols("VecP") %>%
    add_missing_pcols("DFM") %>%
    dplyr::mutate(
      star_VecP_M1 = ifelse(Rel_VecP_M1 < 1, p_to_stars(p_VecP_M1), ""),
      star_VecP_M2 = ifelse(Rel_VecP_M2 < 1, p_to_stars(p_VecP_M2), ""),
      star_VecP_M3 = ifelse(Rel_VecP_M3 < 1, p_to_stars(p_VecP_M3), ""),
      
      star_DFM_M1 = ifelse(Rel_DFM_M1 < 1, p_to_stars(p_DFM_M1), ""),
      star_DFM_M2 = ifelse(Rel_DFM_M2 < 1, p_to_stars(p_DFM_M2), ""),
      star_DFM_M3 = ifelse(Rel_DFM_M3 < 1, p_to_stars(p_DFM_M3), "")
    )
  
  hp_vecc_pre  <- get_hyper_country_table(summary_vector, "pre",  country_order)
  hp_vecc_post <- get_hyper_country_table(summary_vector, "post", country_order)
  
  hp_vecp_pre  <- get_hyper_country_table(summary_vectensor, "pre",  country_order)
  hp_vecp_post <- get_hyper_country_table(summary_vectensor, "post", country_order)
  
  hp_dfm_pre  <- get_hyper_country_table(summary_dfm, "pre",  country_order)
  hp_dfm_post <- get_hyper_country_table(summary_dfm, "post", country_order)
  
  fmt <- function(x) {
    ifelse(is.na(x), "", sprintf("%.3f", x))
  }
  
  fmt_star <- function(x, s = "") {
    if (length(x) == 0L || is.na(x)) {
      return("")
    }
    
    s <- ifelse(
      is.na(s) || s == "",
      "",
      paste0("\\sym{", s, "}")
    )
    
    paste0(sprintf("%.3f", x), s)
  }
  
  block_rows <- function(period_name,
                         period_label,
                         hp_vc = NULL,
                         hp_vp = NULL,
                         hp_dfm = NULL) {
    
    out <- c()
    
    line1 <- paste0(
      "{\\fontsize{5.8}{6.3}\\selectfont\\textbf{",
      period_label["title"],
      "}}"
    )
    
    line2 <- paste0(
      "{\\fontsize{4.5}{4.8}\\selectfont ",
      period_label["date"],
      "}"
    )
    
    ncols <- if (isTRUE(dfm_include_pq)) 21 else 20
    
    out <- c(
      out,
      "\\specialrule{0.08em}{0.15em}{0.08em}",
      paste0("\\multicolumn{", ncols, "}{l}{", line1, "}\\\\[-0.60em]"),
      paste0("\\multicolumn{", ncols, "}{l}{", line2, "}\\\\[-0.32em]"),
      paste0("\\cmidrule{1-", ncols, "}")
    )
    
    dfm_hp_header <- if (isTRUE(dfm_include_pq)) {
      paste0(
        " & \\cellcolor{hypergray}\\hphdr $P_f$",
        " & \\cellcolor{hypergray}\\hphdr $r$",
        " & \\cellcolor{hypergray}\\hphdr $P_q$"
      )
    } else {
      paste0(
        " & \\cellcolor{hypergray}\\hphdr $P_f$",
        " & \\cellcolor{hypergray}\\hphdr $r$"
      )
    }
    
    hp_line <- paste0(
      "\\rowcolor{hypergray}\n",
      "\\textit{\\fontsize{5.0}{5.8}\\selectfont Hyper-params}",
      
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
      
      " & \\multicolumn{3}{c}{}",
      dfm_hp_header,
      "\\\\[0.25em]"
    )
    
    out <- c(out, hp_line)
    
    df_p <- df %>%
      dplyr::filter(period == period_name)
    
    for (cc in country_order) {
      
      row <- df_p %>%
        dplyr::filter(country == cc)
      
      if (nrow(row) == 0L) {
        next
      }
      
      hpvc  <- hp_cells_tprf(hp_vc, cc)
      hpvp  <- hp_cells_tprf(hp_vp, cc)
      hpdfm <- hp_cells_dfm(hp_dfm, cc, include_pq = dfm_include_pq)
      
      dfm_cells <- if (isTRUE(dfm_include_pq)) {
        paste0(
          fmt_star(row$Rel_DFM_M1, row$star_DFM_M1), " & ",
          fmt_star(row$Rel_DFM_M2, row$star_DFM_M2), " & ",
          fmt_star(row$Rel_DFM_M3, row$star_DFM_M3), " & ",
          "\\hpstyle ", hpdfm[1], " & ",
          "\\hpstyle ", hpdfm[2], " & ",
          "\\hpstyle ", hpdfm[3]
        )
      } else {
        paste0(
          fmt_star(row$Rel_DFM_M1, row$star_DFM_M1), " & ",
          fmt_star(row$Rel_DFM_M2, row$star_DFM_M2), " & ",
          fmt_star(row$Rel_DFM_M3, row$star_DFM_M3), " & ",
          "\\hpstyle ", hpdfm[1], " & ",
          "\\hpstyle ", hpdfm[2]
        )
      }
      
      row_line <- paste0(
        cc, " & ",
        
        "\\matcell{", fmt(row$VecC_M1), "} & ",
        "\\matcell{", fmt(row$VecC_M2), "} & ",
        "\\matcell{", fmt(row$VecC_M3), "} & ",
        "\\hpstyle ", hpvc[1], " & ",
        "\\hpstyle ", hpvc[2], " & ",
        "\\hpstyle ", hpvc[3], " & ",
        "\\hpstyle ", hpvc[4], " & ",
        
        fmt_star(row$Rel_VecP_M1, row$star_VecP_M1), " & ",
        fmt_star(row$Rel_VecP_M2, row$star_VecP_M2), " & ",
        fmt_star(row$Rel_VecP_M3, row$star_VecP_M3), " & ",
        "\\hpstyle ", hpvp[1], " & ",
        "\\hpstyle ", hpvp[2], " & ",
        "\\hpstyle ", hpvp[3], " & ",
        "\\hpstyle ", hpvp[4], " & ",
        
        dfm_cells,
        " \\\\"
      )
      
      out <- c(out, row_line)
    }
    
    out
  }
  
  size_caption <- paste0(tolower(Size), " information set")
  
  sel_caption <- ifelse(
    toupper(sel) == "LASSO",
    "LASSO preselection",
    "correlation preselection"
  )
  
  if (isTRUE(dfm_include_pq)) {
    dfm_tabular_spec <- c(
      "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
      "*{3}{>{\\centering\\arraybackslash}p{0.13cm}}"
    )
    dfm_top_header <- " & \\multicolumn{6}{c}{\\textbf{$\\dfrac{\\text{VEC-C MF--TPRF}}{\\text{DFM}}$}}"
    dfm_second_header <- "& M1 & M2 & M3 & \\multicolumn{3}{c}{HP}"
    dfm_cmidrule <- "\\cmidrule(lr){16-21}"
  } else {
    dfm_tabular_spec <- c(
      "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
      "*{2}{>{\\centering\\arraybackslash}p{0.13cm}}"
    )
    dfm_top_header <- " & \\multicolumn{5}{c}{\\textbf{$\\dfrac{\\text{VEC-C MF--TPRF}}{\\text{DFM}}$}}"
    dfm_second_header <- "& M1 & M2 & M3 & \\multicolumn{2}{c}{HP}"
    dfm_cmidrule <- "\\cmidrule(lr){16-20}"
  }
  
  table_note <- paste0(
    "The shaded block reports the VEC-C RMSFE for M1, M2, and M3 nowcasts, ",
    "together with the selected VEC-C hyperparameters. ",
    "The remaining blocks report relative RMSFE with respect to VEC-P and DFM; ",
    "values below one favour VEC-C. ",
    "Asterisks denote one-sided Diebold--Mariano significance in favour of VEC-C ",
    "based on squared forecast errors and a Newey--West HAC variance estimator ",
    "($^{*}\\, p\\!<\\!0.10$, $^{**}\\, p\\!<\\!0.05$, $^{***}\\, p\\!<\\!0.01$). ")
  
  tabular_spec <- c(
    ">{\\centering\\arraybackslash}p{1.05cm}",
    "@{\\hspace{0.16cm}}",
    
    "*{3}{>{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.82cm}}",
    "*{4}{>{\\columncolor{matrixgray}\\centering\\arraybackslash}p{0.13cm}}",
    
    "@{\\hspace{0.16cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{0.72cm}}",
    "*{4}{>{\\centering\\arraybackslash}p{0.13cm}}",
    
    "@{\\hspace{0.16cm}}",
    dfm_tabular_spec
  )
  
  top_header <- paste0(
    "\\multicolumn{1}{>{\\cellcolor{topgray}}c@{\\hspace{0.16cm}}}{\\textbf{", toupper(Size), "}}",
    " & \\multicolumn{7}{@{}>{\\columncolor{matrixgray}}c@{\\hspace{0.16cm}}}{\\textbf{VEC-C MF--TPRF}}",
    " & \\multicolumn{7}{c}{\\textbf{$\\dfrac{\\text{VEC-C MF--TPRF}}{\\text{VEC-P MF--TPRF}}$}}",
    dfm_top_header
  )
  
  cmidrule_line <- paste0(
    "\\cmidrule(lr){2-8}",
    "\\cmidrule(lr){9-15}",
    dfm_cmidrule
  )
  
  second_header <- paste0(
    "\\textbf{Country}",
    " & \\cellcolor{matrixgray}M1",
    " & \\cellcolor{matrixgray}M2",
    " & \\cellcolor{matrixgray}M3",
    " & \\multicolumn{4}{>{\\columncolor{matrixgray}}c}{HP}",
    " & M1 & M2 & M3 & \\multicolumn{4}{c}{HP}",
    dfm_second_header
  )
  
  latex_lines <- c(
    "%=========================================================",
    paste0("% VEC-C BENCHMARK TABLE -- ", toupper(Size), " INFORMATION SET -- ", toupper(sel)),
    "%=========================================================",
    "\\begin{table}[!ht]",
    "\\centering",
    "\\scriptsize",
    "\\renewcommand{\\arraystretch}{1.06}",
    "\\setlength{\\tabcolsep}{1.25pt}",
    "",
    "\\newcommand{\\hpstyle}{\\fontsize{4.2}{4.8}\\selectfont}",
    "\\newcommand{\\hphdr}{\\fontsize{4.3}{5.0}\\selectfont}",
    "\\newcommand{\\matcell}[1]{\\cellcolor{matrixgray}#1}",
    "\\newcommand{\\sym}[1]{\\textsuperscript{\\fontsize{3.8}{3.8}\\selectfont#1}}",
    "",
    "\\definecolor{topgray}{gray}{0.86}",
    "\\definecolor{hypergray}{gray}{0.91}",
    "\\definecolor{matrixgray}{gray}{0.965}",
    "",
    paste0(
      "\\caption{VEC-C MF--TPRF and relative RMSFE by country: ",
      size_caption, " \\& ", sel_caption, "}"
    ),
    paste0("\\label{tab:EA_nowcast_vecc_", tolower(Size), "_", tolower(sel), "}"),
    "",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{",
    paste(tabular_spec, collapse = ""),
    "}",
    "\\toprule",
    top_header,
    "\\\\",
    cmidrule_line,
    second_header,
    "\\\\",
    block_rows(
      "Pre-COVID",
      period_label = period_labels[["Pre-COVID"]],
      hp_vc        = hp_vecc_pre,
      hp_vp        = hp_vecp_pre,
      hp_dfm       = hp_dfm_pre
    ),
    block_rows(
      "COVID period",
      period_label = period_labels[["COVID period"]],
      hp_vc        = hp_vecc_pre,
      hp_vp        = hp_vecp_pre,
      hp_dfm       = hp_dfm_pre
    ),
    block_rows(
      "Post-COVID",
      period_label = period_labels[["Post-COVID"]],
      hp_vc        = hp_vecc_post,
      hp_vp        = hp_vecp_post,
      hp_dfm       = hp_dfm_post
    ),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\vspace{0.15cm}",
    paste0("\\parbox{0.96\\textwidth}{\\footnotesize ", table_note, "}"),
    "\\end{table}"
  )
  
  paste(latex_lines, collapse = "\n")
}

# ==============================================================================
# utils_results_05_plots_and_wins.R
# Country plots, RMSFE-level wins, forecast-level wins, LaTeX generators
# ==============================================================================

extract_country_rt <- function(df_rt, country_code, model_label) {
  df_rt %>%
    dplyr::filter(country == country_code) %>%
    dplyr::mutate(
      model = model_label,
      type  = factor(type, levels = c("M1", "M2", "M3"))
    ) %>%
    dplyr::select(date, country, nowcast, type, model)
}

extract_country_gdp <- function(df_yq, country_code) {
  df_yq %>%
    dplyr::filter(country == country_code) %>%
    dplyr::select(date, country, GDP)
}

extract_country_rt_multi <- function(df_rt, country_codes, model_label) {
  df_rt %>%
    dplyr::filter(country %in% country_codes) %>%
    dplyr::mutate(
      model   = model_label,
      type    = factor(type, levels = c("M1", "M2", "M3")),
      country = factor(country, levels = country_codes)
    ) %>%
    dplyr::select(date, country, nowcast, type, model)
}

extract_country_gdp_multi <- function(df_yq, country_codes) {
  df_yq %>%
    dplyr::filter(country %in% country_codes) %>%
    dplyr::mutate(country = factor(country, levels = country_codes)) %>%
    dplyr::select(date, country, GDP)
}

theme_country_compare <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25"),
      
      strip.background = ggplot2::element_rect(
        fill      = "grey88",
        colour    = "grey65",
        linewidth = 0.4
      ),
      strip.text = ggplot2::element_text(
        face   = "bold",
        colour = "black"
      )
    )
}

get_gdp_ylim <- function(df_gdp, pad_frac = 0.02, min_pad = 0.001) {
  
  rng <- range(df_gdp$GDP, na.rm = TRUE)
  pad <- max(diff(rng) * pad_frac, min_pad)
  
  c(rng[1] - pad, rng[2] + pad)
}
# ==============================================================================
# RMSFE-LEVEL WIN SUMMARY TABLES
# ==============================================================================

reshape_comp_to_long <- function(df_comp,
                                 include_dfm = TRUE,
                                 country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
                                 month_order = c("M1", "M2", "M3")) {
  
  cols <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  if (isTRUE(include_dfm) && all(c("DFM_M1", "DFM_M2", "DFM_M3") %in% names(df_comp))) {
    cols <- c(cols, "DFM_M1", "DFM_M2", "DFM_M3")
  }
  
  df_comp %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(cols),
      names_to  = c("model", "month"),
      names_sep = "_",
      values_to = "rmsfe"
    ) %>%
    tidyr::pivot_wider(
      names_from  = model,
      values_from = rmsfe
    ) %>%
    dplyr::mutate(
      month   = factor(month, levels = month_order),
      period  = factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")),
      country = factor(country, levels = country_order)
    )
}

compute_summary_metrics <- function(df_long) {
  out <- df_long %>%
    dplyr::mutate(
      win_vecp = as.integer(Matrix < VecTensor),
      win_vecc = as.integer(Matrix < Vector)
    )
  
  if ("DFM" %in% names(out)) {
    out <- out %>%
      dplyr::mutate(win_dfm = as.integer(Matrix < DFM))
  }
  
  out
}

build_period_month_table <- function(df_eval,
                                     include_dfm = TRUE,
                                     period_order = c("Pre-COVID", "COVID period", "Post-COVID"),
                                     month_order = c("M1", "M2", "M3")) {
  
  out <- df_eval %>%
    dplyr::filter(period %in% period_order) %>%
    dplyr::group_by(period, month) %>%
    dplyr::summarise(
      share_vecp = 100 * mean(win_vecp, na.rm = TRUE),
      share_vecc = 100 * mean(win_vecc, na.rm = TRUE),
      share_dfm  = if ("win_dfm" %in% names(.)) 100 * mean(win_dfm, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      period = factor(period, levels = period_order),
      month  = factor(month, levels = month_order)
    ) %>%
    dplyr::arrange(period, month) %>%
    tidyr::pivot_wider(
      names_from  = month,
      values_from = c(share_vecp, share_vecc, share_dfm),
      names_glue  = "{.value}_{month}"
    ) %>%
    dplyr::rename(
      Period      = period,
      `VEC-P M1`  = share_vecp_M1,
      `VEC-P M2`  = share_vecp_M2,
      `VEC-P M3`  = share_vecp_M3,
      `VEC-C M1`  = share_vecc_M1,
      `VEC-C M2`  = share_vecc_M2,
      `VEC-C M3`  = share_vecc_M3,
      `DFM M1`    = share_dfm_M1,
      `DFM M2`    = share_dfm_M2,
      `DFM M3`    = share_dfm_M3
    ) %>%
    dplyr::mutate(dplyr::across(-Period, ~ round(.x, 1)))
  
  if (!isTRUE(include_dfm) || !all(c("DFM M1", "DFM M2", "DFM M3") %in% names(out))) {
    out <- out %>%
      dplyr::select(-dplyr::any_of(c("DFM M1", "DFM M2", "DFM M3")))
  }
  
  out
}

# ==============================================================================
# FORECAST-LEVEL WIN TABLES
# ==============================================================================

build_forecast_win_table <- function(df,
                                     group_vars = NULL,
                                     include_dfm = TRUE) {
  
  has_dfm <- isTRUE(include_dfm) && "win_dfm" %in% names(df)
  
  summarise_fun <- function(data) {
    out <- data %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (has_dfm) {
      out <- data %>%
        dplyr::summarise(
          Total             = dplyr::n(),
          `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
          `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
          `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
          `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
          `DFM wins`        = sum(win_dfm, na.rm = TRUE),
          `DFM win share`   = 100 * mean(win_dfm, na.rm = TRUE),
          .groups = "drop"
        )
    }
    
    out
  }
  
  if (is.null(group_vars)) {
    
    out <- summarise_fun(df) %>%
      dplyr::mutate(Group = "Overall", .before = 1)
    
  } else {
    
    out <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      summarise_fun()
  }
  
  out %>%
    dplyr::mutate(
      Total             = as.integer(Total),
      `VEC-P wins`      = as.integer(`VEC-P wins`),
      `VEC-C wins`      = as.integer(`VEC-C wins`),
      `VEC-P win share` = round(`VEC-P win share`, 1),
      `VEC-C win share` = round(`VEC-C win share`, 1),
      dplyr::across(dplyr::any_of("DFM wins"), as.integer),
      dplyr::across(dplyr::any_of("DFM win share"), ~ round(.x, 1))
    )
}

build_country_period_month_table <- function(df,
                                             include_dfm = TRUE,
                                             country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
                                             period_order = c("Pre-COVID", "COVID period", "Post-COVID"),
                                             month_order = c("M1", "M2", "M3")) {
  
  out <- df %>%
    dplyr::filter(period %in% period_order) %>%
    dplyr::group_by(period, country, type) %>%
    dplyr::summarise(
      share_vecp = 100 * mean(win_vecp, na.rm = TRUE),
      share_vecc = 100 * mean(win_vecc, na.rm = TRUE),
      share_dfm  = if ("win_dfm" %in% names(.)) 100 * mean(win_dfm, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      period  = factor(period, levels = period_order),
      country = factor(country, levels = country_order),
      type    = factor(type, levels = month_order)
    ) %>%
    dplyr::arrange(period, country, type) %>%
    tidyr::pivot_wider(
      names_from  = type,
      values_from = c(share_vecp, share_vecc, share_dfm),
      names_glue  = "{.value}_{type}"
    ) %>%
    dplyr::rename(
      Period      = period,
      Country     = country,
      `VEC-P M1`  = share_vecp_M1,
      `VEC-P M2`  = share_vecp_M2,
      `VEC-P M3`  = share_vecp_M3,
      `VEC-C M1`  = share_vecc_M1,
      `VEC-C M2`  = share_vecc_M2,
      `VEC-C M3`  = share_vecc_M3,
      `DFM M1`    = share_dfm_M1,
      `DFM M2`    = share_dfm_M2,
      `DFM M3`    = share_dfm_M3
    ) %>%
    dplyr::mutate(dplyr::across(-c(Period, Country), ~ round(.x, 1)))
  
  if (!isTRUE(include_dfm) || !all(c("DFM M1", "DFM M2", "DFM M3") %in% names(out))) {
    out <- out %>%
      dplyr::select(-dplyr::any_of(c("DFM M1", "DFM M2", "DFM M3")))
  }
  
  out
}

# ==============================================================================
# LATEX HELPERS FOR WIN TABLES
# ==============================================================================

fmt_int <- function(x) {
  ifelse(is.na(x), "--", as.character(as.integer(round(x))))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    "--",
    format(round(x, digits), nsmall = digits, trim = TRUE)
  )
}

bold_if_good_share <- function(x, digits = 1) {
  val <- fmt_pct(x, digits)
  ifelse(!is.na(x) & x > 50, paste0("\\textbf{", val, "}"), val)
}

fmt_small_share <- function(x, digits = 1) {
  val <- fmt_pct(x, digits)
  ifelse(!is.na(x) & x > 50, paste0("\\textbf{", val, "}"), val)
}

format_win_table_for_latex <- function(df) {
  df_tex <- df
  
  for (nm in names(df_tex)) {
    if (grepl("wins$", nm)) {
      df_tex[[nm]] <- fmt_int(df_tex[[nm]])
    }
    
    if (grepl("win share$", nm)) {
      df_tex[[nm]] <- bold_if_good_share(df_tex[[nm]])
    }
  }
  
  if ("Total" %in% names(df_tex)) {
    df_tex$Total <- fmt_int(df_tex$Total)
  }
  
  df_tex
}

# ==============================================================================
# LATEX GENERATORS: SUMMARY WIN TABLES
# ==============================================================================

make_latex_period_month_wins <- function(
    df,
    caption = "Win shares of Matrix MF--TPRF across periods and nowcast vintages",
    label   = "tab:summary_wins_period_month",
    note    = "Entries report the share of cases in which Matrix MF--TPRF has a lower RMSFE than the benchmark."
) {
  
  df_tex <- df
  df_tex$Period <- escape_latex(df_tex$Period)
  
  for (nm in names(df_tex)[-1]) {
    df_tex[[nm]] <- bold_if_good_share(df[[nm]])
  }
  
  has_dfm <- all(c("DFM M1", "DFM M2", "DFM M3") %in% names(df_tex))
  
  if (has_dfm) {
    row_cols <- c(
      "Period",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3",
      "DFM M1", "DFM M2", "DFM M3"
    )
    
    tabular <- "l@{\\hspace{0.55cm}}ccc@{\\hspace{0.50cm}}ccc@{\\hspace{0.50cm}}ccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{3}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Period}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  } else {
    row_cols <- c(
      "Period",
      "VEC-P M1", "VEC-P M2", "VEC-P M3",
      "VEC-C M1", "VEC-C M2", "VEC-C M3"
    )
    
    tabular <- "l@{\\hspace{0.65cm}}ccc@{\\hspace{0.60cm}}ccc"
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Period}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  }
  
  rows <- apply(df_tex[, row_cols, drop = FALSE], 1, function(r) {
    paste(r, collapse = " & ")
  })
  
  paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\footnotesize\n",
    "\\renewcommand{\\arraystretch}{1.12}\n",
    "\\setlength{\\tabcolsep}{3.8pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    paste(rows, collapse = " \\\\\n"), " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.04cm}\n",
    "\\parbox{0.86\\linewidth}{\\scriptsize\\centering ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}


make_latex_forecast_wins_summary <- function(
    overall_df,
    period_df,
    vintage_df,
    caption = "Forecast-level wins of Matrix MF--TPRF against the benchmarks",
    label   = "tab:forecast_wins_summary",
    note    = "Wins count the forecasts for which Matrix MF--TPRF has a lower squared forecast error than the benchmark."
) {
  
  overall_tex <- overall_df
  period_tex  <- period_df
  vintage_tex <- vintage_df
  
  names(overall_tex)[1] <- "Group"
  names(period_tex)[1]  <- "Group"
  names(vintage_tex)[1] <- "Group"
  
  has_dfm <- all(c("DFM wins", "DFM win share") %in% names(overall_tex))
  
  format_block <- function(df_block) {
    df_block[[1]] <- escape_latex(df_block[[1]])
    df_block <- format_win_table_for_latex(df_block)
    
    cols <- if (has_dfm) {
      c(
        "Group", "Total",
        "VEC-P wins", "VEC-P win share",
        "VEC-C wins", "VEC-C win share",
        "DFM wins", "DFM win share"
      )
    } else {
      c(
        "Group", "Total",
        "VEC-P wins", "VEC-P win share",
        "VEC-C wins", "VEC-C win share"
      )
    }
    
    apply(df_block[, cols, drop = FALSE], 1, function(r) {
      paste(r, collapse = " & ")
    })
  }
  
  rows_overall <- format_block(overall_tex)
  rows_period  <- format_block(period_tex)
  rows_vintage <- format_block(vintage_tex)
  
  if (has_dfm) {
    tabular <- "l c@{\\hspace{0.45cm}}cc@{\\hspace{0.45cm}}cc@{\\hspace{0.45cm}}cc"
    ncols <- 8
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{2}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Group} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share} \\\\\n"
    )
  } else {
    tabular <- "l c@{\\hspace{0.50cm}}cc@{\\hspace{0.50cm}}cc"
    ncols <- 6
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Group} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share} \\\\\n"
    )
  }
  
  paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\footnotesize\n",
    "\\renewcommand{\\arraystretch}{1.12}\n",
    "\\setlength{\\tabcolsep}{3.2pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    rows_overall[1], " \\\\\n",
    "\\addlinespace[0.25em]\n",
    "\\multicolumn{", ncols, "}{c}{\\scriptsize\\textbf{By evaluation period}} \\\\\n",
    "\\cmidrule(lr){1-", ncols, "}\n",
    paste(rows_period, collapse = " \\\\\n"), " \\\\\n",
    "\\addlinespace[0.25em]\n",
    "\\multicolumn{", ncols, "}{c}{\\scriptsize\\textbf{By nowcast vintage}} \\\\\n",
    "\\cmidrule(lr){1-", ncols, "}\n",
    paste(rows_vintage, collapse = " \\\\\n"), " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.04cm}\n",
    "\\parbox{0.86\\linewidth}{\\scriptsize\\centering ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}


make_latex_forecast_win_table <- function(
    df,
    caption,
    label,
    first_col_name = NULL,
    note = "Entries report forecast-level wins of Matrix MF--TPRF against each benchmark."
) {
  
  df_tex <- df
  
  if (!is.null(first_col_name)) {
    names(df_tex)[1] <- first_col_name
  }
  
  df_tex[[1]] <- escape_latex(df_tex[[1]])
  df_tex <- format_win_table_for_latex(df_tex)
  
  has_dfm <- all(c("DFM wins", "DFM win share") %in% names(df_tex))
  
  if (has_dfm) {
    row_cols <- c(
      names(df_tex)[1], "Total",
      "VEC-P wins", "VEC-P win share",
      "VEC-C wins", "VEC-C win share",
      "DFM wins", "DFM win share"
    )
    
    tabular <- "l c@{\\hspace{0.45cm}}cc@{\\hspace{0.45cm}}cc@{\\hspace{0.45cm}}cc"
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{2}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{", escape_latex(names(df_tex)[1]), "} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share} \\\\\n"
    )
  } else {
    row_cols <- c(
      names(df_tex)[1], "Total",
      "VEC-P wins", "VEC-P win share",
      "VEC-C wins", "VEC-C win share"
    )
    
    tabular <- "l c@{\\hspace{0.50cm}}cc@{\\hspace{0.50cm}}cc"
    
    header <- paste0(
      " & & \\multicolumn{2}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{2}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{", escape_latex(names(df_tex)[1]), "} & \\textbf{Total}",
      " & \\textbf{Wins} & \\textbf{Share}",
      " & \\textbf{Wins} & \\textbf{Share} \\\\\n"
    )
  }
  
  rows <- apply(df_tex[, row_cols, drop = FALSE], 1, function(r) {
    paste(r, collapse = " & ")
  })
  
  paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\footnotesize\n",
    "\\renewcommand{\\arraystretch}{1.10}\n",
    "\\setlength{\\tabcolsep}{3.2pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    paste(rows, collapse = " \\\\\n"), " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.04cm}\n",
    "\\parbox{0.78\\linewidth}{\\scriptsize\\centering ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}


make_latex_country_period_month_wins <- function(
    df,
    caption = "Forecast-level win shares of Matrix MF--TPRF by country, period, and nowcast vintage",
    label   = "tab:forecast_wins_country_period_month",
    note    = "Entries report the share of forecasts for which Matrix MF--TPRF has a lower squared forecast error than the benchmark.",
    period_date_map = c(
      "Pre-COVID"    = "Jan 2017 -- Feb 2020",
      "COVID period" = "Mar 2020 -- Jul 2021",
      "Post-COVID"   = "Aug 2021 -- Feb 2026"
    )
) {
  
  df_tex <- df
  df_tex$Period  <- as.character(df_tex$Period)
  df_tex$Country <- escape_latex(df_tex$Country)
  
  has_dfm <- all(c("DFM M1", "DFM M2", "DFM M3") %in% names(df_tex))
  
  share_cols <- setdiff(names(df_tex), c("Period", "Country"))
  
  for (nm in share_cols) {
    df_tex[[nm]] <- fmt_small_share(df[[nm]])
  }
  
  make_block <- function(period_name) {
    
    block_df <- df_tex %>%
      dplyr::filter(Period == period_name)
    
    if (nrow(block_df) == 0L) return("")
    
    date_line <- period_date_map[[period_name]]
    if (is.null(date_line) || is.na(date_line)) date_line <- ""
    
    if (has_dfm) {
      row_cols <- c(
        "Country",
        "VEC-P M1", "VEC-P M2", "VEC-P M3",
        "VEC-C M1", "VEC-C M2", "VEC-C M3",
        "DFM M1", "DFM M2", "DFM M3"
      )
      ncols <- 10
    } else {
      row_cols <- c(
        "Country",
        "VEC-P M1", "VEC-P M2", "VEC-P M3",
        "VEC-C M1", "VEC-C M2", "VEC-C M3"
      )
      ncols <- 7
    }
    
    rows <- apply(block_df[, row_cols, drop = FALSE], 1, function(r) {
      paste(r, collapse = " & ")
    })
    
    paste0(
      "\\addlinespace[0.18em]\n",
      "\\specialrule{0.045em}{0.06em}{0.035em}\n",
      "\\multicolumn{", ncols, "}{@{}l}{\\fontsize{6.4}{6.8}\\selectfont\\textbf{",
      escape_latex(period_name),
      "}} \\\\\n",
      "\\multicolumn{", ncols, "}{@{}l}{\\fontsize{5.4}{5.8}\\selectfont ",
      escape_latex(date_line),
      "} \\\\\n",
      "\\cmidrule(lr){1-", ncols, "}\n",
      paste(rows, collapse = " \\\\\n"),
      " \\\\\n"
    )
  }
  
  body <- paste(
    make_block("Pre-COVID"),
    make_block("COVID period"),
    make_block("Post-COVID"),
    sep = "\n"
  )
  
  if (has_dfm) {
    tabular <- paste0(
      ">{\\centering\\arraybackslash}p{0.58cm}",
      "@{\\hspace{0.28cm}}ccc",
      "@{\\hspace{0.28cm}}ccc",
      "@{\\hspace{0.28cm}}ccc"
    )
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}}",
      " & \\multicolumn{3}{c}{\\textbf{DFM}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Country}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  } else {
    tabular <- paste0(
      ">{\\centering\\arraybackslash}p{0.58cm}",
      "@{\\hspace{0.35cm}}ccc",
      "@{\\hspace{0.35cm}}ccc"
    )
    
    header <- paste0(
      " & \\multicolumn{3}{c}{\\textbf{VEC-P}}",
      " & \\multicolumn{3}{c}{\\textbf{VEC-C}} \\\\\n",
      "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\n",
      "\\rowcolor{topgray}\n",
      "\\textbf{Country}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
      " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n"
    )
  }
  
  paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\fontsize{6.6}{7.2}\\selectfont\n",
    "\\renewcommand{\\arraystretch}{1.02}\n",
    "\\setlength{\\tabcolsep}{2.2pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    "\\begin{tabular}{", tabular, "}\n",
    "\\toprule\n",
    header,
    "\\midrule\n",
    body,
    "\\bottomrule\n",
    "\\end{tabular}\n\n",
    "\\vspace{0.04cm}\n",
    "\\parbox{0.86\\linewidth}{\\scriptsize\\centering ",
    escape_latex(note),
    "}\n",
    "\\end{table}"
  )
}

# ==============================================================================
# MATRIX FACTOR DIAGNOSTICS: STATIC FACTOR + ROW/COLUMN LOADINGS
# ==============================================================================

make_simple_latex_table <- function(df, caption, label, digits = 3) {
  
  df_out <- df
  
  num_cols <- sapply(df_out, is.numeric)
  df_out[num_cols] <- lapply(df_out[num_cols], function(x) round(x, digits))
  
  paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{l", paste(rep("c", ncol(df_out) - 1), collapse = ""), "}\n",
    "\\toprule\n",
    paste(colnames(df_out), collapse = " & "), " \\\\\n",
    "\\midrule\n",
    paste(
      apply(df_out, 1, function(x) paste(x, collapse = " & ")),
      collapse = " \\\\\n"
    ),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}


extract_matrix_loadings <- function(summary_matrix) {
  
  fc <- summary_matrix$factor_comparison
  
  possible_R_names <- c(
    "R", "R_hat", "row_loadings", "loadings_row",
    "matrix_R", "R_loadings"
  )
  
  possible_C_names <- c(
    "C", "C_hat", "col_loadings", "column_loadings",
    "loadings_col", "matrix_C", "C_loadings"
  )
  
  R_loadings <- NULL
  C_loadings <- NULL
  
  for (nm in possible_R_names) {
    if (!is.null(fc[[nm]])) R_loadings <- fc[[nm]]
  }
  
  for (nm in possible_C_names) {
    if (!is.null(fc[[nm]])) C_loadings <- fc[[nm]]
  }
  
  if (is.null(R_loadings) && !is.null(summary_matrix$R)) {
    R_loadings <- summary_matrix$R
  }
  
  if (is.null(C_loadings) && !is.null(summary_matrix$C)) {
    C_loadings <- summary_matrix$C
  }
  
  if (is.null(R_loadings) || is.null(C_loadings)) {
    cat("\nAvailable names in summary_matrix$factor_comparison:\n")
    print(names(fc))
    stop("Could not find row/column loadings.")
  }
  
  list(
    R = as.matrix(R_loadings),
    C = as.matrix(C_loadings)
  )
}


make_loading_table <- function(loadings, first_col_name) {
  
  row_names <- rownames(loadings)
  
  if (is.null(row_names)) {
    row_names <- paste0(first_col_name, " ", seq_len(nrow(loadings)))
  }
  
  df <- data.frame(
    name = row_names,
    loadings,
    check.names = FALSE
  )
  
  colnames(df)[1] <- first_col_name
  colnames(df)[-1] <- paste0("Factor ", seq_len(ncol(loadings)))
  
  df
}


plot_static_matrix_factor <- function(summary_matrix, params, Size, sel) {
  
  fc <- summary_matrix$factor_comparison
  
  df_static <- fc$df_q
  
  if (is.null(df_static)) {
    df_static <- data.frame(
      date   = as.Date(fc$dates_q),
      GDP    = as.numeric(fc$gdp_q),
      Matrix = as.numeric(fc$matrix_factor_quarterly)
    )
  }
  
  df_long <- df_static %>%
    dplyr::select(date, GDP, Matrix) %>%
    tidyr::pivot_longer(
      cols      = c(GDP, Matrix),
      names_to  = "series",
      values_to = "value"
    ) %>%
    dplyr::group_by(series) %>%
    dplyr::mutate(value_std = as.numeric(scale(value))) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = date, y = value_std, colour = series)
  ) +
    ggplot2::annotate(
      "rect",
      xmin = params$covid_start,
      xmax = params$covid_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.15
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      colour = "grey70"
    ) +
    ggplot2::geom_line(linewidth = 1.05) +
    ggplot2::scale_x_date(
      breaks = seq(
        lubridate::floor_date(min(df_long$date), unit = "year"),
        lubridate::floor_date(max(df_long$date), unit = "year"),
        by = "1 year"
      ),
      date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title    = "Static Matrix MF-TPRF factor",
      subtitle = paste0("Size = ", Size, ", sel = ", sel),
      x        = "Date",
      y        = "Standardized value",
      colour   = NULL
    ) +
    theme_factor_compare()
}


plot_column_loadings <- function(df_col_loadings, Size, sel) {
  
  df_long <- df_col_loadings %>%
    tidyr::pivot_longer(
      cols      = -Column,
      names_to  = "factor",
      values_to = "loading"
    )
  
  ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = Column, y = loading, fill = factor)
  ) +
    ggplot2::geom_col(position = "dodge", width = 0.75) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.35,
      colour = "grey45"
    ) +
    ggplot2::labs(
      title    = "Matrix MF-TPRF column loadings",
      subtitle = paste0("Size = ", Size, ", sel = ", sel),
      x        = NULL,
      y        = "Loading",
      fill     = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey25")
    )
}


run_matrix_factor_diagnostics <- function(
    summary_matrix,
    params,
    Size,
    sel,
    path_final_results,
    suffix_out
) {
  
  if (!identical(Size, "small") || !(sel %in% c("corr", "LASSO"))) {
    cat("\nSkipping matrix factor diagnostics: available only for small + corr/LASSO.\n")
    return(NULL)
  }
  
  if (is.null(summary_matrix$factor_comparison)) {
    stop("summary_matrix does not contain factor_comparison.")
  }
  
  path_matrix_factor_diag <- file.path(
    path_final_results,
    "matrix_factor_diagnostics"
  )
  
  dir.create(path_matrix_factor_diag, recursive = TRUE, showWarnings = FALSE)
  # ---------------------------------------------------------------------------
  # 1. Static matrix factor plot
  # ---------------------------------------------------------------------------
  
  plot_static_factor_matrix <- plot_static_matrix_factor(
    summary_matrix = summary_matrix,
    params         = params,
    Size           = Size,
    sel            = sel
  )
  
  print(plot_static_factor_matrix)
  
  file_static_factor_matrix <- file.path(
    path_matrix_factor_diag,
    paste0("plot_static_factor_matrix_", suffix_out, ".png")
  )
  
  ggplot2::ggsave(
    filename = file_static_factor_matrix,
    plot     = plot_static_factor_matrix,
    width    = 11,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  
  # ---------------------------------------------------------------------------
  # 2. Row and column loadings
  # ---------------------------------------------------------------------------
  
  loadings <- extract_matrix_loadings(summary_matrix)
  
  df_row_loadings <- make_loading_table(
    loadings       = loadings$R,
    first_col_name = "Row"
  )
  
  df_col_loadings <- make_loading_table(
    loadings       = loadings$C,
    first_col_name = "Column"
  )
  
  cat("\n================ MATRIX ROW LOADINGS TABLE ================\n")
  print(df_row_loadings)
  
  cat("\n================ MATRIX COLUMN LOADINGS TABLE ================\n")
  print(df_col_loadings)
  
  latex_row_loadings <- make_simple_latex_table(
    df      = df_row_loadings,
    caption = paste0("Matrix MF--TPRF row loadings (", Size, ", sel = ", sel, ")"),
    label   = paste0("tab:matrix_row_loadings_", Size, "_", sel),
    digits  = 3
  )
  
  latex_col_loadings <- make_simple_latex_table(
    df      = df_col_loadings,
    caption = paste0("Matrix MF--TPRF column loadings (", Size, ", sel = ", sel, ")"),
    label   = paste0("tab:matrix_column_loadings_", Size, "_", sel),
    digits  = 3
  )
  
  print_latex_block("MATRIX ROW LOADINGS LATEX", latex_row_loadings)
  print_latex_block("MATRIX COLUMN LOADINGS LATEX", latex_col_loadings)
  
  writeLines(
    latex_row_loadings,
    file.path(path_matrix_factor_diag, paste0("matrix_row_loadings_", suffix_out, ".tex"))
  )
  
  writeLines(
    latex_col_loadings,
    file.path(path_matrix_factor_diag, paste0("matrix_column_loadings_", suffix_out, ".tex"))
  )
  
  # ---------------------------------------------------------------------------
  # 3. Column loadings plot
  # ---------------------------------------------------------------------------
  
  plot_col_loadings <- plot_column_loadings(
    df_col_loadings = df_col_loadings,
    Size            = Size,
    sel             = sel
  )
  
  print(plot_col_loadings)
  
  file_col_loadings <- file.path(
    path_matrix_factor_diag,
    paste0("plot_matrix_column_loadings_", suffix_out, ".png")
  )
  
  ggplot2::ggsave(
    filename = file_col_loadings,
    plot     = plot_col_loadings,
    width    = 10.5,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  
  list(
    df_row_loadings           = df_row_loadings,
    df_col_loadings           = df_col_loadings,
    latex_row_loadings        = latex_row_loadings,
    latex_col_loadings        = latex_col_loadings,
    plot_static_factor_matrix = plot_static_factor_matrix,
    plot_col_loadings         = plot_col_loadings,
    file_static_factor_matrix = file_static_factor_matrix,
    file_col_loadings         = file_col_loadings
  )
}

# ==============================================================================
# GDP TARGET-PROXY AVAILABILITY TABLE FOR LATEX
# ==============================================================================
build_missing_matrix_table <- function(
    selection_report,
    var_scope = "union"
) {
  
  size_order     <- c("small", "medium", "large")
  vintage_levels <- c("M1", "M2", "M3")
  
  required_helpers <- c(
    "build_tensor",
    "make_regime_data",
    "unbalancedness_tensor",
    "compute_m_tr"
  )
  
  missing_helpers <- required_helpers[
    !vapply(
      required_helpers,
      exists,
      logical(1),
      mode = "function"
    )
  ]
  
  if (length(missing_helpers) > 0L) {
    stop(
      "Missing required functions: ",
      paste(missing_helpers, collapse = ", "),
      ". Make sure matrix.mf.tprf.now.R is sourced."
    )
  }
  
  params_ref <- selection_report$fit_results$small$params
  
  required_params <- c(
    "start_eval",
    "end_eval"
  )
  
  missing_params <- setdiff(required_params, names(params_ref))
  
  if (length(missing_params) > 0L) {
    stop(
      "Missing parameters: ",
      paste(missing_params, collapse = ", ")
    )
  }
  
  if (is.null(selection_report$fit_results$small$dates_q)) {
    stop(
      "Missing dates_q in selection_report$fit_results$small."
    )
  }
  
  if (is.null(selection_report$dates$updated_regime_start)) {
    stop(
      "Missing updated_regime_start in selection_report$dates."
    )
  }
  
  dates_q <- as.Date(
    selection_report$fit_results$small$dates_q
  )
  
  evaluation_start <- as.Date(params_ref$start_eval)
  evaluation_end   <- as.Date(params_ref$end_eval)
  
  post_start <- as.Date(
    selection_report$dates$updated_regime_start
  )
  
  build_regime_data <- function(
    regime_obj,
    regime_label,
    size_label
  ) {
    
    params_i <- selection_report$fit_results[[size_label]]$params
    
    tensor_i <- build_tensor(
      prep      = regime_obj[[size_label]]$all_countries,
      params    = params_i,
      var_scope = var_scope
    )
    
    selection_end_i <- as.Date(
      regime_obj[[size_label]]$dimensions$selection_end[1L]
    )
    
    regime_data_i <- make_regime_data(
      tensor_obj    = tensor_i,
      label         = regime_label,
      selection_end = selection_end_i,
      params        = params_i
    )
    
    list(
      data    = regime_data_i,
      dates_m = as.Date(tensor_i$dates)
    )
  }
  
  get_block_stats <- function(
    X_cut,
    variable_index
  ) {
    
    X_block <- X_cut[, , variable_index, drop = FALSE]
    
    missing_cells <- sum(is.na(X_block))
    total_cells   <- length(X_block)
    
    c(
      missing_cells = missing_cells,
      total_cells   = total_cells,
      missing_share = 100 * missing_cells / total_cells
    )
  }
  
  build_regime_vintages <- function(
    regime_obj,
    regime_label,
    size_label,
    window_start,
    window_end
  ) {
    
    obj_i <- build_regime_data(
      regime_obj   = regime_obj,
      regime_label = regime_label,
      size_label   = size_label
    )
    
    regime_data_i <- obj_i$data
    dates_m_i     <- obj_i$dates_m
    
    if (length(dates_m_i) != dim(regime_data_i$X_full)[1L]) {
      stop(
        "Incompatible dates_m and X_full dimensions for ",
        regime_label,
        ", ",
        size_label,
        "."
      )
    }
    
    vintage_index <- which(
      dates_m_i >= window_start &
        dates_m_i <= window_end
    )
    
    if (length(vintage_index) == 0L) {
      stop(
        "No evaluation vintages found for ",
        regime_label,
        ", ",
        size_label,
        "."
      )
    }
    
    month_position <- vapply(
      vintage_index,
      function(tt) {
        compute_m_tr(
          date_t  = dates_m_i[tt],
          dates_q = dates_q
        )
      },
      integer(1)
    )
    
    valid_vintage <- !is.na(month_position)
    
    vintage_index  <- vintage_index[valid_vintage]
    month_position <- month_position[valid_vintage]
    
    if (length(vintage_index) == 0L) {
      stop(
        "No valid M1/M2/M3 vintages found for ",
        regime_label,
        ", ",
        size_label,
        "."
      )
    }
    
    idx_M <- seq_len(regime_data_i$N_m)
    
    idx_Q <- seq.int(
      from = regime_data_i$N_m + 1L,
      to   = regime_data_i$N_m + regime_data_i$N_q
    )
    
    vintage_stats <- lapply(
      seq_along(vintage_index),
      function(i) {
        
        tt_i <- vintage_index[i]
        
        X_cut_i <- unbalancedness_tensor(
          X_full    = regime_data_i$X_full,
          Unb       = regime_data_i$Unb,
          current_t = tt_i
        )
        
        stats_M <- get_block_stats(
          X_cut          = X_cut_i,
          variable_index = idx_M
        )
        
        stats_Q <- get_block_stats(
          X_cut          = X_cut_i,
          variable_index = idx_Q
        )
        
        tibble::tibble(
          regime       = regime_label,
          size         = size_label,
          vintage      = paste0("M", month_position[i]),
          window_start = as.Date(window_start),
          window_end   = as.Date(window_end),
          
          frequency = c("Monthly", "Quarterly"),
          
          missing_cells = c(
            stats_M["missing_cells"],
            stats_Q["missing_cells"]
          ),
          
          total_cells = c(
            stats_M["total_cells"],
            stats_Q["total_cells"]
          ),
          
          missing_share = c(
            stats_M["missing_share"],
            stats_Q["missing_share"]
          )
        )
      }
    )
    
    dplyr::bind_rows(vintage_stats)
  }
  
  initial_end <- post_start %m-% lubridate::period(1, "month")
  
  stats_long <- dplyr::bind_rows(
    lapply(
      size_order,
      function(size_i) {
        build_regime_vintages(
          regime_obj   = selection_report$initial,
          regime_label = "Initial regime",
          size_label   = size_i,
          window_start = evaluation_start,
          window_end   = initial_end
        )
      }
    ),
    
    lapply(
      size_order,
      function(size_i) {
        build_regime_vintages(
          regime_obj   = selection_report$updated,
          regime_label = "Post-COVID update",
          size_label   = size_i,
          window_start = post_start,
          window_end   = evaluation_end
        )
      }
    )
  )
  
  frequency_tab <- stats_long %>%
    dplyr::group_by(
      regime,
      size,
      frequency,
      window_start,
      window_end,
      vintage
    ) %>%
    dplyr::summarise(
      mean_missing_share = mean(missing_share),
      n_vintages         = dplyr::n(),
      .groups            = "drop"
    ) %>%
    dplyr::mutate(
      frequency = tolower(as.character(frequency)),
      vintage   = factor(as.character(vintage), levels = vintage_levels)
    ) %>%
    tidyr::pivot_wider(
      id_cols = c(
        regime,
        size,
        window_start,
        window_end
      ),
      names_from  = c(frequency, vintage),
      values_from = c(
        mean_missing_share,
        n_vintages
      ),
      names_glue = "{frequency}_{vintage}_{.value}"
    )
  
  overall_tab <- stats_long %>%
    dplyr::group_by(
      regime,
      size,
      window_start,
      window_end,
      vintage
    ) %>%
    dplyr::summarise(
      overall_share_vintage = 100 *
        sum(missing_cells) /
        sum(total_cells),
      .groups = "drop"
    ) %>%
    dplyr::group_by(
      regime,
      size,
      window_start,
      window_end
    ) %>%
    dplyr::summarise(
      overall_share = mean(overall_share_vintage),
      n_realtime_vintages = dplyr::n(),
      .groups = "drop"
    )
  
  tab <- frequency_tab %>%
    dplyr::left_join(
      overall_tab,
      by = c(
        "regime",
        "size",
        "window_start",
        "window_end"
      )
    )
  
  required_share_cols <- c(
    "monthly_M1_mean_missing_share",
    "monthly_M2_mean_missing_share",
    "monthly_M3_mean_missing_share",
    "quarterly_M1_mean_missing_share",
    "quarterly_M2_mean_missing_share",
    "quarterly_M3_mean_missing_share"
  )
  
  for (column_i in required_share_cols) {
    if (!column_i %in% names(tab)) {
      tab[[column_i]] <- NA_real_
    }
  }
  
  tab %>%
    dplyr::transmute(
      regime,
      size,
      window_start,
      window_end,
      
      monthly_M1 = monthly_M1_mean_missing_share,
      monthly_M2 = monthly_M2_mean_missing_share,
      monthly_M3 = monthly_M3_mean_missing_share,
      
      quarterly_M1 = quarterly_M1_mean_missing_share,
      quarterly_M2 = quarterly_M2_mean_missing_share,
      quarterly_M3 = quarterly_M3_mean_missing_share,
      
      overall_share,
      n_realtime_vintages = as.integer(n_realtime_vintages)
    ) %>%
    dplyr::mutate(
      regime = factor(
        regime,
        levels = c(
          "Initial regime",
          "Post-COVID update"
        )
      ),
      
      size = factor(
        size,
        levels = size_order
      )
    ) %>%
    dplyr::arrange(regime, size)
}

build_latex_missing_matrix_table <- function(
    tbl,
    selection_report,
    caption,
    label,
    note = paste0(
      "\\textit{Notes:} Monthly and quarterly entries report average unavailable ",
      "shares (\\%) in the pseudo-real-time predictor tensors before imputation, ",
      "separately for M1, M2, and M3 vintages. The final column reports the ",
      "mean unavailable share across all real-time M1--M3 datasets in the ",
      "corresponding update regime and information set. Within each vintage, ",
      "this overall share is computed over all predictor cells. Quarterly ",
      "predictors are represented on the monthly grid and include structural ",
      "within-quarter unavailability."
    )
) {
  
  required_cols <- c(
    "regime",
    "size",
    "monthly_M1",
    "monthly_M2",
    "monthly_M3",
    "quarterly_M1",
    "quarterly_M2",
    "quarterly_M3",
    "overall_share"
  )
  
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in real-time availability table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  pre_start <- as.Date(
    selection_report$dates$initial_regime_start
  )
  
  post_start <- as.Date(
    selection_report$dates$updated_regime_start
  )
  
  evaluation_end <- as.Date(
    selection_report$fit_results$small$params$end_eval
  )
  
  pre_end <- post_start %m-% lubridate::period(1, "month")
  
  latex_month_year <- function(x) {
    
    month_names <- c(
      "Jan.", "Feb.", "Mar.", "Apr.",
      "May", "Jun.", "Jul.", "Aug.",
      "Sep.", "Oct.", "Nov.", "Dec."
    )
    
    x <- as.Date(x)
    
    paste0(
      month_names[as.integer(format(x, "%m"))],
      "\\ ",
      format(x, "%Y")
    )
  }
  
  regime_labels <- c(
    pre = paste0(
      "\\textbf{Pre-update:} ",
      latex_month_year(pre_start),
      "--",
      latex_month_year(pre_end)
    ),
    
    post = paste0(
      "\\textbf{Post-update:} ",
      latex_month_year(post_start),
      "--",
      latex_month_year(evaluation_end)
    )
  )
  
  tab <- tbl %>%
    dplyr::mutate(
      regime_id = dplyr::case_when(
        as.character(regime) == "Initial regime" ~ "pre",
        as.character(regime) == "Post-COVID update" ~ "post",
        TRUE ~ NA_character_
      ),
      
      size = factor(
        tolower(as.character(size)),
        levels = c("small", "medium", "large")
      )
    )
  
  if (anyNA(tab$regime_id)) {
    stop("Unknown regime label in tbl$regime.")
  }
  
  if (anyDuplicated(tab[c("regime_id", "size")]) > 0L) {
    stop(
      "More than one row found for at least one regime ",
      "and information-set combination."
    )
  }
  
  fmt_share <- function(x) {
    
    if (is.na(x)) {
      return("--")
    }
    
    sprintf("%.1f", as.numeric(x))
  }
  
  make_regime_block <- function(
    regime_key,
    regime_label
  ) {
    
    block <- tab %>%
      dplyr::filter(.data$regime_id == regime_key) %>%
      dplyr::arrange(size)
    
    if (nrow(block) != 3L) {
      stop(
        "Expected Small, Medium, and Large rows for ",
        regime_key,
        "."
      )
    }
    
    rows <- vapply(
      seq_len(nrow(block)),
      function(i) {
        
        paste0(
          tools::toTitleCase(as.character(block$size[i])), " & ",
          fmt_share(block$monthly_M1[i]), " & ",
          fmt_share(block$monthly_M2[i]), " & ",
          fmt_share(block$monthly_M3[i]), " & ",
          fmt_share(block$quarterly_M1[i]), " & ",
          fmt_share(block$quarterly_M2[i]), " & ",
          fmt_share(block$quarterly_M3[i]), " & ",
          fmt_share(block$overall_share[i]),
          " \\\\"
        )
      },
      character(1)
    )
    
    c(
      paste0(
        "\\rowcolor{regimegray}\n",
        "\\multicolumn{8}{@{}c@{}}{",
        regime_label,
        "} \\\\[-0.18em]"
      ),
      "\\cmidrule(lr){1-8}",
      rows
    )
  }
  
  body <- c(
    make_regime_block(
      regime_key   = "pre",
      regime_label = regime_labels["pre"]
    ),
    
    "\\addlinespace[0.28em]",
    
    make_regime_block(
      regime_key   = "post",
      regime_label = regime_labels["post"]
    )
  )
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{1.06}\n",
    "\\setlength{\\tabcolsep}{2.8pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n",
    "\\definecolor{regimegray}{gray}{0.94}\n\n",
    
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    
    "\\begin{tabularx}{\\textwidth}{@{}",
    ">{\\centering\\arraybackslash}p{2.20cm}",
    "@{\\hspace{0.25cm}}",
    "*{3}{>{\\centering\\arraybackslash}X}",
    "@{\\hspace{0.55cm}}",
    "*{3}{>{\\centering\\arraybackslash}X}",
    "@{\\hspace{0.55cm}}",
    ">{\\centering\\arraybackslash}X",
    "@{}}\n",
    
    "\\toprule\n",
    
    " & \\multicolumn{3}{c}{\\textbf{Monthly predictors}} & ",
    "\\multicolumn{3}{c}{\\textbf{Quarterly predictors}} & ",
    "\\textbf{All predictors} \\\\\n",
    
    "\\cmidrule(lr){2-4}",
    "\\cmidrule(lr){5-7}",
    "\\cmidrule(lr){8-8}\n",
    
    "\\cellcolor{topgray}\\textbf{Information set} & ",
    "\\cellcolor{topgray}\\textbf{M1} & ",
    "\\cellcolor{topgray}\\textbf{M2} & ",
    "\\cellcolor{topgray}\\textbf{M3} & ",
    "\\cellcolor{topgray}\\textbf{M1} & ",
    "\\cellcolor{topgray}\\textbf{M2} & ",
    "\\cellcolor{topgray}\\textbf{M3} & ",
    "\\cellcolor{topgray}\\textbf{Average} \\\\\n",
    
    "\\midrule\n",
    paste(body, collapse = "\n"),
    
    "\n\\bottomrule\n",
    "\\end{tabularx}\n",
    
    "\\vspace{-0.08cm}\n",
    "\\noindent\\parbox{\\textwidth}{\\centering\\scriptsize\n",
    note,
    "}\n",
    
    "\\end{table}\n"
  )
}
# ==============================================================================
# CALIBRATION-VINTAGE PREDICTOR AVAILABILITY TABLE
# ==============================================================================

build_calibration_missing_matrix_table <- function(
    selection_report,
    var_scope = "union"
) {
  
  size_order <- c("small", "medium", "large")
  
  required_helpers <- c(
    "build_tensor",
    "make_regime_data",
    "unbalancedness_tensor"
  )
  
  missing_helpers <- required_helpers[
    !vapply(
      required_helpers,
      exists,
      logical(1),
      mode = "function"
    )
  ]
  
  if (length(missing_helpers) > 0L) {
    stop(
      "Missing required functions: ",
      paste(missing_helpers, collapse = ", "),
      "."
    )
  }
  
  get_missing_stats <- function(X_block) {
    
    n_missing <- sum(is.na(X_block))
    n_total   <- length(X_block)
    
    c(
      missing = n_missing,
      total   = n_total,
      share   = 100 * n_missing / n_total
    )
  }
  
  build_one <- function(
    regime_obj,
    regime_label,
    size_label
  ) {
    
    params_i <- selection_report$fit_results[[size_label]]$params
    
    tensor_i <- build_tensor(
      prep      = regime_obj[[size_label]]$all_countries,
      params    = params_i,
      var_scope = var_scope
    )
    
    selection_end_i <- as.Date(
      regime_obj[[size_label]]$dimensions$selection_end[1L]
    )
    
    regime_data_i <- make_regime_data(
      tensor_obj    = tensor_i,
      label         = regime_label,
      selection_end = selection_end_i,
      params        = params_i
    )
    
    dates_m_i <- as.Date(
      dimnames(regime_data_i$X_full)[[1]]
    )
    
    if (anyNA(dates_m_i)) {
      stop(
        "Could not recover monthly dates for ",
        regime_label,
        ", ",
        size_label,
        "."
      )
    }
    
    calibration_t <- which(
      format(dates_m_i, "%Y-%m") ==
        format(selection_end_i, "%Y-%m")
    )
    
    if (length(calibration_t) != 1L) {
      stop(
        "Could not identify a unique calibration vintage for ",
        regime_label,
        ", ",
        size_label,
        ". Expected month: ",
        format(selection_end_i, "%Y-%m"),
        "."
      )
    }
    
    X_cut_i <- unbalancedness_tensor(
      X_full    = regime_data_i$X_full,
      Unb       = regime_data_i$Unb,
      current_t = calibration_t
    )
    
    idx_M <- seq_len(regime_data_i$N_m)
    
    idx_Q <- seq.int(
      from = regime_data_i$N_m + 1L,
      to   = regime_data_i$N_m + regime_data_i$N_q
    )
    
    stats_M <- get_missing_stats(
      X_cut_i[, , idx_M, drop = FALSE]
    )
    
    stats_Q <- get_missing_stats(
      X_cut_i[, , idx_Q, drop = FALSE]
    )
    
    stats_all <- get_missing_stats(X_cut_i)
    
    tibble::tibble(
      regime = regime_label,
      size   = size_label,
      
      calibration_date = dates_m_i[calibration_t],
      history_start    = dates_m_i[1L],
      history_months   = dim(X_cut_i)[1L],
      
      n_countries = dim(X_cut_i)[2L],
      N_m         = regime_data_i$N_m,
      N_q         = regime_data_i$N_q,
      
      monthly_missing = stats_M["missing"],
      monthly_total   = stats_M["total"],
      monthly_share   = stats_M["share"],
      
      quarterly_missing = stats_Q["missing"],
      quarterly_total   = stats_Q["total"],
      quarterly_share   = stats_Q["share"],
      
      total_missing = stats_all["missing"],
      total_entries = stats_all["total"],
      total_share   = stats_all["share"]
    )
  }
  
  dplyr::bind_rows(
    lapply(
      size_order,
      function(size_i) {
        build_one(
          regime_obj   = selection_report$initial,
          regime_label = "Initial calibration",
          size_label   = size_i
        )
      }
    ),
    
    lapply(
      size_order,
      function(size_i) {
        build_one(
          regime_obj   = selection_report$updated,
          regime_label = "Post-COVID recalibration",
          size_label   = size_i
        )
      }
    )
  ) %>%
    dplyr::mutate(
      regime = factor(
        regime,
        levels = c(
          "Initial calibration",
          "Post-COVID recalibration"
        )
      ),
      
      size = factor(
        size,
        levels = size_order
      )
    ) %>%
    dplyr::arrange(regime, size)
}

build_latex_calibration_missing_matrix_table <- function(
    tbl,
    caption,
    label,
    note = paste0(
      "\\textit{Notes:} Entries report unavailable predictor shares (\\%) ",
      "at the two real-time calibration vintages used to freeze model ",
      "complexity within each update regime. The same publication-delay mask ",
      "used for hyperparameter selection is applied before computing ",
      "availability. GDP is excluded. Quarterly predictors are represented on ",
      "the monthly grid and therefore include structural within-quarter ",
      "unavailability."
    )
) {
  
  required_cols <- c(
    "regime",
    "size",
    "calibration_date",
    "monthly_share",
    "quarterly_share",
    "total_share"
  )
  
  missing_cols <- setdiff(required_cols, names(tbl))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in calibration-availability table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  latex_month_year <- function(x) {
    
    month_names <- c(
      "Jan.", "Feb.", "Mar.", "Apr.",
      "May", "Jun.", "Jul.", "Aug.",
      "Sep.", "Oct.", "Nov.", "Dec."
    )
    
    x <- as.Date(x)
    
    paste0(
      month_names[as.integer(format(x, "%m"))],
      "\\ ",
      format(x, "%Y")
    )
  }
  
  fmt_share <- function(x) {
    
    if (is.na(x)) {
      return("--")
    }
    
    sprintf("%.1f", as.numeric(x))
  }
  
  tab <- tbl %>%
    dplyr::mutate(
      regime_id = dplyr::case_when(
        as.character(regime) == "Initial calibration" ~ "pre",
        as.character(regime) == "Post-COVID recalibration" ~ "post",
        TRUE ~ NA_character_
      ),
      
      size = factor(
        tolower(as.character(size)),
        levels = c("small", "medium", "large")
      ),
      
      calibration_date = as.Date(calibration_date)
    )
  
  if (anyNA(tab$regime_id)) {
    stop("Unknown regime label in tbl$regime.")
  }
  
  make_regime_block <- function(
    regime_key,
    regime_title
  ) {
    
    block <- tab %>%
      dplyr::filter(.data$regime_id == regime_key) %>%
      dplyr::arrange(size)
    
    if (nrow(block) != 3L) {
      stop(
        "Expected Small, Medium, and Large rows for ",
        regime_key,
        "."
      )
    }
    
    calibration_dates <- unique(block$calibration_date)
    
    if (length(calibration_dates) != 1L) {
      stop(
        "More than one calibration date found for ",
        regime_key,
        "."
      )
    }
    
    rows <- vapply(
      seq_len(nrow(block)),
      function(i) {
        
        paste0(
          tools::toTitleCase(as.character(block$size[i])), " & ",
          fmt_share(block$monthly_share[i]), " & ",
          fmt_share(block$quarterly_share[i]), " & ",
          fmt_share(block$total_share[i]),
          " \\\\"
        )
      },
      character(1)
    )
    
    c(
      paste0(
        "\\rowcolor{regimegray}\n",
        "\\multicolumn{4}{@{}c@{}}{\\textbf{",
        regime_title,
        " calibration:} ",
        latex_month_year(calibration_dates),
        "} \\\\[-0.18em]"
      ),
      "\\cmidrule(lr){1-4}",
      rows
    )
  }
  
  body <- c(
    make_regime_block(
      regime_key   = "pre",
      regime_title = "Pre-update"
    ),
    
    "\\addlinespace[0.28em]",
    
    make_regime_block(
      regime_key   = "post",
      regime_title = "Post-update"
    )
  )
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{1.08}\n",
    "\\setlength{\\tabcolsep}{3.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n",
    "\\definecolor{regimegray}{gray}{0.94}\n\n",
    
    "\\caption{", escape_latex(caption), "}\n",
    "\\label{", label, "}\n\n",
    
    "\\begin{tabularx}{\\textwidth}{@{}",
    ">{\\centering\\arraybackslash}p{2.65cm}",
    "@{\\hspace{0.55cm}}",
    "*{3}{>{\\centering\\arraybackslash}X}",
    "@{}}\n",
    
    "\\toprule\n",
    
    "\\cellcolor{topgray}\\textbf{Information set} & ",
    "\\cellcolor{topgray}\\textbf{Monthly (\\%)} & ",
    "\\cellcolor{topgray}\\textbf{Quarterly (\\%)} & ",
    "\\cellcolor{topgray}\\textbf{All predictors (\\%)} \\\\\n",
    
    "\\midrule\n",
    paste(body, collapse = "\n"),
    
    "\n\\bottomrule\n",
    "\\end{tabularx}\n",
    
    "\\vspace{-0.08cm}\n",
    "\\noindent\\parbox{\\textwidth}{\\centering\\scriptsize\n",
    note,
    "}\n",
    
    "\\end{table}\n"
  )
}

# ==============================================================================
# utils_results_05_unit_target_proxy.R
# ==============================================================================
# Unit-target proxy robustness tables:
#   1. Unit-target RMSFEs.
#   2. Unit-target / aggregate-target relative RMSFEs with DM tests.
# ==============================================================================

unit_proxy_sizes <- c("small", "medium", "large")

unit_proxy_period_order <- c(
  "Pre-COVID",
  "COVID period",
  "Post-COVID"
)

unit_proxy_vintage_order <- c("M1", "M2", "M3")


# ==============================================================================
# 1. FILE LOADING
# ==============================================================================

find_proxy_rt_file <- function(path_results,
                               model,
                               size,
                               sel) {
  
  if (!dir.exists(path_results)) {
    stop("`path_results` does not exist: ", path_results)
  }
  
  pattern <- paste0(
    "^rt_", model,
    "_Size-", size,
    "_sel-", sel,
    ".*\\.rds$"
  )
  
  files <- list.files(
    path        = path_results,
    pattern     = pattern,
    full.names  = TRUE,
    ignore.case = TRUE
  )
  
  if (length(files) == 0L) {
    stop(
      "No pseudo-real-time result found for:\n",
      "  model = ", model, "\n",
      "  size  = ", size, "\n",
      "  sel   = ", sel, "\n",
      "  path  = ", path_results
    )
  }
  
  if (length(files) > 1L) {
    
    files <- files[
      order(
        file.info(files)$mtime,
        decreasing = TRUE
      )
    ]
    
    message(
      "Multiple files found for model = ", model,
      ", Size = ", size,
      ", sel = ", sel,
      ". Using:\n",
      files[1L]
    )
  }
  
  files[1L]
}


validate_proxy_rt_result <- function(result,
                                     expected_proxy_mode,
                                     expected_size,
                                     expected_sel) {
  
  required_fields <- c(
    "proxy_mode",
    "stage",
    "Size",
    "sel",
    "params",
    "countries",
    "dates_q",
    "Y_q_all",
    "pseudo_rt_all"
  )
  
  missing_fields <- setdiff(required_fields, names(result))
  
  if (length(missing_fields) > 0L) {
    stop(
      "Saved result is missing: ",
      paste(missing_fields, collapse = ", ")
    )
  }
  
  if (!identical(
    tolower(as.character(result$proxy_mode)[1L]),
    tolower(expected_proxy_mode)
  )) {
    stop(
      "Proxy-mode mismatch. Expected `", expected_proxy_mode,
      "`, found `", result$proxy_mode, "`."
    )
  }
  
  if (!identical(
    tolower(as.character(result$stage)[1L]),
    "pseudo_realtime"
  )) {
    stop(
      "The saved object is not a pseudo-real-time result."
    )
  }
  
  if (!identical(
    tolower(as.character(result$Size)[1L]),
    tolower(expected_size)
  )) {
    stop(
      "Size mismatch. Expected `", expected_size,
      "`, found `", result$Size, "`."
    )
  }
  
  if (!identical(
    tolower(as.character(result$sel)[1L]),
    tolower(expected_sel)
  )) {
    stop(
      "Selection-rule mismatch. Expected `", expected_sel,
      "`, found `", result$sel, "`."
    )
  }
  
  invisible(TRUE)
}


load_proxy_rt_results <- function(path_results,
                                  model,
                                  proxy_mode,
                                  sel,
                                  sizes = unit_proxy_sizes) {
  
  sizes <- tolower(as.character(sizes))
  
  if (!setequal(sizes, unit_proxy_sizes)) {
    stop(
      "`sizes` must contain exactly: ",
      paste(unit_proxy_sizes, collapse = ", ")
    )
  }
  
  sizes <- unit_proxy_sizes
  
  files <- vapply(
    sizes,
    function(size_i) {
      find_proxy_rt_file(
        path_results = path_results,
        model        = model,
        size         = size_i,
        sel          = sel
      )
    },
    character(1)
  )
  
  results <- lapply(files, readRDS)
  names(results) <- sizes
  
  for (size_i in sizes) {
    
    validate_proxy_rt_result(
      result              = results[[size_i]],
      expected_proxy_mode = proxy_mode,
      expected_size       = size_i,
      expected_sel        = sel
    )
  }
  
  list(
    files   = files,
    results = results
  )
}


# ==============================================================================
# 2. EVALUATION DATA FOR RMSFE AND DM TESTS
# ==============================================================================

unit_proxy_period <- function(date,
                              params) {
  
  covid_start <- as.Date(params$covid_start)
  covid_end   <- as.Date(params$covid_end)
  
  dplyr::case_when(
    date < covid_start ~ "Pre-COVID",
    date <= covid_end  ~ "COVID period",
    TRUE               ~ "Post-COVID"
  )
}


unit_proxy_quarter_id <- function(date) {
  paste0(
    lubridate::year(date),
    "Q",
    lubridate::quarter(date)
  )
}


build_unit_proxy_actuals <- function(rt_result,
                                     country_order) {
  
  y_q_all <- as.matrix(rt_result$Y_q_all)
  dates_q <- as.Date(rt_result$dates_q)
  
  if (nrow(y_q_all) != length(dates_q)) {
    stop("`Y_q_all` and `dates_q` have inconsistent dimensions.")
  }
  
  if (is.null(colnames(y_q_all))) {
    stop("`Y_q_all` must have country names.")
  }
  
  missing_countries <- setdiff(country_order, colnames(y_q_all))
  
  if (length(missing_countries) > 0L) {
    stop(
      "The following countries are missing from `Y_q_all`: ",
      paste(missing_countries, collapse = ", ")
    )
  }
  
  as.data.frame(
    y_q_all,
    check.names = FALSE
  ) |>
    tibble::as_tibble(.name_repair = "minimal") |>
    dplyr::mutate(
      date       = dates_q,
      quarter_id = unit_proxy_quarter_id(date),
      period     = unit_proxy_period(date, rt_result$params)
    ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(country_order),
      names_to  = "country",
      values_to = "GDP"
    ) |>
    dplyr::transmute(
      country    = as.character(country),
      quarter_id = as.character(quarter_id),
      GDP        = as.numeric(GDP),
      period     = as.character(period)
    )
}


build_unit_proxy_rt_for_dm <- function(rt_result,
                                       country_order) {
  
  rt_df <- tibble::as_tibble(rt_result$pseudo_rt_all)
  
  required_cols <- c(
    "date",
    "country",
    "nowcast",
    "month_in_quarter"
  )
  
  missing_cols <- setdiff(required_cols, names(rt_df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`pseudo_rt_all` is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- rt_df |>
    dplyr::transmute(
      date    = as.Date(date),
      country = as.character(country),
      nowcast = as.numeric(nowcast),
      type    = as.character(month_in_quarter)
    ) |>
    dplyr::filter(
      country %in% country_order,
      type %in% unit_proxy_vintage_order
    )
  
  duplicated_rows <- out |>
    dplyr::mutate(
      quarter_id = unit_proxy_quarter_id(date)
    ) |>
    dplyr::count(
      country,
      quarter_id,
      type
    ) |>
    dplyr::filter(n > 1L)
  
  if (nrow(duplicated_rows) > 0L) {
    stop(
      "More than one nowcast was found for at least one ",
      "country-quarter-vintage combination."
    )
  }
  
  out
}


build_unit_proxy_evaluation <- function(rt_result,
                                        country_order) {
  
  if (!exists("build_eval_df", mode = "function")) {
    stop(
      "`build_eval_df()` is not available. Source the results utilities ",
      "containing the Diebold--Mariano helpers before running this code."
    )
  }
  
  rt_df <- build_unit_proxy_rt_for_dm(
    rt_result    = rt_result,
    country_order = country_order
  )
  
  actual_df <- build_unit_proxy_actuals(
    rt_result    = rt_result,
    country_order = country_order
  )
  
  build_eval_df(
    df_rt = rt_df,
    df_y  = actual_df
  ) |>
    dplyr::mutate(
      country    = as.character(country),
      quarter_id = as.character(quarter_id),
      type       = as.character(type),
      period     = as.character(period)
    ) |>
    dplyr::filter(
      type %in% unit_proxy_vintage_order,
      period %in% unit_proxy_period_order,
      is.finite(GDP),
      is.finite(nowcast),
      is.finite(se)
    )
}


match_unit_and_aggregate_evaluations <- function(eval_unit,
                                                 eval_aggregate) {
  
  common <- eval_unit |>
    dplyr::transmute(
      country     = as.character(country),
      quarter_id  = as.character(quarter_id),
      type        = as.character(type),
      period_unit = as.character(period),
      GDP_unit    = as.numeric(GDP),
      se_unit     = as.numeric(se)
    ) |>
    dplyr::inner_join(
      eval_aggregate |>
        dplyr::transmute(
          country          = as.character(country),
          quarter_id       = as.character(quarter_id),
          type             = as.character(type),
          period_aggregate = as.character(period),
          GDP_aggregate    = as.numeric(GDP),
          se_aggregate     = as.numeric(se)
        ),
      by = c("country", "quarter_id", "type")
    )
  
  if (nrow(common) == 0L) {
    stop(
      "No common evaluable forecasts were found for the unit-target and ",
      "aggregate-target designs."
    )
  }
  
  period_mismatch <- common |>
    dplyr::filter(period_unit != period_aggregate)
  
  if (nrow(period_mismatch) > 0L) {
    stop(
      "The unit-target and aggregate-target designs assign different ",
      "evaluation periods to common forecast observations."
    )
  }
  
  actual_mismatch <- common |>
    dplyr::filter(
      abs(GDP_unit - GDP_aggregate) > 1e-12
    )
  
  if (nrow(actual_mismatch) > 0L) {
    stop(
      "The realised GDP values differ across the unit-target and ",
      "aggregate-target result files."
    )
  }
  
  common |>
    dplyr::transmute(
      country,
      quarter_id,
      type,
      period = period_unit,
      se_unit,
      se_aggregate
    )
}


# ==============================================================================
# 3. RMSFE RATIOS AND DIEBOLD--MARIANO TESTS
# ==============================================================================

compute_unit_proxy_comparison <- function(unit_result,
                                          aggregate_result,
                                          country_order) {
  
  if (!exists("build_dm_wide", mode = "function")) {
    stop(
      "`build_dm_wide()` is not available. Source the results utilities ",
      "containing the Diebold--Mariano helpers before running this code."
    )
  }
  
  eval_unit <- build_unit_proxy_evaluation(
    rt_result    = unit_result,
    country_order = country_order
  )
  
  eval_aggregate <- build_unit_proxy_evaluation(
    rt_result    = aggregate_result,
    country_order = country_order
  )
  
  common_eval <- match_unit_and_aggregate_evaluations(
    eval_unit      = eval_unit,
    eval_aggregate = eval_aggregate
  )
  
  rmsfe_table <- common_eval |>
    dplyr::group_by(
      country,
      period,
      type
    ) |>
    dplyr::summarise(
      n_forecasts     = dplyr::n(),
      rmsfe_unit      = sqrt(mean(se_unit)),
      rmsfe_aggregate = sqrt(mean(se_aggregate)),
      relative_rmsfe  = rmsfe_unit / rmsfe_aggregate,
      .groups         = "drop"
    )
  
  dm_result <- build_dm_wide(
    eval_matrix    = eval_unit,
    eval_benchmark = eval_aggregate,
    benchmark_tag  = "aggregate",
    alternative    = "less"
  )
  
  dm_wide <- dm_result$wide
  
  p_cols <- paste0(
    "p_aggregate_",
    unit_proxy_vintage_order
  )
  
  for (col_i in p_cols) {
    if (!col_i %in% names(dm_wide)) {
      dm_wide[[col_i]] <- NA_real_
    }
  }
  
  dm_long <- dm_wide |>
    dplyr::select(
      country,
      period,
      dplyr::all_of(p_cols)
    ) |>
    tidyr::pivot_longer(
      cols         = dplyr::all_of(p_cols),
      names_to     = "type",
      names_prefix = "p_aggregate_",
      values_to    = "p_value"
    ) |>
    dplyr::mutate(
      country = as.character(country),
      period  = as.character(period),
      type    = as.character(type)
    )
  
  rmsfe_table <- rmsfe_table |>
    dplyr::left_join(
      dm_long,
      by = c("country", "period", "type")
    )
  
  list(
    evaluation_unit      = eval_unit,
    evaluation_aggregate = eval_aggregate,
    common_evaluation    = common_eval,
    rmsfe                = rmsfe_table,
    dm                   = dm_result
  )
}


combine_unit_proxy_sizes <- function(comparisons_by_size,
                                     country_order) {
  
  if (!setequal(
    names(comparisons_by_size),
    unit_proxy_sizes
  )) {
    stop(
      "`comparisons_by_size` must contain: ",
      paste(unit_proxy_sizes, collapse = ", ")
    )
  }
  
  long_df <- dplyr::bind_rows(
    lapply(
      unit_proxy_sizes,
      function(size_i) {
        comparisons_by_size[[size_i]]$rmsfe |>
          dplyr::mutate(size = size_i)
      }
    )
  ) |>
    dplyr::mutate(
      country = as.character(country),
      period  = as.character(period),
      type    = as.character(type),
      size    = as.character(size)
    )
  
  grid <- tidyr::expand_grid(
    period  = unit_proxy_period_order,
    country = country_order
  )
  
  rmsfe_wide <- long_df |>
    dplyr::select(
      country,
      period,
      size,
      type,
      rmsfe_unit
    ) |>
    tidyr::pivot_wider(
      names_from  = c(size, type),
      values_from = rmsfe_unit,
      names_glue  = "{size}_{type}"
    )
  
  relative_wide <- long_df |>
    dplyr::select(
      country,
      period,
      size,
      type,
      relative_rmsfe,
      p_value
    ) |>
    tidyr::pivot_wider(
      names_from  = c(size, type),
      values_from = c(relative_rmsfe, p_value),
      names_glue  = "{.value}_{size}_{type}"
    )
  
  rmsfe_columns <- unlist(
    lapply(
      unit_proxy_sizes,
      function(size_i) {
        paste0(
          size_i,
          "_",
          unit_proxy_vintage_order
        )
      }
    )
  )
  
  relative_columns <- unlist(
    lapply(
      unit_proxy_sizes,
      function(size_i) {
        paste0(
          "relative_rmsfe_",
          size_i,
          "_",
          unit_proxy_vintage_order
        )
      }
    )
  )
  
  p_value_columns <- unlist(
    lapply(
      unit_proxy_sizes,
      function(size_i) {
        paste0(
          "p_value_",
          size_i,
          "_",
          unit_proxy_vintage_order
        )
      }
    )
  )
  
  for (col_i in rmsfe_columns) {
    if (!col_i %in% names(rmsfe_wide)) {
      rmsfe_wide[[col_i]] <- NA_real_
    }
  }
  
  for (col_i in c(relative_columns, p_value_columns)) {
    if (!col_i %in% names(relative_wide)) {
      relative_wide[[col_i]] <- NA_real_
    }
  }
  
  rmsfe_table <- grid |>
    dplyr::left_join(
      rmsfe_wide,
      by = c("country", "period")
    ) |>
    dplyr::arrange(
      factor(period, levels = unit_proxy_period_order),
      factor(country, levels = country_order)
    )
  
  relative_table <- grid |>
    dplyr::left_join(
      relative_wide,
      by = c("country", "period")
    ) |>
    dplyr::arrange(
      factor(period, levels = unit_proxy_period_order),
      factor(country, levels = country_order)
    )
  
  list(
    rmsfe_table    = rmsfe_table,
    relative_table = relative_table,
    long_table     = long_df
  )
}


# ==============================================================================
# 4. LATEX FORMATTERS
# ==============================================================================

unit_proxy_p_to_stars <- function(p_value) {
  
  dplyr::case_when(
    is.na(p_value)   ~ "",
    p_value < 0.01   ~ "***",
    p_value < 0.05   ~ "**",
    p_value < 0.10   ~ "*",
    TRUE             ~ ""
  )
}


unit_proxy_fmt_number <- function(x,
                                  digits = 3) {
  
  if (length(x) == 0L || !is.finite(x)) {
    return("--")
  }
  
  sprintf(
    paste0("%.", digits, "f"),
    x
  )
}


unit_proxy_fmt_relative <- function(x,
                                    p_value,
                                    digits = 3) {
  
  value <- unit_proxy_fmt_number(
    x      = x,
    digits = digits
  )
  
  stars <- unit_proxy_p_to_stars(p_value)
  
  if (
    is.finite(x) &&
    x < 1 &&
    nzchar(stars)
  ) {
    value <- paste0(
      value,
      "\\textsuperscript{\\fontsize{4.4}{4.8}\\selectfont ",
      stars,
      "}"
    )
  }
  
  value
}


unit_proxy_table_to_latex <- function(table_df,
                                      table_type = c("rmsfe", "relative"),
                                      caption,
                                      label,
                                      country_order,
                                      digits = 3,
                                      note = NULL) {
  
  table_type <- match.arg(table_type)
  
  value_columns <- unlist(
    lapply(
      unit_proxy_sizes,
      function(size_i) {
        if (table_type == "rmsfe") {
          paste0(size_i, "_", unit_proxy_vintage_order)
        } else {
          paste0(
            "relative_rmsfe_",
            size_i,
            "_",
            unit_proxy_vintage_order
          )
        }
      }
    )
  )
  
  required_columns <- c(
    "country",
    "period",
    value_columns
  )
  
  if (table_type == "relative") {
    
    p_columns <- unlist(
      lapply(
        unit_proxy_sizes,
        function(size_i) {
          paste0(
            "p_value_",
            size_i,
            "_",
            unit_proxy_vintage_order
          )
        }
      )
    )
    
    required_columns <- c(
      required_columns,
      p_columns
    )
  }
  
  missing_columns <- setdiff(
    required_columns,
    names(table_df)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "The table is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  if (is.null(note)) {
    
    if (table_type == "rmsfe") {
      note <- paste0(
        "\\textit{Notes:} Entries are unit-target Matrix MF--TPRF RMSFEs ",
        "computed over common evaluable forecasts."
      )
    } else {
      note <- paste0(
        "\\textit{Notes:} Entries are unit-target to aggregate-target RMSFE ",
        "ratios. Values below one favour unit-target. Asterisks denote ",
        "one-sided Diebold--Mariano rejections in favour of unit-target ",
        "based on squared errors and a Newey--West HAC variance estimator ",
        "($^{*}\\,p\\!<\\!0.10$, $^{**}\\,p\\!<\\!0.05$, ",
        "$^{***}\\,p\\!<\\!0.01$)."
      )
    }
  }
  
  table_df <- table_df |>
    dplyr::mutate(
      country = as.character(country),
      period  = as.character(period)
    ) |>
    dplyr::arrange(
      factor(period, levels = unit_proxy_period_order),
      factor(country, levels = country_order)
    )
  
  make_period_block <- function(period_i) {
    
    block_df <- table_df |>
      dplyr::filter(period == period_i) |>
      dplyr::arrange(
        factor(country, levels = country_order)
      )
    
    row_lines <- vapply(
      country_order,
      function(country_i) {
        
        row_i <- block_df |>
          dplyr::filter(country == country_i)
        
        values_i <- vapply(
          unit_proxy_sizes,
          function(size_i) {
            
            vapply(
              unit_proxy_vintage_order,
              function(vintage_i) {
                
                if (nrow(row_i) == 0L) {
                  return("\\utnum{--}")
                }
                
                if (table_type == "rmsfe") {
                  
                  value_col <- paste0(
                    size_i,
                    "_",
                    vintage_i
                  )
                  
                  paste0(
                    "\\utnum{",
                    unit_proxy_fmt_number(
                      x      = row_i[[value_col]][1L],
                      digits = digits
                    ),
                    "}"
                  )
                  
                } else {
                  
                  value_col <- paste0(
                    "relative_rmsfe_",
                    size_i,
                    "_",
                    vintage_i
                  )
                  
                  p_col <- paste0(
                    "p_value_",
                    size_i,
                    "_",
                    vintage_i
                  )
                  
                  paste0(
                    "\\utnum{",
                    unit_proxy_fmt_relative(
                      x       = row_i[[value_col]][1L],
                      p_value = row_i[[p_col]][1L],
                      digits  = digits
                    ),
                    "}"
                  )
                }
              },
              character(1)
            )
          },
          character(3)
        )
        
        paste0(
          paste(
            c(country_i, unlist(values_i, use.names = FALSE)),
            collapse = " & "
          ),
          " \\\\"
        )
      },
      character(1)
    )
    
    c(
      "\\addlinespace[0.18em]",
      "\\specialrule{0.06em}{0.10em}{0.16em}",
      "\\rowcolor{matrixgray}",
      paste0(
        "\\multicolumn{10}{@{}l}{",
        "\\fontsize{7.0}{7.7}\\selectfont\\textbf{",
        period_i,
        "}} \\\\[-0.18em]"
      ),
      "\\cmidrule{1-10}",
      row_lines
    )
  }
  
  body <- paste(
    unlist(
      lapply(
        unit_proxy_period_order,
        make_period_block
      ),
      use.names = FALSE
    ),
    collapse = "\n"
  )
  
  paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begingroup\n",
    "\\fontsize{7.6}{8.4}\\selectfont\n",
    "\\renewcommand{\\arraystretch}{0.86}\n",
    "\\setlength{\\tabcolsep}{1.5pt}\n",
    "\\newcommand{\\utnum}[1]{{\\fontsize{6.4}{7.0}\\selectfont #1}}\n",
    "\\definecolor{topgray}{gray}{0.92}\n",
    "\\definecolor{hypergray}{gray}{0.91}\n",
    "\\definecolor{matrixgray}{gray}{0.965}\n",
    "\\resizebox{0.98\\textwidth}{!}{%\n",
    "\\begin{tabular}{@{}",
    ">{\\centering\\arraybackslash}p{1.35cm}",
    "@{\\hspace{0.18cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{1.40cm}}",
    "@{\\hspace{0.18cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{1.40cm}}",
    "@{\\hspace{0.18cm}}",
    "*{3}{>{\\centering\\arraybackslash}p{1.40cm}}",
    "@{}}\n",
    "\\toprule\n",
    "\\rowcolor{hypergray}\n",
    " & \\multicolumn{3}{c}{\\textbf{Small}}",
    " & \\multicolumn{3}{c}{\\textbf{Medium}}",
    " & \\multicolumn{3}{c}{\\textbf{Large}} \\\\\n",
    "\\cmidrule(lr){2-4}",
    "\\cmidrule(lr){5-7}",
    "\\cmidrule(lr){8-10}\n",
    "\\rowcolor{topgray}\n",
    "\\textbf{Country}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n",
    "\\midrule\n",
    body,
    "\n\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\par\\vspace{0.05cm}\n",
    "\\parbox{0.90\\linewidth}{\\centering\\fontsize{6.3}{6.9}\\selectfont ",
    note,
    "}\n",
    "\\endgroup\n",
    "\\end{table}\n"
  )
}

# ==============================================================================
# 5. MAIN FUNCTION
# ==============================================================================

make_unit_target_proxy_tables_from_files <- function(
    path_results,
    sel,
    country_order = c(
      "DE", "FR", "IT", "ES",
      "NL", "BE", "AT", "PT"
    ),
    model_unit = "matrix_multivariate",
    model_aggregate = "matrix_scalar",
    digits = 3,
    caption_rmsfe = NULL,
    caption_relative = NULL,
    label_rmsfe = NULL,
    label_relative = NULL
) {
  
  unit_loaded <- load_proxy_rt_results(
    path_results = path_results,
    model        = model_unit,
    proxy_mode   = "multivariate",
    sel          = sel
  )
  
  aggregate_loaded <- load_proxy_rt_results(
    path_results = path_results,
    model        = model_aggregate,
    proxy_mode   = "scalar",
    sel          = sel
  )
  
  comparisons_by_size <- lapply(
    unit_proxy_sizes,
    function(size_i) {
      
      compute_unit_proxy_comparison(
        unit_result      = unit_loaded$results[[size_i]],
        aggregate_result = aggregate_loaded$results[[size_i]],
        country_order    = country_order
      )
    }
  )
  
  names(comparisons_by_size) <- unit_proxy_sizes
  
  combined <- combine_unit_proxy_sizes(
    comparisons_by_size = comparisons_by_size,
    country_order       = country_order
  )
  
  sel_label <- if (toupper(sel) == "LASSO") {
    "LASSO-based selection"
  } else {
    "correlation screening"
  }
  
  sel_tag <- if (toupper(sel) == "LASSO") {
    "lasso"
  } else {
    "corr"
  }
  
  if (is.null(caption_rmsfe)) {
    caption_rmsfe <- paste0(
      "Unit-target Matrix MF--TPRF RMSFE: ",
      sel_label
    )
  }
  
  if (is.null(caption_relative)) {
    caption_relative <- paste0(
      "Relative RMSFE of unit-target to aggregate-target Matrix MF--TPRF: ",
      sel_label
    )
  }
  
  if (is.null(label_rmsfe)) {
    label_rmsfe <- paste0(
      "tab:unit_target_rmsfe_",
      sel_tag
    )
  }
  
  if (is.null(label_relative)) {
    label_relative <- paste0(
      "tab:unit_target_relative_rmsfe_",
      sel_tag
    )
  }
  
  latex_rmsfe <- unit_proxy_table_to_latex(
    table_df     = combined$rmsfe_table,
    table_type   = "rmsfe",
    caption      = caption_rmsfe,
    label        = label_rmsfe,
    country_order = country_order,
    digits       = digits
  )
  
  latex_relative <- unit_proxy_table_to_latex(
    table_df     = combined$relative_table,
    table_type   = "relative",
    caption      = caption_relative,
    label        = label_relative,
    country_order = country_order,
    digits       = digits
  )
  
  invisible(
    list(
      files = list(
        unit      = unit_loaded$files,
        aggregate = aggregate_loaded$files
      ),
      comparisons_by_size = comparisons_by_size,
      rmsfe_table         = combined$rmsfe_table,
      relative_table      = combined$relative_table,
      long_table          = combined$long_table,
      latex_rmsfe         = latex_rmsfe,
      latex_relative      = latex_relative
    )
  )
}

# ==============================================================================
# EA NOWCAST COMPARISON HELPERS
# ==============================================================================

print_ea_console_table <- function(x) {
  
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.")
  }
  
  if (nrow(x) == 0L) {
    cat("No rows to print.\n")
    return(invisible(x))
  }
  
  print(
    as.data.frame(x),
    row.names = FALSE
  )
  
  invisible(x)
}


make_ea_quarter_id <- function(date) {
  
  date <- as.Date(date)
  
  paste0(
    lubridate::year(date),
    "Q",
    lubridate::quarter(date)
  )
}


extract_matrix_ea_nowcasts <- function(
    rt_obj,
    month_order = c("M1", "M2", "M3")
) {
  
  if (is.null(rt_obj$pseudo_rt_aggregate)) {
    stop("Matrix RT object does not contain `pseudo_rt_aggregate`.")
  }
  
  required_cols <- c(
    "date",
    "nowcast",
    "month_in_quarter"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(rt_obj$pseudo_rt_aggregate)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Matrix EA nowcast table is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  rt_obj$pseudo_rt_aggregate %>%
    dplyr::transmute(
      date    = as.Date(date),
      type    = as.character(month_in_quarter),
      nowcast = as.numeric(nowcast),
      model   = "Matrix MF-TPRF"
    ) %>%
    dplyr::filter(type %in% month_order)
}


extract_vector_ea_nowcasts <- function(
    rt_obj,
    month_order = c("M1", "M2", "M3")
) {
  
  if (is.null(rt_obj$pseudo_realtime_all)) {
    stop("Vector EA RT object does not contain `pseudo_realtime_all`.")
  }
  
  required_cols <- c(
    "date",
    "nowcast",
    "month_in_quarter"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(rt_obj$pseudo_realtime_all)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Vector EA nowcast table is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  rt_obj$pseudo_realtime_all %>%
    dplyr::transmute(
      date    = as.Date(date),
      type    = as.character(month_in_quarter),
      nowcast = as.numeric(nowcast),
      model   = "VEC-C"
    ) %>%
    dplyr::filter(type %in% month_order)
}


extract_matrix_ea_gdp <- function(rt_obj) {
  
  if (is.null(rt_obj$Y_q_all) || is.null(rt_obj$dates_q)) {
    stop("Matrix RT object does not contain `Y_q_all` and/or `dates_q`.")
  }
  
  y_q_all <- as.matrix(rt_obj$Y_q_all)
  dates_q <- as.Date(rt_obj$dates_q)
  
  if (nrow(y_q_all) != length(dates_q)) {
    stop("Matrix EA GDP and quarterly dates have inconsistent lengths.")
  }
  
  if (is.null(colnames(y_q_all))) {
    stop("Matrix `Y_q_all` has no column names.")
  }
  
  ea_col <- which(
    toupper(colnames(y_q_all)) == "EA"
  )
  
  if (length(ea_col) != 1L) {
    stop("Could not uniquely identify the EA GDP column in Matrix `Y_q_all`.")
  }
  
  data.frame(
    date       = dates_q,
    quarter_id = make_ea_quarter_id(dates_q),
    GDP        = as.numeric(y_q_all[, ea_col]),
    row.names  = NULL
  )
}


extract_vector_ea_gdp <- function(rt_obj) {
  
  if (is.null(rt_obj$y_q) || is.null(rt_obj$dates_q)) {
    stop("Vector EA RT object does not contain `y_q` and/or `dates_q`.")
  }
  
  y_q     <- as.numeric(rt_obj$y_q)
  dates_q <- as.Date(rt_obj$dates_q)
  
  if (length(y_q) != length(dates_q)) {
    stop("Vector EA GDP and quarterly dates have inconsistent lengths.")
  }
  
  data.frame(
    date       = dates_q,
    quarter_id = make_ea_quarter_id(dates_q),
    GDP        = y_q,
    row.names  = NULL
  )
}


check_unique_nowcast_keys <- function(df, label) {
  
  duplicates <- df %>%
    dplyr::count(date, type) %>%
    dplyr::filter(n > 1L)
  
  if (nrow(duplicates) > 0L) {
    stop(
      label,
      " contains duplicate nowcasts for at least one date-vintage pair."
    )
  }
}


check_common_ea_settings <- function(matrix_rt, vector_rt, sel_i) {
  
  matrix_params <- matrix_rt$params
  vector_params <- vector_rt$params
  
  fields_to_match <- c(
    "start_eval",
    "covid_start",
    "covid_end",
    "covid_mask_m",
    "covid_mask_q",
    "n_m",
    "n_q"
  )
  
  for (field_i in fields_to_match) {
    
    if (is.null(matrix_params[[field_i]]) ||
        is.null(vector_params[[field_i]])) {
      stop(
        "Missing parameter `", field_i,
        "` in one of the EA RT objects."
      )
    }
    
    value_matrix <- paste(
      as.character(matrix_params[[field_i]]),
      collapse = "|"
    )
    
    value_vector <- paste(
      as.character(vector_params[[field_i]]),
      collapse = "|"
    )
    
    if (!identical(value_matrix, value_vector)) {
      stop(
        "Matrix and Vector EA runs differ in parameter `",
        field_i,
        "` for selection method ", sel_i, "."
      )
    }
  }
  
  if (!identical(
    tolower(as.character(matrix_rt$sel)),
    tolower(as.character(sel_i))
  )) {
    stop("Loaded Matrix RT object does not match selection method: ", sel_i)
  }
  
  if (!identical(
    tolower(as.character(vector_rt$sel)),
    tolower(as.character(sel_i))
  )) {
    stop("Loaded Vector RT object does not match selection method: ", sel_i)
  }
  
  if (!identical(
    as.character(matrix_params$end_eval),
    as.character(vector_params$end_eval)
  )) {
    cat(
      "\nEA comparison (", sel_i, "): Matrix and Vector have different ",
      "`end_eval` dates. Only common nowcast dates are used.\n",
      sep = ""
    )
  }
}


build_ea_nowcast_evaluation <- function(
    df_nowcasts,
    df_gdp,
    selection
) {
  
  required_nowcasts <- c(
    "date",
    "type",
    "nowcast",
    "model"
  )
  
  required_gdp <- c(
    "date",
    "quarter_id",
    "GDP"
  )
  
  missing_nowcasts <- setdiff(
    required_nowcasts,
    names(df_nowcasts)
  )
  
  missing_gdp <- setdiff(
    required_gdp,
    names(df_gdp)
  )
  
  if (length(missing_nowcasts) > 0L) {
    stop(
      "Nowcast data are missing: ",
      paste(missing_nowcasts, collapse = ", ")
    )
  }
  
  if (length(missing_gdp) > 0L) {
    stop(
      "GDP data are missing: ",
      paste(missing_gdp, collapse = ", ")
    )
  }
  
  df_nowcasts_clean <- df_nowcasts %>%
    dplyr::mutate(
      date       = as.Date(date),
      quarter_id = make_ea_quarter_id(date)
    )
  
  df_gdp_clean <- df_gdp %>%
    dplyr::transmute(
      quarter_id = as.character(quarter_id),
      target_date = as.Date(date),
      GDP = as.numeric(GDP)
    ) %>%
    dplyr::filter(
      is.finite(GDP)
    )
  
  gdp_duplicates <- df_gdp_clean %>%
    dplyr::distinct(
      quarter_id,
      GDP
    ) %>%
    dplyr::count(quarter_id) %>%
    dplyr::filter(n > 1L)
  
  if (nrow(gdp_duplicates) > 0L) {
    stop(
      "More than one realised EA GDP value was found for at least one quarter."
    )
  }
  
  df_gdp_clean <- df_gdp_clean %>%
    dplyr::distinct(
      quarter_id,
      .keep_all = TRUE
    )
  
  df_nowcasts_clean %>%
    dplyr::inner_join(
      df_gdp_clean,
      by = "quarter_id"
    ) %>%
    dplyr::mutate(
      country   = "EA",
      selection = as.character(selection),
      error     = nowcast - GDP,
      abs_error = abs(error),
      se        = error^2
    ) %>%
    dplyr::arrange(
      type,
      model,
      date
    )
}


add_ea_evaluation_periods <- function(
    df_eval,
    params,
    include_full_sample = TRUE
) {
  
  required_cols <- c(
    "target_date",
    "country",
    "quarter_id",
    "type",
    "model",
    "selection",
    "se"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(df_eval)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Evaluation data are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  covid_start <- as.Date(params$covid_start)
  covid_end   <- as.Date(params$covid_end)
  
  df_periods <- df_eval %>%
    dplyr::mutate(
      period = dplyr::case_when(
        target_date < covid_start ~ "Pre-COVID",
        target_date <= covid_end  ~ "COVID period",
        target_date > covid_end   ~ "Post-COVID",
        TRUE                      ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(period)
    )
  
  if (isTRUE(include_full_sample)) {
    
    df_periods <- dplyr::bind_rows(
      df_periods %>%
        dplyr::mutate(
          period = "Full sample"
        ),
      df_periods
    )
  }
  
  df_periods %>%
    dplyr::mutate(
      period = factor(
        period,
        levels = c(
          "Full sample",
          "Pre-COVID",
          "COVID period",
          "Post-COVID"
        )
      )
    ) %>%
    dplyr::arrange(
      selection,
      period,
      type,
      model,
      date
    )
}


report_ea_nowcast_deviations <- function(
    df_eval,
    threshold = 0.05,
    selection_label = NULL
) {
  
  threshold_pp <- 100 * threshold
  
  df_outliers <- df_eval %>%
    dplyr::filter(
      is.finite(GDP),
      is.finite(nowcast),
      is.finite(abs_error),
      abs_error > threshold
    ) %>%
    dplyr::mutate(
      GDP_pct           = 100 * GDP,
      nowcast_pct       = 100 * nowcast,
      forecast_error_pp = 100 * error,
      abs_error_pp      = 100 * abs_error
    ) %>%
    dplyr::select(
      date,
      target_date,
      type,
      model,
      GDP_pct,
      nowcast_pct,
      forecast_error_pp,
      abs_error_pp
    ) %>%
    dplyr::arrange(
      type,
      model,
      date
    )
  
  cat(
    "\n",
    paste(rep("=", 88), collapse = ""),
    "\n",
    "EA NOWCAST DEVIATIONS ABOVE ",
    format(threshold_pp, nsmall = 1),
    " PERCENTAGE POINTS",
    if (!is.null(selection_label)) {
      paste0(" | ", selection_label)
    } else {
      ""
    },
    "\n",
    paste(rep("=", 88), collapse = ""),
    "\n",
    sep = ""
  )
  
  if (nrow(df_outliers) == 0L) {
    
    cat(
      "No Matrix or VEC-C EA nowcast differs from realised EA GDP by more than ",
      format(threshold_pp, nsmall = 1),
      " percentage points.\n",
      sep = ""
    )
    
  } else {
    
    print_ea_console_table(df_outliers)
  }
  
  invisible(df_outliers)
}


get_ea_gdp_ylim <- function(
    df_gdp,
    pad_frac = 0,
    min_pad = 0
) {
  
  if (!is.data.frame(df_gdp) || !"GDP" %in% names(df_gdp)) {
    stop("`df_gdp` must be a data.frame containing `GDP`.")
  }
  
  y_range <- range(
    df_gdp$GDP,
    na.rm = TRUE
  )
  
  if (!all(is.finite(y_range))) {
    stop("Cannot compute EA y-axis limits: GDP is entirely missing.")
  }
  
  if (diff(y_range) == 0) {
    y_pad <- max(
      abs(y_range[1]) * 0.05,
      0.001
    )
  } else {
    y_pad <- max(
      diff(y_range) * pad_frac,
      min_pad
    )
  }
  
  c(
    y_range[1] - y_pad,
    y_range[2] + y_pad
  )
}


# ==============================================================================
# REPLACE EVERYTHING FROM compute_ea_rmsfe_by_period() TO THE END OF
# matrix.mf.tprf.utils.R WITH THIS BLOCK
# ==============================================================================

compute_ea_rmsfe_by_period <- function(df_eval_periods) {
  
  df_eval_periods %>%
    dplyr::filter(
      is.finite(se)
    ) %>%
    dplyr::group_by(
      selection,
      period,
      model,
      type
    ) %>%
    dplyr::summarise(
      N     = dplyr::n(),
      RMSFE = sqrt(mean(se)),
      .groups = "drop"
    )
}


build_ea_rmsfe_table <- function(
    df_rmsfe,
    selection_order = c("LASSO", "corr"),
    period_order = c(
      "Full sample",
      "Pre-COVID",
      "COVID period",
      "Post-COVID"
    ),
    month_order = c("M1", "M2", "M3"),
    digits = 3
) {
  
  required_cols <- c(
    "selection",
    "period",
    "model",
    "type",
    "RMSFE"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(df_rmsfe)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "RMSFE data are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  table_rmsfe <- df_rmsfe %>%
    dplyr::mutate(
      selection = dplyr::case_when(
        tolower(as.character(selection)) == "lasso" ~ "LASSO",
        tolower(as.character(selection)) == "corr"  ~ "corr",
        TRUE ~ as.character(selection)
      ),
      period = as.character(period),
      type   = as.character(type),
      model_key = dplyr::case_when(
        model == "Matrix MF-TPRF" ~ "Matrix",
        model == "VEC-C"          ~ "VecC",
        TRUE                      ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      !is.na(model_key),
      type %in% month_order
    ) %>%
    dplyr::select(
      selection,
      period,
      model_key,
      type,
      RMSFE
    ) %>%
    tidyr::pivot_wider(
      names_from  = c(model_key, type),
      values_from = RMSFE,
      names_glue  = "{model_key}_{type}"
    )
  
  required_value_cols <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "VecC_M1", "VecC_M2", "VecC_M3"
  )
  
  for (col_i in required_value_cols) {
    if (!col_i %in% names(table_rmsfe)) {
      table_rmsfe[[col_i]] <- NA_real_
    }
  }
  
  table_rmsfe %>%
    dplyr::mutate(
      Relative_M1 = Matrix_M1 / VecC_M1,
      Relative_M2 = Matrix_M2 / VecC_M2,
      Relative_M3 = Matrix_M3 / VecC_M3,
      selection   = factor(
        selection,
        levels = selection_order
      ),
      period = factor(
        period,
        levels = period_order
      )
    ) %>%
    dplyr::arrange(
      selection,
      period
    ) %>%
    dplyr::mutate(
      dplyr::across(
        c(
          Matrix_M1,
          Matrix_M2,
          Matrix_M3,
          VecC_M1,
          VecC_M2,
          VecC_M3,
          Relative_M1,
          Relative_M2,
          Relative_M3
        ),
        ~ round(.x, digits)
      ),
      selection = as.character(selection),
      period    = as.character(period)
    )
}


make_latex_ea_rmsfe_table <- function(
    df,
    caption = "Pseudo-real-time EA GDP RMSFE: Matrix MF--TPRF versus VEC-C",
    label = "tab:ea_rmsfe_matrix_vector_small",
    digits = 3
) {
  
  required_cols <- c(
    "selection",
    "period",
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "VecC_M1", "VecC_M2", "VecC_M3",
    "Relative_M1", "Relative_M2", "Relative_M3"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(df)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "RMSFE table is missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  fmt <- function(x) {
    ifelse(
      is.na(x),
      "--",
      sprintf(
        paste0("%.", digits, "f"),
        x
      )
    )
  }
  
  selection_labels <- c(
    "LASSO" = "LASSO screening",
    "corr"  = "Correlation screening"
  )
  
  build_rows <- function(selection_i) {
    
    df_i <- df %>%
      dplyr::filter(
        selection == selection_i
      )
    
    if (nrow(df_i) == 0L) {
      return(character(0))
    }
    
    selection_header <- paste0(
      "\\rowcolor{hypergray}\n",
      "\\multicolumn{10}{c}{\\textbf{",
      selection_labels[[selection_i]],
      "}} \\\\"
    )
    
    row_lines <- vapply(
      seq_len(nrow(df_i)),
      function(i) {
        
        paste(
          c(
            df_i$period[i],
            fmt(df_i$Matrix_M1[i]),
            fmt(df_i$Matrix_M2[i]),
            fmt(df_i$Matrix_M3[i]),
            fmt(df_i$VecC_M1[i]),
            fmt(df_i$VecC_M2[i]),
            fmt(df_i$VecC_M3[i]),
            fmt(df_i$Relative_M1[i]),
            fmt(df_i$Relative_M2[i]),
            fmt(df_i$Relative_M3[i])
          ),
          collapse = " & "
        )
      },
      character(1)
    )
    
    c(
      selection_header,
      "\\cmidrule(lr){1-10}",
      paste0(
        row_lines,
        " \\\\"
      )
    )
  }
  
  body <- c(
    build_rows("LASSO"),
    "\\addlinespace[0.35em]",
    build_rows("corr")
  )
  
  paste0(
    "\\begin{table}[!htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\renewcommand{\\arraystretch}{1.03}\n",
    "\\setlength{\\tabcolsep}{3.0pt}\n",
    "\\definecolor{topgray}{gray}{0.92}\n",
    "\\definecolor{hypergray}{gray}{0.91}\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\resizebox{\\textwidth}{!}{%\n",
    "\\begin{tabular}{lccc@{\\hspace{0.40cm}}ccc@{\\hspace{0.40cm}}ccc}\n",
    "\\toprule\n",
    " & \\multicolumn{3}{c}{\\textbf{Matrix MF--TPRF}}",
    " & \\multicolumn{3}{c}{\\textbf{VEC-C}}",
    " & \\multicolumn{3}{c}{\\textbf{Matrix / VEC-C}} \\\\\n",
    "\\cmidrule(lr){2-4}",
    "\\cmidrule(lr){5-7}",
    "\\cmidrule(lr){8-10}\n",
    "\\rowcolor{topgray}\n",
    "\\textbf{Period}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3}",
    " & \\textbf{M1} & \\textbf{M2} & \\textbf{M3} \\\\\n",
    "\\midrule\n",
    paste(
      body,
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}%\n",
    "}\n",
    "\\vspace{0.05cm}\n",
    "\\parbox{0.90\\linewidth}{\\centering\\footnotesize ",
    "\\textit{Notes:} RMSFEs are computed over common Matrix--VEC-C EA ",
    "pseudo-real-time forecasts. Relative RMSFEs below one favour Matrix MF--TPRF.}",
    "\n",
    "\\end{table}\n"
  )
}