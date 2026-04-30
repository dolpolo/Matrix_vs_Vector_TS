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
                                         country_cols = country_order,
                                         check_symbol = "\\checkmark") {
  
  if (!is.data.frame(df_selection_wide)) {
    stop("df_selection_wide must be a data.frame.")
  }
  
  df <- df_selection_wide
  
  required_cols <- c(
    "base_name",
    "frequency",
    "model_size",
    country_cols
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing columns in df_selection_wide: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
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
    df[, c("base_name", "frequency", "model_size", country_cols)],
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
    "NFCLB.SLN",  "Non-Financial Corporations - Liabilities - Short-Term Loans",            "F", "H", 1, "Q", 45, "Credit Aggregates",
    "NFCLB.LLN",  "Non-Financial Corporations - Liabilities - Long-Term Loans",             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS",      "General Government: Total Financial Assets",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.SLN",  "General Government - Assets: Short-Term Loans",                          "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGASS.LLN",  "General Government - Assets: Long-Term Loans",                           "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB",       "General Government: Total Financial Liabilities",                        "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.SLN",   "General Government - Liabilities: Short-Term Loans",                     "F", "H", 1, "Q", 45, "Credit Aggregates",
    "GGLB.LLN",   "General Government - Liabilities: Long-Term Loans",                      "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB",       "Households: Total Financial Liabilities",                                "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.SLN",   "Households - Liabilities: Short-Term Loans",                             "F", "H", 1, "Q", 45, "Credit Aggregates",
    "HHLB.LLN",   "Households - Liabilities: Long-Term Loans",                              "F", "H", 1, "Q", 45, "Credit Aggregates",
    
    # Labor Costs
    "ULCCON",     "Nominal Unit Labor Costs: Construction",                                 "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCIN",      "Nominal Unit Labor Costs: Industry",                                     "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCMN",      "Nominal Unit Labor Costs: Manufacturing",                                "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCFC",      "Nominal Unit Labor Costs: Financial Activities",                         "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRE",      "Nominal Unit Labor Costs: Real Estate",                                  "N", "H", 1, "Q", 45, "Labor Costs",
    "ULCPR",      "Nominal Unit Labor Costs: Professional, Scientific, Technical activities","N", "H", 1, "Q", 45, "Labor Costs",
    "ULCRT",      "Nominal Unit Labor Costs: Wholesale/Retail Trade, Transport, Food, IT",  "N", "H", 1, "Q", 45, "Labor Costs",
    
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

build_latex_selection_table <- function(sel_small,
                                        sel_medium,
                                        sel_large,
                                        metadata = build_metadata_table(),
                                        country_cols = country_order,
                                        caption = "Macroeconomic predictors: country-specific selection across information-set sizes",
                                        label = "tab:variables_selected_country_size") {
  
  sel_all <- build_selection_size_table(
    sel_small    = sel_small,
    sel_medium   = sel_medium,
    sel_large    = sel_large,
    country_cols = country_cols
  ) %>%
    dplyr::rename(ID = base_name, Fr = frequency) %>%
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
  
  tab <- metadata %>%
    dplyr::filter(ID != "GDP") %>%
    dplyr::inner_join(sel_all, by = c("ID", "Fr")) %>%
    dplyr::mutate(
      dplyr::across(dplyr::all_of(country_cols), to_setcell),
      ID     = escape_latex(ID),
      Series = escape_latex(Series),
      Cl     = escape_latex(Cl),
      Cat    = escape_latex(Cat),
      Fr     = escape_latex(Fr),
      Group  = factor(Group, levels = groups_order)
    ) %>%
    dplyr::arrange(Group, Fr, ID) %>%
    dplyr::mutate(N = dplyr::row_number())
  
  header <- paste0(
    "\\begin{table}[p]\n",
    "\\centering\n",
    "\\tiny\n",
    "\\renewcommand{\\arraystretch}{0.82}\n",
    "\\setlength{\\tabcolsep}{1.6pt}\n",
    "\\newcommand{\\setcell}[1]{{\\fontsize{4.1}{4.4}\\selectfont #1}}\n\n",
    "\\caption{", escape_latex(caption), "}\n",
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
  ncols_total <- 8 + length(country_cols)
  
  for (grp in groups_order) {
    subtab <- tab %>% dplyr::filter(Group == grp)
    
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
    
    rows <- subtab %>%
      dplyr::select(
        N, ID, Series, Cl, Cat, Tr, Fr, Del,
        dplyr::all_of(country_cols)
      )
    
    row_lines <- apply(
      rows,
      1,
      function(r) paste0(paste(r, collapse = " & "), " \\\\")
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
  
  if (is.null(group_vars)) {
    
    out <- df %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        `DFM wins`        = if (isTRUE(include_dfm) && "win_dfm" %in% names(.)) {
          sum(win_dfm, na.rm = TRUE)
        } else {
          NA_real_
        },
        `DFM win share`   = if (isTRUE(include_dfm) && "win_dfm" %in% names(.)) {
          100 * mean(win_dfm, na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      ) %>%
      dplyr::mutate(Group = "Overall", .before = 1)
    
  } else {
    
    out <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        `DFM wins`        = if (isTRUE(include_dfm) && "win_dfm" %in% names(.)) {
          sum(win_dfm, na.rm = TRUE)
        } else {
          NA_real_
        },
        `DFM win share`   = if (isTRUE(include_dfm) && "win_dfm" %in% names(.)) {
          100 * mean(win_dfm, na.rm = TRUE)
        } else {
          NA_real_
        },
        .groups = "drop"
      )
  }
  
  if (!isTRUE(include_dfm) || !"win_dfm" %in% names(df)) {
    out <- out %>%
      dplyr::select(-dplyr::any_of(c("DFM wins", "DFM win share")))
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
  
  summarise_expr <- function(data) {
    out <- data %>%
      dplyr::summarise(
        Total             = dplyr::n(),
        `VEC-P wins`      = sum(win_vecp, na.rm = TRUE),
        `VEC-P win share` = 100 * mean(win_vecp, na.rm = TRUE),
        `VEC-C wins`      = sum(win_vecc, na.rm = TRUE),
        `VEC-C win share` = 100 * mean(win_vecc, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (isTRUE(include_dfm) && "win_dfm" %in% names(data)) {
      out <- out %>%
        dplyr::mutate(
          `DFM wins`      = sum(data$win_dfm, na.rm = TRUE),
          `DFM win share` = 100 * mean(data$win_dfm, na.rm = TRUE)
        )
    }
    
    out
  }
  
  if (is.null(group_vars)) {
    out <- summarise_expr(df) %>%
      dplyr::mutate(Group = "Overall", .before = 1)
  } else {
    out <- df %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
      summarise_expr()
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
