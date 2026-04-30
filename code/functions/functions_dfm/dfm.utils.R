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


# ==============================================================================
# 3. LOCAL HELPERS
# ==============================================================================

build_run_tag_dfm <- function(params) {
  paste0(
    "eval-", format(params$start_eval, "%Y%m"), "_", format(params$end_eval, "%Y%m"),
    "_sel-", params$sel_method,
    "_Nm-", params$n_m,
    "_Nq-", params$n_q,
    "_Kmax-", params$Kmax,
    "_pmax-", params$pmax,
    "_qmax-", params$qmax,
    "_CovidM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-", as.integer(isTRUE(params$covid_mask_q))
  )
}

get_size_tag_dfm <- function(N_m, N_q) {
  paste0("Nm", N_m, "_Nq", N_q)
}

build_result_filename_dfm <- function(path_out,
                                      model,
                                      stage,
                                      Size,
                                      sel,
                                      countries,
                                      N_m,
                                      N_q,
                                      r,
                                      p,
                                      q,
                                      covid_m,
                                      covid_q,
                                      ext = "rds",
                                      timestamp = TRUE) {
  stamp <- if (timestamp) paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")) else ""
  
  file.path(
    path_out,
    paste0(
      toupper(model), "_",
      stage,
      "_Size-", Size,
      "_sel-", sel,
      "_cc-", countries,
      "_Nm-", N_m,
      "_Nq-", N_q,
      "_r-", r,
      "_p-", p,
      "_q-", q,
      "_CovidM-", covid_m,
      "_CovidQ-", covid_q,
      stamp,
      ".",
      ext
    )
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

# ==============================================================================
# 2. HELPERS
# ==============================================================================

make_quarter_avg_from_monthly <- function(x, T_q) {
  x <- as.numeric(x)
  stopifnot(length(x) >= 3 * T_q)
  x_use <- x[seq_len(3 * T_q)]
  as.numeric(tapply(x_use, rep(seq_len(T_q), each = 3), mean))
}

rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

latex_table_periods <- function(mat, caption, label) {
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{lccc}\n",
    "\\toprule\n",
    "Period & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf(
        "%s & %.4f & %.4f & %.4f \\\\",
        rownames(mat),
        mat[, "M1"], mat[, "M2"], mat[, "M3"]
      ),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

destd_mat <- function(X_std, mu, sd) {
  X_std <- as.matrix(X_std)
  stopifnot(length(mu) == ncol(X_std), length(sd) == ncol(X_std))
  sweep(sweep(X_std, 2, sd, `*`), 2, mu, `+`)
}

# ---- file finder like MF-TPRF ----
find_result_file_dfm <- function(path, stage = c("fit", "rt", "summary"),
                                 model = NULL, sel = NULL, country = NULL) {
  stage <- match.arg(stage)
  
  files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0L) stop("No .rds files found in: ", path)
  
  base <- basename(files)
  
  patt_stage <- switch(
    stage,
    fit     = "fit",
    rt      = "rt|RollingNowcast|rolling",
    summary = "summary"
  )
  
  keep <- grepl(patt_stage, base, ignore.case = TRUE)
  
  if (!is.null(model)) {
    keep <- keep & grepl(model, base, ignore.case = TRUE)
  }
  if (!is.null(sel)) {
    keep <- keep & grepl(paste0("sel-", sel), base, ignore.case = TRUE)
  }
  if (!is.null(country)) {
    keep <- keep & grepl(country, dirname(files)) | keep & grepl(country, base, ignore.case = TRUE)
  }
  
  cand <- files[keep]
  if (length(cand) == 0L) {
    stop("No matching file found in ", path, " for stage = ", stage)
  }
  
  cand[which.max(file.info(cand)$mtime)]
}

extract_params_object <- function(obj) {
  if (!is.null(obj$params)) return(obj$params)
  stop("No params object found in RDS.")
}

extract_dates_y_dfm <- function(obj) {
  dates_m <- if (!is.null(obj$dates_m)) as.Date(obj$dates_m) else NULL
  dates_q <- if (!is.null(obj$dates_q)) as.Date(obj$dates_q) else NULL
  y_q     <- if (!is.null(obj$y_q)) as.numeric(obj$y_q) else NULL
  list(dates_m = dates_m, dates_q = dates_q, y_q = y_q)
}

extract_metadata_dfm <- function(obj) {
  if (!is.null(obj$metadata)) return(obj$metadata)
  list()
}

extract_preprocessing_dfm <- function(obj) {
  if (!is.null(obj$preprocessing)) return(obj$preprocessing)
  list()
}

extract_hyper_dfm <- function(obj) {
  if (!is.null(obj$hyper)) return(obj$hyper)
  list(r = NA_real_, p = NA_real_, q = NA_real_)
}

extract_fit_dfm <- function(obj) {
  if (!is.null(obj$fit)) return(obj$fit)
  stop("No fit object found in full-sample DFM RDS.")
}

# ------------------------------------------------------------------------------
# rolling helper: supports both new and old structures
# ------------------------------------------------------------------------------
extract_rt_nowcasts_dfm <- function(rt_obj) {
  
  # --------------------------------------------------------------------------
  # Case 1: rolling object saved with $nowcast
  # --------------------------------------------------------------------------
  if (!is.null(rt_obj$nowcast)) {
    now_obj <- rt_obj$nowcast
    
    return(list(
      r_fix   = now_obj$r_fix %||% now_obj$hyper_pre$r,
      p_fix   = now_obj$p_fix %||% now_obj$hyper_pre$p,
      q_fix   = now_obj$q_fix %||% now_obj$hyper_pre$q,
      M1_std  = now_obj$M1_std,
      M2_std  = now_obj$M2_std,
      M3_std  = now_obj$M3_std,
      M1_orig = now_obj$M1_orig,
      M2_orig = now_obj$M2_orig,
      M3_orig = now_obj$M3_orig
    ))
  }
  
  # --------------------------------------------------------------------------
  # Case 2: new saved structure with $pseudo_realtime_raw
  # --------------------------------------------------------------------------
  if (!is.null(rt_obj$pseudo_realtime_raw)) {
    now_obj <- rt_obj$pseudo_realtime_raw
    
    return(list(
      r_fix   = if (!is.null(now_obj$hyper_pre$r)) now_obj$hyper_pre$r else now_obj$r_fix,
      p_fix   = if (!is.null(now_obj$hyper_pre$p)) now_obj$hyper_pre$p else now_obj$p_fix,
      q_fix   = if (!is.null(now_obj$hyper_pre$q)) now_obj$hyper_pre$q else now_obj$q_fix,
      M1_std  = now_obj$M1_std,
      M2_std  = now_obj$M2_std,
      M3_std  = now_obj$M3_std,
      M1_orig = now_obj$M1_orig,
      M2_orig = now_obj$M2_orig,
      M3_orig = now_obj$M3_orig
    ))
  }
  
  # --------------------------------------------------------------------------
  # Case 3: flat structure with top-level nowcasts
  # --------------------------------------------------------------------------
  if (!is.null(rt_obj$M1_orig) || !is.null(rt_obj$M2_orig) || !is.null(rt_obj$M3_orig)) {
    return(list(
      r_fix   = rt_obj$r_fix,
      p_fix   = rt_obj$p_fix,
      q_fix   = rt_obj$q_fix,
      M1_std  = rt_obj$M1_std,
      M2_std  = rt_obj$M2_std,
      M3_std  = rt_obj$M3_std,
      M1_orig = rt_obj$M1_orig,
      M2_orig = rt_obj$M2_orig,
      M3_orig = rt_obj$M3_orig
    ))
  }
  
  stop("No compatible rolling nowcast object found in DFM rolling RDS.")
}
