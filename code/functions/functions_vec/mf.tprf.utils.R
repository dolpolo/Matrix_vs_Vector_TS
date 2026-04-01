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

get_country_eval_windows <- function(country_inputs, params) {
  
  last_y_date <- max(country_inputs$dates_q, na.rm = TRUE)
  
  extra_m <- if (!is.null(params$extra_months_after_last_y)) {
    as.integer(params$extra_months_after_last_y)
  } else {
    1L
  }
  
  end_eval_rmsfe <- min(params$end_eval, last_y_date)
  end_eval_rt    <- min(params$end_eval, last_y_date %m+% lubridate::months(extra_m))
  
  list(
    last_y_date    = last_y_date,
    end_eval_rmsfe = end_eval_rmsfe,
    end_eval_rt    = end_eval_rt
  )
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
    "_Kmax-", params$Kmax,
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

