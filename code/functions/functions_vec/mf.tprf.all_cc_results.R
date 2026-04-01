# ==============================================================================
# MF-TPRF Cross-Country Utilities
# ------------------------------------------------------------------------------
# Corrected version:
#   1. full-sample estimation aligned to available quarterly target
#   2. pseudo real-time rolling stops country-by-country at last available month
#   3. RMSFE evaluation stops country-by-country at last available target quarter
#   4. corrected RT summary extraction: RT_res$all
#   5. full-sample and rolling hyperparameters propagated to country summaries
#   6. cross-country hyperparameter tables added to cross-country outputs
# ==============================================================================

# ==============================================================================
# 0. HELPERS
# ==============================================================================

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

rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

get_country_eval_windows <- function(country_inputs, params) {
  
  last_y_date <- max(country_inputs$dates_q, na.rm = TRUE)
  last_m_date <- max(country_inputs$dates_m, na.rm = TRUE)
  
  end_eval_rmsfe <- min(params$end_eval, last_y_date)
  end_eval_rt    <- min(params$end_eval, last_m_date)
  
  list(
    last_y_date    = last_y_date,
    last_m_date    = last_m_date,
    end_eval_rmsfe = end_eval_rmsfe,
    end_eval_rt    = end_eval_rt
  )
}

# ==============================================================================
# 1. COUNTRY INPUT PREPARATION
# ==============================================================================

prepare_country_inputs <- function(all_countries, country, params) {
  
  country_data <- all_countries$data[[country]]
  
  data    <- country_data$Data
  dates_m <- as.Date(country_data$Dates)
  series  <- country_data$Series
  
  gdp_col     <- country_data$target_col
  target_name <- paste0(params$target, "_", country)
  
  # ---------------------------------------------------------------------------
  # Target
  # ---------------------------------------------------------------------------
  y <- as.matrix(data[, gdp_col, drop = FALSE])
  
  y_obs_idx <- which(!is.na(y[, 1]))
  y_q       <- as.numeric(y[y_obs_idx, 1])
  dates_q   <- dates_m[y_obs_idx]
  
  # ---------------------------------------------------------------------------
  # Predictors
  # ---------------------------------------------------------------------------
  X <- as.matrix(data[, -gdp_col, drop = FALSE])
  
  N_m <- country_data$nM
  
  N_q_tot      <- country_data$nQ
  q_series_all <- series[(N_m + 1):(N_m + N_q_tot)]
  q_cols       <- which(toupper(q_series_all) != toupper(target_name))
  N_q          <- length(q_cols)
  N            <- N_m + N_q
  
  agg_m <- country_data$agg_m
  agg_q <- country_data$agg_q[q_cols]
  agg   <- c(agg_m, agg_q)
  
  freq_m <- country_data$freq_m
  freq_q <- country_data$freq_q[q_cols]
  freq   <- c(freq_m, freq_q)
  
  unb_m <- country_data$unb_m
  unb_q <- country_data$unb_q[q_cols]
  unb   <- c(unb_m, unb_q)
  
  class_m <- country_data$ClassM
  class_q <- country_data$ClassQ[q_cols]
  class   <- c(class_m, class_q)
  
  list(
    country      = country,
    target_name  = target_name,
    X            = X,
    y_q          = y_q,
    y_obs_idx    = y_obs_idx,
    dates_m      = dates_m,
    dates_q      = dates_q,
    series_names = series,
    N_m          = N_m,
    N_q          = N_q,
    N            = N,
    agg_m        = agg_m,
    agg_q        = agg_q,
    agg          = agg,
    freq         = freq,
    unb          = unb,
    class        = class
  )
}

prepare_country_inputs_from_tensor_vectorized <- function(tensor, data, country, params) {
  
  all_names   <- colnames(data)
  target_name <- paste0(country, "_", params$target)
  gdp_col     <- which(all_names == target_name)
  
  if (length(gdp_col) != 1L) {
    stop("Target column not found uniquely for country: ", country)
  }
  
  # ---------------------------------------------------------------------------
  # Target
  # ---------------------------------------------------------------------------
  y <- as.matrix(data[, gdp_col, drop = FALSE])
  
  dates_m   <- as.Date(dimnames(tensor$Y)[[1]])
  y_obs_idx <- which(!is.na(y[, 1]))
  y_q       <- as.numeric(y[y_obs_idx, 1])
  dates_q   <- dates_m[y_obs_idx]
  
  # ---------------------------------------------------------------------------
  # Predictors
  # ---------------------------------------------------------------------------
  X_vec <- as.matrix(data[, -gdp_col, drop = FALSE])
  predictor_names <- colnames(X_vec)
  
  get_varname <- function(x) sub("^[^_]+_", "", x)
  predictor_var <- get_varname(predictor_names)
  
  keep <- predictor_var != params$target
  X_vec           <- X_vec[, keep, drop = FALSE]
  predictor_names <- predictor_names[keep]
  predictor_var   <- predictor_var[keep]
  
  keep_nonempty <- colSums(!is.na(X_vec)) > 0
  X_vec           <- X_vec[, keep_nonempty, drop = FALSE]
  predictor_names <- predictor_names[keep_nonempty]
  predictor_var   <- predictor_var[keep_nonempty]
  
  vars_m <- tensor$vars[1:tensor$n_M]
  vars_q <- tensor$vars[(tensor$n_M + 1):(tensor$n_M + tensor$n_Q)]
  
  get_base_meta <- function(x) {
    if (is.null(dim(x))) return(as.vector(x))
    if (is.matrix(x) || length(dim(x)) == 2L) return(as.vector(x[1, ]))
    stop("Metadata object has unsupported dimensions.")
  }
  
  agg_m_base  <- get_base_meta(tensor$agg_M)
  agg_q_base  <- get_base_meta(tensor$agg_Q)
  freq_m_base <- get_base_meta(tensor$freq_M)
  freq_q_base <- get_base_meta(tensor$freq_Q)
  unb_m_base  <- get_base_meta(tensor$unb_M)
  unb_q_base  <- get_base_meta(tensor$unb_Q)
  
  agg_m_lookup  <- setNames(agg_m_base,  vars_m)
  agg_q_lookup  <- setNames(agg_q_base,  vars_q)
  freq_m_lookup <- setNames(freq_m_base, vars_m)
  freq_q_lookup <- setNames(freq_q_base, vars_q)
  unb_m_lookup  <- setNames(unb_m_base,  vars_m)
  unb_q_lookup  <- setNames(unb_q_base,  vars_q)
  
  freq <- ifelse(
    predictor_var %in% vars_m, "M",
    ifelse(predictor_var %in% vars_q, "Q", NA)
  )
  
  unb <- ifelse(
    predictor_var %in% vars_m,
    unb_m_lookup[predictor_var],
    unb_q_lookup[predictor_var]
  )
  
  agg_all <- ifelse(
    predictor_var %in% vars_m,
    agg_m_lookup[predictor_var],
    agg_q_lookup[predictor_var]
  )
  
  is_monthly   <- freq == "M"
  is_quarterly <- freq == "Q"
  
  N_m <- sum(is_monthly)
  N_q <- sum(is_quarterly)
  N   <- N_m + N_q
  
  agg_m <- matrix(agg_all[is_monthly], nrow = 1)
  agg_q <- matrix(agg_all[is_quarterly], nrow = 1)
  agg   <- c(agg_all[is_monthly], agg_all[is_quarterly])
  
  list(
    country      = country,
    target_name  = target_name,
    X            = X_vec,
    y_q          = y_q,
    y_obs_idx    = y_obs_idx,
    dates_m      = dates_m,
    dates_q      = dates_q,
    series_names = predictor_names,
    N_m          = N_m,
    N_q          = N_q,
    N            = N,
    agg_m        = agg_m,
    agg_q        = agg_q,
    agg          = agg,
    freq         = freq,
    unb          = unb
  )
}

# ==============================================================================
# 2. FULL-SAMPLE ESTIMATION
# ==============================================================================

estimate_country_mf_tprf <- function(country_inputs, params) {
  
  X   <- country_inputs$X
  y_q <- country_inputs$y_q
  N_m <- country_inputs$N_m
  N_q <- country_inputs$N_q
  
  out_std <- standardize_with_na(X)
  X_std   <- out_std$X_std
  
  imp_xp <- init_XP_ER(X_std, params$Kmax)
  X_xp   <- imp_xp$X_init
  r_hat  <- imp_xp$r
  
  X_m_xp <- X_xp[, 1:N_m, drop = FALSE]
  X_q_xp <- X_xp[, (N_m + 1):(N_m + N_q), drop = FALSE]
  
  X_mq_xp  <- agg_mq(X_m_xp, country_inputs$agg_m)
  X_qq_xp  <- agg_qq(X_q_xp, country_inputs$agg_q)
  X_xp_agg <- cbind(X_mq_xp, X_qq_xp)
  
  # align low-frequency regressors to observed quarterly target
  n_q_y <- length(y_q)
  
  if (n_q_y > nrow(X_xp_agg)) {
    stop(
      "Observed quarterly target length exceeds aggregated regressor rows for country ",
      country_inputs$country,
      ". length(y_q) = ", n_q_y,
      ", nrow(X_xp_agg) = ", nrow(X_xp_agg)
    )
  }
  
  X_xp_agg_est <- X_xp_agg[seq_len(n_q_y), , drop = FALSE]
  
  pls_object <- select_L_autoproxy_3prf(
    X_xp_agg_est,
    y_q,
    Zmax = params$Zmax
  )
  
  Lproxy <- pls_object$L_opt
  
  lag_sel <- choose_UMIDAS_lag(
    X_lf        = X_xp_agg_est,
    X_hf        = X_xp,
    y_q         = y_q,
    Lmax        = params$Lmax,
    Lproxy      = Lproxy,
    p_AR_max    = params$p_AR_max,
    Robust_F    = params$Robust_F,
    alpha       = params$alpha,
    robust_type = params$robust_type,
    nw_lag      = params$nw_lag
  )
  
  L_midas <- lag_sel$best_BIC$L
  p_ar    <- lag_sel$best_BIC$p_AR
  
  fit <- MF_TPRF(
    X_lf        = X_xp_agg_est,
    X_hf        = X_xp,
    y_q         = y_q,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_ar,
    Robust_F    = params$Robust_F,
    alpha       = params$alpha,
    robust_type = params$robust_type,
    nw_lag      = params$nw_lag
  )
  
  list(
    preprocessing = list(
      out_std      = out_std,
      imp_xp       = imp_xp,
      X_xp         = X_xp,
      X_mq_xp      = X_mq_xp,
      X_qq_xp      = X_qq_xp,
      X_xp_agg     = X_xp_agg,
      X_xp_agg_est = X_xp_agg_est
    ),
    hyper = list(
      r_hat   = r_hat,
      Lproxy  = Lproxy,
      L_midas = L_midas,
      p_ar    = p_ar,
      pls     = pls_object,
      lag_sel = lag_sel
    ),
    fit = fit
  )
}

# ==============================================================================
# 3. PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

run_country_mf_tprf_rt <- function(country_inputs, params) {
  
  rt_raw <- pseudo_realtime_MF_TPRF_XP(
    X_full  = country_inputs$X,
    y_q     = country_inputs$y_q,
    params  = params,
    dates   = country_inputs$dates_m,
    dates_q = country_inputs$dates_q,
    Freq    = country_inputs$freq,
    Unb     = country_inputs$unb,
    agg_m   = country_inputs$agg_m,
    agg_q   = country_inputs$agg_q
  )
  
  df_M1 <- list_to_df_nowcast(rt_raw$M1, "M1")
  df_M2 <- list_to_df_nowcast(rt_raw$M2, "M2")
  df_M3 <- list_to_df_nowcast(rt_raw$M3, "M3")
  
  df_rt <- dplyr::bind_rows(df_M1, df_M2, df_M3)
  if (nrow(df_rt) > 0L) {
    df_rt <- dplyr::arrange(df_rt, date, month_in_quarter)
  }
  
  list(
    raw = rt_raw,
    M1  = df_M1,
    M2  = df_M2,
    M3  = df_M3,
    all = df_rt
  )
}

# ==============================================================================
# 4. COUNTRY WRAPPER
# ==============================================================================

run_country_mf_tprf <- function(country_inputs, params, path_results = NULL) {
  
  win <- get_country_eval_windows(country_inputs, params)
  
  params_cc <- params
  params_cc$last_y_date    <- win$last_y_date
  params_cc$last_m_date    <- win$last_m_date
  params_cc$end_eval_rmsfe <- win$end_eval_rmsfe
  params_cc$end_eval_rt    <- win$end_eval_rt
  
  # RT loop should stop at country-specific last available monthly date
  params_rt <- params_cc
  params_rt$end_eval <- params_cc$end_eval_rt
  
  est <- estimate_country_mf_tprf(country_inputs, params_cc)
  rt  <- run_country_mf_tprf_rt(country_inputs, params_rt)
  
  out <- list(
    model_id = "MF_TPRF",
    country  = country_inputs$country,
    params   = params_cc,
    inputs   = country_inputs,
    full_sample     = est,
    pseudo_realtime = rt
  )
  
  if (!is.null(path_results)) {
    dir.create(file.path(path_results, country_inputs$country), recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      out,
      file.path(path_results, country_inputs$country,
                paste0("MF_TPRF_ALL_", country_inputs$country, ".rds"))
    )
  }
  
  out
}

# ==============================================================================
# 5. CROSS-COUNTRY WRAPPERS
# ==============================================================================

run_all_countries_mf_tprf <- function(countries, all_countries, params, path_results = NULL) {
  
  results <- vector("list", length(countries))
  names(results) <- countries
  
  for (cc in countries) {
    cat("\n==============================\n")
    cat("Running country:", cc, "\n")
    cat("==============================\n")
    
    country_inputs <- prepare_country_inputs(
      all_countries = all_countries,
      country       = cc,
      params        = params
    )
    
    results[[cc]] <- run_country_mf_tprf(
      country_inputs = country_inputs,
      params         = params,
      path_results   = path_results
    )
  }
  
  results
}

run_all_countries_mf_tprf_from_tensor <- function(countries, tensor, data, params, path_results = NULL) {
  
  countries_run <- if (!is.null(params$target_cc)) setdiff(countries, params$target_cc) else countries
  
  results <- vector("list", length(countries_run))
  names(results) <- countries_run
  
  for (cc in countries_run) {
    cat("\n==============================\n")
    cat("Running vectorized-tensor country:", cc, "\n")
    cat("==============================\n")
    
    country_inputs <- prepare_country_inputs_from_tensor_vectorized(
      tensor  = tensor,
      data    = data,
      country = cc,
      params  = params
    )
    
    results[[cc]] <- run_country_mf_tprf(
      country_inputs = country_inputs,
      params         = params,
      path_results   = path_results
    )
  }
  
  results
}

# ==============================================================================
# 6. COUNTRY SUMMARY OBJECTS
# ==============================================================================

summarize_mf_tprf_country <- function(
    MF_TPRF_res,
    RT_res = NULL,
    country,
    dates_m,
    dates_q,
    y_q,
    params
) {
  
  model_id <- "MF_TPRF"
  
  # ---------------------------------------------------------------------------
  # A. FULL-SAMPLE NOWCAST OBJECTS
  # ---------------------------------------------------------------------------
  fit_obj <- if (!is.null(MF_TPRF_res$fit)) {
    MF_TPRF_res$fit
  } else if (!is.null(MF_TPRF_res$MF_TPRF)) {
    MF_TPRF_res$MF_TPRF
  } else {
    stop("No valid full-sample fit object found in MF_TPRF_res.")
  }
  
  hyper_full <- if (!is.null(MF_TPRF_res$hyper)) MF_TPRF_res$hyper else NULL
  
  hyper_rt <- NULL
  if (!is.null(RT_res) && !is.null(RT_res$raw)) {
    hyper_rt <- list(
      pre  = RT_res$raw$hyper_pre,
      post = RT_res$raw$hyper_post
    )
  }
  
  y_now_full   <- fit_obj$y_nowcast
  T_m_full     <- length(y_now_full)
  T_q_complete <- length(y_q)
  
  if (T_m_full < 3 * T_q_complete) {
    stop("y_now_full has fewer than 3 * T_q_complete months.")
  }
  
  dates_m_full <- dates_m[1:T_m_full]
  
  n_in  <- 3 * T_q_complete
  n_tot <- T_m_full
  
  period_vec <- rep("in-sample", n_tot)
  if (n_tot > n_in) {
    period_vec[(n_in + 1):n_tot] <- "real-time"
  }
  
  df_now_full <- data.frame(
    date       = dates_m_full,
    y_now_full = y_now_full,
    period     = factor(period_vec, levels = c("in-sample", "real-time"))
  )
  
  df_quarterly <- data.frame(
    date   = dates_q,
    y_true = as.numeric(y_q)
  )
  
  plot_nowcast <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = params$covid_start, xmax = params$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = df_now_full,
      ggplot2::aes(x = date, y = y_now_full, color = period),
      linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_now_full, period == "real-time"),
      ggplot2::aes(x = date, y = y_now_full),
      size = 2.5, color = "black"
    ) +
    ggplot2::geom_line(
      data = df_quarterly,
      ggplot2::aes(x = date, y = y_true, color = "Quarterly GDP"),
      linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "in-sample"     = "#1F77B4",
        "real-time"     = "#2ECC71",
        "Quarterly GDP" = "#D62728"
      ),
      name = ""
    ) +
    ggplot2::labs(
      title = paste0("MF-TPRF Monthly Nowcast (In-Sample + Real-Time) - ", country),
      x = "Date",
      y = "GDP Growth"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  # ---------------------------------------------------------------------------
  # B. IN-SAMPLE RMSFE
  # ---------------------------------------------------------------------------
  start_eval  <- params$start_eval
  end_eval    <- if (!is.null(params$end_eval_rmsfe)) params$end_eval_rmsfe else params$end_eval
  covid_start <- params$covid_start
  covid_end   <- params$covid_end
  
  is_PRE   <- dates_q >= start_eval  & dates_q <  covid_start
  is_COVID <- dates_q >= covid_start & dates_q <= covid_end
  is_POST  <- dates_q >  covid_end   & dates_q <= end_eval
  is_ALL   <- dates_q >= start_eval  & dates_q <= end_eval
  
  y_now_in <- y_now_full[1:n_in]
  
  M1_idx <- seq(1, n_in, by = 3)
  M2_idx <- seq(2, n_in, by = 3)
  M3_idx <- seq(3, n_in, by = 3)
  
  rmsfe_mat_insample <- rbind(
    "Full sample" = c(
      rmsfe_period(is_ALL,   y_q, y_now_in, M1_idx),
      rmsfe_period(is_ALL,   y_q, y_now_in, M2_idx),
      rmsfe_period(is_ALL,   y_q, y_now_in, M3_idx)
    ),
    "Pre-COVID" = c(
      rmsfe_period(is_PRE,   y_q, y_now_in, M1_idx),
      rmsfe_period(is_PRE,   y_q, y_now_in, M2_idx),
      rmsfe_period(is_PRE,   y_q, y_now_in, M3_idx)
    ),
    "COVID period" = c(
      rmsfe_period(is_COVID, y_q, y_now_in, M1_idx),
      rmsfe_period(is_COVID, y_q, y_now_in, M2_idx),
      rmsfe_period(is_COVID, y_q, y_now_in, M3_idx)
    ),
    "Post-COVID" = c(
      rmsfe_period(is_POST,  y_q, y_now_in, M1_idx),
      rmsfe_period(is_POST,  y_q, y_now_in, M2_idx),
      rmsfe_period(is_POST,  y_q, y_now_in, M3_idx)
    )
  )
  colnames(rmsfe_mat_insample) <- c("M1", "M2", "M3")
  
  latex_insample <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    sprintf("\\caption{MF-TPRF RMSFE by period -- %s}\n", country),
    sprintf("\\label{tab:RMSFE_MF_TPRF_%s}\n", country),
    "\\begin{tabular}{lccc}\n",
    "\\toprule\n",
    "Period & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf("%s & %.4f & %.4f & %.4f \\\\",
              rownames(rmsfe_mat_insample),
              rmsfe_mat_insample[, 1],
              rmsfe_mat_insample[, 2],
              rmsfe_mat_insample[, 3]),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  # ---------------------------------------------------------------------------
  # C. ROLLING REAL-TIME
  # ---------------------------------------------------------------------------
  df_rt <- NULL
  df_yq_eval <- NULL
  plot_rt <- NULL
  rmsfe_mat_rt <- NULL
  latex_rt <- NULL
  
  if (!is.null(RT_res) && !is.null(RT_res$all)) {
    
    df_rt <- RT_res$all %>%
      dplyr::rename(
        type = month_in_quarter,
        nowcast = nowcast,
        date = date
      ) %>%
      dplyr::arrange(date)
    
    df_rt$type <- factor(df_rt$type, levels = c("M1", "M2", "M3"))
    
    end_eval_score <- if (!is.null(params$end_eval_rmsfe)) params$end_eval_rmsfe else params$end_eval
    
    idx_eval <- which(dates_q >= params$start_eval & dates_q <= end_eval_score)
    df_yq_eval <- data.frame(
      date = dates_q[idx_eval],
      GDP  = y_q[idx_eval]
    )
    
    plot_rt <- ggplot2::ggplot() +
      ggplot2::annotate(
        "rect",
        xmin = params$covid_start, xmax = params$covid_end,
        ymin = -Inf, ymax = Inf,
        fill = "grey80", alpha = 0.20
      ) +
      ggplot2::geom_line(
        data = df_yq_eval,
        ggplot2::aes(x = date, y = GDP, color = "True GDP"),
        linewidth = 1.2
      ) +
      ggplot2::geom_line(
        data = df_rt,
        ggplot2::aes(x = date, y = nowcast, color = type),
        linewidth = 1
      ) +
      ggplot2::geom_point(
        data = df_rt,
        ggplot2::aes(x = date, y = nowcast, color = type),
        size = 2
      ) +
      ggplot2::scale_color_manual(
        values = c(
          "True GDP" = "#D62728",
          "M1"       = "#1F77B4",
          "M2"       = "#2ECC71",
          "M3"       = "#F1C40F"
        ),
        name = "Series"
      ) +
      ggplot2::labs(
        title    = paste0("MF-TPRF Rolling Real-Time Nowcasts - ", country),
        subtitle = "M1: early • M2: mid-quarter • M3: end-quarter",
        x        = "Date",
        y        = "GDP Growth Rate"
      ) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(
        legend.position  = "bottom",
        plot.title       = ggplot2::element_text(face = "bold"),
        plot.subtitle    = ggplot2::element_text(color = "gray30"),
        panel.grid.minor = ggplot2::element_blank()
      )
    
    df_yq_eval <- df_yq_eval %>%
      dplyr::mutate(
        period = dplyr::case_when(
          date <  params$covid_start                            ~ "Pre-COVID",
          date >= params$covid_start & date <= params$covid_end ~ "COVID period",
          date >  params$covid_end                              ~ "Post-COVID",
          TRUE                                                  ~ NA_character_
        ),
        quarter_id = paste0(lubridate::year(date), "Q", lubridate::quarter(date))
      )
    
    df_rt_eval <- df_rt %>%
      dplyr::mutate(
        quarter_id = paste0(lubridate::year(date), "Q", lubridate::quarter(date))
      )
    
    df_eval <- df_rt_eval %>%
      dplyr::inner_join(
        dplyr::select(df_yq_eval, quarter_id, GDP, period),
        by = "quarter_id"
      ) %>%
      dplyr::arrange(date)
    
    rmsfe_by_period <- df_eval %>%
      dplyr::group_by(period, type) %>%
      dplyr::summarise(
        RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
        .groups = "drop"
      )
    
    rmsfe_full <- df_eval %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(
        RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::mutate(period = "Full sample")
    
    rmsfe_all <- dplyr::bind_rows(rmsfe_full, rmsfe_by_period) %>%
      dplyr::filter(period %in% c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")) %>%
      dplyr::mutate(
        period = factor(
          period,
          levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")
        ),
        type = factor(type, levels = c("M1", "M2", "M3"))
      ) %>%
      dplyr::arrange(period, type)
    
    rmsfe_mat_df <- tidyr::pivot_wider(rmsfe_all, names_from = type, values_from = RMSFE) %>%
      dplyr::arrange(period)
    
    rmsfe_mat_rt <- as.matrix(rmsfe_mat_df[, c("M1", "M2", "M3")])
    rownames(rmsfe_mat_rt) <- rmsfe_mat_df$period
    
    latex_rt <- paste0(
      "\\begin{table}[htbp]\n",
      "\\centering\n",
      sprintf("\\caption{MF-TPRF Rolling RMSFE by period -- %s}\n", country),
      sprintf("\\label{tab:RMSFE_MF_TPRF_Rolling_%s}\n", country),
      "\\begin{tabular}{lccc}\n",
      "\\toprule\n",
      "Period & M1 & M2 & M3 \\\\\n",
      "\\midrule\n",
      paste(
        sprintf("%s & %.4f & %.4f & %.4f \\\\",
                rownames(rmsfe_mat_rt),
                rmsfe_mat_rt[, "M1"],
                rmsfe_mat_rt[, "M2"],
                rmsfe_mat_rt[, "M3"]),
        collapse = "\n"
      ),
      "\n\\bottomrule\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    )
  }
  
  list(
    model_id        = model_id,
    country         = country,
    hyper_full      = hyper_full,
    hyper_rt        = hyper_rt,
    df_now_full     = df_now_full,
    df_quarterly    = df_quarterly,
    plot_nowcast    = plot_nowcast,
    rmsfe_insample  = rmsfe_mat_insample,
    latex_insample  = latex_insample,
    df_rt           = df_rt,
    df_yq_eval      = df_yq_eval,
    plot_rt         = plot_rt,
    rmsfe_rt        = rmsfe_mat_rt,
    latex_rt        = latex_rt
  )
}

# ==============================================================================
# 7. CROSS-COUNTRY OUTPUTS
# ==============================================================================

build_cross_country_outputs <- function(summary_all, params, model_label = "MF-TPRF") {
  
  df_now_full_all <- dplyr::bind_rows(
    lapply(names(summary_all), function(cc) {
      df <- summary_all[[cc]]$df_now_full
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$country <- cc
      df
    })
  )
  
  df_quarterly_all <- dplyr::bind_rows(
    lapply(names(summary_all), function(cc) {
      df <- summary_all[[cc]]$df_quarterly
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$country <- cc
      df
    })
  )
  
  country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
  
  if (!is.null(df_now_full_all) && nrow(df_now_full_all) > 0) {
    df_now_full_all$country <- factor(df_now_full_all$country, levels = country_order)
  }
  if (!is.null(df_quarterly_all) && nrow(df_quarterly_all) > 0) {
    df_quarterly_all$country <- factor(df_quarterly_all$country, levels = country_order)
  }
  
  plot_nowcast_facet <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = params$covid_start, xmax = params$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = df_now_full_all,
      ggplot2::aes(x = date, y = y_now_full, color = period),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = dplyr::filter(df_now_full_all, period == "real-time"),
      ggplot2::aes(x = date, y = y_now_full),
      size = 1.2, color = "black"
    ) +
    ggplot2::geom_line(
      data = df_quarterly_all,
      ggplot2::aes(x = date, y = y_true, color = "Quarterly GDP"),
      linewidth = 1
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "in-sample"     = "#1F77B4",
        "real-time"     = "#2ECC71",
        "Quarterly GDP" = "#D62728"
      ),
      name = ""
    ) +
    ggplot2::facet_wrap(~ country, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = paste0(model_label, " Monthly Nowcast (In-Sample + Real-Time)"),
      x = "Date",
      y = "GDP Growth"
    ) +
    ggplot2::theme_minimal(base_size = 13)
  
  tab_insample_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    mat <- summary_all[[cc]]$rmsfe_insample
    if (is.null(mat)) return(NULL)
    
    data.frame(
      country = cc,
      period  = "Full sample",
      M1 = mat["Full sample", "M1"],
      M2 = mat["Full sample", "M2"],
      M3 = mat["Full sample", "M3"],
      row.names = NULL
    )
  }))
  
  latex_tab_insample_all <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{MF-TPRF in-sample RMSFE across countries}\n",
    "\\label{tab:mf_tprf_insample_all}\n",
    "\\begin{tabular}{lccc}\n",
    "\\toprule\n",
    "Country & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf("%s & %.4f & %.4f & %.4f \\\\",
              tab_insample_all$country,
              tab_insample_all$M1,
              tab_insample_all$M2,
              tab_insample_all$M3),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  df_rt_all <- dplyr::bind_rows(
    lapply(names(summary_all), function(cc) {
      df <- summary_all[[cc]]$df_rt
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$country <- cc
      df
    })
  )
  
  df_yq_eval_all <- dplyr::bind_rows(
    lapply(names(summary_all), function(cc) {
      df <- summary_all[[cc]]$df_yq_eval
      if (is.null(df) || nrow(df) == 0) return(NULL)
      df$country <- cc
      df
    })
  )
  
  if (!is.null(df_rt_all) && nrow(df_rt_all) > 0) {
    df_rt_all$country <- factor(df_rt_all$country, levels = country_order)
  }
  if (!is.null(df_yq_eval_all) && nrow(df_yq_eval_all) > 0) {
    df_yq_eval_all$country <- factor(df_yq_eval_all$country, levels = country_order)
  }
  
  plot_rt_facet <- ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = params$covid_start, xmax = params$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey80", alpha = 0.20
    ) +
    ggplot2::geom_line(
      data = df_yq_eval_all,
      ggplot2::aes(x = date, y = GDP, color = "True GDP"),
      linewidth = 1
    ) +
    ggplot2::geom_line(
      data = df_rt_all,
      ggplot2::aes(x = date, y = nowcast, color = type),
      linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = df_rt_all,
      ggplot2::aes(x = date, y = nowcast, color = type),
      size = 1.2
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "True GDP" = "#D62728",
        "M1"       = "#1F77B4",
        "M2"       = "#2ECC71",
        "M3"       = "#F1C40F"
      ),
      name = "Series"
    ) +
    ggplot2::facet_wrap(~ country, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = paste0(model_label, " Rolling Real-Time Nowcasts"),
      subtitle = "M1: early • M2: mid-quarter • M3: end-quarter",
      x = "Date",
      y = "GDP Growth Rate"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  tab_rt_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    mat <- summary_all[[cc]]$rmsfe_rt
    if (is.null(mat)) return(NULL)
    
    data.frame(
      country = cc,
      period  = rownames(mat),
      M1 = mat[, "M1"],
      M2 = mat[, "M2"],
      M3 = mat[, "M3"],
      row.names = NULL
    )
  }))
  
  latex_tab_rt_all <- paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{MF-TPRF rolling RMSFE across countries and periods}\n",
    "\\label{tab:mf_tprf_rt_all}\n",
    "\\begin{tabular}{llccc}\n",
    "\\toprule\n",
    "Country & Period & M1 & M2 & M3 \\\\\n",
    "\\midrule\n",
    paste(
      sprintf("%s & %s & %.4f & %.4f & %.4f \\\\",
              tab_rt_all$country,
              tab_rt_all$period,
              tab_rt_all$M1,
              tab_rt_all$M2,
              tab_rt_all$M3),
      collapse = "\n"
    ),
    "\n\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  hyper_full_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    hp <- summary_all[[cc]]$hyper_full
    if (is.null(hp)) return(NULL)
    
    data.frame(
      country = cc,
      r_hat   = paste(hp$r_hat, collapse = ","),
      Lproxy  = hp$Lproxy,
      L_midas = hp$L_midas,
      p_ar    = hp$p_ar,
      row.names = NULL
    )
  }))
  
  hyper_rt_pre_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    hp_cc <- summary_all[[cc]]$hyper_rt
    if (is.null(hp_cc) || is.null(hp_cc$pre)) return(NULL)
    hp <- hp_cc$pre
    
    data.frame(
      country = cc,
      Lproxy  = hp$Lproxy,
      L_midas = hp$L_midas,
      p_AR    = hp$p_AR,
      r       = hp$r,
      row.names = NULL
    )
  }))
  
  hyper_rt_post_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    hp_cc <- summary_all[[cc]]$hyper_rt
    if (is.null(hp_cc) || is.null(hp_cc$post)) return(NULL)
    hp <- hp_cc$post
    
    data.frame(
      country   = cc,
      Lproxy    = hp$Lproxy,
      L_midas   = hp$L_midas,
      p_AR      = hp$p_AR,
      r         = hp$r,
      t_recalib = if (!is.null(hp$t_recalib)) as.character(hp$t_recalib) else NA_character_,
      row.names = NULL
    )
  }))
  
  list(
    hyper_full_all         = hyper_full_all,
    hyper_rt_pre_all       = hyper_rt_pre_all,
    hyper_rt_post_all      = hyper_rt_post_all,
    df_now_full_all        = df_now_full_all,
    df_quarterly_all       = df_quarterly_all,
    plot_nowcast_facet     = plot_nowcast_facet,
    tab_insample_all       = tab_insample_all,
    latex_tab_insample_all = latex_tab_insample_all,
    df_rt_all              = df_rt_all,
    df_yq_eval_all         = df_yq_eval_all,
    plot_rt_facet          = plot_rt_facet,
    tab_rt_all             = tab_rt_all,
    latex_tab_rt_all       = latex_tab_rt_all
  )
}

# ==============================================================================
# 8. MODEL COMPARISON HELPERS
# ==============================================================================

extract_country_rmsfe_long <- function(summary_all, model_name, which_rmsfe = c("rt", "insample")) {
  
  which_rmsfe <- match.arg(which_rmsfe)
  
  out <- lapply(names(summary_all), function(cc) {
    
    obj <- summary_all[[cc]]
    mat <- if (which_rmsfe == "rt") obj$rmsfe_rt else obj$rmsfe_insample
    if (is.null(mat)) return(NULL)
    
    data.frame(
      country = cc,
      period  = rownames(mat),
      M1      = mat[, "M1"],
      M2      = mat[, "M2"],
      M3      = mat[, "M3"],
      model   = model_name,
      row.names = NULL
    )
  })
  
  dplyr::bind_rows(out)
}

extract_matrix_rmsfe_long <- function(rmsfe_df, model_name = "Matrix", which_rmsfe = c("rt", "insample")) {
  
  which_rmsfe <- match.arg(which_rmsfe)
  
  df <- rmsfe_df
  df$model <- model_name
  
  df[, c("country", "period", "M1", "M2", "M3", "model")]
}

build_model_comparison_table <- function(
    matrix_rmsfe,
    vector_summary_all,
    vectensor_summary_all,
    which_rmsfe = c("rt", "insample")
) {
  
  which_rmsfe <- match.arg(which_rmsfe)
  
  df_matrix <- extract_matrix_rmsfe_long(
    rmsfe_df    = matrix_rmsfe,
    model_name  = "Matrix",
    which_rmsfe = which_rmsfe
  )
  
  df_vector <- extract_country_rmsfe_long(
    summary_all = vector_summary_all,
    model_name  = "Vector",
    which_rmsfe = which_rmsfe
  )
  
  df_vectensor <- extract_country_rmsfe_long(
    summary_all = vectensor_summary_all,
    model_name  = "VecTensor",
    which_rmsfe = which_rmsfe
  )
  
  df_all <- dplyr::bind_rows(df_matrix, df_vector, df_vectensor)
  
  df_all %>%
    tidyr::pivot_longer(
      cols      = c("M1", "M2", "M3"),
      names_to  = "horizon",
      values_to = "RMSFE"
    ) %>%
    dplyr::mutate(col_id = paste0(model, "_", horizon)) %>%
    dplyr::select(country, period, col_id, RMSFE) %>%
    tidyr::pivot_wider(names_from = col_id, values_from = RMSFE) %>%
    dplyr::arrange(
      country,
      factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID"))
    )
}

comparison_table_to_latex <- function(df_comp, caption, label) {
  
  cols_needed <- c(
    "Matrix_M1", "Matrix_M2", "Matrix_M3",
    "Vector_M1", "Vector_M2", "Vector_M3",
    "VecTensor_M1", "VecTensor_M2", "VecTensor_M3"
  )
  
  for (cc in cols_needed) {
    if (!cc %in% names(df_comp)) df_comp[[cc]] <- NA_real_
  }
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\scriptsize\n",
    sprintf("\\caption{%s}\n", caption),
    sprintf("\\label{%s}\n", label),
    "\\begin{tabular}{llccc ccc ccc}\n",
    "\\toprule\n",
    " & & \\multicolumn{3}{c}{Matrix MF-TPRF} & \\multicolumn{3}{c}{Vector MF-TPRF} & \\multicolumn{3}{c}{Vectorized-from-Tensor MF-TPRF} \\\\\n",
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