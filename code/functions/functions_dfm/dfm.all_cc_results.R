# ==============================================================================
# DFM Cross-Country Utilities
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
  
  list(
    last_y_date    = last_y_date,
    last_m_date    = last_m_date,
    end_eval_rmsfe = min(params$end_eval, last_y_date),
    end_eval_rt    = min(params$end_eval, last_m_date)
  )
}

destd_mat <- function(X_std, mu, sd) {
  X_std <- as.matrix(X_std)
  stopifnot(length(mu) == ncol(X_std), length(sd) == ncol(X_std))
  sweep(sweep(X_std, 2, sd, `*`), 2, mu, `+`)
}

get_idio_spec <- function(params) {
  if (is.null(params$idio_spec)) return("auto")
  params$idio_spec
}

is_q0_spec <- function(params) {
  identical(get_idio_spec(params), "q0")
}

select_or_fix_q <- function(Ehat, NM, params) {
  
  if (is_q0_spec(params)) {
    return(list(
      q          = 0L,
      idio_ar_ic = NULL,
      q_source   = "fixed_q0"
    ))
  }
  
  idio_ar_ic <- select_p_ar_ic_global(
    Ehat[, 1:NM, drop = FALSE],
    pmax          = params$qmax,
    include_const = FALSE,
    min_T         = 30
  )
  
  list(
    q          = idio_ar_ic$p_BIC,
    idio_ar_ic = idio_ar_ic,
    q_source   = "bic_auto"
  )
}

# ==============================================================================
# 1. COUNTRY INPUT PREPARATION
# ==============================================================================

prepare_country_inputs_dfm <- function(all_countries, country, params) {
  
  country_data <- all_countries$data[[country]]
  
  data    <- as.matrix(country_data$Data)
  dates_m <- as.Date(country_data$Dates)
  series  <- country_data$Series
  
  gdp_col     <- country_data$target_col
  target_name <- country_data$target_name
  
  y         <- as.matrix(data[, gdp_col, drop = FALSE])
  y_obs_idx <- which(!is.na(y[, 1]))
  y_q       <- as.numeric(y[y_obs_idx, 1])
  dates_q   <- dates_m[y_obs_idx]
  
  NM <- country_data$nM
  NQ <- country_data$nQ
  
  list(
    country      = country,
    target_name  = target_name,
    X_full       = data,
    y_q          = y_q,
    y_obs_idx    = y_obs_idx,
    gdp_col      = gdp_col,
    dates_m      = dates_m,
    dates_q      = dates_q,
    series_names = series,
    N_m          = NM,
    N_q          = NQ,
    N            = NM + NQ,
    Type         = c(country_data$TypeM,  country_data$TypeQ),
    agg_m        = country_data$agg_m,
    agg_q        = country_data$agg_q,
    agg          = c(country_data$agg_m, country_data$agg_q),
    freq         = c(country_data$freq_m, country_data$freq_q),
    unb          = c(country_data$unb_m,  country_data$unb_q),
    class        = c(country_data$ClassM, country_data$ClassQ)
  )
}

# ==============================================================================
# 2. FULL-SAMPLE ESTIMATION
# ==============================================================================

estimate_country_dfm <- function(country_inputs, params) {
  
  X_full  <- country_inputs$X_full
  gdp_col <- country_inputs$gdp_col
  NM      <- country_inputs$N_m
  NQ      <- country_inputs$N_q
  agg     <- country_inputs$agg
  
  out_std <- standardize_with_na(X_full)
  X_std   <- out_std$X_std
  
  mu_full <- out_std$mean
  sd_full <- out_std$sd
  
  cov_proxy_out <- all_purpose_covariance(X_std)
  
  ER_proxy_out <- select_num_factors_ER(
    Sigma = cov_proxy_out$Sigma_tilde,
    Kmax  = params$Kmax
  )
  
  r   <- ER_proxy_out$r
  eig <- ER_proxy_out$eig
  
  xp <- estimate_factors_XP(X_std, r = r)
  
  L_hat <- xp$Lambda
  F_hat <- xp$F_hat
  C_hat <- xp$C_hat
  
  factor_var_ic <- select_p_var_ic(
    F_hat         = F_hat,
    pmax          = params$pmax,
    include_const = TRUE
  )
  
  p <- factor_var_ic$p_BIC
  
  Ehat <- X_std - C_hat
  
  q_out      <- select_or_fix_q(Ehat = Ehat, NM = NM, params = params)
  q          <- q_out$q
  idio_ar_ic <- q_out$idio_ar_ic
  q_source   <- q_out$q_source
  
  Init <- InitialCond(
    X_std    = X_std,
    r        = r,
    p_factor = p,
    q_idio   = q,
    NM       = NM,
    NQ       = NQ,
    restr    = params$restr,
    agg      = agg,
    kappa    = params$kappa
  )
  
  fit_dfm <- DFM_EM(
    X        = X_std,
    Init     = Init,
    max_iter = params$max_iter,
    tol      = params$tol
  )
  
  n_lags_Q_used <- if (params$restr == "MM") 5L else 3L
  
  fit_obj <- dfm_kalman_fit(
    X_in     = X_std,
    res      = fit_dfm,
    n_lags_Q = n_lags_Q_used
  )
  
  X_fit_std  <- fit_obj$X_fit
  X_fit_orig <- destd_mat(X_fit_std, mu_full, sd_full)
  
  list(
    preprocessing = list(
      X_raw         = X_full,
      out_std       = out_std,
      X_std         = X_std,
      cov_proxy_out = cov_proxy_out,
      xp            = xp,
      Ehat          = Ehat
    ),
    hyper = list(
      r             = r,
      p             = p,
      q             = q,
      q_source      = q_source,
      idio_spec     = get_idio_spec(params),
      eig           = eig,
      ER_proxy_out  = ER_proxy_out,
      factor_var_ic = factor_var_ic,
      idio_ar_ic    = idio_ar_ic
    ),
    fit = list(
      Init          = Init,
      res           = fit_dfm,
      A             = fit_dfm$A,
      C             = fit_dfm$C,
      Q             = fit_dfm$Q,
      R             = fit_dfm$R,
      Z0            = fit_dfm$Z0,
      V0            = fit_dfm$V0,
      loglik        = fit_dfm$loglik,
      crit          = fit_dfm$crit,
      L_hat         = L_hat,
      F_hat         = F_hat,
      C_hat         = C_hat,
      std_map       = list(mean = mu_full, sd = sd_full),
      fit_obj       = fit_obj,
      n_lags_Q_used = n_lags_Q_used,
      X_fit_std     = X_fit_std,
      X_fit_orig    = X_fit_orig,
      gdp_fit_std   = X_fit_std[, gdp_col],
      gdp_fit_orig  = X_fit_orig[, gdp_col]
    )
  )
}

# ==============================================================================
# 3. PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

run_country_dfm_rt <- function(country_inputs, params) {
  
  user_hyper_pre <- if (is_q0_spec(params)) {
    list(r = NULL, p = NULL, q = 0L)
  } else {
    list(r = NULL, p = NULL, q = NULL)
  }
  
  user_hyper_post <- if (is_q0_spec(params)) {
    list(r = NULL, p = NULL, q = 0L)
  } else {
    list(r = NULL, p = NULL, q = NULL)
  }
  
  rt_raw <- pseudo_realtime_DFM_EM_reestimate(
    X_full   = country_inputs$X_full,
    NQ       = country_inputs$N_q,
    params   = params,
    dates_m  = country_inputs$dates_m,
    dates_q  = country_inputs$dates_q,
    Freq     = country_inputs$freq,
    Unb      = country_inputs$unb,
    gdp_col  = country_inputs$gdp_col,
    agg      = country_inputs$agg,
    do_post_covid_recalibration = TRUE,
    user_hyper_pre  = user_hyper_pre,
    user_hyper_post = user_hyper_post,
    max_iter_em = params$max_iter,
    tol_em      = params$tol,
    pmax        = params$pmax,
    qmax        = params$qmax,
    min_est_T   = 24,
    verbose     = TRUE
  )
  
  df_M1 <- list_to_df_nowcast(rt_raw$M1_orig, "M1")
  df_M2 <- list_to_df_nowcast(rt_raw$M2_orig, "M2")
  df_M3 <- list_to_df_nowcast(rt_raw$M3_orig, "M3")
  
  df_rt <- dplyr::bind_rows(df_M1, df_M2, df_M3)
  if (nrow(df_rt) > 0L) {
    df_rt <- dplyr::arrange(df_rt, date, month_in_quarter)
  }
  
  rolling_to_df <- function(now_obj) {
    make_df <- function(x, label, scale) {
      if (length(x) == 0) return(NULL)
      data.frame(
        date    = as.Date(names(x)),
        value   = as.numeric(x),
        vintage = label,
        scale   = scale,
        row.names = NULL
      )
    }
    
    dplyr::bind_rows(
      make_df(now_obj$M1_std,  "M1", "std"),
      make_df(now_obj$M2_std,  "M2", "std"),
      make_df(now_obj$M3_std,  "M3", "std"),
      make_df(now_obj$M1_orig, "M1", "orig"),
      make_df(now_obj$M2_orig, "M2", "orig"),
      make_df(now_obj$M3_orig, "M3", "orig")
    )
  }
  
  list(
    raw  = rt_raw,
    M1   = df_M1,
    M2   = df_M2,
    M3   = df_M3,
    all  = df_rt,
    long = rolling_to_df(rt_raw)
  )
}

# ==============================================================================
# 4. COUNTRY WRAPPER
# ==============================================================================

run_country_dfm <- function(country_inputs, params, path_results = NULL) {
  
  win <- get_country_eval_windows(country_inputs, params)
  
  params_cc <- params
  params_cc$last_y_date    <- win$last_y_date
  params_cc$last_m_date    <- win$last_m_date
  params_cc$end_eval_rmsfe <- win$end_eval_rmsfe
  params_cc$end_eval_rt    <- win$end_eval_rt
  
  params_rt <- params_cc
  params_rt$end_eval <- params_cc$end_eval_rt
  
  est <- estimate_country_dfm(country_inputs, params_cc)
  rt  <- run_country_dfm_rt(country_inputs, params_rt)
  
  out <- list(
    model_id        = "DFM_EM",
    country         = country_inputs$country,
    params          = params_cc,
    inputs          = country_inputs,
    full_sample     = est,
    pseudo_realtime = rt
  )
  
  if (!is.null(path_results)) {
    dir.create(file.path(path_results, country_inputs$country), recursive = TRUE, showWarnings = FALSE)
    
    file_country <- file.path(
      path_results,
      country_inputs$country,
      paste0(
        "DFM_ALL_",
        country_inputs$country,
        "_Idio-",
        get_idio_spec(params_cc),
        ".rds"
      )
    )
    
    saveRDS(out, file_country)
  }
  
  out
}

# ==============================================================================
# 5. CROSS-COUNTRY WRAPPERS
# ==============================================================================

run_all_countries_dfm <- function(countries, all_countries, params, path_results = NULL) {
  
  results <- vector("list", length(countries))
  names(results) <- countries
  
  for (cc in countries) {
    cat("\n==============================\n")
    cat("Running DFM country:", cc, "\n")
    cat("Idiosyncratic specification:", get_idio_spec(params), "\n")
    cat("==============================\n")
    
    country_inputs <- prepare_country_inputs_dfm(
      all_countries = all_countries,
      country       = cc,
      params        = params
    )
    
    results[[cc]] <- run_country_dfm(
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
# Key difference vs MF-TPRF:
#   full-sample "nowcast" is replaced by DFM fitted/smoothed GDP
# ==============================================================================

summarize_dfm_country <- function(
    DFM_res,
    RT_res = NULL,
    country,
    dates_m,
    dates_q,
    y_q,
    params
) {
  
  model_id <- "DFM_EM"
  
  fit_obj   <- DFM_res$fit
  hyper_full <- if (!is.null(DFM_res$hyper)) DFM_res$hyper else NULL
  
  hyper_rt <- NULL
  if (!is.null(RT_res) && !is.null(RT_res$raw)) {
    hyper_rt <- list(
      pre  = RT_res$raw$hyper_pre,
      post = RT_res$raw$hyper_post
    )
  }
  
  # ---------------------------------------------------------------------------
  # A. FULL-SAMPLE FITTED GDP OBJECTS
  # ---------------------------------------------------------------------------
  gdp_fit_full <- fit_obj$gdp_fit_orig
  T_m_full     <- length(gdp_fit_full)
  T_q_complete <- length(y_q)
  
  if (T_m_full < 3 * T_q_complete) {
    stop("gdp_fit_full has fewer than 3 * T_q_complete months.")
  }
  
  dates_m_full <- dates_m[seq_len(T_m_full)]
  
  n_in  <- 3 * T_q_complete
  n_tot <- T_m_full
  
  period_vec <- rep("in-sample", n_tot)
  if (n_tot > n_in) {
    period_vec[(n_in + 1):n_tot] <- "real-time"
  }
  
  df_now_full <- data.frame(
    date       = dates_m_full,
    y_now_full = as.numeric(gdp_fit_full),
    period     = factor(period_vec, levels = c("in-sample", "real-time"))
  )
  
  df_quarterly <- data.frame(
    date   = dates_q,
    y_true = as.numeric(y_q)
  )
  
  year_breaks_full <- seq(
    lubridate::floor_date(min(df_now_full$date, na.rm = TRUE), unit = "year"),
    lubridate::floor_date(max(df_now_full$date, na.rm = TRUE), unit = "year"),
    by = "1 year"
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
    ggplot2::scale_x_date(
      breaks = year_breaks_full,
      date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::labs(
      title = "Monthly fitted GDP and observed GDP",
      subtitle = NULL,
      x = "Date",
      y = "GDP growth"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey25")
    )
  
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
  
  y_now_in <- gdp_fit_full[1:n_in]
  
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
    sprintf("\\caption{DFM RMSFE by period -- %s}\n", country),
    sprintf("\\label{tab:RMSFE_DFM_%s}\n", country),
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
    
    year_breaks_rt <- seq(
      lubridate::floor_date(min(df_rt$date, na.rm = TRUE), unit = "year"),
      lubridate::floor_date(max(df_rt$date, na.rm = TRUE), unit = "year"),
      by = "1 year"
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
      ggplot2::scale_x_date(
        breaks = year_breaks_rt,
        date_labels = "%Y",
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::labs(
        title    = paste("Rolling real-time nowcasts:", country),
        subtitle = NULL,
        x        = "Date",
        y        = "GDP growth"
      ) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::theme(
        legend.position  = "bottom",
        plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(color = "gray30", hjust = 0.5),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
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
      sprintf("\\caption{DFM Rolling RMSFE by period -- %s}\n", country),
      sprintf("\\label{tab:RMSFE_DFM_Rolling_%s}\n", country),
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

build_cross_country_outputs_dfm <- function(summary_all, params, model_label = "DFM") {
  
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
  
  country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  
  if (!is.null(df_now_full_all) && nrow(df_now_full_all) > 0) {
    df_now_full_all$country <- factor(df_now_full_all$country, levels = country_order)
  }
  if (!is.null(df_quarterly_all) && nrow(df_quarterly_all) > 0) {
    df_quarterly_all$country <- factor(df_quarterly_all$country, levels = country_order)
  }
  
  year_breaks_full_all <- seq(
    lubridate::floor_date(min(df_now_full_all$date, na.rm = TRUE), unit = "year"),
    lubridate::floor_date(max(df_now_full_all$date, na.rm = TRUE), unit = "year"),
    by = "1 year"
  )
  
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
    ggplot2::scale_x_date(
      breaks = year_breaks_full_all,
      date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::facet_wrap(~ country, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = "Monthly fitted GDP and observed GDP",
      subtitle = NULL,
      x = "Date",
      y = "GDP growth"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey25")
    )
  
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
    "\\caption{DFM in-sample RMSFE across countries}\n",
    "\\label{tab:dfm_insample_all}\n",
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
  
  year_breaks_rt_all <- seq(
    lubridate::floor_date(min(df_rt_all$date, na.rm = TRUE), unit = "year"),
    lubridate::floor_date(max(df_rt_all$date, na.rm = TRUE), unit = "year"),
    by = "1 year"
  )
  
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
    ggplot2::scale_x_date(
      breaks = year_breaks_rt_all,
      date_labels = "%Y",
      expand = ggplot2::expansion(mult = c(0.01, 0.02))
    ) +
    ggplot2::facet_wrap(~ country, ncol = 2, scales = "free_y") +
    ggplot2::labs(
      title = "Rolling real-time nowcasts",
      subtitle = "Quarterly targets and monthly nowcast updates",
      x = "Date",
      y = "GDP growth"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position  = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25")
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
    "\\caption{DFM rolling RMSFE across countries and periods}\n",
    "\\label{tab:dfm_rt_all}\n",
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
      r       = hp$r,
      p       = hp$p,
      q       = hp$q,
      row.names = NULL
    )
  }))
  
  hyper_rt_pre_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    hp_cc <- summary_all[[cc]]$hyper_rt
    if (is.null(hp_cc) || is.null(hp_cc$pre)) return(NULL)
    hp <- hp_cc$pre
    
    data.frame(
      country = cc,
      r       = hp$r,
      p       = hp$p,
      q       = hp$q,
      row.names = NULL
    )
  }))
  
  hyper_rt_post_all <- do.call(rbind, lapply(names(summary_all), function(cc) {
    hp_cc <- summary_all[[cc]]$hyper_rt
    if (is.null(hp_cc) || is.null(hp_cc$post)) return(NULL)
    hp <- hp_cc$post
    
    data.frame(
      country   = cc,
      r         = hp$r,
      p         = hp$p,
      q         = hp$q,
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
# 7B. REAL-TIME PLOT VARIANTS
# ==============================================================================

build_rt_plot_variants <- function(
    df_rt_all,
    df_yq_eval_all,
    params,
    model_label = "MF-TPRF",
    country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
) {
  
  df_rt_all <- df_rt_all %>%
    dplyr::mutate(country = factor(as.character(country), levels = country_order))
  
  df_yq_eval_all <- df_yq_eval_all %>%
    dplyr::mutate(country = factor(as.character(country), levels = country_order))
  
  make_rt_plot <- function(df_rt, df_yq, title, subtitle = NULL, ncol = 2, shade_covid = TRUE) {
    
    year_breaks <- seq(
      lubridate::floor_date(min(df_rt$date, na.rm = TRUE), unit = "year"),
      lubridate::floor_date(max(df_rt$date, na.rm = TRUE), unit = "year"),
      by = "1 year"
    )
    
    p <- ggplot2::ggplot()
    
    if (isTRUE(shade_covid)) {
      p <- p +
        ggplot2::annotate(
          "rect",
          xmin = params$covid_start, xmax = params$covid_end,
          ymin = -Inf, ymax = Inf,
          fill = "grey80", alpha = 0.20
        )
    }
    
    p +
      ggplot2::geom_line(
        data = df_yq,
        ggplot2::aes(x = date, y = GDP, color = "True GDP"),
        linewidth = 1.0
      ) +
      ggplot2::geom_line(
        data = df_rt,
        ggplot2::aes(x = date, y = nowcast, color = type),
        linewidth = 0.8
      ) +
      ggplot2::geom_point(
        data = df_rt,
        ggplot2::aes(x = date, y = nowcast, color = type),
        size = 1.4
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
      ggplot2::scale_x_date(
        breaks = year_breaks,
        date_labels = "%Y",
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::facet_wrap(~ country, ncol = ncol, scales = "free_y") +
      ggplot2::labs(
        title = title,
        subtitle = subtitle,
        x = "Date",
        y = "GDP growth"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position  = "bottom",
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25"),
        strip.text       = ggplot2::element_text(face = "bold")
      )
  }
  
  countries_big4   <- c("DE", "FR", "IT", "ES")
  countries_other4 <- c("NL", "BE", "AT", "PT")
  countries_post8  <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
  
  plot_all <- make_rt_plot(
    df_rt     = df_rt_all,
    df_yq     = df_yq_eval_all,
    title     = "Rolling real-time nowcasts",
    subtitle  = "Quarterly targets and monthly nowcast updates"
  )
  
  plot_big4 <- make_rt_plot(
    df_rt    = df_rt_all %>% dplyr::filter(country %in% countries_big4),
    df_yq    = df_yq_eval_all %>% dplyr::filter(country %in% countries_big4),
    title    = "Rolling real-time nowcasts: DE, FR, IT, ES",
    subtitle = NULL
  )
  
  plot_other4 <- make_rt_plot(
    df_rt    = df_rt_all %>% dplyr::filter(country %in% countries_other4),
    df_yq    = df_yq_eval_all %>% dplyr::filter(country %in% countries_other4),
    title    = "Rolling real-time nowcasts: NL, BE, AT, PT",
    subtitle = NULL
  )
  
  df_rt_post8 <- df_rt_all %>%
    dplyr::filter(country %in% countries_post8, date > params$covid_end)
  
  df_yq_post8 <- df_yq_eval_all %>%
    dplyr::filter(country %in% countries_post8, date > params$covid_end)
  
  plot_post8 <- make_rt_plot(
    df_rt    = df_rt_post8,
    df_yq    = df_yq_post8,
    title    = "Post-COVID rolling real-time nowcasts",
    subtitle = paste0(
      "Sample: ",
      format(min(df_rt_post8$date, na.rm = TRUE), "%b %Y"),
      "–",
      format(max(df_rt_post8$date, na.rm = TRUE), "%b %Y")
    ),
    shade_covid = FALSE
  )
  
  countries_available <- country_order[country_order %in% unique(as.character(df_rt_all$country))]
  
  plot_by_country <- lapply(countries_available, function(cc) {
    
    df_rt_cc <- df_rt_all %>% dplyr::filter(country == cc)
    df_yq_cc <- df_yq_eval_all %>% dplyr::filter(country == cc)
    
    year_breaks_cc <- seq(
      lubridate::floor_date(min(df_rt_cc$date, na.rm = TRUE), unit = "year"),
      lubridate::floor_date(max(df_rt_cc$date, na.rm = TRUE), unit = "year"),
      by = "1 year"
    )
    
    ggplot2::ggplot() +
      ggplot2::annotate(
        "rect",
        xmin = params$covid_start, xmax = params$covid_end,
        ymin = -Inf, ymax = Inf,
        fill = "grey80", alpha = 0.20
      ) +
      ggplot2::geom_line(
        data = df_yq_cc,
        ggplot2::aes(x = date, y = GDP, color = "True GDP"),
        linewidth = 1.0
      ) +
      ggplot2::geom_line(
        data = df_rt_cc,
        ggplot2::aes(x = date, y = nowcast, color = type),
        linewidth = 0.9
      ) +
      ggplot2::geom_point(
        data = df_rt_cc,
        ggplot2::aes(x = date, y = nowcast, color = type),
        size = 1.5
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
      ggplot2::scale_x_date(
        breaks = year_breaks_cc,
        date_labels = "%Y",
        expand = ggplot2::expansion(mult = c(0.01, 0.02))
      ) +
      ggplot2::labs(
        title = paste("Rolling real-time nowcasts:", cc),
        subtitle = NULL,
        x = "Date",
        y = "GDP growth"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position  = "bottom",
        axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = ggplot2::element_text(hjust = 0.5, colour = "grey25")
      )
  })
  
  names(plot_by_country) <- countries_available
  
  list(
    plot_all        = plot_all,
    plot_big4       = plot_big4,
    plot_other4     = plot_other4,
    plot_post8      = plot_post8,
    plot_by_country = plot_by_country
  )
}

# ==============================================================================
# 7C. SAVE REAL-TIME PLOT VARIANTS
# ==============================================================================

save_rt_plot_variants <- function(
    rt_plots,
    path_graph_rt,
    model_name,
    Size = NULL,
    sel = NULL,
    width_full = 12,
    height_full = 8,
    width_country = 9,
    height_country = 5,
    dpi = 300
) {
  
  path_graph_rt_full   <- file.path(path_graph_rt, "full_oos")
  path_graph_rt_split  <- file.path(path_graph_rt, "split_groups")
  path_graph_rt_post   <- file.path(path_graph_rt, "post_covid")
  path_graph_rt_bycc   <- file.path(path_graph_rt, "by_country")
  
  dir.create(path_graph_rt_full,  recursive = TRUE, showWarnings = FALSE)
  dir.create(path_graph_rt_split, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_graph_rt_post,  recursive = TRUE, showWarnings = FALSE)
  dir.create(path_graph_rt_bycc,  recursive = TRUE, showWarnings = FALSE)
  
  tag_size <- if (!is.null(Size)) paste0("_Size-", Size) else ""
  tag_sel  <- if (!is.null(sel))  paste0("_sel-", sel)  else ""
  
  file_graph_rt_all <- file.path(
    path_graph_rt_full,
    paste0("plot_rt_full_oos_", model_name, tag_size, tag_sel, ".png")
  )
  
  file_graph_rt_big4 <- file.path(
    path_graph_rt_split,
    paste0("plot_rt_big4_", model_name, tag_size, tag_sel, ".png")
  )
  
  file_graph_rt_other4 <- file.path(
    path_graph_rt_split,
    paste0("plot_rt_other4_", model_name, tag_size, tag_sel, ".png")
  )
  
  file_graph_rt_post8 <- file.path(
    path_graph_rt_post,
    paste0("plot_rt_postcovid_8countries_", model_name, tag_size, tag_sel, ".png")
  )
  
  ggplot2::ggsave(file_graph_rt_all,    rt_plots$plot_all,    width = width_full, height = height_full, dpi = dpi)
  ggplot2::ggsave(file_graph_rt_big4,   rt_plots$plot_big4,   width = width_full, height = height_full, dpi = dpi)
  ggplot2::ggsave(file_graph_rt_other4, rt_plots$plot_other4, width = width_full, height = height_full, dpi = dpi)
  ggplot2::ggsave(file_graph_rt_post8,  rt_plots$plot_post8,  width = width_full, height = height_full, dpi = dpi)
  
  file_graph_rt_by_country <- list()
  
  if (!is.null(rt_plots$plot_by_country) && length(rt_plots$plot_by_country) > 0L) {
    for (cc in names(rt_plots$plot_by_country)) {
      file_cc <- file.path(
        path_graph_rt_bycc,
        paste0("plot_rt_", cc, "_", model_name, tag_size, tag_sel, ".png")
      )
      ggplot2::ggsave(
        filename = file_cc,
        plot     = rt_plots$plot_by_country[[cc]],
        width    = width_country,
        height   = height_country,
        dpi      = dpi
      )
      file_graph_rt_by_country[[cc]] <- file_cc
    }
  }
  
  list(
    path_graph_rt_full       = path_graph_rt_full,
    path_graph_rt_split      = path_graph_rt_split,
    path_graph_rt_post       = path_graph_rt_post,
    path_graph_rt_bycc       = path_graph_rt_bycc,
    file_graph_rt_all        = file_graph_rt_all,
    file_graph_rt_big4       = file_graph_rt_big4,
    file_graph_rt_other4     = file_graph_rt_other4,
    file_graph_rt_post8      = file_graph_rt_post8,
    file_graph_rt_by_country = file_graph_rt_by_country
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