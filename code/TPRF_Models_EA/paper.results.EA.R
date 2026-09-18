# ==============================================================================
# EA NOWCAST COMPARISON
# Matrix MF-TPRF vs VEC-C
# Small information set: Corr and LASSO
# ==============================================================================

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ==============================================================================
# 1. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"

path_final_results <- file.path(
  path_main,
  "TPRF_Models_EA/Final_Tab_Graph"
)

path_ea_nowcast_plots <- file.path(
  path_final_results,
  "EA_nowcast_comparison"
)

path_matrix_results <- file.path(
  path_main,
  "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs"
)

path_vector_results_ea <- file.path(
  path_main,
  "TPRF_Models_EA/Vector_MF-TPRF/results/outputs/EA"
)

dir.create(
  path_ea_nowcast_plots,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================================
# 2. FUNCTIONS
# ==============================================================================

path_func <- file.path(
  path_main,
  "functions/functions_mat"
)

source(
  file.path(
    path_func,
    "matrix.mf.tprf.utils.R"
  )
)

# ==============================================================================
# 3. CONFIGURATION
# ==============================================================================

proxy_mode <- "scalar"   # "scalar" | "multivariate"
size_ea    <- "small"

ea_selection_grid <- c(
  "corr",
  "LASSO"
)

ea_month_order <- c(
  "M1",
  "M2",
  "M3"
)

# 0.05 = 5 percentage points.
deviation_threshold <- 0.05

# Margins around the realised-EA-GDP range used for the y-axis.
gdp_axis_pad_frac <- 0.15
gdp_axis_min_pad  <- 0.05

matrix_model_name <- paste0(
  "matrix_",
  proxy_mode
)

# ==============================================================================
# 4. LOAD, ALIGN, AND EVALUATE EA NOWCASTS
# ==============================================================================

ea_nowcast_comparison <- list()
ea_outlier_reports    <- list()
ea_eval_list          <- list()
ea_coverage_list      <- list()

params_plot <- NULL

for (sel_i in ea_selection_grid) {
  
  file_matrix_rt_ea <- find_result_file(
    path  = path_matrix_results,
    model = matrix_model_name,
    stage = "rt",
    Size  = size_ea,
    sel   = sel_i
  )
  
  file_vector_rt_ea <- find_result_file(
    path  = path_vector_results_ea,
    model = "vector",
    stage = "rt",
    Size  = size_ea,
    sel   = sel_i
  )
  
  matrix_rt_ea <- readRDS(file_matrix_rt_ea)
  vector_rt_ea <- readRDS(file_vector_rt_ea)
  
  if (is.null(matrix_rt_ea$proxy_mode)) {
    stop("The Matrix EA RT object does not contain `proxy_mode`.")
  }
  
  if (!identical(
    as.character(matrix_rt_ea$proxy_mode),
    proxy_mode
  )) {
    stop(
      "Matrix proxy mode mismatch for selection method ",
      sel_i,
      ". Requested: ",
      proxy_mode,
      " | Loaded: ",
      matrix_rt_ea$proxy_mode
    )
  }
  
  check_common_ea_settings(
    matrix_rt = matrix_rt_ea,
    vector_rt = vector_rt_ea,
    sel_i     = sel_i
  )
  
  if (is.null(params_plot)) {
    params_plot <- matrix_rt_ea$params
  }
  
  df_matrix_ea <- extract_matrix_ea_nowcasts(
    rt_obj      = matrix_rt_ea,
    month_order = ea_month_order
  )
  
  df_vector_ea <- extract_vector_ea_nowcasts(
    rt_obj      = vector_rt_ea,
    month_order = ea_month_order
  )
  
  check_unique_nowcast_keys(
    df_matrix_ea,
    "Matrix EA nowcasts"
  )
  
  check_unique_nowcast_keys(
    df_vector_ea,
    "Vector EA nowcasts"
  )
  
  common_keys <- df_matrix_ea %>%
    distinct(date, type) %>%
    inner_join(
      df_vector_ea %>%
        distinct(date, type),
      by = c("date", "type")
    ) %>%
    arrange(
      date,
      match(type, ea_month_order)
    )
  
  if (nrow(common_keys) == 0L) {
    stop(
      "No common Matrix--VEC-C EA date-vintage pairs found for ",
      sel_i,
      "."
    )
  }
  
  common_quarters <- common_keys %>%
    mutate(
      quarter_id = make_ea_quarter_id(date)
    ) %>%
    distinct(quarter_id)
  
  df_nowcasts <- bind_rows(
    df_matrix_ea,
    df_vector_ea
  ) %>%
    inner_join(
      common_keys,
      by = c("date", "type")
    ) %>%
    mutate(
      type = factor(
        type,
        levels = ea_month_order
      ),
      model = factor(
        model,
        levels = c(
          "Matrix MF-TPRF",
          "VEC-C"
        )
      )
    ) %>%
    arrange(
      type,
      model,
      date
    )
  
  df_gdp_matrix <- extract_matrix_ea_gdp(
    matrix_rt_ea
  )
  
  df_gdp_vector <- extract_vector_ea_gdp(
    vector_rt_ea
  )
  
  df_gdp_check <- df_gdp_matrix %>%
    select(
      quarter_id,
      GDP_matrix = GDP
    ) %>%
    inner_join(
      df_gdp_vector %>%
        select(
          quarter_id,
          GDP_vector = GDP
        ),
      by = "quarter_id"
    )
  
  if (nrow(df_gdp_check) == 0L) {
    stop(
      "No common realised EA GDP observations found for ",
      sel_i,
      "."
    )
  }
  
  gdp_difference <- max(
    abs(
      df_gdp_check$GDP_matrix -
        df_gdp_check$GDP_vector
    ),
    na.rm = TRUE
  )
  
  if (
    is.finite(gdp_difference) &&
    gdp_difference > 1e-10
  ) {
    warning(
      "Matrix and Vector EA GDP series differ for ",
      sel_i,
      ". The graph and evaluation use the Matrix GDP series."
    )
  }
  
  df_gdp_common <- df_gdp_matrix %>%
    semi_join(
      common_quarters,
      by = "quarter_id"
    ) %>%
    arrange(date)
  
  df_eval_i <- build_ea_nowcast_evaluation(
    df_nowcasts = df_nowcasts,
    df_gdp      = df_gdp_common,
    selection   = sel_i
  )
  
  df_coverage_i <- df_eval_i %>%
    distinct(
      date,
      type,
      model
    ) %>%
    count(
      model,
      type,
      name = "N"
    ) %>%
    arrange(
      model,
      match(type, ea_month_order)
    )
  
  cat(
    "\n",
    paste(rep("=", 88), collapse = ""),
    "\nEA EVALUATION COVERAGE | ",
    sel_i,
    "\n",
    paste(rep("=", 88), collapse = ""),
    "\n",
    sep = ""
  )
  
  print_ea_console_table(
    df_coverage_i
  )
  
  selection_label <- if (identical(sel_i, "corr")) {
    "Correlation screening"
  } else {
    "LASSO screening"
  }
  
  ea_outlier_reports[[sel_i]] <- report_ea_nowcast_deviations(
    df_eval         = df_eval_i,
    threshold       = deviation_threshold,
    selection_label = selection_label
  )
  
  df_gdp_facet <- df_eval_i %>%
    distinct(
      date,
      type,
      GDP
    ) %>%
    transmute(
      date,
      type = factor(
        type,
        levels = ea_month_order
      ),
      value  = GDP,
      series = "Observed EA GDP"
    )
  
  df_plot <- bind_rows(
    df_gdp_facet,
    df_nowcasts %>%
      transmute(
        date,
        type,
        value  = nowcast,
        series = as.character(model)
      )
  ) %>%
    mutate(
      series = factor(
        series,
        levels = c(
          "Observed EA GDP",
          "Matrix MF-TPRF",
          "VEC-C"
        )
      )
    ) %>%
    arrange(
      type,
      series,
      date
    )
  
  ea_nowcast_comparison[[sel_i]] <- list(
    selection       = sel_i,
    file_matrix_rt  = file_matrix_rt_ea,
    file_vector_rt  = file_vector_rt_ea,
    common_keys     = common_keys,
    common_quarters = common_quarters,
    nowcasts        = df_nowcasts,
    gdp             = df_gdp_common,
    evaluation      = df_eval_i,
    coverage        = df_coverage_i,
    outliers        = ea_outlier_reports[[sel_i]],
    plot_data       = df_plot
  )
  
  ea_eval_list[[sel_i]]     <- df_eval_i
  ea_coverage_list[[sel_i]] <- df_coverage_i
}

# ==============================================================================
# 5. FIGURES
# Corr and LASSO are separate figures.
# Each figure contains M1, M2, and M3.
# The y-axis is based only on realised EA GDP.
# Nowcasts outside the axis range are squished at the graph border.
# ==============================================================================

plot_ea_matrix_vector <- list()
file_ea_matrix_vector <- list()
ylim_ea_by_selection  <- list()

for (sel_i in ea_selection_grid) {
  
  df_plot_i <- ea_nowcast_comparison[[sel_i]]$plot_data
  
  ylim_i <- get_ea_gdp_ylim(
    df_gdp   = ea_nowcast_comparison[[sel_i]]$gdp,
    pad_frac = gdp_axis_pad_frac,
    min_pad  = gdp_axis_min_pad
  )
  
  ylim_ea_by_selection[[sel_i]] <- ylim_i
  
  cat(
    "\nEA graph limits (", sel_i, "): [",
    sprintf("%.4f", ylim_i[1]),
    ", ",
    sprintf("%.4f", ylim_i[2]),
    "] based only on realised EA GDP.\n",
    sep = ""
  )
  
  p_ea <- ggplot(
    df_plot_i,
    aes(
      x        = date,
      y        = value,
      colour   = series,
      linetype = series,
      group    = interaction(type, series)
    )
  ) +
    annotate(
      "rect",
      xmin  = params_plot$covid_start,
      xmax  = params_plot$covid_end,
      ymin  = -Inf,
      ymax  = Inf,
      fill  = "grey72",
      alpha = 0.11
    ) +
    geom_hline(
      yintercept = 0,
      linewidth  = 0.35,
      colour     = "grey55"
    ) +
    geom_line(
      linewidth = 1.00,
      alpha     = 0.98,
      lineend   = "round"
    ) +
    facet_wrap(
      ~ type,
      ncol = 3
    ) +
    scale_colour_manual(
      values = c(
        "Observed EA GDP" = "#C00000",
        "Matrix MF-TPRF"  = "#228B22",
        "VEC-C"           = "#1F77B4"
      ),
      breaks = c(
        "Observed EA GDP",
        "Matrix MF-TPRF",
        "VEC-C"
      ),
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(
        "Observed EA GDP" = "solid",
        "Matrix MF-TPRF"  = "solid",
        "VEC-C"           = "33"
      ),
      breaks = c(
        "Observed EA GDP",
        "Matrix MF-TPRF",
        "VEC-C"
      ),
      name = NULL
    ) +
    scale_x_date(
      breaks = seq(
        floor_date(
          min(df_plot_i$date, na.rm = TRUE),
          unit = "year"
        ),
        floor_date(
          max(df_plot_i$date, na.rm = TRUE),
          unit = "year"
        ),
        by = "2 years"
      ),
      date_labels = "%Y",
      expand = expansion(
        mult = c(0.01, 0.02)
      )
    ) +
    scale_y_continuous(
      limits = ylim_i,
      oob    = scales::squish,
      breaks = pretty(
        ylim_i,
        n = 6
      ),
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    coord_cartesian(
      expand = FALSE,
      clip   = "on"
    ) +
    labs(
      title = "Expanding pseudo-real-time EA nowcasts",
      x = NULL,
      y = "EA GDP growth"
    ) +
    theme_country_compare(
      base_size = 13
    )
  
  print(p_ea)
  
  file_plot <- file.path(
    path_ea_nowcast_plots,
    paste0(
      "plot_EA_Matrix_vs_Vector_proxy-",
      proxy_mode,
      "_Size-small_sel-",
      sel_i,
      ".png"
    )
  )
  
  ggsave(
    filename = file_plot,
    plot     = p_ea,
    width    = 12.5,
    height   = 4.8,
    dpi      = 500,
    bg       = "white"
  )
  
  plot_ea_matrix_vector[[sel_i]] <- p_ea
  file_ea_matrix_vector[[sel_i]] <- file_plot
  
  ea_nowcast_comparison[[sel_i]]$plot      <- p_ea
  ea_nowcast_comparison[[sel_i]]$file_plot <- file_plot
  ea_nowcast_comparison[[sel_i]]$ylim      <- ylim_i
}

# ==============================================================================
# 6. RMSFE TABLE
# ==============================================================================

ea_eval_all <- bind_rows(
  ea_eval_list
)

ea_eval_periods <- add_ea_evaluation_periods(
  df_eval             = ea_eval_all,
  params              = params_plot,
  include_full_sample = TRUE
)

ea_rmsfe_long <- compute_ea_rmsfe_by_period(
  ea_eval_periods
)

ea_rmsfe_table <- build_ea_rmsfe_table(
  df_rmsfe = ea_rmsfe_long,
  selection_order = c(
    "LASSO",
    "corr"
  ),
  month_order = ea_month_order,
  digits = 3
)

latex_ea_rmsfe <- make_latex_ea_rmsfe_table(
  df = ea_rmsfe_table,
  caption = paste0(
    "Pseudo-real-time EA GDP RMSFE: Matrix MF--TPRF versus VEC-C ",
    "(small information set, ",
    proxy_mode,
    " proxy)"
  ),
  label = paste0(
    "tab:ea_rmsfe_matrix_vector_small_",
    proxy_mode
  ),
  digits = 3
)

cat(
  "\n\n",
  paste(rep("=", 88), collapse = ""),
  "\nEA RMSFE TABLE\n",
  paste(rep("=", 88), collapse = ""),
  "\n",
  sep = ""
)

print_ea_console_table(
  ea_rmsfe_table
)

cat(
  "\n\n",
  paste(rep("=", 88), collapse = ""),
  "\nEA RMSFE LATEX\n",
  paste(rep("=", 88), collapse = ""),
  "\n",
  latex_ea_rmsfe,
  "\n",
  paste(rep("=", 88), collapse = ""),
  "\n",
  sep = ""
)

writeLines(
  latex_ea_rmsfe,
  con = file.path(
    path_ea_nowcast_plots,
    paste0(
      "table_EA_RMSFE_Matrix_vs_Vector_proxy-",
      proxy_mode,
      "_Size-small.tex"
    )
  )
)

# ==============================================================================
# 7. SAVE OUTPUT
# ==============================================================================

saveRDS(
  list(
    proxy_mode          = proxy_mode,
    size                = size_ea,
    params              = params_plot,
    deviation_threshold = deviation_threshold,
    comparison          = ea_nowcast_comparison,
    outlier_reports     = ea_outlier_reports,
    coverage            = ea_coverage_list,
    evaluation          = ea_eval_all,
    evaluation_periods  = ea_eval_periods,
    rmsfe_long          = ea_rmsfe_long,
    rmsfe_table         = ea_rmsfe_table,
    latex_rmsfe         = latex_ea_rmsfe,
    ylim_by_selection   = ylim_ea_by_selection,
    plots               = plot_ea_matrix_vector,
    files               = file_ea_matrix_vector
  ),
  file.path(
    path_ea_nowcast_plots,
    paste0(
      "EA_nowcast_comparison_Matrix_vs_Vector_proxy-",
      proxy_mode,
      "_Size-small.rds"
    )
  )
)

cat(
  "\nEA Matrix-vs-Vector graphs and EA RMSFE-DM output saved to:\n",
  path_ea_nowcast_plots,
  "\n"
)

# ==============================================================================
# 8. IN-SAMPLE FACTOR INTERPRETATION
# QoQ vs YoY EA GDP growth
# ==============================================================================

# IMPORTANT:
# The in-sample factor interpretation exercise reported in the slides
# uses the SMALL + CORRELATION specification.

sel_factor <- "LASSO"


# ------------------------------------------------------------------------------
# 8.1 LOAD FULL-SAMPLE FIT OBJECTS
# ------------------------------------------------------------------------------

file_matrix_fit_ea <- find_result_file(
  path  = path_matrix_results,
  model = matrix_model_name,
  stage = "fit",
  Size  = size_ea,
  sel   = sel_factor
)

file_vector_fit_ea <- find_result_file(
  path  = path_vector_results_ea,
  model = "vector",
  stage = "fit",
  Size  = size_ea,
  sel   = sel_factor
)

matrix_fit_ea <- readRDS(file_matrix_fit_ea)
vector_fit_ea <- readRDS(file_vector_fit_ea)

cat(
  "\nLoaded Matrix fit:\n",
  file_matrix_fit_ea,
  "\n"
)

cat(
  "\nLoaded Vector fit:\n",
  file_vector_fit_ea,
  "\n"
)


# ------------------------------------------------------------------------------
# 8.2 EA GDP: QoQ AND YoY LOG GROWTH
# ------------------------------------------------------------------------------

# Use the GDP target stored in the SAME Matrix full-sample fit object.
# This avoids relying on `matrix_rt_ea`, which is overwritten inside the
# previous Corr/LASSO real-time loop.

dates_q <- as.Date(matrix_fit_ea$dates_q)

Y_q_all_fit <- matrix_fit_ea$target$y_q_all

if (!is.null(colnames(Y_q_all_fit)) &&
    "EA" %in% colnames(Y_q_all_fit)) {
  
  y_q_qoq <- as.numeric(
    Y_q_all_fit[, "EA"]
  )
  
} else {
  
  gdp_col_fit <- matrix_fit_ea$target$gdp_col
  
  y_q_qoq <- as.numeric(
    Y_q_all_fit[, gdp_col_fit]
  )
}

if (length(y_q_qoq) != length(dates_q)) {
  stop(
    "GDP target and quarterly-date vectors have different lengths."
  )
}

# Current target:
#
# GDP_QoQ_t = log(GDP_t) - log(GDP_{t-1})
#
# Therefore:
#
# GDP_YoY_t = log(GDP_t) - log(GDP_{t-4})
#           = GDP_QoQ_t
#             + GDP_QoQ_{t-1}
#             + GDP_QoQ_{t-2}
#             + GDP_QoQ_{t-3}

y_q_yoy <- y_q_qoq +
  dplyr::lag(y_q_qoq, 1) +
  dplyr::lag(y_q_qoq, 2) +
  dplyr::lag(y_q_qoq, 3)

df_gdp_growth <- data.frame(
  date    = dates_q,
  GDP_QoQ = y_q_qoq,
  GDP_YoY = y_q_yoy
)


# ------------------------------------------------------------------------------
# 8.3 BUILD MATRIX FACTOR DATA
# ------------------------------------------------------------------------------

df_factor_matrix <- data.frame(
  date = as.Date(
    matrix_fit_ea$fit$factors$dates_q
  ),
  M1 = as.numeric(
    matrix_fit_ea$fit$factors$F1[, 1]
  ),
  M2 = as.numeric(
    matrix_fit_ea$fit$factors$F2[, 1]
  ),
  M3 = as.numeric(
    matrix_fit_ea$fit$factors$F3[, 1]
  )
) %>%
  dplyr::mutate(
    F_quarter = (M1 + M2 + M3) / 3
  )


# ------------------------------------------------------------------------------
# 8.4 BUILD VECTOR FACTOR DATA
# ------------------------------------------------------------------------------

df_factor_vector <- data.frame(
  date = as.Date(
    vector_fit_ea$dates_q
  ),
  M1 = as.numeric(
    vector_fit_ea$fit$F1[, 1]
  ),
  M2 = as.numeric(
    vector_fit_ea$fit$F2[, 1]
  ),
  M3 = as.numeric(
    vector_fit_ea$fit$F3[, 1]
  )
) %>%
  dplyr::mutate(
    F_quarter = (M1 + M2 + M3) / 3
  )


# ------------------------------------------------------------------------------
# 8.5 CHECK FACTOR COVERAGE
# ------------------------------------------------------------------------------

cat(
  "\nMatrix factor sample:",
  format(
    min(df_factor_matrix$date),
    "%Y-%m-%d"
  ),
  "to",
  format(
    max(df_factor_matrix$date),
    "%Y-%m-%d"
  ),
  "| N =",
  nrow(df_factor_matrix),
  "\n"
)

cat(
  "Vector factor sample:",
  format(
    min(df_factor_vector$date),
    "%Y-%m-%d"
  ),
  "to",
  format(
    max(df_factor_vector$date),
    "%Y-%m-%d"
  ),
  "| N =",
  nrow(df_factor_vector),
  "\n"
)


# ------------------------------------------------------------------------------
# 8.6 MERGE GDP AND FACTORS
# ------------------------------------------------------------------------------

df_matrix_fit <- df_gdp_growth %>%
  dplyr::inner_join(
    df_factor_matrix,
    by = "date"
  ) %>%
  dplyr::arrange(date)

df_vector_fit <- df_gdp_growth %>%
  dplyr::inner_join(
    df_factor_vector,
    by = "date"
  ) %>%
  dplyr::arrange(date)


# ------------------------------------------------------------------------------
# 8.7 FUNCTION: IN-SAMPLE FACTOR FIT
# ------------------------------------------------------------------------------

estimate_factor_fit <- function(
    data,
    y_var,
    model_name
) {
  
  df <- data %>%
    dplyr::select(
      y = dplyr::all_of(y_var),
      F_quarter,
      M1,
      M2,
      M3
    ) %>%
    tidyr::drop_na()
  
  # ------------------------------------------------------------
  # Specification 1:
  # GDP growth on quarterly-average factor
  # ------------------------------------------------------------
  
  fit_q <- stats::lm(
    y ~ F_quarter,
    data = df
  )
  
  # ------------------------------------------------------------
  # Specification 2:
  # GDP growth on M1, M2, M3 separately
  # ------------------------------------------------------------
  
  fit_m <- stats::lm(
    y ~ M1 + M2 + M3,
    data = df
  )
  
  # With one regressor + intercept:
  #
  # R2 = Corr(y, F_quarter)^2
  
  corr_q <- stats::cor(
    df$y,
    df$F_quarter
  )
  
  # With M1/M2/M3 jointly, sqrt(R2) is the multiple correlation:
  #
  # Corr(y, fitted(y))
  
  multiple_corr_m <- stats::cor(
    df$y,
    stats::fitted(fit_m)
  )
  
  dplyr::bind_rows(
    
    data.frame(
      dependent       = y_var,
      model           = model_name,
      specification   = "Quarterly average",
      N               = stats::nobs(fit_q),
      R2              = summary(fit_q)$r.squared,
      Adj_R2          = summary(fit_q)$adj.r.squared,
      Correlation     = corr_q,
      Multiple_Corr   = abs(corr_q)
    ),
    
    data.frame(
      dependent       = y_var,
      model           = model_name,
      specification   = "M1 + M2 + M3",
      N               = stats::nobs(fit_m),
      R2              = summary(fit_m)$r.squared,
      Adj_R2          = summary(fit_m)$adj.r.squared,
      Correlation     = NA_real_,
      Multiple_Corr   = multiple_corr_m
    )
  )
}


# ==============================================================================
# 8.8 REPLICATION CHECK:
# ORIGINAL QoQ EXERCISE
# ==============================================================================

fit_qoq_original <- dplyr::bind_rows(
  
  estimate_factor_fit(
    data       = df_matrix_fit,
    y_var      = "GDP_QoQ",
    model_name = "Matrix MF-TPRF"
  ),
  
  estimate_factor_fit(
    data       = df_vector_fit,
    y_var      = "GDP_QoQ",
    model_name = "Vector MF-TPRF"
  )
)


cat(
  "\n",
  paste(rep("=", 90), collapse = ""),
  "\nORIGINAL IN-SAMPLE EXERCISE: QoQ GDP GROWTH | CORR SCREENING\n",
  paste(rep("=", 90), collapse = ""),
  "\n",
  sep = ""
)

print(
  fit_qoq_original %>%
    dplyr::mutate(
      R2_pct        = round(100 * R2, 1),
      Adj_R2_pct    = round(100 * Adj_R2, 1),
      Correlation   = round(Correlation, 3),
      Multiple_Corr = round(Multiple_Corr, 3)
    ) %>%
    dplyr::select(
      model,
      specification,
      N,
      R2_pct,
      Adj_R2_pct,
      Correlation,
      Multiple_Corr
    )
)


# ==============================================================================
# 8.9 BUILD ONE IDENTICAL COMMON SAMPLE
# Matrix vs Vector AND QoQ vs YoY
# ==============================================================================

df_common_factor <- df_gdp_growth %>%
  dplyr::inner_join(
    df_factor_matrix %>%
      dplyr::rename(
        Matrix_M1 = M1,
        Matrix_M2 = M2,
        Matrix_M3 = M3,
        Matrix_F  = F_quarter
      ),
    by = "date"
  ) %>%
  dplyr::inner_join(
    df_factor_vector %>%
      dplyr::rename(
        Vector_M1 = M1,
        Vector_M2 = M2,
        Vector_M3 = M3,
        Vector_F  = F_quarter
      ),
    by = "date"
  ) %>%
  dplyr::filter(
    stats::complete.cases(
      GDP_QoQ,
      GDP_YoY,
      Matrix_M1,
      Matrix_M2,
      Matrix_M3,
      Matrix_F,
      Vector_M1,
      Vector_M2,
      Vector_M3,
      Vector_F
    )
  ) %>%
  dplyr::arrange(date)


cat(
  "\nCommon QoQ--YoY / Matrix--Vector sample:",
  format(
    min(df_common_factor$date),
    "%Y-%m-%d"
  ),
  "to",
  format(
    max(df_common_factor$date),
    "%Y-%m-%d"
  ),
  "| N =",
  nrow(df_common_factor),
  "\n"
)


# ------------------------------------------------------------------------------
# 8.10 MODEL-SPECIFIC DATA ON THE COMMON SAMPLE
# ------------------------------------------------------------------------------

df_matrix_common <- df_common_factor %>%
  dplyr::transmute(
    date,
    GDP_QoQ,
    GDP_YoY,
    M1        = Matrix_M1,
    M2        = Matrix_M2,
    M3        = Matrix_M3,
    F_quarter = Matrix_F
  )

df_vector_common <- df_common_factor %>%
  dplyr::transmute(
    date,
    GDP_QoQ,
    GDP_YoY,
    M1        = Vector_M1,
    M2        = Vector_M2,
    M3        = Vector_M3,
    F_quarter = Vector_F
  )


# ==============================================================================
# 8.11 QoQ vs YoY ON EXACTLY THE SAME SAMPLE
# ==============================================================================

factor_fit_comparison <- dplyr::bind_rows(
  
  # QoQ -- Matrix
  estimate_factor_fit(
    data       = df_matrix_common,
    y_var      = "GDP_QoQ",
    model_name = "Matrix MF-TPRF"
  ),
  
  # QoQ -- Vector
  estimate_factor_fit(
    data       = df_vector_common,
    y_var      = "GDP_QoQ",
    model_name = "Vector MF-TPRF"
  ),
  
  # YoY -- Matrix
  estimate_factor_fit(
    data       = df_matrix_common,
    y_var      = "GDP_YoY",
    model_name = "Matrix MF-TPRF"
  ),
  
  # YoY -- Vector
  estimate_factor_fit(
    data       = df_vector_common,
    y_var      = "GDP_YoY",
    model_name = "Vector MF-TPRF"
  )
)


# ------------------------------------------------------------------------------
# 8.12 FINAL TABLE
# ------------------------------------------------------------------------------

factor_fit_table <- factor_fit_comparison %>%
  dplyr::mutate(
    
    # Explicit namespace avoids car::recode vs dplyr::recode conflict
    Growth = dplyr::recode(
      dependent,
      GDP_QoQ = "QoQ",
      GDP_YoY = "YoY"
    ),
    
    R2_pct        = 100 * R2,
    Adj_R2_pct    = 100 * Adj_R2,
    Abs_Corr      = abs(Correlation)
    
  ) %>%
  dplyr::select(
    Growth,
    model,
    specification,
    N,
    R2_pct,
    Adj_R2_pct,
    Correlation,
    Abs_Corr,
    Multiple_Corr
  )


cat(
  "\n",
  paste(rep("=", 90), collapse = ""),
  "\nIN-SAMPLE FACTOR INTERPRETATION: QoQ vs YoY GDP GROWTH\n",
  paste(rep("=", 90), collapse = ""),
  "\n",
  sep = ""
)

print(
  factor_fit_table %>%
    dplyr::mutate(
      R2_pct        = round(R2_pct, 1),
      Adj_R2_pct    = round(Adj_R2_pct, 1),
      Correlation   = round(Correlation, 3),
      Abs_Corr      = round(Abs_Corr, 3),
      Multiple_Corr = round(Multiple_Corr, 3)
    )
)


# ------------------------------------------------------------------------------
# 8.13 COMPACT R2 TABLE
# ------------------------------------------------------------------------------

factor_r2_compact <- factor_fit_table %>%
  dplyr::select(
    Growth,
    model,
    specification,
    R2_pct
  ) %>%
  dplyr::mutate(
    R2_pct = round(R2_pct, 1)
  ) %>%
  tidyr::pivot_wider(
    names_from  = model,
    values_from = R2_pct
  ) %>%
  dplyr::arrange(
    factor(
      Growth,
      levels = c("QoQ", "YoY")
    ),
    factor(
      specification,
      levels = c(
        "Quarterly average",
        "M1 + M2 + M3"
      )
    )
  )


cat(
  "\n",
  paste(rep("-", 70), collapse = ""),
  "\nCOMPACT R2 COMPARISON (%)\n",
  paste(rep("-", 70), collapse = ""),
  "\n",
  sep = ""
)

print(
  factor_r2_compact
)


# ------------------------------------------------------------------------------
# 8.14 SAVE RESULTS
# ------------------------------------------------------------------------------

utils::write.csv(
  factor_fit_table,
  file.path(
    path_ea_nowcast_plots,
    paste0(
      "EA_factor_interpretation_QoQ_vs_YoY_sel-",
      sel_factor,
      ".csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  list(
    selection         = sel_factor,
    gdp_growth        = df_gdp_growth,
    matrix_factors    = df_factor_matrix,
    vector_factors    = df_factor_vector,
    original_qoq      = fit_qoq_original,
    common_sample     = df_common_factor,
    fit_comparison    = factor_fit_comparison,
    summary_table     = factor_fit_table,
    compact_r2_table  = factor_r2_compact
  ),
  file.path(
    path_ea_nowcast_plots,
    paste0(
      "EA_factor_interpretation_QoQ_vs_YoY_sel-",
      sel_factor,
      ".rds"
    )
  )
)

# ==============================================================================
# 8.15 FIGURE: R2 COMPARISON
# QoQ vs YoY | Matrix vs Vector
# ==============================================================================

df_r2_plot <- factor_fit_table %>%
  dplyr::mutate(
    Growth = factor(
      Growth,
      levels = c("QoQ", "YoY")
    ),
    specification = factor(
      specification,
      levels = c(
        "Quarterly average",
        "M1 + M2 + M3"
      )
    ),
    model = factor(
      model,
      levels = c(
        "Matrix MF-TPRF",
        "Vector MF-TPRF"
      )
    )
  )

p_factor_r2 <- ggplot(
  df_r2_plot,
  aes(
    x = Growth,
    y = R2_pct,
    fill = model
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.64
  ) +
  geom_text(
    aes(
      label = sprintf("%.1f", R2_pct)
    ),
    position = position_dodge(width = 0.72),
    vjust = -0.35,
    size = 3.6
  ) +
  facet_wrap(
    ~ specification,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Matrix MF-TPRF" = "#228B22",
      "Vector MF-TPRF" = "#1F77B4"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(df_r2_plot$R2_pct, na.rm = TRUE) * 1.12
    ),
    labels = function(x) paste0(x, "%"),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    title = "In-sample EA GDP fit",
    subtitle = paste0(
      "Small information set — ",
      ifelse(
        sel_factor == "corr",
        "Correlation screening",
        "LASSO screening"
      )
    ),
    x = NULL,
    y = expression(R^2)
  ) +
  theme_country_compare(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom"
  )

print(p_factor_r2)

file_factor_r2 <- file.path(
  path_ea_nowcast_plots,
  paste0(
    "plot_EA_factor_R2_QoQ_vs_YoY_sel-",
    sel_factor,
    ".png"
  )
)

ggsave(
  filename = file_factor_r2,
  plot = p_factor_r2,
  width = 9.2,
  height = 4.8,
  dpi = 500,
  bg = "white"
)

# ==============================================================================
# 8.16 FIGURE: YoY GDP AND ESTIMATED FACTORS
# ==============================================================================

sign_matrix <- sign(
  cor(
    df_common_factor$GDP_YoY,
    df_common_factor$Matrix_F,
    use = "complete.obs"
  )
)

sign_vector <- sign(
  cor(
    df_common_factor$GDP_YoY,
    df_common_factor$Vector_F,
    use = "complete.obs"
  )
)

df_factor_yoy_plot <- df_common_factor %>%
  dplyr::transmute(
    date,
    
    `EA GDP YoY` =
      as.numeric(scale(GDP_YoY)),
    
    `Matrix MF-TPRF` =
      as.numeric(
        scale(sign_matrix * Matrix_F)
      ),
    
    `Vector MF-TPRF` =
      as.numeric(
        scale(sign_vector * Vector_F)
      )
  ) %>%
  tidyr::pivot_longer(
    cols = -date,
    names_to = "series",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    series = factor(
      series,
      levels = c(
        "EA GDP YoY",
        "Matrix MF-TPRF",
        "Vector MF-TPRF"
      )
    )
  )

p_factor_yoy <- ggplot(
  df_factor_yoy_plot,
  aes(
    x = date,
    y = value,
    colour = series,
    linetype = series
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "grey60"
  ) +
  geom_line(
    linewidth = 0.95
  ) +
  scale_colour_manual(
    values = c(
      "EA GDP YoY"     = "#C00000",
      "Matrix MF-TPRF" = "#228B22",
      "Vector MF-TPRF" = "#1F77B4"
    ),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(
      "EA GDP YoY"     = "solid",
      "Matrix MF-TPRF" = "solid",
      "Vector MF-TPRF" = "33"
    ),
    name = NULL
  ) +
  scale_x_date(
    date_breaks = "3 years",
    date_labels = "%Y",
    expand = expansion(
      mult = c(0.01, 0.02)
    )
  ) +
  labs(
    title = "EA year-on-year GDP growth and estimated factors",
    subtitle = paste0(
      "Standardised series — ",
      ifelse(
        sel_factor == "corr",
        "Correlation screening",
        "LASSO screening"
      )
    ),
    x = NULL,
    y = "Standardised units"
  ) +
  theme_country_compare(
    base_size = 13
  ) +
  theme(
    legend.position = "bottom"
  )

print(p_factor_yoy)

file_factor_yoy <- file.path(
  path_ea_nowcast_plots,
  paste0(
    "plot_EA_factor_vs_YoY_GDP_sel-",
    sel_factor,
    ".png"
  )
)

ggsave(
  filename = file_factor_yoy,
  plot = p_factor_yoy,
  width = 11.5,
  height = 5.0,
  dpi = 500,
  bg = "white"
)


# ==============================================================================
# 5.1 SINGLE-VINTAGE EA NOWCAST PLOTS
# Full sample + Pre-COVID / COVID / Post-COVID
#
# For each screening rule:
#   LASSO / corr
#
# and for each vintage:
#   M1 / M2 / M3
#
# save:
#   1. full-sample single-vintage figure
#   2. three-panel regime figure
# ==============================================================================

path_ea_regime_split <- file.path(
  path_ea_nowcast_plots,
  "regime_split"
)

dir.create(
  path_ea_regime_split,
  recursive = TRUE,
  showWarnings = FALSE
)

ea_vintage_regime_plots <- list()
ea_vintage_regime_files <- list()


# ------------------------------------------------------------------------------
# Display labels
# ------------------------------------------------------------------------------

ea_regime_series_levels <- c(
  "Observed EA GDP",
  "Matrix MF-TPRF",
  "EA-vector MF-TPRF"
)

ea_regime_colors <- c(
  "Observed EA GDP"      = "#C00000",
  "Matrix MF-TPRF"       = "#228B22",
  "EA-vector MF-TPRF"    = "#1F77B4"
)

ea_regime_linetypes <- c(
  "Observed EA GDP"      = "solid",
  "Matrix MF-TPRF"       = "solid",
  "EA-vector MF-TPRF"    = "33"
)


# ==============================================================================
# LOOP OVER SCREENING RULE AND VINTAGE
# ==============================================================================

for (sel_i in ea_selection_grid) {
  
  selection_label <- if (
    identical(
      sel_i,
      "corr"
    )
  ) {
    "Correlation screening"
  } else {
    "LASSO screening"
  }
  
  
  # ---------------------------------------------------------------------------
  # Use already aligned plot data from Section 4.
  #
  # Only the display name VEC-C -> EA-vector MF-TPRF is changed here.
  # ---------------------------------------------------------------------------
  
  df_plot_selection <- ea_nowcast_comparison[[sel_i]]$plot_data %>%
    dplyr::mutate(
      series = dplyr::recode(
        as.character(series),
        "VEC-C" = "EA-vector MF-TPRF"
      )
    )
  
  
  ea_vintage_regime_plots[[sel_i]] <- list()
  ea_vintage_regime_files[[sel_i]] <- list()
  
  
  for (vintage_i in ea_month_order) {
    
    # -------------------------------------------------------------------------
    # OUTPUT DIRECTORY
    #
    # EA_nowcast_comparison/
    #   regime_split/
    #     LASSO/
    #       M1/
    #       M2/
    #       M3/
    #     corr/
    #       M1/
    #       M2/
    #       M3/
    # -------------------------------------------------------------------------
    
    path_vintage <- file.path(
      path_ea_regime_split,
      sel_i,
      vintage_i
    )
    
    dir.create(
      path_vintage,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    
    # -------------------------------------------------------------------------
    # BUILD PLOTS
    # -------------------------------------------------------------------------
    
    plot_out <- make_nowcast_regime_plots(
      df_plot          = df_plot_selection,
      vintage          = vintage_i,
      params           = params_plot,
      observed_series  = "Observed EA GDP",
      series_levels    = ea_regime_series_levels,
      series_colors    = ea_regime_colors,
      series_linetypes = ea_regime_linetypes,
      title_prefix     = "Expanding pseudo-real-time EA nowcasts",
      subtitle         = paste0(
        "Small information set — ",
        selection_label
      ),
      y_label          = "EA GDP growth",
      pad_frac         = gdp_axis_pad_frac,
      min_pad          = gdp_axis_min_pad
    )
    
    
    # -------------------------------------------------------------------------
    # FILE NAMES
    # -------------------------------------------------------------------------
    
    file_full <- file.path(
      path_vintage,
      paste0(
        "plot_EA_",
        vintage_i,
        "_full_sample",
        "_proxy-",
        proxy_mode,
        "_Size-",
        size_ea,
        "_sel-",
        sel_i,
        ".png"
      )
    )
    
    file_regimes <- file.path(
      path_vintage,
      paste0(
        "plot_EA_",
        vintage_i,
        "_regimes",
        "_proxy-",
        proxy_mode,
        "_Size-",
        size_ea,
        "_sel-",
        sel_i,
        ".png"
      )
    )
    
    
    # -------------------------------------------------------------------------
    # SAVE FULL SAMPLE
    # -------------------------------------------------------------------------
    
    ggplot2::ggsave(
      filename = file_full,
      plot     = plot_out$full_sample,
      width    = 10.8,
      height   = 5.0,
      dpi      = 500,
      bg       = "white"
    )
    
    
    # -------------------------------------------------------------------------
    # SAVE THREE-PANEL REGIME FIGURE
    # -------------------------------------------------------------------------
    
    ggplot2::ggsave(
      filename = file_regimes,
      plot     = plot_out$regimes,
      width    = 13.2,
      height   = 4.9,
      dpi      = 500,
      bg       = "white"
    )
    
    
    # -------------------------------------------------------------------------
    # PRINT TO DEVICE
    # -------------------------------------------------------------------------
    
    print(
      plot_out$full_sample
    )
    
    print(
      plot_out$regimes
    )
    
    
    # -------------------------------------------------------------------------
    # STORE
    # -------------------------------------------------------------------------
    
    ea_vintage_regime_plots[[sel_i]][[vintage_i]] <- plot_out
    
    ea_vintage_regime_files[[sel_i]][[vintage_i]] <- list(
      full_sample = file_full,
      regimes     = file_regimes
    )
    
    
    cat(
      "\nSaved EA ",
      vintage_i,
      " plots | ",
      sel_i,
      "\n  Full sample: ",
      file_full,
      "\n  Regimes:     ",
      file_regimes,
      "\n",
      sep = ""
    )
  }
}

