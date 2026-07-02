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