# ==============================================================================
# Final Results Script
# Matrix MF-TPRF vs VEC-P vs VEC-C vs DFM
# ==============================================================================

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(sandwich)

# ==============================================================================
# 1. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"

path_final_results <- file.path(path_main, "TPRF_Models_EA/Final_Tab_Graph")
dir.create(path_final_results, recursive = TRUE, showWarnings = FALSE)

path_final_factor_graphs <- file.path(path_final_results, "factor_comparison")
dir.create(path_final_factor_graphs, recursive = TRUE, showWarnings = FALSE)

path_matrix_results    <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_vector_results    <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs")
path_vector_results_ea <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs/EA")
path_vectensor_results <- file.path(path_main, "TPRF_Models_EA/VecTensor_MF-TPRF/results/outputs_vec")
path_dfm_results       <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/outputs/q0")
path_dfm_results_auto  <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/outputs")

path_func <- file.path(path_main, "functions/functions_mat")
source(file.path(path_func, "matrix.mf.tprf.utils.R"))

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

params <- list(
  start_est   = as.Date("2000-04-01"),
  start_eval  = as.Date("2017-01-01"),
  end_eval    = as.Date("2026-02-01"),
  covid_start = as.Date("2020-03-01"),
  covid_end   = as.Date("2021-07-01")
)

Size <- "medium"   # "small" | "medium" | "large"
sel  <- "corr"   # "corr"  | "LASSO"

model_matrix    <- "matrix"
model_vector    <- "vector"
model_vectensor <- "vectensor"

country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
period_order  <- c("Pre-COVID", "COVID period", "Post-COVID")
month_order   <- c("M1", "M2", "M3")

period_date_map <- c(
  "Pre-COVID"    = "Jan 2017 -- Feb 2020",
  "COVID period" = "Mar 2020 -- Jul 2021",
  "Post-COVID"   = "Aug 2021 -- Feb 2026"
)

suffix_out <- paste0("Size-", Size, "_sel-", sel)

run_factor_analysis <- is_factor_case(Size, sel)

# ==============================================================================
# 3. SMALL OUTPUT HELPERS
# ==============================================================================

print_latex_block <- function(title, latex_code) {
  cat("\n\n")
  cat("================================================================================\n")
  cat(title, "\n")
  cat("================================================================================\n")
  cat(latex_code)
  cat("\n================================================================================\n\n")
}

safe_print_plot <- function(x, plot_name = "plot") {
  if (!is.null(x)) {
    print(x)
  } else {
    cat("\nSkipping", plot_name, ": object is NULL.\n")
  }
}

safe_write_lines <- function(x, file) {
  if (!is.null(x)) {
    writeLines(x, con = file)
  }
}

# ==============================================================================
# 4. LOAD SUMMARY OBJECTS
# ==============================================================================

file_matrix_summary <- find_result_file(
  path  = path_matrix_results,
  model = model_matrix,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

file_vector_summary <- find_result_file(
  path  = path_vector_results,
  model = model_vector,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

file_vectensor_summary <- find_result_file(
  path  = path_vectensor_results,
  model = model_vectensor,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

file_dfm_summary <- find_dfm_summary_file(
  path = path_dfm_results,
  Size = Size,
  sel  = sel
)

file_dfm_summary_auto <- find_dfm_summary_file(
  path = path_dfm_results_auto,
  Size = Size,
  sel  = sel
)

summary_matrix    <- readRDS(file_matrix_summary)
summary_vector    <- readRDS(file_vector_summary)
summary_vectensor <- readRDS(file_vectensor_summary)
summary_dfm      <- readRDS(file_dfm_summary)
summary_dfm_auto <- readRDS(file_dfm_summary_auto)

cat("\nLoaded matrix summary from:\n", file_matrix_summary, "\n")
cat("\nLoaded vector summary from:\n", file_vector_summary, "\n")
cat("\nLoaded vectensor summary from:\n", file_vectensor_summary, "\n")
cat("\nLoaded DFM summary from:\n", file_dfm_summary, "\n")
cat("\nLoaded DFM summary for q auto from:\n", file_dfm_summary_auto, "\n")


Size <- if (!is.null(summary_matrix$Size)) summary_matrix$Size else Size
sel  <- if (!is.null(summary_matrix$sel))  summary_matrix$sel  else sel

suffix_out <- paste0("Size-", Size, "_sel-", sel)
run_factor_analysis <- is_factor_case(Size, sel)

file_vector_summary_ea <- NULL
summary_vector_ea      <- NULL

if (run_factor_analysis) {
  file_vector_summary_ea <- find_result_file(
    path  = path_vector_results_ea,
    model = model_vector,
    stage = "summary",
    Size  = Size,
    sel   = sel
  )
  
  summary_vector_ea <- readRDS(file_vector_summary_ea)
  cat("\nLoaded vector EA summary from:\n", file_vector_summary_ea, "\n")
} else {
  cat("\nSkipping vector EA summary load: required only for Size = 'small' and sel = 'corr'.\n")
}

# ==============================================================================
# 5. OPTIONAL STATIC FACTOR COMPARISON
# ==============================================================================

df_factor_compare         <- NULL
df_factor_compare_long    <- NULL
plot_factor_compare       <- NULL
file_graph_factor_compare <- NULL

fit_static_factors_long  <- NULL
fit_static_factors_R2    <- NULL
fit_static_factors_AdjR2 <- NULL
latex_factor_fit         <- NULL

if (run_factor_analysis) {
  
  cat("\nRunning factor comparison and static-factor fit.\n")
  
  if (is.null(summary_matrix$factor_comparison)) {
    stop("summary_matrix does not contain factor_comparison.")
  }
  
  if (is.null(summary_vector_ea$factor_comparison)) {
    stop("summary_vector_ea does not contain factor_comparison.")
  }
  
  df_factor_matrix_q <- summary_matrix$factor_comparison$df_q
  df_factor_vector_q <- summary_vector_ea$factor_comparison$df_q
  
  if (is.null(df_factor_matrix_q) || is.null(df_factor_vector_q)) {
    stop("Missing df_q inside factor_comparison.")
  }
  
  df_factor_compare <- df_factor_matrix_q %>%
    dplyr::select(date, GDP, Matrix) %>%
    dplyr::inner_join(
      df_factor_vector_q %>%
        dplyr::select(date, Vector),
      by = "date"
    )
  
  df_factor_compare_long <- df_factor_compare %>%
    tidyr::pivot_longer(
      cols      = -date,
      names_to  = "series",
      values_to = "value"
    ) %>%
    standardize_by_series("value")
  
  plot_factor_compare <- ggplot(
    df_factor_compare_long,
    aes(x = date, y = value_std, colour = series)
  ) +
    annotate(
      "rect",
      xmin = params$covid_start,
      xmax = params$covid_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.15
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 1.05) +
    scale_x_date(
      breaks = seq(
        floor_date(min(df_factor_compare_long$date), unit = "year"),
        floor_date(max(df_factor_compare_long$date), unit = "year"),
        by = "1 year"
      ),
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title  = "Static TPRF factor comparison",
      x      = "Date",
      y      = "Standardized value",
      colour = NULL
    ) +
    theme_factor_compare()
  
  print(plot_factor_compare)
  
  file_graph_factor_compare <- file.path(
    path_final_factor_graphs,
    paste0("plot_factor_comparison_EA_", suffix_out, ".png")
  )
  
  ggsave(
    filename = file_graph_factor_compare,
    plot     = plot_factor_compare,
    width    = 11,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  
  fit_rows_list <- list()
  
  df_mat_q <- make_quarter_design(
    f_q     = summary_matrix$factor_comparison$matrix_factor_quarterly,
    dates_q = summary_matrix$factor_comparison$dates_q,
    y_q     = summary_matrix$factor_comparison$gdp_q
  )
  
  fit_rows_list[["Matrix_Q"]] <- build_fit_stats(
    df         = df_mat_q,
    y_col      = "GDP",
    x_cols     = "Factor",
    spec_name  = "Quarter average",
    model_name = "Matrix"
  )
  
  df_mat_m <- make_monthly_design(
    f_m     = summary_matrix$factor_comparison$matrix_factor_monthly,
    dates_q = summary_matrix$factor_comparison$dates_q,
    y_q     = summary_matrix$factor_comparison$gdp_q
  )
  
  fit_rows_list[["Matrix_M"]] <- build_fit_stats(
    df         = df_mat_m,
    y_col      = "GDP",
    x_cols     = c("M1", "M2", "M3"),
    spec_name  = "Monthly: M1, M2, M3",
    model_name = "Matrix"
  )
  
  df_vec_q <- make_quarter_design(
    f_q     = summary_vector_ea$factor_comparison$vector_factor_quarterly,
    dates_q = summary_vector_ea$factor_comparison$dates_q,
    y_q     = summary_vector_ea$factor_comparison$gdp_q
  )
  
  fit_rows_list[["Vector_Q"]] <- build_fit_stats(
    df         = df_vec_q,
    y_col      = "GDP",
    x_cols     = "Factor",
    spec_name  = "Quarter average",
    model_name = "Vector"
  )
  
  df_vec_m <- make_monthly_design(
    f_m     = summary_vector_ea$factor_comparison$vector_factor_monthly,
    dates_q = summary_vector_ea$factor_comparison$dates_q,
    y_q     = summary_vector_ea$factor_comparison$gdp_q
  )
  
  fit_rows_list[["Vector_M"]] <- build_fit_stats(
    df         = df_vec_m,
    y_col      = "GDP",
    x_cols     = c("M1", "M2", "M3"),
    spec_name  = "Monthly: M1, M2, M3",
    model_name = "Vector"
  )
  
  fit_static_factors_long <- dplyr::bind_rows(fit_rows_list)
  
  fit_static_factors_R2 <- fit_static_factors_long %>%
    dplyr::mutate(R2 = round(100 * R2, 1)) %>%
    dplyr::select(model, specification, R2) %>%
    tidyr::pivot_wider(names_from = model, values_from = R2)
  
  fit_static_factors_AdjR2 <- fit_static_factors_long %>%
    dplyr::mutate(Adj_R2 = round(100 * Adj_R2, 1)) %>%
    dplyr::select(model, specification, Adj_R2) %>%
    tidyr::pivot_wider(names_from = model, values_from = Adj_R2)
  
  latex_factor_fit <- make_factor_fit_latex(
    df_long = fit_static_factors_long,
    caption = paste0(
      "In-sample fit of the extracted static factors (",
      Size, ", sel = ", sel, ")"
    ),
    label = paste0("tab:static_factor_fit_", Size, "_", sel)
  )
  
  print_latex_block("STATIC FACTOR FIT LATEX", latex_factor_fit)
  
} else {
  cat("\nSkipping factor comparison and static-factor fit.\n")
}


# ==============================================================================
# 5.5 MATRIX FACTOR DIAGNOSTICS
# Factor over time + row/column loadings
# Available for small + corr/LASSO
# ==============================================================================

matrix_factor_diag <- NULL

run_matrix_factor_diagnostics <- identical(Size, "small") && sel %in% c("corr", "LASSO")

if (run_matrix_factor_diagnostics) {
  
  cat("\nRunning matrix factor diagnostics: factor over time + row/column loadings.\n")
  
  if (is.null(summary_matrix$factor_comparison)) {
    stop("summary_matrix$factor_comparison is missing.")
  }
  
  if (is.null(summary_matrix$file_fit) || !file.exists(summary_matrix$file_fit)) {
    stop("summary_matrix$file_fit is missing or does not exist.")
  }
  
  path_matrix_factor_diag <- file.path(path_final_results, "matrix_factor_diagnostics")
  dir.create(path_matrix_factor_diag, recursive = TRUE, showWarnings = FALSE)
  
  fit_matrix <- readRDS(summary_matrix$file_fit)
  
  if (is.null(fit_matrix$fit$R) || is.null(fit_matrix$fit$C)) {
    cat("\nAvailable names in fit_matrix$fit:\n")
    print(names(fit_matrix$fit))
    stop("R and/or C loadings not found in fit_matrix$fit.")
  }
  
  # Same sign convention as in the old script
  R_loadings <- -as.matrix(fit_matrix$fit$R)
  C_loadings <- -as.matrix(fit_matrix$fit$C)
  
  # ---------------------------------------------------------------------------
  # 0. Estimated latent factor over time
  # ---------------------------------------------------------------------------
  
  F_hf <- fit_matrix$fit$factors$F_hf
  
  if (is.null(F_hf)) {
    stop("fit_matrix$fit$factors$F_hf is missing.")
  }
  
  dates_f <- as.Date(fit_matrix$fit$factors$dates_m)
  
  df_latent_factor <- expand.grid(
    date = dates_f,
    row_factor = seq_len(dim(F_hf)[2]),
    col_factor = seq_len(dim(F_hf)[3])
  )
  
  df_latent_factor$value <- as.numeric(F_hf)
  
  df_latent_factor <- df_latent_factor %>%
    dplyr::mutate(
      row_factor = paste0("Row ", row_factor),
      col_factor = paste0("Col ", col_factor)
    )
  
  plot_estimated_latent_factor <- ggplot(
    df_latent_factor,
    aes(x = date, y = value)
  ) +
    annotate(
      "rect",
      xmin = params$covid_start,
      xmax = params$covid_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.15
    ) +
    geom_hline(yintercept = 0, linewidth = 0.30, colour = "grey70") +
    geom_line(linewidth = 0.85, colour = "black") +
    facet_grid(row_factor ~ col_factor) +
    scale_x_date(
      breaks = seq(
        floor_date(min(df_latent_factor$date), unit = "year"),
        floor_date(max(df_latent_factor$date), unit = "year"),
        by = "1 year"
      ),
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = "Estimated latent factor",
      x     = "Date",
      y     = "Factor value"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      strip.text       = element_text(face = "bold")
    )
  
  print(plot_estimated_latent_factor)
  
  file_estimated_latent_factor <- file.path(
    path_matrix_factor_diag,
    paste0("plot_estimated_latent_factor_", suffix_out, ".png")
  )
  
  ggsave(
    filename = file_estimated_latent_factor,
    plot     = plot_estimated_latent_factor,
    width    = 11,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  
  # ---------------------------------------------------------------------------
  # 1. Matrix static factor over time
  # ---------------------------------------------------------------------------
  
  fc <- summary_matrix$factor_comparison
  
  df_factor_matrix_q <- fc$df_q
  
  if (is.null(df_factor_matrix_q)) {
    df_factor_matrix_q <- data.frame(
      date   = as.Date(fc$dates_q),
      GDP    = as.numeric(fc$gdp_q),
      Matrix = as.numeric(fc$matrix_factor_quarterly)
    )
  }
  
  df_factor_matrix_long <- df_factor_matrix_q %>%
    dplyr::select(date, GDP, Matrix) %>%
    tidyr::pivot_longer(
      cols      = c(GDP, Matrix),
      names_to  = "series",
      values_to = "value"
    ) %>%
    standardize_by_series("value")
  
  plot_matrix_factor_over_time <- ggplot(
    df_factor_matrix_long,
    aes(x = date, y = value_std, colour = series)
  ) +
    annotate(
      "rect",
      xmin = params$covid_start,
      xmax = params$covid_end,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey70",
      alpha = 0.15
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_line(linewidth = 1.05) +
    scale_x_date(
      breaks = seq(
        floor_date(min(df_factor_matrix_long$date), unit = "year"),
        floor_date(max(df_factor_matrix_long$date), unit = "year"),
        by = "1 year"
      ),
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title    = "Static Matrix MF-TPRF factor over time",
      subtitle = paste0("Size = ", Size, ", sel = ", sel),
      x        = "Date",
      y        = "Standardized value",
      colour   = NULL
    ) +
    theme_factor_compare()
  
  print(plot_matrix_factor_over_time)
  
  file_matrix_factor_over_time <- file.path(
    path_matrix_factor_diag,
    paste0("plot_matrix_factor_over_time_", suffix_out, ".png")
  )
  
  ggsave(
    filename = file_matrix_factor_over_time,
    plot     = plot_matrix_factor_over_time,
    width    = 11,
    height   = 6,
    dpi      = 300,
    bg       = "white"
  )
  
  # ---------------------------------------------------------------------------
  # 2. Row loadings table
  # ---------------------------------------------------------------------------
  
  row_names <- rownames(R_loadings)
  if (is.null(row_names)) {
    row_names <- paste0("Row ", seq_len(nrow(R_loadings)))
  }
  
  df_row_loadings <- data.frame(
    Row = row_names,
    R_loadings,
    check.names = FALSE
  )
  
  colnames(df_row_loadings)[-1] <- paste0("Factor ", seq_len(ncol(R_loadings)))
  
  cat("\n================ MATRIX ROW LOADINGS TABLE ================\n")
  print(df_row_loadings)
  
  latex_row_loadings <- paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{Matrix MF--TPRF row loadings (", Size, ", sel = ", sel, ")}\n",
    "\\label{tab:matrix_row_loadings_", Size, "_", sel, "}\n",
    "\\begin{tabular}{l", paste(rep("c", ncol(R_loadings)), collapse = ""), "}\n",
    "\\toprule\n",
    paste(colnames(df_row_loadings), collapse = " & "), " \\\\\n",
    "\\midrule\n",
    paste(
      apply(df_row_loadings, 1, function(x) paste(x, collapse = " & ")),
      collapse = " \\\\\n"
    ),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  print_latex_block("MATRIX ROW LOADINGS LATEX", latex_row_loadings)
  
  file_row_loadings_tex <- file.path(
    path_matrix_factor_diag,
    paste0("matrix_row_loadings_", suffix_out, ".tex")
  )
  
  writeLines(latex_row_loadings, file_row_loadings_tex)
  
  # ---------------------------------------------------------------------------
  # 3. Column loadings table
  # ---------------------------------------------------------------------------
  
  col_names <- rownames(C_loadings)
  if (is.null(col_names)) {
    col_names <- paste0("Column ", seq_len(nrow(C_loadings)))
  }
  
  df_col_loadings <- data.frame(
    Column = col_names,
    C_loadings,
    check.names = FALSE
  )
  
  colnames(df_col_loadings)[-1] <- paste0("Factor ", seq_len(ncol(C_loadings)))
  
  df_col_loadings <- df_col_loadings %>%
    dplyr::arrange(dplyr::desc(`Factor 1`))
  
  cat("\n================ MATRIX COLUMN LOADINGS TABLE ================\n")
  print(df_col_loadings)
  
  latex_col_loadings <- paste0(
    "\\begin{table}[!ht]\n",
    "\\centering\n",
    "\\scriptsize\n",
    "\\caption{Matrix MF--TPRF column loadings (", Size, ", sel = ", sel, ")}\n",
    "\\label{tab:matrix_column_loadings_", Size, "_", sel, "}\n",
    "\\begin{tabular}{l", paste(rep("c", ncol(C_loadings)), collapse = ""), "}\n",
    "\\toprule\n",
    paste(colnames(df_col_loadings), collapse = " & "), " \\\\\n",
    "\\midrule\n",
    paste(
      apply(df_col_loadings, 1, function(x) paste(x, collapse = " & ")),
      collapse = " \\\\\n"
    ),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  
  print_latex_block("MATRIX COLUMN LOADINGS LATEX", latex_col_loadings)
  
  file_col_loadings_tex <- file.path(
    path_matrix_factor_diag,
    paste0("matrix_column_loadings_", suffix_out, ".tex")
  )
  
  writeLines(latex_col_loadings, file_col_loadings_tex)
  
  # ---------------------------------------------------------------------------
  # 4. Column loadings barplot, ordered decreasing
  # ---------------------------------------------------------------------------
  
  col_order <- df_col_loadings$Column
  
  df_col_bar <- df_col_loadings %>%
    tidyr::pivot_longer(
      cols      = -Column,
      names_to  = "factor",
      values_to = "loading"
    ) %>%
    dplyr::mutate(
      Column = factor(Column, levels = col_order),
      factor = factor(factor, levels = paste0("Factor ", seq_len(ncol(C_loadings))))
    )
  
  plot_col_loadings <- ggplot(
    df_col_bar,
    aes(x = Column, y = loading, fill = factor)
  ) +
    geom_col(position = "dodge", width = 0.75) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey45") +
    coord_flip() +
    labs(
      title    = "Matrix MF-TPRF column loadings",
      subtitle = paste0("Size = ", Size, ", sel = ", sel),
      x        = NULL,
      y        = "Loading",
      fill     = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey25")
    )
  
  print(plot_col_loadings)
  
  file_col_loadings_png <- file.path(
    path_matrix_factor_diag,
    paste0("plot_matrix_column_loadings_", suffix_out, ".png")
  )
  
  ggsave(
    filename = file_col_loadings_png,
    plot     = plot_col_loadings,
    width    = 8,
    height   = 10.5,
    dpi      = 300,
    bg       = "white"
  )
  
  # ---------------------------------------------------------------------------
  # 5. Column loadings heatmap, ordered decreasing from top to bottom
  # ---------------------------------------------------------------------------
  
  df_col_heatmap <- df_col_loadings %>%
    tidyr::pivot_longer(
      cols      = -Column,
      names_to  = "factor",
      values_to = "loading"
    ) %>%
    dplyr::mutate(
      Column = factor(Column, levels = rev(col_order)),
      factor = factor(factor, levels = paste0("Factor ", seq_len(ncol(C_loadings))))
    )
  
  plot_col_loadings_heatmap <- ggplot(
    df_col_heatmap,
    aes(x = factor, y = Column, fill = loading)
  ) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_gradient2(
      low      = "#2166ac",
      mid      = "white",
      high     = "#b2182b",
      midpoint = 0,
      name     = "Loading"
    ) +
    labs(
      title    = "Matrix MF-TPRF column loadings heatmap",
      subtitle = paste0("Size = ", Size, ", sel = ", sel),
      x        = NULL,
      y        = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid      = element_blank(),
      axis.text.x     = element_text(face = "bold"),
      axis.text.y     = element_text(size = 7),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      plot.subtitle   = element_text(hjust = 0.5, colour = "grey25")
    )
  
  print(plot_col_loadings_heatmap)
  
  file_col_loadings_heatmap <- file.path(
    path_matrix_factor_diag,
    paste0("heatmap_matrix_column_loadings_", suffix_out, ".png")
  )
  
  ggsave(
    filename = file_col_loadings_heatmap,
    plot     = plot_col_loadings_heatmap,
    width    = 7.5,
    height   = 10.5,
    dpi      = 300,
    bg       = "white"
  )
  
  matrix_factor_diag <- list(
    fit_matrix_file              = summary_matrix$file_fit,
    R_loadings                   = R_loadings,
    C_loadings                   = C_loadings,
    df_row_loadings              = df_row_loadings,
    df_col_loadings              = df_col_loadings,
    latex_row_loadings           = latex_row_loadings,
    latex_col_loadings           = latex_col_loadings,
    plot_matrix_factor_over_time = plot_matrix_factor_over_time,
    plot_col_loadings            = plot_col_loadings,
    plot_col_loadings_heatmap    = plot_col_loadings_heatmap,
    file_matrix_factor_over_time = file_matrix_factor_over_time,
    file_row_loadings_tex        = file_row_loadings_tex,
    file_col_loadings_tex        = file_col_loadings_tex,
    file_col_loadings_png        = file_col_loadings_png,
    file_col_loadings_heatmap    = file_col_loadings_heatmap
  )
  
} else {
  cat("\nSkipping matrix factor diagnostics: available only for small + corr/LASSO.\n")
}

# ==============================================================================
# 6. NORMALIZE MODEL-SPECIFIC TABLES
# ==============================================================================

df_matrix_insample <- normalize_summary_table(summary_matrix$tab_insample_all, "Matrix")
df_matrix_rt       <- normalize_summary_table(summary_matrix$tab_rt_all, "Matrix")

df_vector_insample <- normalize_summary_table(summary_vector$tab_insample_all, "Vector")
df_vector_rt       <- normalize_summary_table(summary_vector$tab_rt_all, "Vector")

df_vectensor_insample <- normalize_summary_table(summary_vectensor$tab_insample_all, "VecTensor")
df_vectensor_rt       <- normalize_summary_table(summary_vectensor$tab_rt_all, "VecTensor")

df_dfm_insample <- normalize_summary_table(summary_dfm$tab_insample_all, "DFM")
df_dfm_rt       <- normalize_summary_table(summary_dfm$tab_rt_all, "DFM")

common_countries_insample <- Reduce(
  intersect,
  list(
    unique(df_matrix_insample$country),
    unique(df_vector_insample$country),
    unique(df_vectensor_insample$country),
    unique(df_dfm_insample$country)
  )
)

common_countries_rt <- Reduce(
  intersect,
  list(
    unique(df_matrix_rt$country),
    unique(df_vector_rt$country),
    unique(df_vectensor_rt$country),
    unique(df_dfm_rt$country)
  )
)

df_matrix_insample    <- df_matrix_insample    %>% dplyr::filter(country %in% common_countries_insample)
df_vector_insample    <- df_vector_insample    %>% dplyr::filter(country %in% common_countries_insample)
df_vectensor_insample <- df_vectensor_insample %>% dplyr::filter(country %in% common_countries_insample)
df_dfm_insample       <- df_dfm_insample       %>% dplyr::filter(country %in% common_countries_insample)

df_matrix_rt    <- df_matrix_rt    %>% dplyr::filter(country %in% common_countries_rt)
df_vector_rt    <- df_vector_rt    %>% dplyr::filter(country %in% common_countries_rt)
df_vectensor_rt <- df_vectensor_rt %>% dplyr::filter(country %in% common_countries_rt)
df_dfm_rt       <- df_dfm_rt       %>% dplyr::filter(country %in% common_countries_rt)

# ==============================================================================
# 7. MODEL-SPECIFIC PLOTS
# ==============================================================================

safe_print_plot(summary_matrix$plot_nowcast_facet, "matrix nowcast facet")
safe_print_plot(summary_matrix$plot_rt_facet,      "matrix rolling facet")

safe_print_plot(summary_vector$plot_nowcast_facet, "vector nowcast facet")
safe_print_plot(summary_vector$plot_rt_facet,      "vector rolling facet")

safe_print_plot(summary_vectensor$plot_nowcast_facet, "vectensor nowcast facet")
safe_print_plot(summary_vectensor$plot_rt_facet,      "vectensor rolling facet")

safe_print_plot(summary_dfm$plot_nowcast_facet, "DFM nowcast facet")
safe_print_plot(summary_dfm$plot_rt_facet,      "DFM rolling facet")

# ==============================================================================
# 8. COUNTRY COMPARISON PLOTS
# ==============================================================================

label_matrix <- "Matrix MF-TPRF"
label_vecp   <- "VEC-P"
label_vecc   <- "VEC-C"
label_dfm    <- "DFM"

model_levels <- c(label_matrix, label_vecp, label_vecc, label_dfm)

model_colors <- c(
  "Matrix MF-TPRF" = "#1b6ca8",
  "VEC-P"          = "#55a868",
  "VEC-C"          = "#c44e52",
  "DFM"            = "grey35"
)

model_linetypes <- c(
  "Matrix MF-TPRF" = "solid",
  "VEC-P"          = "11",
  "VEC-C"          = "33",
  "DFM"            = "44"
)

country_labels <- c(
  "NL" = "Netherlands",
  "BE" = "Belgium",
  "AT" = "Austria",
  "PT" = "Portugal",
  "DE" = "Germany",
  "FR" = "France",
  "IT" = "Italy",
  "ES" = "Spain"
)

# ------------------------------------------------------------------------------
# 8.1 Single-country plot
# ------------------------------------------------------------------------------

country_focus <- "PT"
country_focus_label <- unname(country_labels[country_focus])

df_single_all <- dplyr::bind_rows(
  extract_country_rt(summary_matrix$df_rt_all,    country_focus, label_matrix),
  extract_country_rt(summary_vectensor$df_rt_all, country_focus, label_vecp),
  extract_country_rt(summary_vector$df_rt_all,    country_focus, label_vecc),
  extract_country_rt(summary_dfm$df_rt_all,       country_focus, label_dfm)
) %>%
  dplyr::mutate(
    model = factor(model, levels = model_levels),
    type  = factor(type, levels = month_order)
  )

df_single_gdp <- extract_country_gdp(
  summary_matrix$df_yq_eval_all,
  country_focus
)

df_single_gdp_facet <- tidyr::crossing(
  df_single_gdp,
  type = factor(month_order, levels = month_order)
)

plot_single_models_facet <- ggplot() +
  annotate(
    "rect",
    xmin = params$covid_start,
    xmax = params$covid_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey72",
    alpha = 0.11
  ) +
  geom_line(
    data = df_single_gdp_facet,
    aes(x = date, y = GDP, group = type),
    colour = "black",
    linewidth = 1.25,
    linetype = "solid",
    lineend = "round"
  ) +
  geom_line(
    data = df_single_all,
    aes(x = date, y = nowcast, colour = model, linetype = model),
    linewidth = 1.00,
    alpha = 0.98,
    lineend = "round"
  ) +
  facet_wrap(~ type, ncol = 3) +
  scale_color_manual(values = model_colors, breaks = model_levels, name = NULL) +
  scale_linetype_manual(values = model_linetypes, breaks = model_levels, name = NULL) +
  scale_x_date(
    breaks = seq(
      floor_date(min(df_single_all$date, na.rm = TRUE), unit = "year"),
      floor_date(max(df_single_all$date, na.rm = TRUE), unit = "year"),
      by = "2 years"
    ),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = paste0("Expanding pseudo-real-time nowcasts for ", country_focus_label),
    x = NULL,
    y = "GDP growth"
  ) +
  theme_country_compare(base_size = 13)

print(plot_single_models_facet)

file_single_png <- file.path(
  path_final_results,
  paste0("plot_", country_focus, "_models_by_month_paper_", suffix_out, ".png")
)

ggsave(
  filename = file_single_png,
  plot     = plot_single_models_facet,
  width    = 12.5,
  height   = 4.8,
  dpi      = 500,
  bg       = "white"
)

# ------------------------------------------------------------------------------
# 8.2 Multi-country plot
# ------------------------------------------------------------------------------

countries_focus <- c("DE", "FR", "IT", "ES")

df_multi_all <- dplyr::bind_rows(
  extract_country_rt_multi(summary_matrix$df_rt_all,    countries_focus, label_matrix),
  extract_country_rt_multi(summary_vectensor$df_rt_all, countries_focus, label_vecp),
  extract_country_rt_multi(summary_vector$df_rt_all,    countries_focus, label_vecc),
  extract_country_rt_multi(summary_dfm$df_rt_all,       countries_focus, label_dfm)
) %>%
  dplyr::mutate(
    model   = factor(model, levels = model_levels),
    type    = factor(type, levels = month_order),
    country = factor(country, levels = countries_focus)
  )

df_multi_gdp <- extract_country_gdp_multi(
  summary_matrix$df_yq_eval_all,
  countries_focus
)

df_multi_gdp_facet <- tidyr::crossing(
  df_multi_gdp,
  type = factor(month_order, levels = month_order)
)

countries_tag  <- paste(countries_focus, collapse = "-")
countries_text <- paste(countries_focus, collapse = ", ")

plot_multi_models_facet <- ggplot() +
  annotate(
    "rect",
    xmin = params$covid_start,
    xmax = params$covid_end,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey72",
    alpha = 0.11
  ) +
  geom_line(
    data = df_multi_gdp_facet,
    aes(x = date, y = GDP, group = interaction(country, type)),
    colour = "black",
    linewidth = 1.20,
    linetype = "solid",
    lineend = "round"
  ) +
  geom_line(
    data = df_multi_all,
    aes(x = date, y = nowcast, colour = model, linetype = model),
    linewidth = 0.95,
    alpha = 0.99,
    lineend = "round"
  ) +
  facet_grid(country ~ type, scales = "fixed") +
  scale_color_manual(values = model_colors, breaks = model_levels, name = NULL) +
  scale_linetype_manual(values = model_linetypes, breaks = model_levels, name = NULL) +
  scale_x_date(
    breaks = seq(
      floor_date(min(df_multi_all$date, na.rm = TRUE), unit = "year"),
      floor_date(max(df_multi_all$date, na.rm = TRUE), unit = "year"),
      by = "2 years"
    ),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = "Expanding pseudo-real-time nowcasts",
    subtitle = paste0("Size = ", Size, ", sel = ", sel, " | Countries: ", countries_text),
    x = NULL,
    y = "GDP growth"
  ) +
  theme_country_compare(base_size = 10.4)

print(plot_multi_models_facet)

file_multi_png <- file.path(
  path_final_results,
  paste0("mcplot_", suffix_out, "_cty-", countries_tag, ".png")
)

ggsave(
  filename = file_multi_png,
  plot     = plot_multi_models_facet,
  width    = 12.4,
  height   = 2.5 + 1.80 * length(countries_focus),
  dpi      = 600,
  bg       = "white"
)

# ==============================================================================
# 9. MODEL-SPECIFIC LATEX TABLES
# ==============================================================================

print_latex_block("MATRIX IN-SAMPLE", summary_matrix$latex_tab_insample_all)
print_latex_block("MATRIX ROLLING",   summary_matrix$latex_tab_rt_all)

print_latex_block("VECTOR IN-SAMPLE", summary_vector$latex_tab_insample_all)
print_latex_block("VECTOR ROLLING",   summary_vector$latex_tab_rt_all)

print_latex_block("VECTENSOR IN-SAMPLE", summary_vectensor$latex_tab_insample_all)
print_latex_block("VECTENSOR ROLLING",   summary_vectensor$latex_tab_rt_all)

print_latex_block("DFM IN-SAMPLE", summary_dfm$latex_tab_insample_all)
print_latex_block("DFM ROLLING",   summary_dfm$latex_tab_rt_all)

# ==============================================================================
# 10. COMPARISON AND RELATIVE RMSFE TABLES
# ==============================================================================

comp_insample <- build_comparison_table(
  df_matrix_insample,
  df_vector_insample,
  df_vectensor_insample,
  df_dfm_insample
)

comp_rt <- build_comparison_table(
  df_matrix_rt,
  df_vector_rt,
  df_vectensor_rt,
  df_dfm_rt
)

latex_comp_insample <- comparison_to_latex(
  comp_insample,
  caption = paste0("In-sample RMSFE comparison across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:comparison_insample_models_", Size, "_", sel)
)

latex_comp_rt <- comparison_to_latex(
  comp_rt,
  caption = paste0("Rolling pseudo-real-time RMSFE comparison across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:comparison_rt_models_", Size, "_", sel)
)

vecc_benchmark_rt <- build_vecc_benchmark_table(
  df_comp     = comp_rt,
  include_dfm = TRUE
)

rel_insample <- build_relative_table(comp_insample)
rel_rt       <- build_relative_table(comp_rt)

latex_rel_insample <- relative_to_latex(
  rel_insample,
  caption = paste0("Relative in-sample RMSFE across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:relative_insample_models_", Size, "_", sel)
)

latex_rel_rt <- relative_to_latex(
  rel_rt,
  caption = paste0("Relative rolling pseudo-real-time RMSFE across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:relative_rt_models_", Size, "_", sel)
)

latex_vecc_benchmark_rt <- vecc_benchmark_to_latex(
  df      = vecc_benchmark_rt,
  caption = paste0("VEC-C RMSFE and relative RMSFE against VEC-P and DFM (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:vecc_benchmark_rt_", Size, "_", sel)
)

print_latex_block("COMPARISON IN-SAMPLE LATEX", latex_comp_insample)
print_latex_block("COMPARISON ROLLING LATEX",   latex_comp_rt)

print_latex_block("RELATIVE IN-SAMPLE LATEX", latex_rel_insample)
print_latex_block("RELATIVE ROLLING LATEX",   latex_rel_rt)

print_latex_block("VEC-C BENCHMARK ROLLING LATEX", latex_vecc_benchmark_rt)


# ==============================================================================
# 11. DIEBOLD-MARIANO TESTS
# ==============================================================================

eval_matrix <- build_eval_df(summary_matrix$df_rt_all, summary_matrix$df_yq_eval_all)
eval_vecp   <- build_eval_df(summary_vectensor$df_rt_all, summary_matrix$df_yq_eval_all)
eval_vecc   <- build_eval_df(summary_vector$df_rt_all, summary_matrix$df_yq_eval_all)
eval_dfm    <- build_eval_df(summary_dfm$df_rt_all, summary_matrix$df_yq_eval_all)

dm_vecp_obj <- build_dm_wide(
  eval_matrix    = eval_matrix,
  eval_benchmark = eval_vecp,
  benchmark_tag  = "VecP",
  alternative    = "less"
)

dm_vecc_obj <- build_dm_wide(
  eval_matrix    = eval_matrix,
  eval_benchmark = eval_vecc,
  benchmark_tag  = "VecC",
  alternative    = "less"
)

dm_dfm_obj <- build_dm_wide(
  eval_matrix    = eval_matrix,
  eval_benchmark = eval_dfm,
  benchmark_tag  = "DFM",
  alternative    = "less"
)

dm_vecp <- dm_vecp_obj$long
dm_vecc <- dm_vecc_obj$long
dm_dfm  <- dm_dfm_obj$long

dm_vecp_wide <- dm_vecp_obj$wide
dm_vecc_wide <- dm_vecc_obj$wide
dm_dfm_wide  <- dm_dfm_obj$wide


dm_vecc_vecp_obj <- build_dm_wide(
  eval_matrix    = eval_vecc,
  eval_benchmark = eval_vecp,
  benchmark_tag  = "VecP",
  alternative    = "less"
)

dm_vecc_dfm_obj <- build_dm_wide(
  eval_matrix    = eval_vecc,
  eval_benchmark = eval_dfm,
  benchmark_tag  = "DFM",
  alternative    = "less"
)

dm_vecc_vecp_wide <- dm_vecc_vecp_obj$wide
dm_vecc_dfm_wide  <- dm_vecc_dfm_obj$wide

# ==============================================================================
# 12. FINAL CUSTOM LATEX TABLE
# ==============================================================================

latex_final_custom <- build_final_large_style_latex(
  summary_matrix    = summary_matrix,
  summary_vector    = summary_vector,
  summary_vectensor = summary_vectensor,
  summary_dfm       = summary_dfm,
  params            = params,
  dm_vecp_wide      = dm_vecp_wide,
  dm_vecc_wide      = dm_vecc_wide,
  dm_dfm_wide       = dm_dfm_wide,
  Size              = Size,
  sel               = sel
)

print_latex_block("FINAL CUSTOM TABLE LATEX", latex_final_custom)


latex_final_vecc_custom <- build_final_vecc_style_latex(
  summary_vector    = summary_vector,
  summary_vectensor = summary_vectensor,
  summary_dfm       = summary_dfm,
  params            = params,
  dm_vecp_wide      = dm_vecc_vecp_wide,
  dm_dfm_wide       = dm_vecc_dfm_wide,
  Size              = Size,
  sel               = sel
)

print_latex_block("FINAL VEC-C CUSTOM TABLE LATEX", latex_final_vecc_custom)

# ==============================================================================
# 13. SYNTHETIC SUMMARY + FORECAST-LEVEL WIN ANALYSIS
# ==============================================================================

comp_long <- reshape_comp_to_long(comp_rt, include_dfm = TRUE)
df_eval   <- compute_summary_metrics(comp_long)

tab_period_month <- build_period_month_table(df_eval, include_dfm = TRUE)

latex_period_month <- make_latex_period_month_wins(
  df = tab_period_month,
  note = "Share of cases in which Matrix MF--TPRF outperforms each benchmark."
)

wins_forecast <- eval_matrix %>%
  dplyr::select(country, quarter_id, type, period, se_matrix = se) %>%
  dplyr::inner_join(
    eval_vecp %>%
      dplyr::select(country, quarter_id, type, se_vecp = se),
    by = c("country", "quarter_id", "type")
  ) %>%
  dplyr::inner_join(
    eval_vecc %>%
      dplyr::select(country, quarter_id, type, se_vecc = se),
    by = c("country", "quarter_id", "type")
  ) %>%
  dplyr::inner_join(
    eval_dfm %>%
      dplyr::select(country, quarter_id, type, se_dfm = se),
    by = c("country", "quarter_id", "type")
  ) %>%
  dplyr::mutate(
    country  = factor(country, levels = country_order),
    type     = factor(type, levels = month_order),
    period   = factor(period, levels = c("Full sample", period_order)),
    win_vecp = as.integer(se_matrix < se_vecp),
    win_vecc = as.integer(se_matrix < se_vecc),
    win_dfm  = as.integer(se_matrix < se_dfm)
  ) %>%
  dplyr::filter(!is.na(period))

tab_forecast_overall <- build_forecast_win_table(
  wins_forecast,
  include_dfm = TRUE
)

tab_forecast_period <- wins_forecast %>%
  dplyr::filter(period %in% period_order) %>%
  build_forecast_win_table(group_vars = "period", include_dfm = TRUE) %>%
  dplyr::mutate(period = factor(period, levels = period_order)) %>%
  dplyr::arrange(period)

tab_forecast_vintage <- wins_forecast %>%
  build_forecast_win_table(group_vars = "type", include_dfm = TRUE) %>%
  dplyr::mutate(type = factor(type, levels = month_order)) %>%
  dplyr::arrange(type)

tab_forecast_country_period_month <- build_country_period_month_table(
  wins_forecast,
  include_dfm = TRUE
)

latex_forecast_summary <- make_latex_forecast_wins_summary(
  overall_df = tab_forecast_overall,
  period_df  = tab_forecast_period,
  vintage_df = tab_forecast_vintage
)

latex_forecast_country_period_month <- make_latex_country_period_month_wins(
  df = tab_forecast_country_period_month
)

print_latex_block("SUMMARY WINS BY PERIOD x MONTH LATEX", latex_period_month)
print_latex_block("FORECAST-LEVEL SUMMARY LATEX",         latex_forecast_summary)
print_latex_block("FORECAST WINS COUNTRY x PERIOD x MONTH LATEX", latex_forecast_country_period_month)

# ==============================================================================
# 14. SAVE TEX OUTPUTS
# ==============================================================================

safe_write_lines(
  summary_matrix$latex_tab_insample_all,
  file.path(path_final_results, paste0("matrix_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_matrix$latex_tab_rt_all,
  file.path(path_final_results, paste0("matrix_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_vector$latex_tab_insample_all,
  file.path(path_final_results, paste0("vector_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_vector$latex_tab_rt_all,
  file.path(path_final_results, paste0("vector_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_vectensor$latex_tab_insample_all,
  file.path(path_final_results, paste0("vectensor_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_vectensor$latex_tab_rt_all,
  file.path(path_final_results, paste0("vectensor_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_dfm$latex_tab_insample_all,
  file.path(path_final_results, paste0("dfm_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  summary_dfm$latex_tab_rt_all,
  file.path(path_final_results, paste0("dfm_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_comp_insample,
  file.path(path_final_results, paste0("comparison_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_comp_rt,
  file.path(path_final_results, paste0("comparison_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_rel_insample,
  file.path(path_final_results, paste0("relative_insample_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_rel_rt,
  file.path(path_final_results, paste0("relative_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_vecc_benchmark_rt,
  file.path(path_final_results, paste0("vecc_benchmark_rt_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_final_custom,
  file.path(path_final_results, paste0("final_custom_table_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_period_month,
  file.path(path_final_results, paste0("summary_wins_period_month_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_forecast_summary,
  file.path(path_final_results, paste0("forecast_wins_summary_", suffix_out, ".tex"))
)

safe_write_lines(
  latex_forecast_country_period_month,
  file.path(path_final_results, paste0("forecast_wins_country_period_month_", suffix_out, ".tex"))
)

if (!is.null(latex_factor_fit)) {
  safe_write_lines(
    latex_factor_fit,
    file.path(path_final_results, paste0("static_factor_fit_", suffix_out, ".tex"))
  )
}

# ==============================================================================
# 15. SAVE RDS OUTPUT
# ==============================================================================

saveRDS(
  list(
    Size                = Size,
    sel                 = sel,
    params              = params,
    run_factor_analysis = run_factor_analysis,
    
    run_matrix_factor_diagnostics  = run_matrix_factor_diagnostics,
    matrix_factor_diag             = matrix_factor_diag,
    
    file_matrix_summary    = file_matrix_summary,
    file_vector_summary    = file_vector_summary,
    file_vector_summary_ea = file_vector_summary_ea,
    file_vectensor_summary = file_vectensor_summary,
    file_dfm_summary       = file_dfm_summary,
    
    summary_matrix    = summary_matrix,
    summary_vector    = summary_vector,
    summary_vector_ea = summary_vector_ea,
    summary_vectensor = summary_vectensor,
    summary_dfm       = summary_dfm,
    
    df_matrix_insample    = df_matrix_insample,
    df_matrix_rt          = df_matrix_rt,
    df_vector_insample    = df_vector_insample,
    df_vector_rt          = df_vector_rt,
    df_vectensor_insample = df_vectensor_insample,
    df_vectensor_rt       = df_vectensor_rt,
    df_dfm_insample       = df_dfm_insample,
    df_dfm_rt             = df_dfm_rt,
    
    comp_insample = comp_insample,
    comp_rt       = comp_rt,
    rel_insample  = rel_insample,
    rel_rt        = rel_rt,
    
    vecc_benchmark_rt       = vecc_benchmark_rt,
    latex_vecc_benchmark_rt = latex_vecc_benchmark_rt,
    
    eval_matrix = eval_matrix,
    eval_vecp   = eval_vecp,
    eval_vecc   = eval_vecc,
    eval_dfm    = eval_dfm,
    
    dm_vecp      = dm_vecp,
    dm_vecc      = dm_vecc,
    dm_dfm       = dm_dfm,
    dm_vecp_wide = dm_vecp_wide,
    dm_vecc_wide = dm_vecc_wide,
    dm_dfm_wide  = dm_dfm_wide,
    
    wins_forecast = wins_forecast,
    
    tab_period_month                  = tab_period_month,
    tab_forecast_overall              = tab_forecast_overall,
    tab_forecast_period               = tab_forecast_period,
    tab_forecast_vintage              = tab_forecast_vintage,
    tab_forecast_country_period_month = tab_forecast_country_period_month,
    
    fit_static_factors_long  = fit_static_factors_long,
    fit_static_factors_R2    = fit_static_factors_R2,
    fit_static_factors_AdjR2 = fit_static_factors_AdjR2,
    
    df_factor_compare         = df_factor_compare,
    df_factor_compare_long    = df_factor_compare_long,
    plot_factor_compare       = plot_factor_compare,
    file_graph_factor_compare = file_graph_factor_compare,
    
    plot_single_models_facet = plot_single_models_facet,
    file_single_png          = file_single_png,
    plot_multi_models_facet  = plot_multi_models_facet,
    file_multi_png           = file_multi_png,
    
    latex_factor_fit                    = latex_factor_fit,
    latex_comp_insample                 = latex_comp_insample,
    latex_comp_rt                       = latex_comp_rt,
    latex_rel_insample                  = latex_rel_insample,
    latex_rel_rt                        = latex_rel_rt,
    latex_final_custom                  = latex_final_custom,
    latex_period_month                  = latex_period_month,
    latex_forecast_summary              = latex_forecast_summary,
    latex_forecast_country_period_month = latex_forecast_country_period_month
  ),
  file.path(path_final_results, paste0("final_model_comparison_", suffix_out, ".rds"))
)

cat("\nAll final outputs saved to:\n", path_final_results, "\n")


