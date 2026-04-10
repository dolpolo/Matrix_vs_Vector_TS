# ==============================================================================
# Final Results Script
# Matrix MF-TPRF vs Vector MF-TPRF vs Vectorized-from-Tensor MF-TPRF
# ------------------------------------------------------------------------------
# This script:
#   1. Loads cross-country summary objects for a chosen Size and sel
#   2. Prints model-specific plots and LaTeX tables
#   3. Builds comparison tables across models
#   4. Builds relative RMSFE tables
#   5. Builds optional static-factor comparison only for:
#        Size = "small" and sel = "corr"
#   6. Saves final comparison outputs
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

path_final_results <- file.path(path_main, "TPRF_Models_EA/Final_Tab_Graph")
dir.create(path_final_results, recursive = TRUE, showWarnings = FALSE)

path_matrix_results <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_vector_results <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs")
path_vector_results_ea <- file.path(
  path_main,
  "TPRF_Models_EA/Vector_MF-TPRF/results/outputs",
  "EA"
)
path_vectensor_results <- file.path(
  path_main,
  "TPRF_Models_EA/VecTensor_MF-TPRF/results/outputs_vec"
)

path_final_factor_graphs <- file.path(path_final_results, "factor_comparison")
dir.create(path_final_factor_graphs, recursive = TRUE, showWarnings = FALSE)

path_func <- file.path(path_main, "functions/functions_mat")
source(file.path(path_func, "matrix.mf.tprf.utils.R"))

# ==============================================================================
# 2. HELPERS
# ==============================================================================

theme_factor_compare <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey25")
    )
}

standardize_by_series <- function(df, value_col = "value") {
  df %>%
    group_by(series) %>%
    mutate(value_std = as.numeric(scale(.data[[value_col]]))) %>%
    ungroup()
}

is_factor_case <- function(Size, sel) {
  identical(Size, "small") && identical(sel, "corr")
}

build_fit_stats <- function(df, y_col, x_cols, spec_name, model_name) {
  form <- as.formula(
    paste(y_col, "~", paste(x_cols, collapse = " + "))
  )
  
  fit <- lm(form, data = df)
  s   <- summary(fit)
  
  data.frame(
    model              = model_name,
    specification      = spec_name,
    R2                 = unname(s$r.squared),
    Adj_R2             = unname(s$adj.r.squared),
    stringsAsFactors   = FALSE
  )
}

make_quarter_design <- function(f_q, dates_q, y_q) {
  data.frame(
    date   = as.Date(dates_q),
    GDP    = as.numeric(y_q),
    Factor = as.numeric(f_q)
  )
}

make_monthly_design <- function(f_m, dates_q, y_q) {
  T_q <- length(y_q)
  stopifnot(length(f_m) >= 3 * T_q)
  
  f_use <- as.numeric(f_m[seq_len(3 * T_q)])
  
  data.frame(
    date = as.Date(dates_q),
    GDP  = as.numeric(y_q),
    M1   = f_use[seq(1, 3 * T_q, by = 3)],
    M2   = f_use[seq(2, 3 * T_q, by = 3)],
    M3   = f_use[seq(3, 3 * T_q, by = 3)]
  )
}

make_factor_fit_latex <- function(df_long, caption, label) {
  
  df_r2 <- df_long %>%
    mutate(R2 = round(100 * R2, 1)) %>%
    select(model, specification, R2) %>%
    pivot_wider(names_from = model, values_from = R2)
  
  df_adj <- df_long %>%
    mutate(Adj_R2 = round(100 * Adj_R2, 1)) %>%
    select(model, specification, Adj_R2) %>%
    pivot_wider(names_from = model, values_from = Adj_R2)
  
  model_cols    <- setdiff(colnames(df_r2), "specification")
  header_models <- paste(model_cols, collapse = " & ")
  
  rows_r2 <- paste(
    apply(df_r2, 1, function(x) {
      paste(c(x["specification"], x[model_cols]), collapse = " & ")
    }),
    collapse = " \\\\\n"
  )
  
  rows_adj <- paste(
    apply(df_adj, 1, function(x) {
      paste(c(x["specification"], x[model_cols]), collapse = " & ")
    }),
    collapse = " \\\\\n"
  )
  
  paste0(
    "\\begin{table}[htbp]\n",
    "\\centering\n",
    "\\caption{", caption, "}\n",
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

extract_country_rt <- function(df_rt, country_code, model_label) {
  df_rt %>%
    filter(country == country_code) %>%
    mutate(
      model = model_label,
      type  = factor(type, levels = c("M1", "M2", "M3"))
    ) %>%
    select(date, country, nowcast, type, model)
}

extract_country_gdp <- function(df_yq, country_code) {
  df_yq %>%
    filter(country == country_code) %>%
    select(date, country, GDP)
}

theme_country_compare <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey25")
    )
}

# ==============================================================================
# 3. CHOOSE CONFIGURATION
# ==============================================================================

params <- list(
  start_est    = as.Date("2000-04-01"),
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2026-02-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01")
)

Size <- "small"   # "small" | "medium" | "large"
sel  <- "LASSO"    # "corr"  | "LASSO"

model_matrix    <- "matrix"
model_vector    <- "vector"
model_vectensor <- "vectensor"

run_factor_analysis <- is_factor_case(Size, sel)

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

summary_matrix    <- readRDS(file_matrix_summary)
summary_vector    <- readRDS(file_vector_summary)
summary_vectensor <- readRDS(file_vectensor_summary)

file_vector_summary_ea <- NULL
summary_vector_ea      <- NULL

cat("\nLoaded matrix summary from:\n", file_matrix_summary, "\n")
cat("\nLoaded vector summary from:\n", file_vector_summary, "\n")
cat("\nLoaded vectensor summary from:\n", file_vectensor_summary, "\n")

# Safety overwrite from loaded objects if available
Size <- if (!is.null(summary_matrix$Size)) summary_matrix$Size else Size
sel  <- if (!is.null(summary_matrix$sel))  summary_matrix$sel  else sel

run_factor_analysis <- is_factor_case(Size, sel)

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
  cat(
    "\nSkipping vector EA summary load: required only for Size = 'small' and sel = 'corr'.\n"
  )
}
# ==============================================================================
# 5. OPTIONAL FACTOR COMPARISON (ONLY FOR small + corr)
# ==============================================================================

df_factor_compare         <- NULL
df_factor_compare_long    <- NULL
plot_factor_compare       <- NULL
file_graph_factor_compare <- NULL

fit_static_factors_long   <- NULL
fit_static_factors_R2     <- NULL
fit_static_factors_AdjR2  <- NULL
latex_factor_fit          <- NULL

if (run_factor_analysis) {
  
  cat("\nRunning factor comparison and static-factor fit (small + corr).\n")
  
  if (is.null(summary_matrix$factor_comparison)) {
    stop("summary_matrix does not contain factor_comparison.")
  }
  if (is.null(summary_vector_ea$factor_comparison)) {
    stop("summary_vector_ea does not contain factor_comparison.")
  }
  
  df_factor_matrix_q <- summary_matrix$factor_comparison$df_q
  df_factor_vector_q <- summary_vector_ea$factor_comparison$df_q
  
  if (is.null(df_factor_matrix_q) || is.null(df_factor_vector_q)) {
    stop("Missing df_q inside factor_comparison for matrix or vector EA summary.")
  }
  
  df_factor_compare <- df_factor_matrix_q %>%
    select(date, GDP, Matrix) %>%
    inner_join(
      df_factor_vector_q %>%
        select(date, Vector),
      by = "date"
    )
  
  df_factor_compare_long <- df_factor_compare %>%
    pivot_longer(
      cols      = -date,
      names_to  = "series",
      values_to = "value"
    ) %>%
    standardize_by_series("value")
  
  plot_factor_compare <- ggplot(
    df_factor_compare_long,
    aes(x = date, y = value_std, color = series)
  ) +
    annotate(
      "rect",
      xmin = params$covid_start, xmax = params$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.15
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
      title    = "Static TPRF factor comparison",
      subtitle = NULL,
      x        = "Date",
      y        = "Factor value",
      color    = ""
    ) +
    theme_factor_compare()
  
  print(plot_factor_compare)
  
  file_graph_factor_compare <- file.path(
    path_final_factor_graphs,
    paste0("plot_factor_comparison_EA_Size-", Size, "_sel-", sel, ".png")
  )
  
  ggsave(
    filename = file_graph_factor_compare,
    plot     = plot_factor_compare,
    width    = 11,
    height   = 6,
    dpi      = 300
  )
  
  cat("\nSaved factor comparison graph to:\n", file_graph_factor_compare, "\n")
  
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
  
  fit_static_factors_long <- bind_rows(fit_rows_list)
  
  fit_static_factors_R2 <- fit_static_factors_long %>%
    mutate(R2 = round(100 * R2, 1)) %>%
    select(model, specification, R2) %>%
    pivot_wider(names_from = model, values_from = R2)
  
  fit_static_factors_AdjR2 <- fit_static_factors_long %>%
    mutate(Adj_R2 = round(100 * Adj_R2, 1)) %>%
    select(model, specification, Adj_R2) %>%
    pivot_wider(names_from = model, values_from = Adj_R2)
  
  cat("\n================ STATIC FACTOR FIT: R2 ================\n")
  print(fit_static_factors_R2)
  
  cat("\n================ STATIC FACTOR FIT: ADJUSTED R2 ================\n")
  print(fit_static_factors_AdjR2)
  
  latex_factor_fit <- make_factor_fit_latex(
    df_long = fit_static_factors_long,
    caption = paste0(
      "In-sample fit of the extracted static factors (",
      Size, ", sel = ", sel, ")"
    ),
    label = paste0("tab:static_factor_fit_", Size, "_", sel)
  )
  
  cat("\n================ STATIC FACTOR FIT LATEX ================\n")
  cat(latex_factor_fit, "\n")
  
} else {
  cat(
    "\nSkipping factor comparison and static-factor fit: available only for Size = 'small' and sel = 'corr'.\n"
  )
}

# ==============================================================================
# 6. MODEL-SPECIFIC TABLES
# ==============================================================================

df_matrix_insample <- normalize_summary_table(
  summary_matrix$tab_insample_all,
  "Matrix"
)

df_matrix_rt <- normalize_summary_table(
  summary_matrix$tab_rt_all,
  "Matrix"
)

df_vector_insample <- normalize_summary_table(
  summary_vector$tab_insample_all,
  "Vector"
)

df_vector_rt <- normalize_summary_table(
  summary_vector$tab_rt_all,
  "Vector"
)

df_vectensor_insample <- normalize_summary_table(
  summary_vectensor$tab_insample_all,
  "VecTensor"
)

df_vectensor_rt <- normalize_summary_table(
  summary_vectensor$tab_rt_all,
  "VecTensor"
)

common_countries_insample <- Reduce(
  intersect,
  list(
    unique(df_matrix_insample$country),
    unique(df_vector_insample$country),
    unique(df_vectensor_insample$country)
  )
)

common_countries_rt <- Reduce(
  intersect,
  list(
    unique(df_matrix_rt$country),
    unique(df_vector_rt$country),
    unique(df_vectensor_rt$country)
  )
)

df_matrix_insample    <- df_matrix_insample    %>% filter(country %in% common_countries_insample)
df_vector_insample    <- df_vector_insample    %>% filter(country %in% common_countries_insample)
df_vectensor_insample <- df_vectensor_insample %>% filter(country %in% common_countries_insample)

df_matrix_rt    <- df_matrix_rt    %>% filter(country %in% common_countries_rt)
df_vector_rt    <- df_vector_rt    %>% filter(country %in% common_countries_rt)
df_vectensor_rt <- df_vectensor_rt %>% filter(country %in% common_countries_rt)

# ==============================================================================
# 7. SHOW MODEL-SPECIFIC PLOTS
# ==============================================================================

print(summary_matrix$plot_nowcast_facet)
print(summary_matrix$plot_rt_facet)

print(summary_vector$plot_nowcast_facet)
print(summary_vector$plot_rt_facet)

print(summary_vectensor$plot_nowcast_facet)
print(summary_vectensor$plot_rt_facet)

# ==============================================================================
# COUNTRY COMPARISON PLOTS
# Single-country and multi-country versions
# ==============================================================================

label_matrix <- "Matrix MF-TPRF"
label_vecp   <- "VEC-P"
label_vecc   <- "VEC-C"

model_levels <- c(label_matrix, label_vecp, label_vecc)

model_colors <- c(
  "Matrix MF-TPRF" = "#1b6ca8",
  "VEC-P"          = "#55a868",
  "VEC-C"          = "#c44e52"
)

# solid = Matrix
# "11"   = dots very close to each other
# "33"   = short dashes, close but still distinguishable
model_linetypes <- c(
  "Matrix MF-TPRF" = "solid",
  "VEC-P"          = "11",
  "VEC-C"          = "33"
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

# ==============================================================================
# SINGLE-COUNTRY PLOT
# ==============================================================================

country_focus <- "PT"
country_focus_label <- unname(country_labels[country_focus])

df_at_matrix <- extract_country_rt(
  summary_matrix$df_rt_all,
  country_code = country_focus,
  model_label  = label_matrix
)

df_at_vector <- extract_country_rt(
  summary_vector$df_rt_all,
  country_code = country_focus,
  model_label  = label_vecc
)

df_at_vectensor <- extract_country_rt(
  summary_vectensor$df_rt_all,
  country_code = country_focus,
  model_label  = label_vecp
)

df_at_all <- bind_rows(
  df_at_matrix,
  df_at_vectensor,
  df_at_vector
) %>%
  mutate(
    model = factor(model, levels = model_levels),
    type  = factor(type, levels = c("M1", "M2", "M3"))
  )

df_at_gdp <- extract_country_gdp(
  summary_matrix$df_yq_eval_all,
  country_code = country_focus
)

df_at_gdp_facet <- tidyr::crossing(
  df_at_gdp,
  type = factor(c("M1", "M2", "M3"), levels = c("M1", "M2", "M3"))
)

plot_at_models_facet <- ggplot() +
  annotate(
    "rect",
    xmin = params$covid_start, xmax = params$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey72", alpha = 0.11
  ) +
  geom_line(
    data = df_at_gdp_facet,
    aes(x = date, y = GDP, group = type),
    colour    = "black",
    linewidth = 1.25,
    linetype  = "solid",
    lineend   = "round"
  ) +
  geom_line(
    data = df_at_all,
    aes(x = date, y = nowcast, colour = model, linetype = model),
    linewidth = 1.00,
    alpha     = 0.98,
    lineend   = "round"
  ) +
  facet_wrap(~ type, ncol = 3) +
  scale_color_manual(
    values = model_colors,
    breaks = model_levels,
    name   = NULL
  ) +
  scale_linetype_manual(
    values = model_linetypes,
    breaks = model_levels,
    name   = NULL
  ) +
  scale_x_date(
    breaks = seq(
      floor_date(min(df_at_all$date, na.rm = TRUE), unit = "year"),
      floor_date(max(df_at_all$date, na.rm = TRUE), unit = "year"),
      by = "2 years"
    ),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title = paste0("Rolling pseudo-real-time nowcasts for ", country_focus_label),
    x     = NULL,
    y     = "GDP growth"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.28),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.28),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y        = element_text(size = 10),
    strip.text         = element_text(face = "bold", size = 11),
    strip.background   = element_rect(fill = "grey94", colour = "grey80", linewidth = 0.3),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 15),
    axis.title.y       = element_text(size = 12),
    plot.margin        = margin(8, 10, 8, 8)
  )

print(plot_at_models_facet)

ggsave(
  filename = file.path(
    path_final_results,
    paste0(
      "plot_", country_focus,
      "_models_by_month_paper_Size-", Size,
      "_sel-", sel, ".png"
    )
  ),
  plot   = plot_at_models_facet,
  width  = 12.5,
  height = 4.8,
  dpi    = 500,
  bg     = "white"
)

# ==============================================================================
# MULTI-COUNTRY PLOT
# ==============================================================================

countries_focus <- c("DE", "FR", "IT", "ES")
# countries_focus <- c("NL", "BE", "AT", "PT")

extract_country_rt_multi <- function(df_rt, country_codes, model_label) {
  df_rt %>%
    filter(country %in% country_codes) %>%
    mutate(
      model   = model_label,
      type    = factor(type, levels = c("M1", "M2", "M3")),
      country = factor(country, levels = country_codes)
    ) %>%
    select(date, country, nowcast, type, model)
}

extract_country_gdp_multi <- function(df_yq, country_codes) {
  df_yq %>%
    filter(country %in% country_codes) %>%
    mutate(country = factor(country, levels = country_codes)) %>%
    select(date, country, GDP)
}

df_multi_matrix <- extract_country_rt_multi(
  summary_matrix$df_rt_all,
  country_codes = countries_focus,
  model_label   = label_matrix
)

df_multi_vector <- extract_country_rt_multi(
  summary_vector$df_rt_all,
  country_codes = countries_focus,
  model_label   = label_vecc
)

df_multi_vectensor <- extract_country_rt_multi(
  summary_vectensor$df_rt_all,
  country_codes = countries_focus,
  model_label   = label_vecp
)

df_multi_all <- bind_rows(
  df_multi_matrix,
  df_multi_vectensor,
  df_multi_vector
) %>%
  mutate(
    model   = factor(model, levels = model_levels),
    type    = factor(type, levels = c("M1", "M2", "M3")),
    country = factor(country, levels = countries_focus)
  )

df_multi_gdp <- extract_country_gdp_multi(
  summary_matrix$df_yq_eval_all,
  country_codes = countries_focus
)

df_multi_gdp_facet <- tidyr::crossing(
  df_multi_gdp,
  type = factor(c("M1", "M2", "M3"), levels = c("M1", "M2", "M3"))
)

countries_tag  <- paste(countries_focus, collapse = "-")
countries_text <- paste(countries_focus, collapse = ", ")

plot_multi_models_facet <- ggplot() +
  annotate(
    "rect",
    xmin = params$covid_start, xmax = params$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey72", alpha = 0.11
  ) +
  geom_line(
    data = df_multi_gdp_facet,
    aes(x = date, y = GDP, group = interaction(country, type)),
    colour    = "black",
    linewidth = 1.20,
    linetype  = "solid",
    lineend   = "round"
  ) +
  geom_line(
    data = df_multi_all,
    aes(x = date, y = nowcast, colour = model, linetype = model),
    linewidth = 0.95,
    alpha     = 0.99,
    lineend   = "round"
  ) +
  facet_grid(
    country ~ type,
    scales = "fixed"
  ) +
  scale_color_manual(
    values = model_colors,
    breaks = model_levels,
    name   = NULL
  ) +
  scale_linetype_manual(
    values = model_linetypes,
    breaks = model_levels,
    name   = NULL
  ) +
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
    title    = "Rolling pseudo-real-time nowcasts",
    subtitle = paste0("Size = ", Size, ", sel = ", sel, " | Countries: ", countries_text),
    x        = NULL,
    y        = "GDP growth"
  ) +
  theme_minimal(base_size = 10.4) +
  theme(
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.text        = element_text(size = 8.8),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.22),
    panel.grid.major.y = element_line(colour = "grey91", linewidth = 0.22),
    panel.spacing.x    = unit(0.50, "lines"),
    panel.spacing.y    = unit(0.35, "lines"),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.6),
    axis.text.y        = element_text(size = 7.6),
    strip.text.x       = element_text(face = "bold", size = 9.4),
    strip.text.y       = element_text(
      face   = "bold",
      size   = 9.2,
      angle  = 270,
      margin = margin(1.5, 2.8, 1.5, 2.8)
    ),
    strip.background   = element_rect(fill = "grey91", colour = "grey68", linewidth = 0.40),
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle      = element_text(hjust = 0.5, size = 9.5, colour = "grey20"),
    axis.title.y       = element_text(size = 9.6),
    plot.margin        = margin(6, 8, 6, 6)
  )

print(plot_multi_models_facet)

file_multi_png <- file.path(
  path_final_results,
  paste0(
    "mcplot_",
    "Size-", Size,
    "_sel-", sel,
    "_cty-", countries_tag,
    ".png"
  )
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
# 8. MODEL-SPECIFIC LATEX TABLES
# ==============================================================================

cat("\n================ MATRIX IN-SAMPLE ================\n")
cat(summary_matrix$latex_tab_insample_all)

cat("\n================ MATRIX ROLLING ================\n")
cat(summary_matrix$latex_tab_rt_all)

cat("\n================ VECTOR IN-SAMPLE ================\n")
cat(summary_vector$latex_tab_insample_all)

cat("\n================ VECTOR ROLLING ================\n")
cat(summary_vector$latex_tab_rt_all)

cat("\n================ VECTENSOR IN-SAMPLE ================\n")
cat(summary_vectensor$latex_tab_insample_all)

cat("\n================ VECTENSOR ROLLING ================\n")
cat(summary_vectensor$latex_tab_rt_all)

# ==============================================================================
# 9. COMPARISON TABLES
# ==============================================================================

comp_insample <- build_comparison_table(
  df_matrix_insample,
  df_vector_insample,
  df_vectensor_insample
)

comp_rt <- build_comparison_table(
  df_matrix_rt,
  df_vector_rt,
  df_vectensor_rt
)

latex_comp_insample <- comparison_to_latex(
  comp_insample,
  caption = paste0(
    "In-sample RMSFE comparison across models (",
    Size, ", sel = ", sel, ")"
  ),
  label = paste0("tab:comparison_insample_models_", Size, "_", sel)
)

latex_comp_rt <- comparison_to_latex(
  comp_rt,
  caption = paste0(
    "Rolling pseudo-real-time RMSFE comparison across models (",
    Size, ", sel = ", sel, ")"
  ),
  label = paste0("tab:comparison_rt_models_", Size, "_", sel)
)

cat("\n================ COMPARISON IN-SAMPLE ================\n")
print(comp_insample)
cat("\n", latex_comp_insample, "\n")

cat("\n================ COMPARISON ROLLING ================\n")
print(comp_rt)
cat("\n", latex_comp_rt, "\n")

# ==============================================================================
# 10. RELATIVE RMSFE TABLES
# ==============================================================================

rel_insample <- build_relative_table(comp_insample)
rel_rt       <- build_relative_table(comp_rt)

latex_rel_insample <- relative_to_latex(
  rel_insample,
  caption = paste0(
    "Relative in-sample RMSFE across models (",
    Size, ", sel = ", sel, ")"
  ),
  label = paste0("tab:relative_insample_models_", Size, "_", sel)
)

latex_rel_rt <- relative_to_latex(
  rel_rt,
  caption = paste0(
    "Relative rolling pseudo-real-time RMSFE across models (",
    Size, ", sel = ", sel, ")"
  ),
  label = paste0("tab:relative_rt_models_", Size, "_", sel)
)

cat("\n================ RELATIVE IN-SAMPLE ================\n")
print(rel_insample)
cat("\n", latex_rel_insample, "\n")

cat("\n================ RELATIVE ROLLING ================\n")
print(rel_rt)
cat("\n", latex_rel_rt, "\n")

# ==============================================================================
# 11. DM TEST
# ==============================================================================

eval_matrix <- build_eval_df(summary_matrix$df_rt_all, summary_matrix$df_yq_eval_all)
eval_vecp   <- build_eval_df(summary_vectensor$df_rt_all, summary_matrix$df_yq_eval_all)
eval_vecc   <- build_eval_df(summary_vector$df_rt_all, summary_matrix$df_yq_eval_all)

dm_base_vecp <- eval_matrix %>%
  select(country, quarter_id, type, period, se_matrix = se) %>%
  inner_join(
    eval_vecp %>%
      select(country, quarter_id, type, se_vecp = se),
    by = c("country", "quarter_id", "type")
  ) %>%
  filter(!is.na(se_matrix), !is.na(se_vecp)) %>%
  mutate(d = se_matrix - se_vecp)

dm_base_vecc <- eval_matrix %>%
  select(country, quarter_id, type, period, se_matrix = se) %>%
  inner_join(
    eval_vecc %>%
      select(country, quarter_id, type, se_vecc = se),
    by = c("country", "quarter_id", "type")
  ) %>%
  filter(!is.na(se_matrix), !is.na(se_vecc)) %>%
  mutate(d = se_matrix - se_vecc)

dm_vecp <- dm_base_vecp %>%
  group_by(country, period, type) %>%
  group_modify(~ dm_test_hac(.x$d, alternative = "less")) %>%
  ungroup()

dm_vecc <- dm_base_vecc %>%
  group_by(country, period, type) %>%
  group_modify(~ dm_test_hac(.x$d, alternative = "less")) %>%
  ungroup()

dm_vecp_wide <- dm_vecp %>%
  select(country, period, type, p_value) %>%
  pivot_wider(
    names_from   = type,
    values_from  = p_value,
    names_prefix = "p_VecP_"
  )

dm_vecc_wide <- dm_vecc %>%
  select(country, period, type, p_value) %>%
  pivot_wider(
    names_from   = type,
    values_from  = p_value,
    names_prefix = "p_VecC_"
  )

# ==============================================================================
# 12. FINAL LATEX TABLE
# ==============================================================================

latex_final_custom <- build_final_large_style_latex(
  summary_matrix    = summary_matrix,
  summary_vector    = summary_vector,
  summary_vectensor = summary_vectensor,
  params            = params,
  dm_vecp_wide      = dm_vecp_wide,
  dm_vecc_wide      = dm_vecc_wide,
  Size              = Size,
  sel               = sel
)

cat(latex_final_custom)

# ==============================================================================
# 13. SAVE FINAL OUTPUTS
# ==============================================================================

suffix_out <- paste0("Size-", Size, "_sel-", sel)

writeLines(
  summary_matrix$latex_tab_insample_all,
  file.path(path_final_results, paste0("matrix_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_matrix$latex_tab_rt_all,
  file.path(path_final_results, paste0("matrix_rt_", suffix_out, ".tex"))
)

writeLines(
  summary_vector$latex_tab_insample_all,
  file.path(path_final_results, paste0("vector_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_vector$latex_tab_rt_all,
  file.path(path_final_results, paste0("vector_rt_", suffix_out, ".tex"))
)

writeLines(
  summary_vectensor$latex_tab_insample_all,
  file.path(path_final_results, paste0("vectensor_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_vectensor$latex_tab_rt_all,
  file.path(path_final_results, paste0("vectensor_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_comp_insample,
  file.path(path_final_results, paste0("comparison_insample_", suffix_out, ".tex"))
)
writeLines(
  latex_comp_rt,
  file.path(path_final_results, paste0("comparison_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_rel_insample,
  file.path(path_final_results, paste0("relative_insample_", suffix_out, ".tex"))
)
writeLines(
  latex_rel_rt,
  file.path(path_final_results, paste0("relative_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_final_custom,
  file.path(path_final_results, paste0("final_custom_table_", suffix_out, ".tex"))
)

if (!is.null(latex_factor_fit)) {
  writeLines(
    latex_factor_fit,
    file.path(path_final_results, paste0("static_factor_fit_", suffix_out, ".tex"))
  )
}

saveRDS(
  list(
    Size                   = Size,
    sel                    = sel,
    run_factor_analysis    = run_factor_analysis,
    
    file_matrix_summary    = file_matrix_summary,
    file_vector_summary    = file_vector_summary,
    file_vector_summary_ea = file_vector_summary_ea,
    file_vectensor_summary = file_vectensor_summary,
    
    summary_matrix         = summary_matrix,
    summary_vector         = summary_vector,
    summary_vector_ea      = summary_vector_ea,
    summary_vectensor      = summary_vectensor,
    
    fit_static_factors_long   = fit_static_factors_long,
    fit_static_factors_R2     = fit_static_factors_R2,
    fit_static_factors_AdjR2  = fit_static_factors_AdjR2,
    latex_factor_fit          = latex_factor_fit,
    
    df_factor_compare         = df_factor_compare,
    df_factor_compare_long    = df_factor_compare_long,
    plot_factor_compare       = plot_factor_compare,
    file_graph_factor_compare = file_graph_factor_compare,
    
    df_matrix_insample     = df_matrix_insample,
    df_matrix_rt           = df_matrix_rt,
    df_vector_insample     = df_vector_insample,
    df_vector_rt           = df_vector_rt,
    df_vectensor_insample  = df_vectensor_insample,
    df_vectensor_rt        = df_vectensor_rt,
    
    comp_insample          = comp_insample,
    comp_rt                = comp_rt,
    rel_insample           = rel_insample,
    rel_rt                 = rel_rt,
    
    dm_vecp                = dm_vecp,
    dm_vecc                = dm_vecc,
    dm_vecp_wide           = dm_vecp_wide,
    dm_vecc_wide           = dm_vecc_wide,
    
    latex_comp_insample    = latex_comp_insample,
    latex_comp_rt          = latex_comp_rt,
    latex_rel_insample     = latex_rel_insample,
    latex_rel_rt           = latex_rel_rt,
    latex_final_custom     = latex_final_custom
  ),
  file.path(path_final_results, paste0("final_model_comparison_", suffix_out, ".rds"))
)

cat("\nAll final outputs saved to:\n", path_final_results, "\n")

