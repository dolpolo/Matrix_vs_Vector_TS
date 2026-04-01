# ==============================================================================
# Matrix MF-TPRF - Results Script
# ------------------------------------------------------------------------------
# This script:
#   1. Loads full-sample and pseudo real-time results
#   2. Extracts true GDP and model nowcasts
#   3. Produces full-sample and rolling plots
#   4. Computes RMSFE by country and period
#   5. Saves graphs and a compact summary object
# ==============================================================================

# ==============================================================================
# 0. PATHS AND PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

path_main     <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
file_out      <- file.path(path_main, "TPRF_Models_EA/Final_Tab_Graph")
path_results  <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_graph    <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/graph_tensor")
path_graph_rt <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/graph_tensor_RT")
path_func     <- file.path(path_main, "functions/functions_mat")

source(file.path(path_func, "matrix.mf.tprf.utils.R"))

dir.create(path_graph,    recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_rt, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. LOAD RESULTS
# ==============================================================================

# Tags of the run you want to analyze
model_name <- "matrix"
Size       <- "large"   # "small" | "medium" | "large"
sel        <- "LASSO"    # "corr" | "LASSO"

file_fit <- find_result_file(
  path  = path_results,
  model = model_name,
  stage = "fit",
  Size  = Size,
  sel   = sel
)

file_rt <- find_result_file(
  path  = path_results,
  model = model_name,
  stage = "rt",
  Size  = Size,
  sel   = sel
)

res_fit <- readRDS(file_fit)
res_rt  <- readRDS(file_rt)

cat("\nLoaded full-sample results from:\n", file_fit, "\n")
cat("\nLoaded pseudo real-time results from:\n", file_rt, "\n")

# Overwrite tags from saved objects if available
model_name <- if (!is.null(res_fit$model)) res_fit$model else model_name
Size       <- if (!is.null(res_fit$Size))  res_fit$Size  else Size
sel        <- if (!is.null(res_fit$sel))   res_fit$sel   else sel

# ==============================================================================
# 2. EXTRACT CORE OBJECTS
# ==============================================================================

params <- res_fit$params
tensor <- res_fit$tensor

df_selection <- res_fit$selections
df_selection_wide <- res_fit$selections_wide
sel_raw <- res_fit$sel_raw
na_pct <- res_fit$na_pct

fit_obj <- if (!is.null(res_fit$fit)) {
  res_fit$fit
} else if (!is.null(res_fit$Tensor_MF_TPRF)) {
  res_fit$Tensor_MF_TPRF
} else {
  stop("No full-sample fit object found in res_fit.")
}

countries_fit <- names(fit_obj$by_country)
country_order <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
countries_eval <- intersect(country_order, countries_fit)

hyper_fit <- if (!is.null(res_fit$hyper)) {
  res_fit$hyper
} else {
  list(
    Lproxy     = res_fit$Lproxy,
    L_midas    = res_fit$L_midas,
    p_ar       = res_fit$p_ar,
    r_selected = res_fit$r_selected
  )
}

Lproxy  <- hyper_fit$Lproxy
L_midas <- hyper_fit$L_midas
p_ar    <- hyper_fit$p_ar
r1      <- if (!is.null(hyper_fit$r_selected)) hyper_fit$r_selected[1] else NA
r2      <- if (!is.null(hyper_fit$r_selected)) hyper_fit$r_selected[2] else NA

N_m <- if (!is.null(res_fit$metadata$N_m)) res_fit$metadata$N_m else NA
N_q <- if (!is.null(res_fit$metadata$N_q)) res_fit$metadata$N_q else NA

# ==============================================================================
# 3. TRUE QUARTERLY GDP
# ==============================================================================

Y_tens  <- tensor$Y
gdp_col <- tensor$target_col
dates_m <- as.Date(dimnames(Y_tens)[[1]])

gdp_tens <- Y_tens[, , gdp_col, drop = FALSE]
gdp_mat  <- matrix(
  gdp_tens[, , 1, drop = TRUE],
  nrow = dim(gdp_tens)[1],
  ncol = dim(gdp_tens)[2]
)
colnames(gdp_mat) <- dimnames(Y_tens)[[2]]

idx_q   <- which(rowSums(!is.na(gdp_mat)) > 0)
dates_q <- dates_m[idx_q]

y_true_q_all <- gdp_mat[idx_q, countries_eval, drop = FALSE]
T_q_complete <- nrow(y_true_q_all)

# ==============================================================================
# 4. FULL-SAMPLE MONTHLY NOWCAST
# ==============================================================================

stopifnot(!is.null(fit_obj$by_country))

y_now_by_country <- lapply(countries_eval, function(cc) {
  out_cc <- fit_obj$by_country[[cc]]
  if (is.null(out_cc$y_nowcast)) stop("Missing y_nowcast for country: ", cc)
  as.numeric(out_cc$y_nowcast)
})
names(y_now_by_country) <- countries_eval

T_m_full <- length(y_now_by_country[[1]])
if (any(vapply(y_now_by_country, length, 1L) != T_m_full)) {
  stop("Monthly nowcast lengths differ across countries.")
}

dates_m_full <- dates_m[seq_len(T_m_full)]

n_in <- 3 * T_q_complete
if (T_m_full < n_in) stop("Not enough months: T_m_full < 3 * T_q_complete.")

period_vec <- rep("in-sample", T_m_full)
if (T_m_full > n_in) period_vec[(n_in + 1):T_m_full] <- "real-time"
period_vec <- factor(period_vec, levels = c("in-sample", "real-time"))

M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

# ==============================================================================
# 5. PERIOD MASKS
# ==============================================================================

start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

is_PRE   <- dates_q >= start_eval  & dates_q <  covid_start
is_COVID <- dates_q >= covid_start & dates_q <= covid_end
is_POST  <- dates_q >  covid_end   & dates_q <= end_eval
is_ALL   <- dates_q >= start_eval  & dates_q <= end_eval

# ==============================================================================
# 6. DATA FRAMES FOR FULL-SAMPLE ANALYSIS
# ==============================================================================

df_now_full <- bind_rows(lapply(countries_eval, function(cc) {
  data.frame(
    date       = dates_m_full,
    country    = cc,
    y_now_full = y_now_by_country[[cc]],
    period     = period_vec
  )
}))

df_quarterly <- data.frame(
  date    = rep(dates_q, times = length(countries_eval)),
  country = rep(countries_eval, each = T_q_complete),
  y_true  = as.vector(y_true_q_all)
)

df_now_full$country  <- factor(df_now_full$country,  levels = countries_eval)
df_quarterly$country <- factor(df_quarterly$country, levels = countries_eval)

# ==============================================================================
# 7. FULL-SAMPLE PLOT
# ==============================================================================

plot_nowcast <- ggplot() +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.15
  ) +
  geom_line(
    data = df_now_full,
    aes(x = date, y = y_now_full, color = period),
    linewidth = 0.9
  ) +
  geom_point(
    data = df_now_full %>% filter(period == "real-time"),
    aes(x = date, y = y_now_full),
    size = 1.5, color = "black"
  ) +
  geom_line(
    data = df_quarterly,
    aes(x = date, y = y_true, color = "Quarterly GDP"),
    linewidth = 1.0
  ) +
  scale_color_manual(
    values = c(
      "in-sample"     = "#1F77B4",
      "real-time"     = "#2ECC71",
      "Quarterly GDP" = "#D62728"
    ),
    name = ""
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = "Matrix MF-TPRF Monthly Nowcast",
    x = "Date",
    y = "GDP Growth"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(plot_nowcast)

# ==============================================================================
# 8. IN-SAMPLE RMSFE
# ==============================================================================

rmsfe_list <- list()

for (cc in countries_eval) {
  y_true_cc <- y_true_q_all[, cc]
  y_now_cc  <- y_now_by_country[[cc]]
  y_now_in  <- y_now_cc[1:n_in]
  
  rmsfe_list[[cc]] <- data.frame(
    country = cc,
    period  = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID"),
    M1 = c(
      rmsfe_period(is_ALL,   y_true_cc, y_now_in, M1_idx),
      rmsfe_period(is_PRE,   y_true_cc, y_now_in, M1_idx),
      rmsfe_period(is_COVID, y_true_cc, y_now_in, M1_idx),
      rmsfe_period(is_POST,  y_true_cc, y_now_in, M1_idx)
    ),
    M2 = c(
      rmsfe_period(is_ALL,   y_true_cc, y_now_in, M2_idx),
      rmsfe_period(is_PRE,   y_true_cc, y_now_in, M2_idx),
      rmsfe_period(is_COVID, y_true_cc, y_now_in, M2_idx),
      rmsfe_period(is_POST,  y_true_cc, y_now_in, M2_idx)
    ),
    M3 = c(
      rmsfe_period(is_ALL,   y_true_cc, y_now_in, M3_idx),
      rmsfe_period(is_PRE,   y_true_cc, y_now_in, M3_idx),
      rmsfe_period(is_COVID, y_true_cc, y_now_in, M3_idx),
      rmsfe_period(is_POST,  y_true_cc, y_now_in, M3_idx)
    )
  )
}

rmsfe_insample <- bind_rows(rmsfe_list) %>%
  mutate(country = factor(country, levels = countries_eval)) %>%
  arrange(country, period)

latex_tab_insample <- list_to_latex_table(
  df      = rmsfe_insample,
  caption = paste0("Matrix MF-TPRF in-sample RMSFE by country and period (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:matrix_mf_tprf_insample_", Size, "_", sel)
)

cat("\n", latex_tab_insample, "\n")

dir.create(file_out, recursive = TRUE, showWarnings = FALSE)

latex_selection <- save_selection_wide_to_latex(
  df_selection_wide = df_selection_wide,
  file    = file.path(file_out, "selection_wide_table.tex"),
  caption = "Country-specific variable selection",
  label   = "tab:selection_wide"
)

cat(latex_selection)
# ==============================================================================
# 9. SAVE FULL-SAMPLE GRAPH
# ==============================================================================

file_graph_now <- file.path(
  path_graph,
  paste0(
    "plot_fit_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-", p_ar,
    "_r1-", r1,
    "_r2-", r2,
    ".png"
  )
)

ggsave(
  filename = file_graph_now,
  plot     = plot_nowcast,
  width    = 12,
  height   = 8,
  dpi      = 300
)

cat("\nSaved full-sample graph to:\n", file_graph_now, "\n")

# ==============================================================================
# 10. BUILD PSEUDO REAL-TIME DATA
# ==============================================================================

params_rt <- res_rt$params

df_rt <- res_rt$pseudo_rt_all %>%
  mutate(
    date    = as.Date(date),
    type    = factor(month_in_quarter, levels = c("M1", "M2", "M3")),
    country = factor(country)
  ) %>%
  select(date, country, nowcast, type) %>%
  arrange(country, date, type)

countries_rt <- unique(as.character(df_rt$country))
countries_eval_rt <- intersect(country_order, countries_rt)

Y_q_all_rt <- as.matrix(res_rt$Y_q_all)
dates_q_rt <- as.Date(res_rt$dates_q)

idx_eval <- which(dates_q_rt >= params_rt$start_eval & dates_q_rt <= params_rt$end_eval)

df_yq <- data.frame(
  date    = rep(dates_q_rt[idx_eval], times = length(countries_eval_rt)),
  country = rep(countries_eval_rt, each = length(idx_eval)),
  GDP     = as.vector(Y_q_all_rt[idx_eval, countries_eval_rt, drop = FALSE])
)

df_rt <- df_rt %>%
  filter(country %in% countries_eval_rt) %>%
  mutate(country = factor(as.character(country), levels = countries_eval_rt))

df_yq$country <- factor(df_yq$country, levels = countries_eval_rt)

# ==============================================================================
# 11. PSEUDO REAL-TIME PLOT
# ==============================================================================

plot_rt <- ggplot() +
  annotate(
    "rect",
    xmin = params_rt$covid_start, xmax = params_rt$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  geom_line(
    data = df_yq,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.0
  ) +
  geom_line(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.8
  ) +
  geom_point(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    size = 1.6
  ) +
  scale_color_manual(
    values = c(
      "True GDP" = "#D62728",
      "M1"       = "#1F77B4",
      "M2"       = "#2ECC71",
      "M3"       = "#F1C40F"
    ),
    name = "Series"
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title    = "Matrix MF-TPRF Rolling Real-Time Nowcasts",
    subtitle = "M1: early quarter, M2: mid quarter, M3: end quarter",
    x        = "Date",
    y        = "GDP Growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(plot_rt)

# ==============================================================================
# 12. ROLLING RMSFE
# ==============================================================================

df_yq_eval <- df_yq %>%
  mutate(
    period = case_when(
      date <  params_rt$covid_start                               ~ "Pre-COVID",
      date >= params_rt$covid_start & date <= params_rt$covid_end ~ "COVID period",
      date >  params_rt$covid_end                                 ~ "Post-COVID",
      TRUE                                                        ~ NA_character_
    ),
    quarter_id = paste0(year(date), "Q", quarter(date))
  )

df_rt_eval <- df_rt %>%
  mutate(quarter_id = paste0(year(date), "Q", quarter(date)))

df_eval <- df_rt_eval %>%
  inner_join(
    df_yq_eval %>% select(country, quarter_id, GDP, period),
    by = c("country", "quarter_id")
  ) %>%
  arrange(country, date)

rmsfe_by_period <- df_eval %>%
  filter(!is.na(period)) %>%
  group_by(country, period, type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  )

rmsfe_full <- df_eval %>%
  group_by(country, type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(period = "Full sample")

rmsfe_rt <- bind_rows(rmsfe_full, rmsfe_by_period) %>%
  mutate(
    period  = factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")),
    type    = factor(type, levels = c("M1", "M2", "M3")),
    country = factor(country, levels = countries_eval_rt)
  ) %>%
  arrange(country, period, type)

rmsfe_rt_wide <- rmsfe_rt %>%
  pivot_wider(names_from = type, values_from = RMSFE) %>%
  arrange(country, period)

latex_tab_rt <- list_to_latex_table(
  df      = rmsfe_rt_wide,
  caption = paste0("Matrix MF-TPRF rolling RMSFE by country and period (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:matrix_mf_tprf_rt_", Size, "_", sel)
)

cat("\n", latex_tab_rt, "\n")

# ==============================================================================
# 13. SAVE ROLLING GRAPH
# ==============================================================================

hyper_rt  <- if (!is.null(res_rt$hyper)) res_rt$hyper else res_rt$pseudo_rt_raw
Lproxy_rt <- hyper_rt$pre$Lproxy
Lmidas_rt <- hyper_rt$pre$L_midas
pAR_rt    <- hyper_rt$pre$p_AR

file_graph_rt <- file.path(
  path_graph_rt,
  paste0(
    "plot_rt_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    "_Lproxy-", Lproxy_rt,
    "_Lmidas-", Lmidas_rt,
    "_pAR-", pAR_rt,
    "_", format(params_rt$start_eval, "%Y-%m"),
    "_to_", format(params_rt$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(
  filename = file_graph_rt,
  plot     = plot_rt,
  width    = 12,
  height   = 8,
  dpi      = 300
)

cat("\nSaved rolling graph to:\n", file_graph_rt, "\n")

# ==============================================================================
# 14. SAVE SUMMARY OBJECT
# ==============================================================================

tab_insample_all <- rmsfe_insample %>%
  mutate(
    country = as.character(country),
    period  = as.character(period)
  ) %>%
  select(country, period, M1, M2, M3)

tab_rt_all <- rmsfe_rt_wide %>%
  mutate(
    country = as.character(country),
    period  = as.character(period)
  ) %>%
  select(country, period, M1, M2, M3)

hyper_summary <- if (!is.null(res_rt$hyper)) {
  res_rt$hyper
} else {
  list(pre = NULL, post = NULL)
}

summary_tensor_out <- list(
  model_id               = "T_MF_TPRF",
  model                  = model_name,
  stage                  = "cross_country",
  Size                   = Size,
  sel                    = sel,
  params                 = params_rt,
  
  hyper = hyper_summary,
  
  countries              = countries_eval_rt,
  
  df_now_full_all        = df_now_full,
  df_quarterly_all       = df_quarterly,
  plot_nowcast_facet     = plot_nowcast,
  tab_insample_all       = tab_insample_all,
  latex_tab_insample_all = latex_tab_insample,
  
  df_rt_all              = df_rt,
  df_yq_eval_all         = df_yq_eval,
  plot_rt_facet          = plot_rt,
  tab_rt_all             = tab_rt_all,
  latex_tab_rt_all       = latex_tab_rt,
  
  selection_details = df_selection,
  selection_wide    = df_selection_wide,
  selection_raw     = sel_raw,
  na_pct            = na_pct,
  
  file_fit               = file_fit,
  file_rt                = file_rt,
  file_graph_fit         = file_graph_now,
  file_graph_rt          = file_graph_rt
)

file_summary_matrix <- build_result_filename(
  path_out         = path_results,
  model            = model_name,
  stage            = "summary",
  Size             = Size,
  sel              = sel,
  countries        = countries_eval_rt,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = Lproxy_rt,
  L_midas          = Lmidas_rt,
  p_ar             = pAR_rt,
  r1               = r1,
  r2               = r2,
  robust_f         = as.integer(isTRUE(params_rt$Robust_F)),
  covid_m          = as.integer(isTRUE(params_rt$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params_rt$covid_mask_q)),
  ext              = "rds",
  timestamp        = FALSE,
  include_details  = TRUE
)

if (file.exists(file_summary_matrix)) file.remove(file_summary_matrix)

saveRDS(summary_tensor_out, file_summary_matrix)

cat("\nSaved summary object to:\n", file_summary_matrix, "\n")
