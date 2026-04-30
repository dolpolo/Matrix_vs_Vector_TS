# ==============================================================================
# Vector DFM - Results Script
# ------------------------------------------------------------------------------
# This script:
#   1. Loads full-sample and pseudo real-time DFM results
#   2. Extracts true GDP and nowcasts
#   3. Produces full-sample diagnostics and rolling plots
#   4. Computes RMSFE by period
#   5. Saves graphs and a compact summary object
# ==============================================================================

# ==============================================================================
# 0. PATHS AND PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
path_results <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/outputs")
path_graph   <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/graph")

country <- "FR"

path_country       <- file.path(path_results, country)
path_graph_country <- file.path(path_graph, country)

dir.create(path_graph_country, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. TAGS OF THE RUN TO LOAD
# ==============================================================================

model_name <- "vector_dfm"
Size       <- "small"  
sel        <- "LASSO"

# ==============================================================================
# 3. LOAD FULL-SAMPLE RESULTS
# ==============================================================================

file_fit <- find_result_file_dfm(
  path    = path_country,
  stage   = "fit",
  model   = model_name,
  sel     = sel,
  country = country
)

fit_res <- readRDS(file_fit)

cat("\nLoaded full-sample DFM results from:\n", file_fit, "\n")

model_name <- if (!is.null(fit_res$model)) fit_res$model else model_name
Size       <- if (!is.null(fit_res$Size))  fit_res$Size  else Size
sel        <- if (!is.null(fit_res$sel))   fit_res$sel   else sel
country    <- if (!is.null(fit_res$country)) fit_res$country else country

# ==============================================================================
# 4. EXTRACT FULL-SAMPLE OBJECTS
# ==============================================================================

params     <- extract_params_object(fit_res)
meta_obj   <- extract_dates_y_dfm(fit_res)
meta_data  <- extract_metadata_dfm(fit_res)
pre_obj    <- extract_preprocessing_dfm(fit_res)
hyper_fit  <- extract_hyper_dfm(fit_res)
fit_obj    <- extract_fit_dfm(fit_res)

dates_m  <- meta_obj$dates_m
dates_q  <- meta_obj$dates_q
y_true_q <- meta_obj$y_q

if (is.null(dates_m) || is.null(dates_q) || is.null(y_true_q)) {
  stop("The full-sample DFM RDS must contain dates_m, dates_q, and y_q.")
}

N_m <- if (!is.null(meta_data$N_m)) meta_data$N_m else NA
N_q <- if (!is.null(meta_data$N_q)) meta_data$N_q else NA
N   <- if (!is.null(meta_data$N))   meta_data$N   else NA

gdp_col     <- if (!is.null(meta_data$gdp_col)) meta_data$gdp_col else NA_integer_
target_name <- if (!is.null(meta_data$target_name)) meta_data$target_name else NA_character_

r_sel <- hyper_fit$r
p_sel <- hyper_fit$p
q_sel <- hyper_fit$q

start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

X_raw   <- pre_obj$X_raw
out_std <- pre_obj$out_std
X_std   <- pre_obj$X_std

if (is.null(X_raw) || is.null(X_std) || is.null(out_std)) {
  stop("Full-sample DFM file must contain preprocessing$X_raw, preprocessing$X_std and preprocessing$out_std.")
}

mu_full <- out_std$mean
sd_full <- out_std$sd

# ==============================================================================
# 5. FULL-SAMPLE IN-SAMPLE FITTED GDP SERIES
# ==============================================================================
# We recover the smoother fit saved in fit$X_fit_orig if present.
# Otherwise we reconstruct it from fit_obj$fit_obj / fit_obj$res.
# ==============================================================================

X_fit_std  <- fit_obj$X_fit_std
X_fit_orig <- fit_obj$X_fit_orig
fit_kalman <- fit_obj$fit_obj

if (is.null(X_fit_std) || is.null(X_fit_orig)) {
  if (!exists("dfm_kalman_fit")) stop("dfm_kalman_fit() is required to reconstruct fitted values.")
  
  res_full <- if (!is.null(fit_obj$res)) fit_obj$res else NULL
  if (is.null(res_full)) {
    res_full <- list(
      A  = fit_obj$A,
      C  = fit_obj$C,
      Q  = fit_obj$Q,
      R  = fit_obj$R,
      Z0 = fit_obj$Z0,
      V0 = fit_obj$V0
    )
  }
  
  n_lags_Q_used <- if (!is.null(fit_obj$n_lags_Q_used)) fit_obj$n_lags_Q_used else 3L
  
  fit_kalman <- dfm_kalman_fit(
    X_in      = X_std,
    res       = res_full,
    n_lags_Q  = as.integer(n_lags_Q_used)
  )
  
  X_fit_std  <- fit_kalman$X_fit
  X_fit_orig <- destd_mat(X_fit_std, mu_full, sd_full)
}

t_idx <- if (!is.null(fit_kalman$t_index)) fit_kalman$t_index else seq_len(nrow(X_fit_std))
dates_fit <- as.Date(dates_m[t_idx])

# observed aligned blocks
Xf_std  <- X_std[t_idx, , drop = FALSE]
Xf_orig <- X_raw[t_idx, , drop = FALSE]

gdp_fit_std  <- if (!is.null(fit_obj$gdp_fit_std))  fit_obj$gdp_fit_std  else X_fit_std[,  gdp_col]
gdp_fit_orig <- if (!is.null(fit_obj$gdp_fit_orig)) fit_obj$gdp_fit_orig else X_fit_orig[, gdp_col]

# quarterly observed GDP
df_quarterly <- data.frame(
  date   = as.Date(dates_q),
  y_true = as.numeric(y_true_q)
)

# ==============================================================================
# 6. BUILD FULL-SAMPLE MONTHLY GDP FIT SERIES
# ==============================================================================
# To mimic MF-TPRF results, we construct a monthly series:
#   - fitted GDP from the smoother
#   - tag months up to last observed GDP quarter as "in-sample"
# ==============================================================================

df_now_full <- data.frame(
  date       = as.Date(dates_fit),
  y_now_full = as.numeric(gdp_fit_orig),
  period     = factor("in-sample", levels = c("in-sample", "real-time"))
)

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
    linewidth = 1.0
  ) +
  geom_line(
    data = df_quarterly,
    aes(x = date, y = y_true, color = "Quarterly GDP"),
    linewidth = 1.1
  ) +
  scale_color_manual(
    values = c(
      "in-sample"     = "#1F77B4",
      "real-time"     = "#2ECC71",
      "Quarterly GDP" = "#D62728"
    ),
    name = ""
  ) +
  labs(
    title = paste0("Vector DFM Monthly Fitted GDP - ", country),
    x = "Date",
    y = "GDP Growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(plot_nowcast)

# ==============================================================================
# 8. IN-SAMPLE RMSFE
# ==============================================================================
# We compute quarterly M1/M2/M3 from the monthly fitted GDP series exactly
# like in MF-TPRF: month 1, 2, 3 of each quarter.
# ==============================================================================

n_in <- 3 * length(y_true_q)

if (nrow(df_now_full) < n_in) {
  stop("The monthly fitted GDP series has fewer than 3*T_q observations.")
}

y_now_in <- df_now_full$y_now_full[seq_len(n_in)]

M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

is_PRE   <- dates_q >= start_eval  & dates_q <  covid_start
is_COVID <- dates_q >= covid_start & dates_q <= covid_end
is_POST  <- dates_q >  covid_end   & dates_q <= end_eval
is_ALL   <- dates_q >= start_eval  & dates_q <= end_eval

rmsfe_insample <- rbind(
  "Full sample" = c(
    rmsfe_period(is_ALL,   y_true_q, y_now_in, M1_idx),
    rmsfe_period(is_ALL,   y_true_q, y_now_in, M2_idx),
    rmsfe_period(is_ALL,   y_true_q, y_now_in, M3_idx)
  ),
  "Pre-COVID" = c(
    rmsfe_period(is_PRE,   y_true_q, y_now_in, M1_idx),
    rmsfe_period(is_PRE,   y_true_q, y_now_in, M2_idx),
    rmsfe_period(is_PRE,   y_true_q, y_now_in, M3_idx)
  ),
  "COVID period" = c(
    rmsfe_period(is_COVID, y_true_q, y_now_in, M1_idx),
    rmsfe_period(is_COVID, y_true_q, y_now_in, M2_idx),
    rmsfe_period(is_COVID, y_true_q, y_now_in, M3_idx)
  ),
  "Post-COVID" = c(
    rmsfe_period(is_POST,  y_true_q, y_now_in, M1_idx),
    rmsfe_period(is_POST,  y_true_q, y_now_in, M2_idx),
    rmsfe_period(is_POST,  y_true_q, y_now_in, M3_idx)
  )
)

colnames(rmsfe_insample) <- c("M1", "M2", "M3")

latex_insample <- latex_table_periods(
  mat     = rmsfe_insample,
  caption = paste0("Vector DFM RMSFE by period -- ", country, " (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:vector_dfm_", country, "_", Size, "_", sel)
)

cat("\n", latex_insample, "\n")

tab_insample_all <- data.frame(
  country = country,
  period  = rownames(rmsfe_insample),
  M1      = rmsfe_insample[, "M1"],
  M2      = rmsfe_insample[, "M2"],
  M3      = rmsfe_insample[, "M3"],
  row.names = NULL
)

# ==============================================================================
# 9. SAVE FULL-SAMPLE GRAPH
# ==============================================================================

file_graph_fit <- file.path(
  path_graph_country,
  paste0(
    "plot_fit_", model_name, "_", country,
    "_Size-", Size,
    "_sel-", sel,
    "_r-", r_sel,
    "_p-", p_sel,
    "_q-", q_sel,
    ".png"
  )
)

ggsave(
  filename = file_graph_fit,
  plot     = plot_nowcast,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\nSaved full-sample graph to:\n", file_graph_fit, "\n")

# ==============================================================================
# 10. LOAD PSEUDO REAL-TIME RESULTS
# ==============================================================================

file_rt <- find_result_file_dfm(
  path    = path_country,
  stage   = "rt",
  model   = model_name,
  sel     = sel,
  country = country
)

rt_res <- readRDS(file_rt)

cat("\nLoaded pseudo real-time DFM results from:\n", file_rt, "\n")

# ==============================================================================
# 11. EXTRACT PSEUDO REAL-TIME OBJECTS
# ==============================================================================

params_rt <- extract_params_object(rt_res)

meta_rt_dates_q <- if (!is.null(rt_res$dates_q)) as.Date(rt_res$dates_q) else dates_q
meta_rt_y_q     <- if (!is.null(rt_res$y_q)) as.numeric(rt_res$y_q) else y_true_q

rt_now <- extract_rt_nowcasts_dfm(rt_res)

r_fix <- rt_now$r_fix
p_fix <- rt_now$p_fix
q_fix <- rt_now$q_fix

now_dfm_orig <- list(
  M1 = rt_now$M1_orig,
  M2 = rt_now$M2_orig,
  M3 = rt_now$M3_orig
)

# build df_rt exactly like MF-TPRF
df_rt <- bind_rows(
  data.frame(date = as.Date(names(now_dfm_orig$M1)), nowcast = as.numeric(now_dfm_orig$M1), type = "M1"),
  data.frame(date = as.Date(names(now_dfm_orig$M2)), nowcast = as.numeric(now_dfm_orig$M2), type = "M2"),
  data.frame(date = as.Date(names(now_dfm_orig$M3)), nowcast = as.numeric(now_dfm_orig$M3), type = "M3")
) %>%
  mutate(type = factor(type, levels = c("M1", "M2", "M3"))) %>%
  arrange(date, type)

idx_eval <- which(meta_rt_dates_q >= params_rt$start_eval & meta_rt_dates_q <= params_rt$end_eval)

df_yq <- data.frame(
  date = meta_rt_dates_q[idx_eval],
  GDP  = meta_rt_y_q[idx_eval]
)

# ==============================================================================
# 12. PSEUDO REAL-TIME PLOT
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
    linewidth = 1.1
  ) +
  geom_line(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.9
  ) +
  geom_point(
    data = df_rt,
    aes(x = date, y = nowcast, color = type),
    size = 1.8
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
  labs(
    title    = paste0("Vector DFM Rolling Real-Time Nowcasts - ", country),
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
# 13. ROLLING RMSFE
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
    df_yq_eval %>% select(quarter_id, GDP, period),
    by = "quarter_id"
  ) %>%
  arrange(date)

rmsfe_by_period <- df_eval %>%
  group_by(period, type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  )

rmsfe_full <- df_eval %>%
  group_by(type) %>%
  summarise(
    RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(period = "Full sample")

rmsfe_rt_long <- bind_rows(rmsfe_full, rmsfe_by_period) %>%
  filter(period %in% c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")) %>%
  mutate(
    period = factor(period, levels = c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")),
    type   = factor(type, levels = c("M1", "M2", "M3"))
  ) %>%
  arrange(period, type)

rmsfe_rt_df <- rmsfe_rt_long %>%
  pivot_wider(names_from = type, values_from = RMSFE) %>%
  arrange(period)

rmsfe_rt <- as.matrix(rmsfe_rt_df[, c("M1", "M2", "M3")])
rownames(rmsfe_rt) <- rmsfe_rt_df$period

latex_rt <- latex_table_periods(
  mat     = rmsfe_rt,
  caption = paste0("Vector DFM Rolling RMSFE by period -- ", country, " (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:vector_dfm_rt_", country, "_", Size, "_", sel)
)

cat("\n", latex_rt, "\n")

tab_rt_all <- data.frame(
  country = country,
  period  = rownames(rmsfe_rt),
  M1      = rmsfe_rt[, "M1"],
  M2      = rmsfe_rt[, "M2"],
  M3      = rmsfe_rt[, "M3"],
  row.names = NULL
)

# ==============================================================================
# 14. SAVE PSEUDO REAL-TIME GRAPH
# ==============================================================================

file_graph_rt <- file.path(
  path_graph_country,
  paste0(
    "plot_rt_", model_name, "_", country,
    "_Size-", Size,
    "_sel-", sel,
    "_rfix-", r_fix,
    "_pfix-", p_fix,
    "_qfix-", q_fix,
    ".png"
  )
)

ggsave(
  filename = file_graph_rt,
  plot     = plot_rt,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\nSaved pseudo real-time graph to:\n", file_graph_rt, "\n")

# ==============================================================================
# 15. EXTRA DIAGNOSTICS: EM CONVERGENCE + GDP SMOOTHER FIT
# ==============================================================================

plot_dfm_fit_gdp <- function(dates, y_obs, y_fit, title, ylab) {
  df_fit <- tibble(
    date   = as.Date(dates),
    value  = as.numeric(y_fit),
    series = "Fitted"
  ) %>% filter(!is.na(value))
  
  df_obs_q <- tibble(
    date   = as.Date(dates),
    value  = as.numeric(y_obs),
    series = "Observed"
  ) %>%
    filter(!is.na(value)) %>%
    arrange(date)
  
  ggplot() +
    geom_line(data = df_fit,   aes(x = date, y = value, color = series), linewidth = 0.9) +
    geom_line(data = df_obs_q, aes(x = date, y = value, color = series), linewidth = 0.9) +
    theme_minimal(base_size = 14) +
    labs(title = title, y = ylab, x = NULL, color = NULL)
}

p_gdp_orig <- plot_dfm_fit_gdp(
  dates = dates_fit,
  y_obs = Xf_orig[, gdp_col],
  y_fit = X_fit_orig[, gdp_col],
  title = paste0("DFM – In-sample smoother fit (GDP original scale) | ", country),
  ylab  = "GDP"
)

file_graph_gdp_orig <- file.path(
  path_graph_country,
  paste0("DFM_SmootherFit_GDP_ORIG_", country, "_sel-", sel, "_r-", r_sel, "_p-", p_sel, "_q-", q_sel, ".png")
)

ggsave(file_graph_gdp_orig, p_gdp_orig, width = 12, height = 6, dpi = 300)

# ==============================================================================
# 16. SAVE SUMMARY OBJECT
# ==============================================================================

summary_vector_out <- list(
  model_id        = "DFM_EM",
  model           = model_name,
  stage           = "country_summary",
  Size            = Size,
  sel             = sel,
  country         = country,
  params_fit      = params,
  params_rt       = params_rt,
  
  df_now_full     = df_now_full,
  df_quarterly    = df_quarterly,
  plot_nowcast    = plot_nowcast,
  rmsfe_insample  = rmsfe_insample,
  latex_insample  = latex_insample,
  
  df_rt           = df_rt,
  df_yq_eval      = df_yq_eval,
  plot_rt         = plot_rt,
  rmsfe_rt        = rmsfe_rt,
  latex_rt        = latex_rt,
  
  file_fit        = file_fit,
  file_rt         = file_rt,
  file_graph_fit  = file_graph_fit,
  file_graph_rt   = file_graph_rt,
  
  diagnostics = list(
    X_raw        = X_raw,
    X_std        = X_std,
    X_fit_std    = X_fit_std,
    X_fit_orig   = X_fit_orig,
    dates_fit    = dates_fit,
    gdp_fit_orig = gdp_fit_orig
  ),
  
  # aliases for compatibility with final script
  plot_nowcast_facet     = plot_nowcast,
  plot_rt_facet          = plot_rt,
  tab_insample_all       = tab_insample_all,
  tab_rt_all             = tab_rt_all,
  df_rt_all              = df_rt,
  df_yq_eval_all         = df_yq_eval,
  latex_tab_insample_all = latex_insample,
  latex_tab_rt_all       = latex_rt
)

file_summary <- file.path(
  path_country,
  paste0(
    "VECTOR_DFM_summary",
    "_Size-", Size,
    "_sel-", sel,
    "_cc-", country,
    "_Nm-", N_m,
    "_Nq-", N_q,
    "_r-", r_sel,
    "_p-", p_sel,
    "_q-", q_sel,
    ".rds"
  )
)

if (file.exists(file_summary)) file.remove(file_summary)

saveRDS(summary_vector_out, file_summary)

cat("\nSaved summary object to:\n", file_summary, "\n")
