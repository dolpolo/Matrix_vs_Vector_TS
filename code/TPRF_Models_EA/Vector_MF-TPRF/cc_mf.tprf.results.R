# ==============================================================================
# Vector MF-TPRF - Results Script
# ------------------------------------------------------------------------------
# This script:
#   1. Loads full-sample and pseudo real-time results
#   2. Extracts true GDP and nowcasts
#   3. Produces full-sample and rolling plots
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
path_results <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/graph")

country <- "EA"

path_country       <- file.path(path_results, country)
path_graph_country <- file.path(path_graph, country)

path_graph_factors <- file.path(path_graph_country, "factor_interpretation")
dir.create(path_graph_country, recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_factors, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. TAGS OF THE RUN TO LOAD
# ==============================================================================

model_name <- "vector"
Size       <- "small"
sel        <- "LASSO"

# ==============================================================================
# 2. HELPERS
# ==============================================================================

make_quarter_avg_from_monthly <- function(x, T_q) {
  x <- as.numeric(x)
  stopifnot(length(x) >= 3 * T_q)
  x_use <- x[seq_len(3 * T_q)]
  as.numeric(tapply(x_use, rep(seq_len(T_q), each = 3), mean))
}

orient_factor_to_target <- function(f_m, y_q) {
  T_q <- length(y_q)
  f_q <- make_quarter_avg_from_monthly(f_m, T_q)
  
  if (cor(f_q, y_q, use = "pairwise.complete.obs") < 0) {
    f_m <- -f_m
    f_q <- -f_q
  }
  
  list(
    monthly   = as.numeric(f_m),
    quarterly = as.numeric(f_q)
  )
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

extract_fit_object <- function(obj) {
  if (!is.null(obj$fit)) return(obj$fit)
  if (!is.null(obj$MF_TPRF)) return(obj$MF_TPRF)
  if (!is.null(obj$full_sample$fit)) return(obj$full_sample$fit)
  if (!is.null(obj$full_sample$MF_TPRF)) return(obj$full_sample$MF_TPRF)
  stop("No fit object found in full-sample RDS.")
}

extract_hyper_object <- function(obj) {
  if (!is.null(obj$hyper)) return(obj$hyper)
  if (!is.null(obj$full_sample$hyper)) return(obj$full_sample$hyper)
  list(Lproxy = NA_real_, L_midas = NA_real_, p_ar = NA_real_, r_hat = NA_real_)
}

extract_dates_y <- function(obj) {
  dates_m <- NULL
  dates_q <- NULL
  y_q     <- NULL
  
  if (!is.null(obj$dates_m)) dates_m <- as.Date(obj$dates_m)
  if (!is.null(obj$dates_q)) dates_q <- as.Date(obj$dates_q)
  if (!is.null(obj$y_q))     y_q     <- as.numeric(obj$y_q)
  
  if ((is.null(dates_m) || is.null(dates_q) || is.null(y_q)) && !is.null(obj$inputs)) {
    if (is.null(dates_m) && !is.null(obj$inputs$dates_m)) dates_m <- as.Date(obj$inputs$dates_m)
    if (is.null(dates_q) && !is.null(obj$inputs$dates_q)) dates_q <- as.Date(obj$inputs$dates_q)
    if (is.null(y_q)     && !is.null(obj$inputs$y_q))     y_q     <- as.numeric(obj$inputs$y_q)
  }
  
  list(dates_m = dates_m, dates_q = dates_q, y_q = y_q)
}

extract_params_object <- function(obj) {
  if (!is.null(obj$params)) return(obj$params)
  stop("No params object found in RDS.")
}

# ==============================================================================
# 2B. FACTOR HELPERS
# ==============================================================================

theme_factor_plot <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey25"),
      axis.text.x      = element_text(angle = 45, hjust = 1)
    )
}

extract_vector_factors <- function(fit_obj) {
  cand <- list(
    fit_obj$F_hat,
    fit_obj$factors$F,
    fit_obj$F_hf,
    fit_obj$F,
    fit_obj$factor,
    fit_obj$factors
  )
  
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (length(cand) == 0L) {
    stop("No factor object found inside fit_obj.")
  }
  
  Fobj <- cand[[1]]
  
  if (is.vector(Fobj)) {
    Fmat <- matrix(as.numeric(Fobj), ncol = 1)
    colnames(Fmat) <- "F1"
    return(Fmat)
  }
  
  if (is.matrix(Fobj)) {
    Fmat <- Fobj
    if (nrow(Fmat) < ncol(Fmat)) {
      Fmat <- t(Fmat)
    }
    Fmat <- as.matrix(Fmat)
    if (is.null(colnames(Fmat))) {
      colnames(Fmat) <- paste0("F", seq_len(ncol(Fmat)))
    }
    return(Fmat)
  }
  
  stop("Factor object found, but its structure is not supported.")
}

# ==============================================================================
# 3. LOAD FULL-SAMPLE RESULTS
# ==============================================================================

file_fit <- find_result_file(
  path  = path_country,
  model = model_name,
  stage = "fit",
  Size  = Size,
  sel   = sel
)

fit_res <- readRDS(file_fit)

cat("\nLoaded full-sample results from:\n", file_fit, "\n")

model_name <- if (!is.null(fit_res$model)) fit_res$model else model_name
Size       <- if (!is.null(fit_res$Size))  fit_res$Size  else Size
sel        <- if (!is.null(fit_res$sel))   fit_res$sel   else sel
country    <- if (!is.null(fit_res$country)) fit_res$country else country

# ==============================================================================
# 4. EXTRACT FULL-SAMPLE OBJECTS
# ==============================================================================

params    <- extract_params_object(fit_res)
fit_obj   <- extract_fit_object(fit_res)
hyper_fit <- extract_hyper_object(fit_res)

meta_obj <- extract_dates_y(fit_res)
dates_m  <- meta_obj$dates_m
dates_q  <- meta_obj$dates_q
y_true_q <- meta_obj$y_q

if (is.null(dates_m) || is.null(dates_q) || is.null(y_true_q)) {
  stop("The full-sample RDS must contain dates_m, dates_q, and y_q.")
}

Lproxy  <- hyper_fit$Lproxy
L_midas <- hyper_fit$L_midas
p_ar    <- hyper_fit$p_ar
r_hat   <- hyper_fit$r_hat

N_m <- if (!is.null(fit_res$metadata$N_m)) fit_res$metadata$N_m else NA
N_q <- if (!is.null(fit_res$metadata$N_q)) fit_res$metadata$N_q else NA

T_q_complete <- length(y_true_q)

start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

# ==============================================================================
# 5. EXTRACT FULL-SAMPLE NOWCAST
# ==============================================================================

y_now_full <- fit_obj$y_nowcast
T_m_full   <- length(y_now_full)

if (is.null(y_now_full)) stop("The fit object does not contain y_nowcast.")
if (T_m_full < 3 * T_q_complete) stop("y_now_full has fewer than 3 * T_q_complete months.")

dates_m_full <- dates_m[seq_len(T_m_full)]

n_in  <- 3 * T_q_complete
n_tot <- T_m_full

period_vec <- rep("in-sample", n_tot)
if (n_tot > n_in) period_vec[(n_in + 1):n_tot] <- "real-time"

df_now_full <- data.frame(
  date       = dates_m_full,
  y_now_full = as.numeric(y_now_full),
  period     = factor(period_vec, levels = c("in-sample", "real-time"))
)

df_quarterly <- data.frame(
  date   = dates_q,
  y_true = y_true_q
)

# ==============================================================================
# 5B. VECTOR FACTOR INTERPRETATION
# ==============================================================================

title_factors_vec    <- "Estimated latent factor"
subtitle_factors_vec <- "MF-TPRF on EA data"

F_vec <- extract_vector_factors(fit_obj)

T_fac <- nrow(F_vec)
dates_f_vec <- dates_m[seq_len(T_fac)]

df_factors_vec <- as.data.frame(F_vec)
df_factors_vec$date <- as.Date(dates_f_vec)

df_factors_vec <- df_factors_vec %>%
  pivot_longer(
    cols = -date,
    names_to = "factor",
    values_to = "value"
  )

plot_factors_vec <- ggplot(df_factors_vec, aes(x = date, y = value)) +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.16
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.8, colour = "black") +
  facet_wrap(~ factor, ncol = 1, scales = "free_y") +
  labs(
    title = title_factors_vec,
    subtitle = subtitle_factors_vec,
    x = "Date",
    y = "Factor value"
  ) +
  theme_factor_plot()

print(plot_factors_vec)

# ==============================================================================
# 5C. PCA FACTOR ON EA DATA
# ==============================================================================

title_pca1    <- "First principal component"
subtitle_pca1 <- "PCA on imputed EA monthly predictors"

if (is.null(fit_res$preprocessing$X_xp)) {
  stop("X_xp not found in fit_res$preprocessing.")
}

X_xp <- as.matrix(fit_res$preprocessing$X_xp)

pca_fit <- prcomp(X_xp, center = TRUE, scale. = FALSE)

pca1 <- as.numeric(pca_fit$x[, 1])

df_pca1 <- data.frame(
  date  = as.Date(dates_m[seq_len(length(pca1))]),
  value = pca1
)

plot_pca1 <- ggplot(df_pca1, aes(x = date, y = value)) +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.16
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.8, colour = "black") +
  labs(
    title = title_pca1,
    subtitle = subtitle_pca1,
    x = "Date",
    y = "PC1 value"
  ) +
  theme_factor_plot()

print(plot_pca1)

# ==============================================================================
# 5D. COMPARISON: VECTOR FACTOR VS PCA1
# ==============================================================================

title_compare    <- "Estimated factor vs first principal component"
subtitle_compare <- "Standardized series, sign-aligned to EA GDP"

f_vec_1 <- as.numeric(F_vec[, 1])

f_vec_q  <- make_quarter_avg_from_monthly(f_vec_1, length(y_true_q))
pca1_q   <- make_quarter_avg_from_monthly(pca1,    length(y_true_q))

if (cor(f_vec_q, y_true_q, use = "pairwise.complete.obs") < 0) {
  f_vec_1 <- -f_vec_1
  f_vec_q <- -f_vec_q
}

if (cor(pca1_q, y_true_q, use = "pairwise.complete.obs") < 0) {
  pca1 <- -pca1
  pca1_q <- -pca1_q
}

df_compare <- data.frame(
  date    = as.Date(dates_m[seq_len(length(f_vec_1))]),
  MF_TPRF = as.numeric(scale(f_vec_1)),
  PCA1    = as.numeric(scale(pca1))
) %>%
  pivot_longer(
    cols = c("MF_TPRF", "PCA1"),
    names_to = "series",
    values_to = "value"
  )

plot_compare <- ggplot(df_compare, aes(x = date, y = value, colour = series)) +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.16
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.9) +
  labs(
    title = title_compare,
    subtitle = subtitle_compare,
    x = "Date",
    y = "Standardized value",
    colour = ""
  ) +
  theme_factor_plot()

print(plot_compare)

# ==============================================================================
# 5E. FACTORS TO SAVE FOR FINAL CROSS-MODEL COMPARISON
# ==============================================================================

factor_vector_monthly_raw <- as.numeric(F_vec[, 1])

factor_vector_oriented <- orient_factor_to_target(
  f_m = factor_vector_monthly_raw,
  y_q = y_true_q
)

factor_vector_monthly   <- factor_vector_oriented$monthly
factor_vector_quarterly <- factor_vector_oriented$quarterly

pca1_oriented <- orient_factor_to_target(
  f_m = pca1,
  y_q = y_true_q
)

pca1_monthly_oriented <- pca1_oriented$monthly
pca1_quarterly        <- pca1_oriented$quarterly

dates_factor_m <- as.Date(dates_m[seq_len(length(factor_vector_monthly))])
dates_factor_q <- as.Date(dates_q)

df_factor_compare_vector_q <- data.frame(
  date   = dates_factor_q,
  GDP    = as.numeric(y_true_q),
  Vector = factor_vector_quarterly,
  PCA1   = pca1_quarterly
)

file_graph_factor_vec <- file.path(
  path_graph_factors,
  paste0("plot_factors_vector_mf_tprf_", country,
         "_Size-", Size, "_sel-", sel, ".png")
)

file_graph_pca1 <- file.path(
  path_graph_factors,
  paste0("plot_pca1_", country,
         "_Size-", Size, "_sel-", sel, ".png")
)

file_graph_compare <- file.path(
  path_graph_factors,
  paste0("plot_compare_factor_vs_pca1_", country,
         "_Size-", Size, "_sel-", sel, ".png")
)

ggsave(file_graph_factor_vec, plot_factors_vec, width = 10, height = 6, dpi = 300)
ggsave(file_graph_pca1,       plot_pca1,        width = 10, height = 4.8, dpi = 300)
ggsave(file_graph_compare,    plot_compare,     width = 10, height = 5.2, dpi = 300)

# ==============================================================================
# 6. FULL-SAMPLE PLOT
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
  geom_point(
    data = df_now_full %>% filter(period == "real-time"),
    aes(x = date, y = y_now_full),
    size = 2.0, color = "black"
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
    title = paste0("Vector MF-TPRF Monthly Nowcast - ", country),
    x = "Date",
    y = "GDP Growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

print(plot_nowcast)

# ==============================================================================
# 7. IN-SAMPLE RMSFE
# ==============================================================================

is_PRE   <- dates_q >= start_eval  & dates_q <  covid_start
is_COVID <- dates_q >= covid_start & dates_q <= covid_end
is_POST  <- dates_q >  covid_end   & dates_q <= end_eval
is_ALL   <- dates_q >= start_eval  & dates_q <= end_eval

y_now_in <- y_now_full[1:n_in]

M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

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
  caption = paste0("Vector MF-TPRF RMSFE by period -- ", country, " (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:vector_mf_tprf_", country, "_", Size, "_", sel)
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
# 8. SAVE FULL-SAMPLE GRAPH
# ==============================================================================

file_graph_fit <- file.path(
  path_graph_country,
  paste0(
    "plot_fit_", model_name, "_", country,
    "_Size-", Size,
    "_sel-", sel,
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-", p_ar,
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
# 9. LOAD PSEUDO REAL-TIME RESULTS
# ==============================================================================

file_rt <- find_result_file(
  path  = path_country,
  model = model_name,
  stage = "rt",
  Size  = Size,
  sel   = sel
)

rt_res <- readRDS(file_rt)

cat("\nLoaded pseudo real-time results from:\n", file_rt, "\n")

# ==============================================================================
# 10. EXTRACT PSEUDO REAL-TIME OBJECTS
# ==============================================================================

params_rt <- extract_params_object(rt_res)
meta_rt   <- extract_dates_y(rt_res)

dates_q_rt <- meta_rt$dates_q
y_q_rt     <- meta_rt$y_q

if (is.null(dates_q_rt) || is.null(y_q_rt)) {
  stop("The pseudo real-time RDS must contain dates_q and y_q.")
}

df_rt <- rt_res$pseudo_realtime_all %>%
  mutate(
    date = as.Date(date),
    type = factor(month_in_quarter, levels = c("M1", "M2", "M3"))
  ) %>%
  select(date, nowcast, type) %>%
  arrange(date, type)

idx_eval <- which(dates_q_rt >= params_rt$start_eval & dates_q_rt <= params_rt$end_eval)

df_yq <- data.frame(
  date = dates_q_rt[idx_eval],
  GDP  = y_q_rt[idx_eval]
)

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
    title    = paste0("Vector MF-TPRF Rolling Real-Time Nowcasts - ", country),
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
  caption = paste0("Vector MF-TPRF Rolling RMSFE by period -- ", country, " (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:vector_mf_tprf_rt_", country, "_", Size, "_", sel)
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
# 13. SAVE PSEUDO REAL-TIME GRAPH
# ==============================================================================

hyper_rt <- if (!is.null(rt_res$pseudo_realtime_raw$hyper_pre)) {
  rt_res$pseudo_realtime_raw$hyper_pre
} else {
  list(Lproxy = NA, L_midas = NA, p_AR = NA)
}

file_graph_rt <- file.path(
  path_graph_country,
  paste0(
    "plot_rt_", model_name, "_", country,
    "_Size-", Size,
    "_sel-", sel,
    "_Lproxy-", hyper_rt$Lproxy,
    "_Lmidas-", hyper_rt$L_midas,
    "_pAR-", hyper_rt$p_AR,
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
# 14. SAVE SUMMARY OBJECT
# ==============================================================================

summary_vector_out <- list(
  model_id        = "MF_TPRF",
  model           = model_name,
  stage           = "country_summary",
  Size            = Size,
  sel             = sel,
  country         = country,
  params_fit      = params,
  params_rt       = params_rt,
  
  factor_interpretation = list(
    df_factors_vec        = df_factors_vec,
    df_pca1               = df_pca1,
    df_compare            = df_compare,
    plot_factors_vec      = plot_factors_vec,
    plot_pca1             = plot_pca1,
    plot_compare          = plot_compare,
    file_graph_factor_vec = file_graph_factor_vec,
    file_graph_pca1       = file_graph_pca1,
    file_graph_compare    = file_graph_compare
  ),
  
  factor_comparison = list(
    dates_m = dates_factor_m,
    dates_q = dates_factor_q,
    gdp_q   = as.numeric(y_true_q),
    
    vector_factor_monthly   = factor_vector_monthly,
    vector_factor_quarterly = factor_vector_quarterly,
    
    pca1_monthly   = pca1_monthly_oriented,
    pca1_quarterly = pca1_quarterly,
    
    vector_factor_name = "Vector MF-TPRF",
    pca1_name          = "PCA1",
    
    df_q = df_factor_compare_vector_q
  ),
  
  # original fields
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

file_summary <- build_result_filename(
  path_out         = path_country,
  model            = model_name,
  stage            = "summary",
  Size             = Size,
  sel              = sel,
  countries        = country,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = Lproxy,
  L_midas          = L_midas,
  p_ar             = p_ar,
  r1               = r_hat,
  r2               = NA,
  robust_f         = as.integer(isTRUE(params$Robust_F)),
  covid_m          = as.integer(isTRUE(params$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params$covid_mask_q)),
  ext              = "rds",
  timestamp        = FALSE,
  include_details  = TRUE
)

if (file.exists(file_summary)) file.remove(file_summary)

saveRDS(summary_vector_out, file_summary)

cat("\nSaved summary object to:\n", file_summary, "\n")