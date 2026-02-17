# ==============================================================================
# DFM – DIAGNOSTICS + VISUALIZATION SECTION 
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)


build_dfm_full_filename <- function(path_results, country, params, NM, NQ, r, p) {
  
  path_country <- file.path(path_results, country)
  
  file_out_dfm <- file.path(
    path_country,
    paste0(
      "DFM_EM_", country,
      "_sel-",   params$sel_method,
      "_Covid-", params$covid_mask,
      "_Nm-",    NM,
      "_Nq-",    NQ,
      "_r-",     r,
      "_p-",     p,
      "_", format(params$start_est, "%Y-%m"),
      "_to_", format(params$end_eval, "%Y-%m"),
      ".rds"
    )
  )
  
  return(file_out_dfm)
}


build_dfm_rt_filename <- function(path_results, country, params, Freq) {
  
  path_country <- file.path(path_results, country)
  
  sel_method <- if (!is.null(params$sel_method)) params$sel_method else "none"
  covid_mask <- if (!is.null(params$covid_mask)) params$covid_mask else "TRUE"
  
  start_est_str <- format(params$start_est, "%Y-%m")
  end_eval_str  <- format(params$end_eval, "%Y-%m")
  
  file.path(
    path_country,
    paste0(
      "DFM_RollingNowcast_", country,
      "_sel-", sel_method,
      "_Covid-", covid_mask,
      "_Nm-", sum(Freq == "M"),
      "_Nq-", sum(Freq == "Q"),
      "_", start_est_str,
      "_to_", end_eval_str,
      ".rds"
    )
  )
}



# ----------------------------------------------------------------------
# 1. LOAD RESULTS (FULL SAMPLE + ROLLING NOWCAST)
# ----------------------------------------------------------------------

file_dfm_full <- build_dfm_full_filename(
  path_results = path_results,
  country      = country,
  params       = params,
  NM           = NM,
  NQ           = NQ,
  r            = r,
  p            = p
)

dfm_full_sample <- readRDS(file_dfm_full)    # full-sample EM DFM

# Extract full-sample objects
X_raw     <- dfm_full_sample$X_raw
X_std     <- dfm_full_sample$X_std
X_in      <- dfm_full_sample$X_in
dates_m   <- dfm_full_sample$dates_m
dates_q   <- dfm_full_sample$dates_q
Init_full <- dfm_full_sample$Init
res_full  <- list(
  A = dfm_full_sample$A,
  C = dfm_full_sample$C,
  Q = dfm_full_sample$Q,
  R = dfm_full_sample$R,
  Z0 = dfm_full_sample$Z0,
  V0 = dfm_full_sample$V0
)

file_dfm_rt <- build_dfm_rt_filename(
  path_results = path_results,
  country      = country,
  params       = params,
  Freq         = Freq
)

dfm_realtime    <- readRDS(file_dfm_rt)      # rolling pseudo-real-time nowcast

now_dfm_re <- dfm_realtime$nowcast
gdp_col    <- dfm_realtime$gdp_col
Freq       <- dfm_realtime$Freq
Unb        <- dfm_realtime$Unb

# ----------------------------------------------------------------------
# 2. BASIC DIAGNOSTICS ON RAW DATA
# ----------------------------------------------------------------------

cat("Missing values:", sum(is.na(data)), "\n")
cat("Share missing:", round(100 * sum(is.na(data)) / prod(dim(data)), 2), "%\n")

# ----------------------------------------------------------------------
# 3. KALMAN SMOOTHER (IN-SAMPLE RECONSTRUCTION)
# ----------------------------------------------------------------------

fit_res_dfm <- dfm_kalman_fit(
  X_in     = X_in,                 # standardized predictors (complete)
  res      = res_full,             # EM result
  n_lags_Q = Init_full$n_lags_Q
)

X_fit <- fit_res_dfm$X_fit         # T_eff × N fitted signals
t_idx <- fit_res_dfm$t_index       # time indices used

# ----------------------------------------------------------------------
# 4. METRICS FUNCTION
# ----------------------------------------------------------------------

dfm_fit_metrics <- function(X_in, X_fit, t_idx, var_names = NULL) {
  
  X  <- as.matrix(X_in)
  Xf <- X[t_idx, , drop = FALSE]     # T_eff × N
  
  N  <- ncol(Xf)
  
  if (is.null(var_names)) {
    var_names <- colnames(X_in)
    if (is.null(var_names)) var_names <- paste0("Var", seq_len(N))
  }
  
  resid <- Xf - X_fit
  
  var_y <- apply(Xf, 2, var)
  mse   <- apply(resid^2, 2, mean)
  rmse  <- sqrt(mse)
  corr  <- sapply(1:N, function(j)
    cor(Xf[, j], X_fit[, j], use = "pairwise.complete.obs")
  )
  
  R2 <- 1 - mse / var_y
  
  data.frame(
    variable = var_names,
    var_y    = var_y,
    mse      = mse,
    rmse     = rmse,
    R2       = R2,
    corr_fit = corr
  )
}

metrics <- dfm_fit_metrics(
  X_in     = X_in,
  X_fit    = X_fit,
  t_idx    = t_idx,
  var_names = colnames(X_in)
)

# ----------------------------------------------------------------------
# 5. SINGLE GDP FIT GRAPH
# ----------------------------------------------------------------------

plot_dfm_fit_gdp <- function(X_in, X_fit, t_idx,
                             dates, gdp_idx,
                             var_names = NULL,
                             title = "DFM – In-sample fit (standardized GDP)") {
  
  X  <- as.matrix(X_in)
  Xf <- X[t_idx, , drop = FALSE]
  
  if (is.null(var_names)) var_names <- colnames(X_in)
  
  dates_fit <- dates[t_idx]
  
  df_plot <- tibble(
    date  = rep(dates_fit, 2),
    value = c(Xf[, gdp_idx], X_fit[, gdp_idx]),
    type  = rep(c("Observed", "Fitted"), each = length(dates_fit))
  )
  
  ggplot(df_plot, aes(x = date, y = value, color = type)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = c("Observed" = "black", "Fitted" = "red")) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      y     = "Standardized GDP",
      x     = NULL,
      color = NULL
    ) +
    coord_cartesian(ylim = c(-2, 2))
}

plot_dfm_gdp_fit <- plot_dfm_fit_gdp(
  X_in     = X_in,
  X_fit    = X_fit,
  t_idx    = t_idx,
  dates    = dates_m,
  gdp_idx  = gdp_col,
  var_names = colnames(X_in)
)

print(plot_dfm_gdp_fit)

# ==============================================================================
# SAVE GRAPH: DFM smoother fit for GDP
# ==============================================================================

# Directory: /graphs/<country>/
path_graph_country <- file.path(path_graph, country)

if (!dir.exists(path_graph_country)) {
  dir.create(path_graph_country, recursive = TRUE)
}

# Build filename
file_graph_dfm_smoother <- file.path(
  path_graph_country,
  paste0(
    "DFM_SmootherFit_GDP_", country,
    "_sel-", params$sel_method,
    "_Covid-", params$covid_mask,
    ".png"
  )
)

# Save the plot
ggsave(
  filename = file_graph_dfm_smoother,
  plot     = plot_dfm_gdp_fit,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\n*** SAVED DFM GDP SMOOTHER FIT TO ***\n", file_graph_dfm_smoother, "\n")


# ----------------------------------------------------------------------
# 6. NOWCAST VISUALIZATION – M1/M2/M3 vs TRUE GDP
# ----------------------------------------------------------------------

# destandardize GDP if needed
gdp_mean <- out_std$mean[gdp_col]
gdp_sd   <- out_std$sd[gdp_col]

now_dfm_re_raw <- list(
  M1 = now_dfm_re$M1 * gdp_sd + gdp_mean,
  M2 = now_dfm_re$M2 * gdp_sd + gdp_mean,
  M3 = now_dfm_re$M3 * gdp_sd + gdp_mean
)

# Build quarterly DF
build_nowcast_quarter_df <- function(now_dfm_re, y_q, dates_q) {
  q_M1 <- as.Date(names(now_dfm_re$M1))
  q_M2 <- as.Date(names(now_dfm_re$M2))
  q_M3 <- as.Date(names(now_dfm_re$M3))
  
  q_all <- Reduce(intersect, list(dates_q, q_M1, q_M2, q_M3))
  
  tibble(
    quarter_end = q_all,
    y_true      = y_q[match(q_all, dates_q)],
    M1          = now_dfm_re$M1[match(q_all, q_M1)],
    M2          = now_dfm_re$M2[match(q_all, q_M2)],
    M3          = now_dfm_re$M3[match(q_all, q_M3)]
  )
}

df_q <- build_nowcast_quarter_df(
  now_dfm_re = now_dfm_re_raw,
  y_q        = y_q,
  dates_q    = dates_q
) %>%
  filter(quarter_end >= params$start_eval)

df_q_long <- df_q %>%
  pivot_longer(
    cols = c(y_true, M1, M2, M3),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(series = factor(series,
                         levels = c("y_true", "M1", "M2", "M3"),
                         labels = c("True GDP", "Nowcast M1", "Nowcast M2", "Nowcast M3")))

plot_dfm_nowcast_M123 <- ggplot(df_q_long,
                                aes(x = quarter_end, y = value, color = series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  annotate("rect",
           xmin = params$covid_start,
           xmax = params$covid_end,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.15) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – Nowcast M1/M2/M3 vs True GDP (", country, ")"),
    x     = "Quarter end",
    y     = "GDP",
    color = NULL
  ) +
  scale_color_manual(values = c(
    "True GDP" = "black",
    "Nowcast M1" = "#1f77b4",
    "Nowcast M2" = "#ff7f0e",
    "Nowcast M3" = "#2ca02c"
  ))

print(plot_dfm_nowcast_M123)

# ==============================================================================
# SAVE GRAPH: DFM Nowcast M1/M2/M3 vs True GDP
# ==============================================================================

# Directory: /graphs/<country>/
path_graph_country <- file.path(path_graph, country)

if (!dir.exists(path_graph_country)) {
  dir.create(path_graph_country, recursive = TRUE)
}

# File name
file_graph_dfm_M123 <- file.path(
  path_graph_country,
  paste0(
    "DFM_Nowcast_M123_", country,
    "_sel-", params$sel_method,
    "_Covid-", params$covid_mask,
    ".png"
  )
)

# Save the plot
ggsave(
  filename = file_graph_dfm_M123,
  plot     = plot_dfm_nowcast_M123,
  width    = 12,
  height   = 6,
  dpi      = 300
)

cat("\n*** SAVED DFM NOWCAST M1/M2/M3 TO ***\n", file_graph_dfm_M123, "\n")

