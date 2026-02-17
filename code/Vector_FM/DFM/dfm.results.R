# ==============================================================================
# DFM – DIAGNOSTICS + VISUALIZATION SECTION (CORRECTED)
# What this fixes vs your draft:
#   1) Uses the objects you actually saved in the FULL-SAMPLE rds:
#        - X_raw, X_std, dates_m, dates_q, Init, A,C,Q,R,Z0,V0, loglik, crit
#        - (and optionally: out_std if you saved it; otherwise we reconstruct it)
#   2) Correct GDP column handling: use saved gdp_col when available.
#   3) Smooth/Fitted of GDP shown in:
#        - standardized scale (as before)
#        - ORIGINAL SCALE (using the FULL-SAMPLE standardization map)
#   4) EM convergence plot: loglik + criterion c_j from the saved res
#   5) Rolling nowcast reading consistent with the NEW rolling save format:
#        dfm_realtime$nowcast with M1/M2/M3 std & orig.
#   6) No "X_in" dependency: we use X_std for fit metrics and smoother.
#   7) Keeps quarterly alignment robust.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ------------------------------------------------------------------------------
# Helpers: COVID tag + filenames
# ------------------------------------------------------------------------------

covid_tag_from_params <- function(params) {
  paste0(
    "M", as.integer(isTRUE(params$covid_mask_m)),
    "_Q", as.integer(isTRUE(params$covid_mask_q))
  )
}

build_dfm_full_filename <- function(path_results, country, params, NM, NQ, r, p) {
  path_country <- file.path(path_results, country)
  covid_tag <- covid_tag_from_params(params)
  
  file.path(
    path_country,
    paste0(
      "DFM_EM_", country,
      "_sel-",   params$sel_method,
      "_Covid-", covid_tag,
      "_Nm-",    NM,
      "_Nq-",    NQ,
      "_r-",     r,
      "_p-",     p,
      "_", format(params$start_est, "%Y-%m"),
      "_to_", format(params$end_eval, "%Y-%m"),
      ".rds"
    )
  )
}

build_dfm_rt_filename <- function(path_results, country, params, NM, NQ,
                                  r_fix = NULL, p_fix = NULL) {
  path_country <- file.path(path_results, country)
  
  sel_method <- if (!is.null(params$sel_method)) params$sel_method else "none"
  covid_tag  <- covid_tag_from_params(params)
  
  start_est_str <- format(params$start_est, "%Y-%m")
  end_eval_str  <- format(params$end_eval, "%Y-%m")
  
  rp_tag <- ""
  if (!is.null(r_fix) && !is.null(p_fix)) rp_tag <- paste0("_rfix-", r_fix, "_pfix-", p_fix)
  
  file.path(
    path_country,
    paste0(
      "DFM_RollingNowcast_", country,
      "_sel-", sel_method,
      "_Covid-", covid_tag,
      "_Nm-", NM,
      "_Nq-", NQ,
      rp_tag,
      "_", start_est_str,
      "_to_", end_eval_str,
      ".rds"
    )
  )
}

# ------------------------------------------------------------------------------
# 0. REQUIREMENTS (must exist in workspace)
#   path_results, path_graph, country, params, NM, NQ, r, p
#   and your functions: dfm_kalman_fit(), kalman.na()
#   (dfm_kalman_fit uses X_in standardized; we'll pass X_std)
# ------------------------------------------------------------------------------

stopifnot(exists("path_results"), exists("path_graph"), exists("country"), exists("params"))
stopifnot(exists("NM"), exists("NQ"), exists("r"), exists("p"))
stopifnot(exists("dfm_kalman_fit"))

# ------------------------------------------------------------------------------
# 1. LOAD RESULTS (FULL SAMPLE + ROLLING)
# ------------------------------------------------------------------------------

# ---- full sample ----
file_dfm_full <- build_dfm_full_filename(
  path_results = path_results,
  country      = country,
  params       = params,
  NM           = NM,
  NQ           = NQ,
  r            = r,
  p            = p
)
if (!file.exists(file_dfm_full)) stop("Full-sample DFM file not found: ", file_dfm_full)
dfm_full_sample <- readRDS(file_dfm_full)

# Objects you SAVED in the full-sample file
X_raw   <- dfm_full_sample$X_raw
X_std   <- dfm_full_sample$X_std
dates_m <- dfm_full_sample$dates_m
dates_q <- dfm_full_sample$dates_q
Init_full <- dfm_full_sample$Init

res_full <- list(
  A  = dfm_full_sample$A,
  C  = dfm_full_sample$C,
  Q  = dfm_full_sample$Q,
  R  = dfm_full_sample$R,
  Z0 = dfm_full_sample$Z0,
  V0 = dfm_full_sample$V0
)

loglik_full <- dfm_full_sample$loglik
crit_full   <- dfm_full_sample$crit

# GDP column index (prefer saved gdp_col; otherwise fallback)
gdp_col <- dfm_full_sample$gdp_col
if (is.null(gdp_col) || length(gdp_col) != 1) {
  # try infer from params$target if you have column names with GDP_...
  target_name <- dfm_full_sample$params$target
  if (!is.null(target_name) && !is.null(colnames(X_raw))) {
    cand <- which(startsWith(colnames(X_raw), paste0(target_name, "_")))
    if (length(cand) == 1) gdp_col <- cand
  }
}
if (is.null(gdp_col) || length(gdp_col) != 1) {
  gdp_col <- get0("gdp_col", ifnotfound = NA_integer_)
}
if (!is.finite(gdp_col)) stop("gdp_col not found (not in full-sample RDS and not in workspace).")

# ---- rolling nowcast ----
file_dfm_rt <- build_dfm_rt_filename(
  path_results = path_results,
  country      = country,
  params       = params,
  NM           = NM,
  NQ           = NQ,
  r_fix        = NULL,
  p_fix        = NULL
)

if (!file.exists(file_dfm_rt)) {
  # autodiscover (pick most recent among matching pattern)
  path_country <- file.path(path_results, country)
  
  sel_method <- if (!is.null(params$sel_method)) params$sel_method else "none"
  covid_tag  <- covid_tag_from_params(params)
  
  patt <- paste0(
    "^DFM_RollingNowcast_", country,
    "_sel-", sel_method,
    "_Covid-", covid_tag,
    "_Nm-", NM,
    "_Nq-", NQ
  )
  
  cand <- list.files(path_country, pattern = patt, full.names = TRUE)
  if (length(cand) == 0) stop("Rolling nowcast file not found (check naming/tag).")
  file_dfm_rt <- cand[which.max(file.info(cand)$mtime)]
  message("Rolling nowcast file auto-selected: ", basename(file_dfm_rt))
}

dfm_realtime <- readRDS(file_dfm_rt)

# NEW rolling structure:
# dfm_realtime$nowcast contains r_fix, p_fix, M1_std/M2_std/M3_std, M1_orig/M2_orig/M3_orig
if (!is.null(dfm_rolling$nowcast)) {
  now_obj <- dfm_rolling$nowcast
  r_fix <- now_obj$r_fix
  p_fix <- now_obj$p_fix
  
  now_dfm_std <- list(M1 = now_obj$M1_std, M2 = now_obj$M2_std, M3 = now_obj$M3_std)
  now_dfm_orig <- list(M1 = now_obj$M1_orig, M2 = now_obj$M2_orig, M3 = now_obj$M3_orig)
} else {
  # backward compat
  r_fix <- dfm_realtime$r_fix
  p_fix <- dfm_realtime$p_fix
  now_dfm_std <- list(M1 = dfm_realtime$M1_std, M2 = dfm_realtime$M2_std, M3 = dfm_realtime$M3_std)
  now_dfm_orig <- list(M1 = dfm_realtime$M1_orig, M2 = dfm_realtime$M2_orig, M3 = dfm_realtime$M3_orig)
}

# metadata fallback (prefer rolling if present)
if (!is.null(dfm_realtime$gdp_col)) gdp_col <- dfm_realtime$gdp_col
if (!is.null(dfm_realtime$dates_m)) dates_m <- dfm_realtime$dates_m
if (!is.null(dfm_realtime$dates_q)) dates_q <- dfm_realtime$dates_q

# ------------------------------------------------------------------------------
# 1b. TRUE QUARTERLY GDP SERIES from raw data
# ------------------------------------------------------------------------------

idx_q_obs <- which(!is.na(X_raw[, gdp_col]))
y_q <- as.numeric(X_raw[idx_q_obs, gdp_col])
dates_q_from_data <- as.Date(dates_m[idx_q_obs])

# ensure dates_q coherence
if (length(dates_q) != length(dates_q_from_data) ||
    any(as.Date(dates_q) != as.Date(dates_q_from_data))) {
  message("WARNING: dates_q in file does not match GDP-observed rows; using dates_q derived from dates_m.")
  dates_q <- dates_q_from_data
}
stopifnot(length(y_q) == length(dates_q))

# ------------------------------------------------------------------------------
# 2. BASIC DIAGNOSTICS ON RAW DATA
# ------------------------------------------------------------------------------

cat("Missing values (raw):", sum(is.na(X_raw)), "\n")
cat("Share missing (raw):", round(100 * sum(is.na(X_raw)) / prod(dim(X_raw)), 2), "%\n")

# ------------------------------------------------------------------------------
# 3. FULL-SAMPLE STANDARDIZATION MAP (to destandardize smoother/fits)
#   If you saved out_std in the full-sample file, use it.
#   Otherwise, reconstruct it from X_raw the same way you standardized.
# ------------------------------------------------------------------------------

# recommended: in your estimation script, save out_std as dfm_full_sample$out_std
out_std_full <- dfm_full_sample$out_std
if (is.null(out_std_full)) {
  if (!exists("standardize_with_na")) stop("standardize_with_na() not found to reconstruct out_std.")
  out_std_full <- standardize_with_na(X_raw)  # MUST match your original standardize step
}

mu_full <- out_std_full$mean
sd_full <- out_std_full$sd

# helpers: vectorized de-standardization (keeps NA)
destd_vec <- function(x_std, mu, sd) mu + sd * x_std
destd_mat <- function(X_std_mat, mu, sd) {
  stopifnot(ncol(X_std_mat) == length(mu), ncol(X_std_mat) == length(sd))
  sweep(sweep(X_std_mat, 2, sd, `*`), 2, mu, `+`)
}

# ------------------------------------------------------------------------------
# 4. KALMAN SMOOTHER FIT (IN-SAMPLE) on standardized panel
#   IMPORTANT: dfm_kalman_fit returns X_fit for rows t_idx (starting at n_lags_Q)
# ------------------------------------------------------------------------------

# n_lags_Q must be consistent with your estimated model:
# stock_flow -> 3, MM -> 5. If you stored it in Init, prefer that.
n_lags_Q_full <- dfm_full_sample$Init$meta$L
if (is.null(n_lags_Q_full) || !is.finite(n_lags_Q_full)) n_lags_Q_full <- 3L

fit_res_dfm <- dfm_kalman_fit(
  X_in     = X_std,          # standardized input
  res      = res_full,
  n_lags_Q = as.integer(n_lags_Q_full)
)

X_fit_std <- fit_res_dfm$X_fit           # (T_eff x N) fitted standardized
t_idx     <- fit_res_dfm$t_index         # indices in original monthly timeline

# Build aligned standardized observed block:
Xf_std <- X_std[t_idx, , drop = FALSE]
dates_fit <- as.Date(dates_m[t_idx])

# De-standardize fitted + observed (full-sample mapping)
X_fit_orig <- destd_mat(X_fit_std, mu_full, sd_full)
Xf_orig    <- destd_mat(Xf_std,    mu_full, sd_full)

# ------------------------------------------------------------------------------
# 5. FIT METRICS (STANDARDIZED + ORIGINAL) and SAVE
# ------------------------------------------------------------------------------

dfm_fit_metrics <- function(X_obs, X_fit, var_names = NULL) {
  X_obs <- as.matrix(X_obs)
  X_fit <- as.matrix(X_fit)
  stopifnot(all(dim(X_obs) == dim(X_fit)))
  
  N <- ncol(X_obs)
  if (is.null(var_names)) {
    var_names <- colnames(X_obs)
    if (is.null(var_names)) var_names <- paste0("Var", seq_len(N))
  }
  
  resid <- X_obs - X_fit
  
  var_y <- apply(X_obs, 2, var, na.rm = TRUE)
  mse   <- apply(resid^2, 2, mean, na.rm = TRUE)
  rmse  <- sqrt(mse)
  corr  <- sapply(seq_len(N), function(j) cor(X_obs[, j], X_fit[, j], use = "pairwise.complete.obs"))
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

metrics_std  <- dfm_fit_metrics(X_obs = Xf_std,  X_fit = X_fit_std,  var_names = colnames(X_std))
metrics_orig <- dfm_fit_metrics(X_obs = Xf_orig, X_fit = X_fit_orig, var_names = colnames(X_raw))

# save metrics
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) dir.create(path_country, recursive = TRUE)

covid_tag <- covid_tag_from_params(params)

metrics_base <- paste0(
  "DFM_FitMetrics_", country,
  "_sel-", params$sel_method,
  "_Covid-", covid_tag,
  "_r-", r,
  "_p-", p,
  "_", format(params$start_est, "%Y-%m"),
  "_to_", format(params$end_eval, "%Y-%m")
)

saveRDS(list(std = metrics_std, orig = metrics_orig),
        file.path(path_country, paste0(metrics_base, ".rds")))
write.csv(metrics_std,  file.path(path_country, paste0(metrics_base, "_std.csv")),  row.names = FALSE)
write.csv(metrics_orig, file.path(path_country, paste0(metrics_base, "_orig.csv")), row.names = FALSE)

cat("\n*** SAVED DFM FIT METRICS (std+orig) TO ***\n", file.path(path_country, paste0(metrics_base, ".rds")), "\n")

# ------------------------------------------------------------------------------
# 6. EM CONVERGENCE PLOT (FULL SAMPLE) + SAVE
# ------------------------------------------------------------------------------

df_em <- tibble(
  iter = seq_along(loglik_full) - 1L,
  loglik = as.numeric(loglik_full)
)

# crit has length (iters_run), loglik length = iters_run + 1
df_crit <- tibble(
  iter = seq_along(crit_full),
  c_j  = as.numeric(crit_full)
)

p_em_loglik <- ggplot(df_em, aes(x = iter, y = loglik)) +
  geom_line(linewidth = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – EM convergence (log-likelihood) | ", country),
    x = "EM iteration",
    y = "Log-likelihood"
  )

p_em_crit <- ggplot(df_crit, aes(x = iter, y = c_j)) +
  geom_line(linewidth = 0.9) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – EM convergence (c_j) | ", country),
    x = "EM iteration",
    y = "Convergence criterion c_j"
  )

path_graph_country <- file.path(path_graph, country)
if (!dir.exists(path_graph_country)) dir.create(path_graph_country, recursive = TRUE)

file_em_loglik <- file.path(
  path_graph_country,
  paste0("DFM_EM_LogLik_", country, "_sel-", params$sel_method, "_Covid-", covid_tag,
         "_r-", r, "_p-", p, ".png")
)
file_em_crit <- file.path(
  path_graph_country,
  paste0("DFM_EM_Crit_", country, "_sel-", params$sel_method, "_Covid-", covid_tag,
         "_r-", r, "_p-", p, ".png")
)

ggsave(file_em_loglik, p_em_loglik, width = 12, height = 6, dpi = 300)
ggsave(file_em_crit,   p_em_crit,   width = 12, height = 6, dpi = 300)

cat("\n*** SAVED EM CONVERGENCE PLOTS TO ***\n", file_em_loglik, "\n", file_em_crit, "\n")

# ------------------------------------------------------------------------------
# 7. GDP SMOOTHER FIT PLOT (STANDARDIZED + ORIGINAL) + SAVE
# ------------------------------------------------------------------------------

plot_dfm_fit_gdp2 <- function(dates, y_obs, y_fit,
                              title = "DFM – In-sample smoother fit (GDP)",
                              ylab  = "GDP") {
  
  dates <- as.Date(dates)
  
  # Fitted: linea continua (mensile)
  df_fit <- tibble(
    date   = dates,
    value  = as.numeric(y_fit),
    series = "Fitted"
  ) %>% filter(!is.na(value))
  
  # Observed: linea TRIMESTRALE (solo quarter-end, quindi senza NA)
  df_obs_q <- tibble(
    date   = dates,
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

# -----------------------------
# STANDARDIZED GDP series aligned
# -----------------------------
p_gdp_std <- plot_dfm_fit_gdp2(
  dates = dates_fit,
  y_obs = Xf_std[, gdp_col],
  y_fit = X_fit_std[, gdp_col],
  title = paste0("DFM – In-sample smoother fit (standardized GDP) | ", country),
  ylab  = "Standardized GDP"
)

# -----------------------------
# ORIGINAL-SCALE GDP series aligned
# -----------------------------
p_gdp_orig <- plot_dfm_fit_gdp2(
  dates = dates_fit,
  y_obs = Xf_orig[, gdp_col],
  y_fit = X_fit_orig[, gdp_col],
  title = paste0("DFM – In-sample smoother fit (GDP original scale) | ", country),
  ylab  = "GDP"
)

print(p_gdp_std)
print(p_gdp_orig)



file_graph_gdp_std <- file.path(
  path_graph_country,
  paste0("DFM_SmootherFit_GDP_STD_", country, "_sel-", params$sel_method,
         "_Covid-", covid_tag, "_r-", r, "_p-", p, ".png")
)
file_graph_gdp_orig <- file.path(
  path_graph_country,
  paste0("DFM_SmootherFit_GDP_ORIG_", country, "_sel-", params$sel_method,
         "_Covid-", covid_tag, "_r-", r, "_p-", p, ".png")
)

ggsave(file_graph_gdp_std,  p_gdp_std,  width = 12, height = 6, dpi = 300)
ggsave(file_graph_gdp_orig, p_gdp_orig, width = 12, height = 6, dpi = 300)

cat("\n*** SAVED GDP SMOOTHER FIT (std+orig) TO ***\n",
    file_graph_gdp_std, "\n", file_graph_gdp_orig, "\n")

# ------------------------------------------------------------------------------
# 8. NOWCAST VISUALIZATION – M1/M2/M3 vs TRUE GDP (ORIGINAL SCALE)
#   Nowcasts are stored BY MONTH (names = monthly dates).
#   Map quarter_end to M1/M2/M3 keys:
#     M1 = q - 2m,  M2 = q - 1m,  M3 = q
# ------------------------------------------------------------------------------

build_nowcast_quarter_df_from_monthly <- function(now_monthly, y_q, dates_q) {
  vM1 <- now_monthly$M1; vM2 <- now_monthly$M2; vM3 <- now_monthly$M3
  
  dM1 <- as.Date(names(vM1))
  dM2 <- as.Date(names(vM2))
  dM3 <- as.Date(names(vM3))
  
  month_key <- function(d) as.Date(floor_date(as.Date(d), unit = "month"))
  
  q_all  <- sort(unique(as.Date(dates_q)))
  m3_key <- month_key(q_all)
  m2_key <- month_key(q_all %m-% months(1))
  m1_key <- month_key(q_all %m-% months(2))
  
  tibble(
    quarter_end = q_all,
    y_true      = y_q[match(q_all, as.Date(dates_q))],
    M1          = vM1[match(m1_key, dM1)],
    M2          = vM2[match(m2_key, dM2)],
    M3          = vM3[match(m3_key, dM3)]
  )
}

df_q <- build_nowcast_quarter_df_from_monthly(
  now_monthly = now_dfm_orig,
  y_q         = y_q,
  dates_q     = dates_q
) %>%
  filter(quarter_end >= as.Date(params$start_eval))

df_q_long <- df_q %>%
  pivot_longer(
    cols      = c(y_true, M1, M2, M3),
    names_to  = "series",
    values_to = "value"
  ) %>%
  mutate(series = factor(
    series,
    levels = c("y_true", "M1", "M2", "M3"),
    labels = c("True GDP", "Nowcast M1", "Nowcast M2", "Nowcast M3")
  ))

plot_dfm_nowcast_M123 <- ggplot(df_q_long, aes(x = quarter_end, y = value, color = series)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  annotate("rect",
           xmin = as.Date(params$covid_start),
           xmax = as.Date(params$covid_end),
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.15) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – Nowcast M1/M2/M3 vs True GDP (", country, ")",
                   " | r_fix=", r_fix, ", p_fix=", p_fix),
    x     = "Quarter",
    y     = "GDP",
    color = NULL
  )

print(plot_dfm_nowcast_M123)

file_graph_dfm_M123 <- file.path(
  path_graph_country,
  paste0(
    "DFM_Nowcast_M123_", country,
    "_sel-", params$sel_method,
    "_Covid-", covid_tag,
    "_rfix-", r_fix,
    "_pfix-", p_fix,
    ".png"
  )
)

ggsave(file_graph_dfm_M123, plot_dfm_nowcast_M123, width = 12, height = 6, dpi = 300)
cat("\n*** SAVED DFM NOWCAST M1/M2/M3 TO ***\n", file_graph_dfm_M123, "\n")

# ------------------------------------------------------------------------------
# 9. NOWCAST EVALUATION: MSE/RMSFE total + pre/post COVID + boxplot
# ------------------------------------------------------------------------------

score_nowcast <- function(df_q, covid_start, covid_end) {
  df_s <- df_q %>%
    mutate(
      err_M1 = M1 - y_true,
      err_M2 = M2 - y_true,
      err_M3 = M3 - y_true,
      period = case_when(
        quarter_end < as.Date(covid_start) ~ "pre-COVID",
        quarter_end > as.Date(covid_end)   ~ "post-COVID",
        TRUE                                ~ "COVID"
      )
    )
  
  summ <- df_s %>%
    group_by(period) %>%
    summarise(
      n_M1 = sum(!is.na(err_M1) & !is.na(y_true)),
      n_M2 = sum(!is.na(err_M2) & !is.na(y_true)),
      n_M3 = sum(!is.na(err_M3) & !is.na(y_true)),
      MSE_M1   = mean(err_M1^2, na.rm = TRUE),
      MSE_M2   = mean(err_M2^2, na.rm = TRUE),
      MSE_M3   = mean(err_M3^2, na.rm = TRUE),
      RMSFE_M1 = sqrt(MSE_M1),
      RMSFE_M2 = sqrt(MSE_M2),
      RMSFE_M3 = sqrt(MSE_M3),
      .groups  = "drop"
    )
  
  total <- df_s %>%
    summarise(
      period  = "total",
      n_M1 = sum(!is.na(err_M1) & !is.na(y_true)),
      n_M2 = sum(!is.na(err_M2) & !is.na(y_true)),
      n_M3 = sum(!is.na(err_M3) & !is.na(y_true)),
      MSE_M1   = mean(err_M1^2, na.rm = TRUE),
      MSE_M2   = mean(err_M2^2, na.rm = TRUE),
      MSE_M3   = mean(err_M3^2, na.rm = TRUE),
      RMSFE_M1 = sqrt(MSE_M1),
      RMSFE_M2 = sqrt(MSE_M2),
      RMSFE_M3 = sqrt(MSE_M3)
    )
  
  list(df_err = df_s, table = bind_rows(total, summ))
}

sc <- score_nowcast(df_q, params$covid_start, params$covid_end)
df_err <- sc$df_err
tab_score <- sc$table

cat("\nNowcast evaluation table (MSE / RMSFE):\n")
print(tab_score)

score_base <- paste0(
  "DFM_NowcastScores_", country,
  "_sel-", params$sel_method,
  "_Covid-", covid_tag,
  "_rfix-", r_fix,
  "_pfix-", p_fix,
  "_", format(params$start_eval, "%Y-%m"),
  "_to_", format(params$end_eval, "%Y-%m")
)

file_score_rds <- file.path(path_country, paste0(score_base, ".rds"))
file_score_csv <- file.path(path_country, paste0(score_base, ".csv"))
saveRDS(tab_score, file_score_rds)
write.csv(tab_score, file_score_csv, row.names = FALSE)
cat("\n*** SAVED NOWCAST SCORES TO ***\n", file_score_rds, "\n", file_score_csv, "\n")

# boxplot of squared errors by horizon and period
df_err_long <- df_err %>%
  select(quarter_end, period, err_M1, err_M2, err_M3) %>%
  pivot_longer(cols = starts_with("err_"), names_to = "horizon", values_to = "err") %>%
  mutate(
    horizon = recode(horizon, err_M1 = "M1", err_M2 = "M2", err_M3 = "M3"),
    se = err^2
  )

plot_dfm_se_box <- ggplot(df_err_long, aes(x = horizon, y = se)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ period, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – Squared forecast errors by horizon (", country, ")",
                   " | r_fix=", r_fix, ", p_fix=", p_fix),
    x = NULL,
    y = "Squared Error"
  )

print(plot_dfm_se_box)

file_graph_dfm_se_box <- file.path(
  path_graph_country,
  paste0(
    "DFM_SE_Boxplot_", country,
    "_sel-", params$sel_method,
    "_Covid-", covid_tag,
    "_rfix-", r_fix,
    "_pfix-", p_fix,
    ".png"
  )
)

ggsave(file_graph_dfm_se_box, plot_dfm_se_box, width = 12, height = 6, dpi = 300)
cat("\n*** SAVED NOWCAST SE BOXPLOT TO ***\n", file_graph_dfm_se_box, "\n")

# ------------------------------------------------------------------------------
# 10. FINAL NOWCAST at quarter-end (M3) vs True GDP
# ------------------------------------------------------------------------------

df_concat_q <- df_q %>%
  transmute(
    quarter_end,
    y_true,
    nowcast_final = M3
  )

plot_dfm_final_nowcast <- ggplot(df_concat_q, aes(x = quarter_end)) +
  geom_line(aes(y = y_true), linewidth = 1) +
  geom_line(aes(y = nowcast_final), linewidth = 1, linetype = 2) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("DFM – Final nowcast at quarter-end (M3) vs True GDP (", country, ")",
                   " | r_fix=", r_fix, ", p_fix=", p_fix),
    x = "Quarter",
    y = "GDP"
  )

print(plot_dfm_final_nowcast)

file_graph_dfm_final <- file.path(
  path_graph_country,
  paste0(
    "DFM_FinalNowcast_vs_True_", country,
    "_sel-", params$sel_method,
    "_Covid-", covid_tag,
    "_rfix-", r_fix,
    "_pfix-", p_fix,
    ".png"
  )
)

ggsave(file_graph_dfm_final, plot_dfm_final_nowcast, width = 12, height = 6, dpi = 300)
cat("\n*** SAVED FINAL NOWCAST PLOT TO ***\n", file_graph_dfm_final, "\n")

