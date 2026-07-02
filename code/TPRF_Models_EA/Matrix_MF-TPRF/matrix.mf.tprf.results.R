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
# 0. PACKAGES, PATHS, AND GRAPH FOLDERS
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(grid)

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
file_out     <- file.path(path_main, "TPRF_Models_EA/Final_Tab_Graph")
path_results <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_func    <- file.path(path_main, "functions/functions_mat")

# main figure folders
path_fig_main        <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/figures_main")
path_fig_realtime    <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/figures_realtime")
path_fig_factors     <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/figures_factor_interpretation")
path_fig_appendix    <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/figures_fixed_rank_appendix")

# real-time subfolders
path_fig_rt_full     <- file.path(path_fig_realtime, "full_sample_realtime")
path_fig_rt_groups   <- file.path(path_fig_realtime, "country_groups")
path_fig_rt_post     <- file.path(path_fig_realtime, "post_covid")
path_fig_rt_bycc     <- file.path(path_fig_realtime, "single_country")

dir.create(path_fig_main,     recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_realtime, recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_factors,  recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_appendix, recursive = TRUE, showWarnings = FALSE)

dir.create(path_fig_rt_full,   recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_rt_groups, recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_rt_post,   recursive = TRUE, showWarnings = FALSE)
dir.create(path_fig_rt_bycc,   recursive = TRUE, showWarnings = FALSE)

source(file.path(path_func, "matrix.mf.tprf.utils.R"))

# ==============================================================================
# 1. RUN IDENTIFIERS
# ==============================================================================

proxy_mode <- "scalar"   # "scalar" | "multivariate"
Size       <- "small"    # "small" | "medium" | "large"
sel        <- "LASSO"    # "corr" | "LASSO"

model_name <- paste0("matrix_", proxy_mode)

# Old fixed-rank appendix result already available on disk.
model_name_fixed <- "matrix"

run_id <- paste(proxy_mode, Size, sel, sep = "_")

proxy_label <- switch(
  proxy_mode,
  scalar = "Scalar EA-residual proxy",
  multivariate = "Multivariate national-GDP proxy"
)

country_order <- c(
  "DE", "FR", "IT", "ES",
  "NL", "BE", "AT", "PT", "EA"
)
# ==============================================================================
# 2. GRAPH HELPERS
# ==============================================================================

theme_paper_plot <- function(base_size = 13) {
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

theme_factor_plot <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(hjust = 0.5, colour = "grey25"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 7)
    )
}

make_year_breaks <- function(x) {
  seq(
    floor_date(min(as.Date(x), na.rm = TRUE), unit = "year"),
    floor_date(max(as.Date(x), na.rm = TRUE), unit = "year"),
    by = "1 year"
  )
}

scale_x_yearly <- function(date_vec) {
  scale_x_date(
    breaks = make_year_breaks(date_vec),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  )
}

# ==============================================================================
# SAVE FIGURES FOR LATEX
# ==============================================================================

save_plot_latex <- function(
    filename,
    plot,
    width = 7,
    height = 4.5,
    dpi = 180,
    bg = "white"
) {
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = bg,
    units    = "in",
    limitsize = FALSE
  )
}

# ==============================================================================
# 3. LOAD RESULTS
# ==============================================================================

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

# Fixed-rank appendix remains the pre-existing benchmark.
file_fit_fixed <- find_result_file(
  path  = path_results,
  model_name_fixed <- model_name,
  stage = "fit_fixed",
  Size  = Size,
  sel   = sel
)

res_fit       <- readRDS(file_fit)
res_rt        <- readRDS(file_rt)
res_fit_fixed <- readRDS(file_fit_fixed)

fit_obj_fixed <- if (!is.null(res_fit_fixed$fit)) {
  res_fit_fixed$fit
} else {
  res_fit_fixed
}

if (is.null(res_fit$proxy_mode) || is.null(res_rt$proxy_mode)) {
  stop("The selected fit or real-time result does not contain proxy_mode.")
}

if (!identical(as.character(res_fit$proxy_mode), proxy_mode)) {
  stop(
    "Full-sample result proxy mode does not match the selected proxy mode. ",
    "Requested: ", proxy_mode,
    " | Loaded: ", res_fit$proxy_mode
  )
}

if (!identical(as.character(res_rt$proxy_mode), proxy_mode)) {
  stop(
    "Real-time result proxy mode does not match the selected proxy mode. ",
    "Requested: ", proxy_mode,
    " | Loaded: ", res_rt$proxy_mode
  )
}

cat("\nLoaded full-sample results from:\n", file_fit, "\n")
cat("\nLoaded pseudo real-time results from:\n", file_rt, "\n")
cat("\nLoaded fixed-rank appendix results from:\n", file_fit_fixed, "\n")

model_name <- as.character(res_fit$model)
proxy_mode <- as.character(res_fit$proxy_mode)
Size       <- as.character(res_fit$Size)
sel        <- as.character(res_fit$sel)

run_id <- paste(proxy_mode, Size, sel, sep = "_")

proxy_label <- switch(
  proxy_mode,
  scalar = "Scalar EA-residual proxy",
  multivariate = "Multivariate national-GDP proxy"
)


hyper_fit_fixed <- list(
  fixed_r = fit_obj_fixed$r_selected,
  Lproxy = if (!is.null(fit_obj_fixed$meta$Lproxy)) {
    fit_obj_fixed$meta$Lproxy
  } else {
    NA_integer_
  },
  L_midas = if (!is.null(fit_obj_fixed$meta$L_midas)) {
    fit_obj_fixed$meta$L_midas
  } else {
    NA_integer_
  },
  p_ar = if (!is.null(fit_obj_fixed$meta$p_AR)) {
    fit_obj_fixed$meta$p_AR
  } else {
    NA_integer_
  },
  r_selected = fit_obj_fixed$r_selected
)

# ==============================================================================
# 4. EXTRACT CORE OBJECTS
# ==============================================================================

required_fit_fields <- c(
  "params",
  "dates_m",
  "dates_q",
  "target",
  "metadata",
  "selections",
  "hyper",
  "fit"
)

missing_fit_fields <- setdiff(required_fit_fields, names(res_fit))

if (length(missing_fit_fields) > 0L) {
  stop(
    "The full-sample result is missing: ",
    paste(missing_fit_fields, collapse = ", ")
  )
}

params <- res_fit$params

df_selection      <- res_fit$selections$table
df_selection_wide <- res_fit$selections$wide
sel_raw           <- res_fit$selections$raw
na_pct            <- res_fit$na_pct

fit_obj <- res_fit$fit

if (is.null(fit_obj$by_country)) {
  stop("fit_obj$by_country is missing.")
}

countries_fit  <- names(fit_obj$by_country)
countries_eval <- intersect(country_order, countries_fit)

if (length(countries_eval) == 0L) {
  stop("No country-specific forecasts found in fit_obj$by_country.")
}

hyper_fit <- res_fit$hyper

Lproxy  <- hyper_fit$Lproxy
L_midas <- hyper_fit$L_midas
p_ar    <- hyper_fit$p_ar

r_selected_fit <- hyper_fit$r_selected_fit

if (is.null(r_selected_fit)) {
  r_selected_fit <- fit_obj$r_selected
}

r1 <- r_selected_fit[1]
r2 <- r_selected_fit[2]

N_m <- res_fit$metadata$N_m
N_q <- res_fit$metadata$N_q
N   <- res_fit$metadata$N

proxy_columns <- res_fit$proxy$columns
# ==============================================================================
# 5. TRUE QUARTERLY GDP
# ==============================================================================

dates_m <- as.Date(res_fit$dates_m)
dates_q <- as.Date(res_fit$dates_q)

y_q_all_saved <- as.matrix(res_fit$target$y_q_all)

if (!all(countries_eval %in% colnames(y_q_all_saved))) {
  stop("Some country GDP series are missing from res_fit$target$y_q_all.")
}

y_true_q_all <- y_q_all_saved[
  ,
  countries_eval,
  drop = FALSE
]

T_q_complete <- nrow(y_true_q_all)

if (length(dates_q) != T_q_complete) {
  stop("dates_q and y_true_q_all have incompatible dimensions.")
}
# ==============================================================================
# 6. FULL-SAMPLE MONTHLY NOWCAST
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
# 6B. FULL-SAMPLE AGGREGATE EA FIT
# ==============================================================================

if (!is.null(fit_obj$aggregate)) {
  
  if (is.null(fit_obj$aggregate$y_nowcast)) {
    stop("fit_obj$aggregate exists but y_nowcast is missing.")
  }
  
  y_now_EA_full <- as.numeric(fit_obj$aggregate$y_nowcast)
  T_m_EA_full   <- length(y_now_EA_full)
  
  dates_m_EA_full <- dates_m[seq_len(T_m_EA_full)]
  
  # True EA quarterly target used in the saved object
  y_EA_true_q <- as.numeric(res_fit$target$y_q_all[, "EA"])
  
  if (length(y_EA_true_q) != length(dates_q)) {
    stop("Length mismatch between EA quarterly target and dates_q.")
  }
  
  n_in_EA <- 3 * length(y_EA_true_q)
  if (T_m_EA_full < n_in_EA) {
    stop("Not enough months in aggregate EA fit: T_m_EA_full < 3*T_q.")
  }
  
  period_EA_vec <- rep("in-sample", T_m_EA_full)
  if (T_m_EA_full > n_in_EA) {
    period_EA_vec[(n_in_EA + 1):T_m_EA_full] <- "real-time"
  }
  period_EA_vec <- factor(period_EA_vec, levels = c("in-sample", "real-time"))
  
  df_now_EA_full <- data.frame(
    date       = dates_m_EA_full,
    country    = "EA",
    y_now_full = y_now_EA_full,
    period     = period_EA_vec
  )
  
  df_quarterly_EA <- data.frame(
    date   = dates_q,
    country = "EA",
    y_true = y_EA_true_q
  )
  
} else {
  
  warning(
    "fit_obj$aggregate is NULL. To plot aggregate EA fit, rerun Tensor_MF_TPRF with forecast_mode = 'aggregate' or 'both'."
  )
  
  df_now_EA_full <- NULL
  df_quarterly_EA <- NULL
}

# ==============================================================================
# 7. PERIOD MASKS
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
# 8. DATA FRAMES FOR FULL-SAMPLE ANALYSIS
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
# 9. FACTOR INTERPRETATION HELPERS
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
  
  list(monthly = as.numeric(f_m), quarterly = as.numeric(f_q))
}

make_factor_df <- function(F_hf, dates_m) {
  r1_loc <- dim(F_hf)[2]
  r2_loc <- dim(F_hf)[3]
  
  bind_rows(lapply(1:r1_loc, function(i) {
    bind_rows(lapply(1:r2_loc, function(j) {
      data.frame(
        date       = as.Date(dates_m),
        row_factor = factor(paste0("Row ", i), levels = rev(paste0("Row ", 1:r1_loc))),
        col_factor = factor(paste0("Col ", j), levels = paste0("Col ", 1:r2_loc)),
        factor     = paste0("F(", i, ",", j, ")"),
        value      = as.numeric(F_hf[, i, j]),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

make_loading_df <- function(Cmat) {
  df_C <- as.data.frame(Cmat)
  df_C$variable <- rownames(Cmat)
  
  df_C %>%
    pivot_longer(
      cols = -variable,
      names_to = "loading",
      values_to = "value"
    )
}

make_corr_heatmap_df <- function(F_q_mat, Y_q_mat, month_label) {
  stopifnot(nrow(F_q_mat) == nrow(Y_q_mat))
  
  out <- expand.grid(
    factor = colnames(F_q_mat),
    series = colnames(Y_q_mat),
    stringsAsFactors = FALSE
  )
  out$corr  <- NA_real_
  out$month <- month_label
  
  for (k in seq_len(ncol(F_q_mat))) {
    for (j in seq_len(ncol(Y_q_mat))) {
      out$corr[out$factor == colnames(F_q_mat)[k] &
                 out$series == colnames(Y_q_mat)[j]] <-
        cor(F_q_mat[, k], Y_q_mat[, j], use = "pairwise.complete.obs")
    }
  }
  
  out$series <- factor(out$series, levels = country_order)
  out
}

orient_RC_signs <- function(Rmat, Cmat, ref_country = "EA") {
  R_out <- Rmat
  C_out <- Cmat
  
  if (!ref_country %in% rownames(R_out)) {
    ref_country <- rownames(R_out)[1]
  }
  
  for (j in seq_len(ncol(R_out))) {
    if (R_out[ref_country, j] > 0) {
      R_out[, j] <- -R_out[, j]
    }
  }
  
  for (j in seq_len(ncol(C_out))) {
    C_out[, j] <- -C_out[, j]
  }
  
  list(R = R_out, C = C_out)
}

# quarterly GDP matrix in desired order
Y_q_corr_all <- as.matrix(res_fit$target$y_q_all)
colnames(Y_q_corr_all)[1] <- "EA"
Y_q_corr_all <- Y_q_corr_all[, intersect(country_order, colnames(Y_q_corr_all)), drop = FALSE]

# restrict factor-GDP correlation analysis to the evaluation window
Y_q_corr_oos <- Y_q_corr_all[is_ALL, , drop = FALSE]

# ==============================================================================
# 10. FACTOR INTERPRETATION: BASELINE
# ==============================================================================

title_factors_base    <- "Estimated latent factor"
subtitle_factors_base <- NULL
title_C_base          <- "Column loadings"
subtitle_C_base       <- NULL
title_corr_base       <- "Factor-GDP correlations"
subtitle_corr_base    <- NULL

F_hf_base     <- fit_obj$factors$F_hf
Rhat_base_raw <- fit_obj$R
Chat_base_raw <- fit_obj$C

sign_base <- orient_RC_signs(Rhat_base_raw, Chat_base_raw, ref_country = "EA")
Rhat_base <- sign_base$R
Chat_base <- sign_base$C

dates_f_m_base <- as.Date(fit_obj$factors$dates_m)

r1_base <- dim(F_hf_base)[2]
r2_base <- dim(F_hf_base)[3]

df_factors_base <- make_factor_df(F_hf_base, dates_f_m_base)

df_factors_base$row_factor <- factor(
  df_factors_base$row_factor,
  levels = paste0("Row ", 1:r1_base)
)

plot_factors_base <- ggplot(df_factors_base, aes(x = date, y = value)) +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.16
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.8, colour = "black") +
  facet_grid(row_factor ~ col_factor, scales = "free_y") +
  scale_x_yearly(dates_f_m_base) +
  labs(
    title = title_factors_base,
    subtitle = subtitle_factors_base,
    x = "Date",
    y = "Factor value"
  ) +
  theme_factor_plot()

print(plot_factors_base)

# ------------------------------------------------------------------------------
# 10.2 Baseline row loadings table
# ------------------------------------------------------------------------------

Rhat_base_tab <- Rhat_base[
  country_order[country_order %in% rownames(Rhat_base)],
  ,
  drop = FALSE
]

tab_R_base <- as.data.frame(round(Rhat_base_tab, 3))
tab_R_base$country <- rownames(tab_R_base)
tab_R_base <- tab_R_base %>%
  select(country, everything())

print(tab_R_base)

# ------------------------------------------------------------------------------
# 10.3 Baseline column loadings
# ------------------------------------------------------------------------------

df_C_base_all <- make_loading_df(Chat_base)

loading_levels_base <- paste0("Column factor ", seq_len(ncol(Chat_base)))
names(loading_levels_base) <- colnames(Chat_base)

df_C_base_all <- df_C_base_all %>%
  mutate(
    loading = dplyr::recode(loading, !!!as.list(loading_levels_base))
  )

var_order_base <- df_C_base_all %>%
  group_by(variable) %>%
  summarise(score = max(abs(value), na.rm = TRUE), .groups = "drop") %>%
  arrange(score)

df_C_base_all$variable <- factor(df_C_base_all$variable, levels = var_order_base$variable)

plot_C_base <- ggplot(df_C_base_all, aes(x = 1, y = variable, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.10) +
  facet_wrap(~ loading, nrow = 1) +
  scale_x_continuous(breaks = NULL) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    name = "Loading",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(55, "pt"),
      barwidth  = unit(10, "pt"),
      frame.colour = "grey70",
      ticks.colour = "grey40"
    )
  ) +
  labs(
    title = title_C_base,
    subtitle = subtitle_C_base,
    x = "",
    y = ""
  ) +
  theme_factor_plot() +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.spacing.x  = unit(0.8, "lines"),
    axis.text.y      = element_text(size = 5.4),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    strip.text       = element_text(face = "bold", size = 11)
  )

print(plot_C_base)

# ------------------------------------------------------------------------------
# 10.4 Baseline factor-GDP correlations
# ------------------------------------------------------------------------------

F1_base <- fit_obj$factors$F1
F2_base <- fit_obj$factors$F2
F3_base <- fit_obj$factors$F3

if (is.null(colnames(F1_base))) {
  factor_names_base <- paste0(
    "F(",
    rep(1:r1_base, each = r2_base), ",",
    rep(1:r2_base, times = r1_base), ")"
  )
  colnames(F1_base) <- factor_names_base
  colnames(F2_base) <- factor_names_base
  colnames(F3_base) <- factor_names_base
}

df_corr_base <- bind_rows(
  make_corr_heatmap_df(F1_base[is_ALL, , drop = FALSE], Y_q_corr_oos, "M1"),
  make_corr_heatmap_df(F2_base[is_ALL, , drop = FALSE], Y_q_corr_oos, "M2"),
  make_corr_heatmap_df(F3_base[is_ALL, , drop = FALSE], Y_q_corr_oos, "M3")
)

df_corr_base$month <- factor(df_corr_base$month, levels = c("M1", "M2", "M3"))

plot_corr_base <- ggplot(df_corr_base, aes(x = series, y = factor, fill = corr)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 3.0) +
  facet_wrap(~ month, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), name = "Corr."
  ) +
  labs(
    title = title_corr_base,
    subtitle = subtitle_corr_base,
    x = "",
    y = ""
  ) +
  theme_factor_plot()

print(plot_corr_base)

# ------------------------------------------------------------------------------
# 10.4B Baseline EA factor to use in final cross-model comparison
# ------------------------------------------------------------------------------

# GDP EA quarterly
y_EA_q <- as.numeric(Y_q_corr_all[, "EA"])

# baseline choice: first matrix factor
factor_matrix_monthly_raw <- as.numeric(F_hf_base[, 1, 1])

# orient sign with respect to EA GDP
factor_matrix_oriented <- orient_factor_to_target(
  f_m = factor_matrix_monthly_raw,
  y_q = y_EA_q
)

factor_matrix_monthly <- factor_matrix_oriented$monthly
factor_matrix_quarterly <- factor_matrix_oriented$quarterly

dates_factor_matrix_m <- as.Date(dates_f_m_base)
dates_factor_matrix_q <- as.Date(dates_q)

df_factor_compare_matrix_q <- data.frame(
  date   = dates_factor_matrix_q,
  GDP    = y_EA_q,
  Matrix = factor_matrix_quarterly
)

# ------------------------------------------------------------------------------
# 10.5 Save baseline factor figures
# ------------------------------------------------------------------------------
file_graph_factors_base <- file.path(
  path_fig_factors,
  paste0("plot_factors_", run_id, ".png")
)

file_graph_C_base <- file.path(
  path_fig_factors,
  paste0("plot_column_loadings_", run_id, ".png")
)

file_graph_corr_base <- file.path(
  path_fig_factors,
  paste0("plot_factor_gdp_corr_", run_id, ".png")
)

save_plot_latex(file_graph_factors_base, plot_factors_base, width = 8, height = 4.5, dpi = 180)
save_plot_latex(file_graph_C_base,       plot_C_base,       width = 7, height = 8.0, dpi = 200)
save_plot_latex(file_graph_corr_base,    plot_corr_base,    width = 8, height = 4.0, dpi = 180)

# ==============================================================================
# 11. FIXED-RANK APPENDIX: r1 = 2, r2 = 3
# ==============================================================================

if (is.null(fit_obj_fixed$factors$F_hf)) stop("Missing fixed model factors$F_hf")
if (is.null(fit_obj_fixed$R)) stop("Missing fixed model R")
if (is.null(fit_obj_fixed$C)) stop("Missing fixed model C")

title_factors_fixed    <- "Estimated latent factors"
subtitle_factors_fixed <- "Fixed-rank specification: r1 = 2, r2 = 3"
title_C_fixed          <- "Column loadings"
subtitle_C_fixed       <- NULL
title_corr_fixed       <- "Factor-GDP correlations"
subtitle_corr_fixed    <- NULL

F_hf_fix     <- fit_obj_fixed$factors$F_hf
Rhat_fix_raw <- fit_obj_fixed$R
Chat_fix_raw <- fit_obj_fixed$C

sign_fix <- orient_RC_signs(Rhat_fix_raw, Chat_fix_raw, ref_country = "EA")
Rhat_fix <- sign_fix$R
Chat_fix <- sign_fix$C

dates_f_m_fix <- as.Date(fit_obj_fixed$factors$dates_m)

r1_fix <- dim(F_hf_fix)[2]
r2_fix <- dim(F_hf_fix)[3]

df_factors_fixed <- make_factor_df(F_hf_fix, dates_f_m_fix)

# ------------------------------------------------------------------------------
# 11.1 Fixed-rank factor paths
# ------------------------------------------------------------------------------

df_factors_fixed$row_factor <- factor(
  df_factors_fixed$row_factor,
  levels = c("Row 1", "Row 2")
)

plot_factors_fixed <- ggplot(df_factors_fixed, aes(x = date, y = value)) +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.16
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_line(linewidth = 0.8, colour = "black") +
  facet_grid(row_factor ~ col_factor, scales = "free_y") +
  scale_x_yearly(dates_f_m_fix) +
  labs(
    title = title_factors_fixed,
    subtitle = subtitle_factors_fixed,
    x = "Date",
    y = "Factor value"
  ) +
  theme_factor_plot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )

print(plot_factors_fixed)

# ------------------------------------------------------------------------------
# 11.2 Fixed-rank row loadings table
# ------------------------------------------------------------------------------

Rhat_fix_tab <- Rhat_fix[
  country_order[country_order %in% rownames(Rhat_fix)],
  ,
  drop = FALSE
]

tab_R_fix <- as.data.frame(round(Rhat_fix_tab, 3))
tab_R_fix$country <- rownames(tab_R_fix)
tab_R_fix <- tab_R_fix %>%
  select(country, everything())

print(tab_R_fix)

# ------------------------------------------------------------------------------
# 11.3 Fixed-rank column loadings
# ------------------------------------------------------------------------------

df_C_fix_all <- make_loading_df(Chat_fix)

var_order_fix <- df_C_fix_all %>%
  group_by(variable) %>%
  summarise(score = max(abs(value), na.rm = TRUE), .groups = "drop") %>%
  arrange(score)

df_C_fix_all$variable <- factor(df_C_fix_all$variable, levels = var_order_fix$variable)

plot_C_fixed <- ggplot(df_C_fix_all, aes(x = loading, y = variable, fill = value)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, name = "Loading"
  ) +
  labs(
    title = title_C_fixed,
    subtitle = subtitle_C_fixed,
    x = "",
    y = ""
  ) +
  theme_factor_plot() +
  theme(axis.text.y = element_text(size = 5.1))

print(plot_C_fixed)

# ------------------------------------------------------------------------------
# 11.4 Fixed-rank factor-GDP correlations
# ------------------------------------------------------------------------------

F1_fix <- fit_obj_fixed$factors$F1
F2_fix <- fit_obj_fixed$factors$F2
F3_fix <- fit_obj_fixed$factors$F3

if (is.null(colnames(F1_fix))) {
  factor_names_fix <- paste0(
    "F(",
    rep(1:r1_fix, each = r2_fix), ",",
    rep(1:r2_fix, times = r1_fix), ")"
  )
  colnames(F1_fix) <- factor_names_fix
  colnames(F2_fix) <- factor_names_fix
  colnames(F3_fix) <- factor_names_fix
}

factor_levels_fix <- c(
  paste0("F(1,", 1:r2_fix, ")"),
  paste0("F(2,", 1:r2_fix, ")")
)

df_corr_fixed <- bind_rows(
  make_corr_heatmap_df(F1_fix[is_ALL, , drop = FALSE], Y_q_corr_oos, "M1"),
  make_corr_heatmap_df(F2_fix[is_ALL, , drop = FALSE], Y_q_corr_oos, "M2"),
  make_corr_heatmap_df(F3_fix[is_ALL, , drop = FALSE], Y_q_corr_oos, "M3")
)

df_corr_fixed$month  <- factor(df_corr_fixed$month, levels = c("M1", "M2", "M3"))
df_corr_fixed$factor <- factor(df_corr_fixed$factor, levels = rev(factor_levels_fix))

plot_corr_fixed <- ggplot(df_corr_fixed, aes(x = series, y = factor, fill = corr)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 2.8) +
  facet_wrap(~ month, ncol = 3) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1, 1), name = "Corr."
  ) +
  labs(
    title = title_corr_fixed,
    subtitle = subtitle_corr_fixed,
    x = "",
    y = ""
  ) +
  theme_factor_plot()

print(plot_corr_fixed)

# ------------------------------------------------------------------------------
# 11.5 Save fixed-rank figures
# ------------------------------------------------------------------------------

file_graph_factors_fixed <- file.path(
  path_fig_appendix,
  paste0("plot_factors_fixed_r1-2_r2-3_Size-", Size, "_sel-", sel, ".png")
)

file_graph_C_fixed <- file.path(
  path_fig_appendix,
  paste0("plot_column_loadings_fixed_allvars_r1-2_r2-3_Size-", Size, "_sel-", sel, ".png")
)

file_graph_corr_fixed <- file.path(
  path_fig_appendix,
  paste0("plot_factor_gdp_corr_fixed_oos_M123_r1-2_r2-3_Size-", Size, "_sel-", sel, ".png")
)

save_plot_latex(file_graph_factors_fixed, plot_factors_fixed, width = 8, height = 5.0, dpi = 180)
save_plot_latex(file_graph_C_fixed,       plot_C_fixed,       width = 7, height = 8.0, dpi = 200)
save_plot_latex(file_graph_corr_fixed,    plot_corr_fixed,    width = 8, height = 4.0, dpi = 180)

cat("\nSaved factor interpretation graphs.\n")
cat("Baseline factor figures saved in:\n", path_fig_factors, "\n")
cat("Fixed-rank appendix figures saved in:\n", path_fig_appendix, "\n")

# ==============================================================================
# 12. FULL-SAMPLE PLOT
# ==============================================================================

title_nowcast_full    <- "Monthly nowcasts and observed GDP"
subtitle_nowcast_full <- NULL

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
  scale_x_yearly(df_now_full$date) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = title_nowcast_full,
    subtitle = subtitle_nowcast_full,
    x = "Date",
    y = "GDP growth"
  ) +
  theme_paper_plot(base_size = 13)

print(plot_nowcast)

# ==============================================================================
# 12B. FULL-SAMPLE PLOT: AGGREGATE EA
# ==============================================================================

if (!is.null(df_now_EA_full) && !is.null(df_quarterly_EA)) {
  
  title_nowcast_EA_full    <- "Aggregate EA fit"
  subtitle_nowcast_EA_full <- "Matrix MF-TPRF aggregate-target specification"
  
  plot_nowcast_EA <- ggplot() +
    annotate(
      "rect",
      xmin = covid_start, xmax = covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey70", alpha = 0.15
    ) +
    geom_line(
      data = df_now_EA_full,
      aes(x = date, y = y_now_full, color = period),
      linewidth = 0.9
    ) +
    geom_line(
      data = df_quarterly_EA,
      aes(x = date, y = y_true, color = "Quarterly EA GDP"),
      linewidth = 1.0
    ) +
    scale_color_manual(
      values = c(
        "in-sample"        = "#1F77B4",
        "real-time"        = "#2ECC71",
        "Quarterly EA GDP" = "#D62728"
      ),
      name = ""
    ) +
    scale_x_yearly(df_now_EA_full$date) +
    labs(
      title = title_nowcast_EA_full,
      subtitle = subtitle_nowcast_EA_full,
      x = "Date",
      y = "EA GDP growth"
    ) +
    theme_paper_plot(base_size = 13)
  
  print(plot_nowcast_EA)
  
  file_graph_now_EA <- file.path(
    path_fig_main,
    paste0(
      "plot_fit_EA_aggregate_",
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
  
  save_plot_latex(
    filename = file_graph_now_EA,
    plot     = plot_nowcast_EA,
    width    = 8,
    height   = 6,
    dpi      = 180
  )
  
  cat("\nSaved aggregate EA full-sample graph to:\n", file_graph_now_EA, "\n")
}

# ==============================================================================
# 13. IN-SAMPLE RMSFE
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
  df = rmsfe_insample,
  caption = paste0(
    "Matrix MF--TPRF in-sample RMSFE by country and period (",
    proxy_label,
    ", ",
    Size,
    ", sel = ",
    sel,
    ")"
  ),
  label = paste0(
    "tab:matrix_mf_tprf_insample_",
    run_id
  )
)

cat("\n", latex_tab_insample, "\n")

dir.create(file_out, recursive = TRUE, showWarnings = FALSE)

latex_selection <- save_selection_wide_to_latex(
  df_selection_wide = df_selection_wide,
  file = file.path(
    file_out,
    paste0("selection_wide_table_", run_id, ".tex")
  ),
  caption = paste0(
    "Country-specific variable selection (",
    proxy_label,
    ", ",
    Size,
    ", sel = ",
    sel,
    ")"
  ),
  label = paste0(
    "tab:selection_wide_",
    run_id
  )
)

cat(latex_selection)

# ==============================================================================
# 14. SAVE FULL-SAMPLE GRAPH
# ==============================================================================

file_graph_now <- file.path(
  path_fig_main,
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

save_plot_latex(
  filename = file_graph_now,
  plot     = plot_nowcast,
  width    = 8,
  height   = 6,
  dpi      = 180
)

cat("\nSaved full-sample graph to:\n", file_graph_now, "\n")

# ==============================================================================
# 15. BUILD PSEUDO REAL-TIME DATA
# ==============================================================================

params_rt <- res_rt$params

df_rt <- res_rt$pseudo_rt_all %>%
  mutate(
    date    = as.Date(date),
    type    = factor(month_in_quarter, levels = c("M1", "M2", "M3")),
    country = as.character(country)
  ) %>%
  select(date, country, nowcast, type) %>%
  arrange(country, date, type)

# ------------------------------------------------------------
# Add aggregate EA real-time nowcasts, if available
# ------------------------------------------------------------

df_rt_EA_aggregate <- NULL

if (!is.null(res_rt$pseudo_rt_aggregate)) {
  
  df_rt_EA_aggregate <- res_rt$pseudo_rt_aggregate %>%
    mutate(
      date    = as.Date(date),
      type    = factor(month_in_quarter, levels = c("M1", "M2", "M3")),
      country = "EA"
    ) %>%
    select(date, country, nowcast, type) %>%
    arrange(date, type)
  
} else if (!is.null(res_rt$pseudo_rt_raw$aggregate)) {
  
  A <- res_rt$pseudo_rt_raw$aggregate
  
  df_A_M1 <- data.frame(
    date = as.Date(names(A$M1)),
    country = "EA",
    nowcast = as.numeric(A$M1),
    type = factor("M1", levels = c("M1", "M2", "M3"))
  )
  
  df_A_M2 <- data.frame(
    date = as.Date(names(A$M2)),
    country = "EA",
    nowcast = as.numeric(A$M2),
    type = factor("M2", levels = c("M1", "M2", "M3"))
  )
  
  df_A_M3 <- data.frame(
    date = as.Date(names(A$M3)),
    country = "EA",
    nowcast = as.numeric(A$M3),
    type = factor("M3", levels = c("M1", "M2", "M3"))
  )
  
  df_rt_EA_aggregate <- bind_rows(df_A_M1, df_A_M2, df_A_M3) %>%
    filter(!is.na(date), is.finite(nowcast)) %>%
    arrange(date, type)
}

if (!is.null(df_rt_EA_aggregate) && nrow(df_rt_EA_aggregate) > 0) {
  
  # avoid duplicate EA rows if already present in pseudo_rt_all
  df_rt <- df_rt %>%
    filter(country != "EA") %>%
    bind_rows(df_rt_EA_aggregate)
}

df_rt <- df_rt %>%
  mutate(country = factor(country)) %>%
  arrange(country, date, type)

countries_rt      <- unique(as.character(df_rt$country))
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
# 16. PSEUDO REAL-TIME PLOT
# ==============================================================================

title_rt_full    <- "Rolling real-time nowcasts"
subtitle_rt_full <- "Quarterly targets and monthly nowcast updates"

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
  scale_x_yearly(df_rt$date) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = title_rt_full,
    subtitle = subtitle_rt_full,
    x = "Date",
    y = "GDP growth"
  ) +
  theme_paper_plot(base_size = 14)

print(plot_rt)

# ==============================================================================
# 16B. PSEUDO REAL-TIME PLOT: AGGREGATE EA
# ==============================================================================

if ("EA" %in% as.character(df_rt$country)) {
  
  df_rt_EA <- df_rt %>%
    filter(country == "EA")
  
  df_yq_EA <- df_yq %>%
    filter(country == "EA")
  
  title_rt_EA    <- "Rolling real-time nowcasts: aggregate EA"
  subtitle_rt_EA <- "Matrix MF-TPRF aggregate-target specification"
  
  plot_rt_EA <- ggplot() +
    annotate(
      "rect",
      xmin = params_rt$covid_start, xmax = params_rt$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey80", alpha = 0.20
    ) +
    geom_line(
      data = df_yq_EA,
      aes(x = date, y = GDP, color = "True EA GDP"),
      linewidth = 1.0
    ) +
    geom_line(
      data = df_rt_EA,
      aes(x = date, y = nowcast, color = type),
      linewidth = 0.8
    ) +
    geom_point(
      data = df_rt_EA,
      aes(x = date, y = nowcast, color = type),
      size = 1.6
    ) +
    scale_color_manual(
      values = c(
        "True EA GDP" = "#D62728",
        "M1"          = "#1F77B4",
        "M2"          = "#2ECC71",
        "M3"          = "#F1C40F"
      ),
      name = "Series"
    ) +
    scale_x_yearly(df_rt_EA$date) +
    labs(
      title = title_rt_EA,
      subtitle = subtitle_rt_EA,
      x = "Date",
      y = "EA GDP growth"
    ) +
    theme_paper_plot(base_size = 14)
  
  print(plot_rt_EA)
  
  file_graph_rt_EA <- file.path(
    path_fig_rt_bycc,
    paste0(
      "plot_rt_EA_aggregate_",
      model_name,
      "_Size-", Size,
      "_sel-", sel,
      ".png"
    )
  )
  
  save_plot_latex(
    filename = file_graph_rt_EA,
    plot     = plot_rt_EA,
    width    = 7,
    height   = 4,
    dpi      = 180
  )
  
  cat("\nSaved aggregate EA real-time graph to:\n", file_graph_rt_EA, "\n")
}

# ==============================================================================
# 17. SPLIT REAL-TIME PLOTS
# ==============================================================================

title_rt_big4    <- "Rolling real-time nowcasts: DE, FR, IT, ES"
subtitle_rt_big4 <- NULL

title_rt_other    <- "Rolling real-time nowcasts: NL, BE, AT, PT, EA"
subtitle_rt_other <- NULL

countries_rt_big4  <- c("DE", "FR", "IT", "ES")
countries_rt_other <- c("NL", "BE", "AT", "PT", "EA")

df_rt_big4 <- df_rt %>% filter(country %in% countries_rt_big4)
df_yq_big4 <- df_yq %>% filter(country %in% countries_rt_big4)

df_rt_other <- df_rt %>% filter(country %in% countries_rt_other)
df_yq_other <- df_yq %>% filter(country %in% countries_rt_other)

plot_rt_big4 <- ggplot() +
  annotate(
    "rect",
    xmin = params_rt$covid_start, xmax = params_rt$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  geom_line(
    data = df_yq_big4,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.0
  ) +
  geom_line(
    data = df_rt_big4,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.8
  ) +
  geom_point(
    data = df_rt_big4,
    aes(x = date, y = nowcast, color = type),
    size = 1.4
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
  scale_x_yearly(df_rt_big4$date) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = title_rt_big4,
    subtitle = subtitle_rt_big4,
    x = "Date",
    y = "GDP growth"
  ) +
  theme_paper_plot(base_size = 14)

print(plot_rt_big4)

plot_rt_other <- ggplot() +
  annotate(
    "rect",
    xmin = params_rt$covid_start, xmax = params_rt$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  geom_line(
    data = df_yq_other,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.0
  ) +
  geom_line(
    data = df_rt_other,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.8
  ) +
  geom_point(
    data = df_rt_other,
    aes(x = date, y = nowcast, color = type),
    size = 1.4
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
  scale_x_yearly(df_rt_other$date) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = title_rt_other,
    subtitle = subtitle_rt_other,
    x = "Date",
    y = "GDP growth"
  ) +
  theme_paper_plot(base_size = 14)

print(plot_rt_other)

file_graph_rt_big4 <- file.path(
  path_fig_rt_groups,
  paste0(
    "plot_rt_big4_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    ".png"
  )
)

file_graph_rt_other <- file.path(
  path_fig_rt_groups,
  paste0(
    "plot_rt_other_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    ".png"
  )
)

save_plot_latex(
  filename = file_graph_rt_big4,
  plot     = plot_rt_big4,
  width    = 8,
  height   = 6,
  dpi      = 180
)

save_plot_latex(
  filename = file_graph_rt_other,
  plot     = plot_rt_other,
  width    = 8,
  height   = 6,
  dpi      = 180
)

# ==============================================================================
# 18. COUNTRY-BY-COUNTRY REAL-TIME PLOTS
# ==============================================================================

plot_rt_country_list <- list()
graph_titles_country_rt <- list()

for (cc in countries_eval_rt) {
  
  df_rt_cc <- dplyr::filter(df_rt, country == cc)
  df_yq_cc <- dplyr::filter(df_yq, country == cc)
  
  title_rt_cc    <- paste("Rolling real-time nowcasts:", cc)
  subtitle_rt_cc <- NULL
  
  p_cc <- ggplot() +
    annotate(
      "rect",
      xmin = params_rt$covid_start, xmax = params_rt$covid_end,
      ymin = -Inf, ymax = Inf,
      fill = "grey80", alpha = 0.20
    ) +
    geom_line(
      data = df_yq_cc,
      aes(x = date, y = GDP, color = "True GDP"),
      linewidth = 1.0
    ) +
    geom_line(
      data = df_rt_cc,
      aes(x = date, y = nowcast, color = type),
      linewidth = 0.9
    ) +
    geom_point(
      data = df_rt_cc,
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
    scale_x_yearly(df_rt_cc$date) +
    labs(
      title = title_rt_cc,
      subtitle = subtitle_rt_cc,
      x = "Date",
      y = "GDP growth"
    ) +
    theme_paper_plot(base_size = 14)
  
  print(p_cc)
  
  plot_rt_country_list[[cc]] <- p_cc
  graph_titles_country_rt[[cc]] <- list(
    title    = title_rt_cc,
    subtitle = subtitle_rt_cc
  )
  
  save_plot_latex(
    filename = file.path(
      path_fig_rt_bycc,
      paste0(
        "plot_rt_",
        cc,
        "_", model_name,
        "_Size-", Size,
        "_sel-", sel,
        ".png"
      )
    ),
    plot   = p_cc,
    width  = 7,
    height = 4,
    dpi    = 180
  )
}

# ==============================================================================
# 19. POST-COVID REAL-TIME PLOT FOR THE 8 COUNTRIES
# ==============================================================================

countries_post8 <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")

df_rt_post8 <- df_rt %>%
  filter(country %in% countries_post8, date > params_rt$covid_end) %>%
  mutate(country = factor(as.character(country), levels = countries_post8))

df_yq_post8 <- df_yq %>%
  filter(country %in% countries_post8, date > params_rt$covid_end) %>%
  mutate(country = factor(as.character(country), levels = countries_post8))

title_rt_post8 <- "Post-COVID rolling real-time nowcasts"
subtitle_rt_post8 <- paste0(
  "Sample: ",
  format(min(df_rt_post8$date, na.rm = TRUE), "%b %Y"),
  "–",
  format(max(df_rt_post8$date, na.rm = TRUE), "%b %Y")
)

plot_rt_post8 <- ggplot() +
  geom_line(
    data = df_yq_post8,
    aes(x = date, y = GDP, color = "True GDP"),
    linewidth = 1.0
  ) +
  geom_line(
    data = df_rt_post8,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.8
  ) +
  geom_point(
    data = df_rt_post8,
    aes(x = date, y = nowcast, color = type),
    size = 1.4
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
  scale_x_date(
    breaks = seq(
      floor_date(min(df_rt_post8$date, na.rm = TRUE), unit = "year"),
      floor_date(max(df_rt_post8$date, na.rm = TRUE), unit = "year"),
      by = "1 year"
    ),
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = title_rt_post8,
    subtitle = subtitle_rt_post8,
    x = "Date",
    y = "GDP growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey25")
  )

print(plot_rt_post8)

file_graph_rt_post8 <- file.path(
  path_fig_rt_post,
  paste0(
    "plot_rt_postcovid_8countries_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    ".png"
  )
)

save_plot_latex(
  filename = file_graph_rt_post8,
  plot     = plot_rt_post8,
  width    = 12,
  height   = 8,
  dpi      = 300
)


cat("\nSaved post-COVID 8-country rolling graph to:\n", file_graph_rt_post8, "\n")

# ==============================================================================
# 20. ROLLING RMSFE
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
  df = rmsfe_rt_wide,
  caption = paste0(
    "Matrix MF--TPRF rolling RMSFE by country and period (",
    proxy_label,
    ", ",
    Size,
    ", sel = ",
    sel,
    ")"
  ),
  label = paste0(
    "tab:matrix_mf_tprf_rt_",
    run_id
  )
)

cat("\n", latex_tab_rt, "\n")

# ==============================================================================
# 20B. ROLLING RMSFE: AGGREGATE EA
# ==============================================================================

if ("EA" %in% as.character(df_rt$country)) {
  
  rmsfe_rt_EA <- rmsfe_rt_wide %>%
    filter(country == "EA")
  
  print(rmsfe_rt_EA)
  
  latex_tab_rt_EA <- list_to_latex_table(
    df = rmsfe_rt_EA,
    caption = paste0(
      "Matrix MF--TPRF aggregate EA rolling RMSFE (",
      proxy_label,
      ", ",
      Size,
      ", sel = ",
      sel,
      ")"
    ),
    label = paste0(
      "tab:matrix_mf_tprf_EA_aggregate_rt_",
      run_id
    )
  )
  
  cat("\n", latex_tab_rt_EA, "\n")
}

# ==============================================================================
# 21. SAVE ROLLING GRAPH
# ==============================================================================

hyper_rt  <- if (!is.null(res_rt$hyper)) res_rt$hyper else res_rt$pseudo_rt_raw
Lproxy_rt <- hyper_rt$pre$Lproxy
Lmidas_rt <- hyper_rt$pre$L_midas
pAR_rt    <- hyper_rt$pre$p_AR

file_graph_rt_full <- file.path(
  path_fig_rt_full,
  paste0(
    "plot_rt_full_oos_",
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

save_plot_latex(
  filename = file_graph_rt_full,
  plot     = plot_rt,
  width    = 8,
  height   = 6,
  dpi      = 180
)

cat("\nSaved rolling graph to:\n", file_graph_rt_full, "\n")

# ==============================================================================
# 22. GRAPH TITLES / SUBTITLES REGISTER
# ==============================================================================

graph_titles <- list(
  baseline = list(
    factors = list(
      title    = title_factors_base,
      subtitle = subtitle_factors_base
    ),
    column_loadings = list(
      title    = title_C_base,
      subtitle = subtitle_C_base
    ),
    factor_gdp_correlations = list(
      title    = title_corr_base,
      subtitle = subtitle_corr_base
    )
  ),
  fixed_rank_appendix = list(
    factors = list(
      title    = title_factors_fixed,
      subtitle = subtitle_factors_fixed
    ),
    column_loadings = list(
      title    = title_C_fixed,
      subtitle = subtitle_C_fixed
    ),
    factor_gdp_correlations = list(
      title    = title_corr_fixed,
      subtitle = subtitle_corr_fixed
    )
  ),
  main_full_sample = list(
    monthly_nowcasts_vs_gdp = list(
      title    = title_nowcast_full,
      subtitle = subtitle_nowcast_full
    )
  ),
  realtime = list(
    full_facet = list(
      title    = title_rt_full,
      subtitle = subtitle_rt_full
    ),
    big4 = list(
      title    = title_rt_big4,
      subtitle = subtitle_rt_big4
    ),
    other = list(
      title    = title_rt_other,
      subtitle = subtitle_rt_other
    ),
    by_country = graph_titles_country_rt,
    post_covid_8 = list(
      title    = title_rt_post8,
      subtitle = subtitle_rt_post8
    )
  )
)

# ==============================================================================
# 23. SAVE SUMMARY OBJECT
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
  model_id = "T_MF_TPRF",
  model    = model_name,
  stage    = "cross_country",
  Size     = Size,
  sel      = sel,
  
  proxy_mode  = proxy_mode,
  proxy_label = proxy_label,
  run_id      = run_id,
  
  params = params_rt,
  
  hyper = hyper_summary,
  graph_titles = graph_titles,
  
  factor_interpretation = list(
    baseline = list(
      df_factors_base         = df_factors_base,
      tab_R_base              = tab_R_base,
      df_C_base_all           = df_C_base_all,
      df_corr_base            = df_corr_base,
      plot_factors_base       = plot_factors_base,
      plot_C_base             = plot_C_base,
      plot_corr_base          = plot_corr_base,
      file_graph_factors_base = file_graph_factors_base,
      file_graph_C_base       = file_graph_C_base,
      file_graph_corr_base    = file_graph_corr_base
    ),
    
    fixed_rank_appendix = list(
      file_fit_fixed           = file_fit_fixed,
      r_selected_fixed         = fit_obj_fixed$r_selected,
      hyper_fixed              = hyper_fit_fixed,
      df_factors_fixed         = df_factors_fixed,
      tab_R_fix                = tab_R_fix,
      df_C_fix_all             = df_C_fix_all,
      df_corr_fixed            = df_corr_fixed,
      plot_factors_fixed       = plot_factors_fixed,
      plot_C_fixed             = plot_C_fixed,
      plot_corr_fixed          = plot_corr_fixed,
      file_graph_factors_fixed = file_graph_factors_fixed,
      file_graph_C_fixed       = file_graph_C_fixed,
      file_graph_corr_fixed    = file_graph_corr_fixed
    )
  ),
  
  factor_comparison = list(
    dates_m = dates_factor_matrix_m,
    dates_q = dates_factor_matrix_q,
    gdp_q   = y_EA_q,
    
    matrix_factor_monthly   = factor_matrix_monthly,
    matrix_factor_quarterly = factor_matrix_quarterly,
    
    factor_name = "Matrix MF-TPRF",
    
    df_q = df_factor_compare_matrix_q
  ),
  
  countries = countries_eval_rt,
  
  file_graph_rt_big4  = file_graph_rt_big4,
  file_graph_rt_other = file_graph_rt_other,
  plot_rt_big4        = plot_rt_big4,
  plot_rt_other       = plot_rt_other,
  plot_rt_by_country  = plot_rt_country_list,
  
  plot_rt_post8       = plot_rt_post8,
  file_graph_rt_post8 = file_graph_rt_post8,
  
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
  
  file_fit          = file_fit,
  file_rt           = file_rt,
  file_fit_fixed    = file_fit_fixed,
  file_graph_fit    = file_graph_now,
  file_graph_rt     = file_graph_rt_full
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
  include_details  = FALSE
)

if (file.exists(file_summary_matrix)) file.remove(file_summary_matrix)

saveRDS(summary_tensor_out, file_summary_matrix)

cat("\nSaved summary object to:\n", file_summary_matrix, "\n")

