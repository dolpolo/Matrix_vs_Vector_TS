# ==============================================================================
# 0. SETUP: PATHS, LIBRARIES, LOAD TENSOR T-MF-TPRF RESULTS (VECTOR-STYLE)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
path_results <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph_tensor")

if (!dir.exists(path_graph)) dir.create(path_graph, recursive = TRUE)

# ------------------------------------------------------------------------------
# 0.1 Build filename EXACTLY like in saving step (simple)
# ------------------------------------------------------------------------------

cc_lab <- paste(countries, collapse = "-")   # same 'countries' used when saving
r1 <- r_hat[1]; r2 <- r_hat[2]

file_results <- file.path(
  path_results,
  paste0(
    "T_MF_TPRF_tensor",
    "_cc-",      cc_lab,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    "_Lproxy-",  Lproxy,
    "_Lmidas-",  L_midas,
    "_pAR-",     p_ar,
    "_r1-",      r1,
    "_r2-",      r2,
    "_RobustF-", as.integer(isTRUE(params$Robust_F)),
    "_CovidM-",  as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-",  as.integer(isTRUE(params$covid_mask_q)),
    ".rds"
  )
)

res_tensor <- readRDS(file_results)

# ------------------------------------------------------------------------------
# 0.2 Extract objects (same habit as vector)
# ------------------------------------------------------------------------------

params             <- res_tensor$params
tensor             <- res_tensor$tensor
Tensor_MF_TPRF_out <- res_tensor$Tensor_MF_TPRF

countries_eval <- setdiff(res_tensor$countries, params$target_cc)
country_order  <- c("DE","FR","IT","ES")   # optional ordering
countries_eval <- intersect(country_order, countries_eval)

Lproxy <- res_tensor$Lproxy
L_midas <- res_tensor$L_midas
p_ar <- res_tensor$p_ar

cat("\n*** LOADED TENSOR RESULTS FROM ***\n", file_results, "\n")


# ==============================================================================
# 1. TRUE QUARTERLY GDP FOR ALL COUNTRIES (FROM TENSOR)  [vector-style]
# ==============================================================================

Y_tens     <- tensor$Y
gdp_col    <- tensor$target_col
dates_m    <- as.Date(dimnames(Y_tens)[[1]])

gdp_tens   <- Y_tens[, , gdp_col, drop = FALSE]     # [T_m x P x 1]
mask_gdp   <- !is.na(gdp_tens[, , 1])
idx_q      <- which(rowSums(mask_gdp) > 0)          # months where quarterly GDP observed
dates_q    <- dates_m[idx_q]                        # quarterly dates (common)

y_true_q_all <- gdp_tens[idx_q, , 1]                # [T_q x P]
T_q_complete <- nrow(y_true_q_all)

# country labels & reorder
colnames(y_true_q_all) <- countries_eval
y_true_q_all <- y_true_q_all[, countries_eval, drop = FALSE]


# ==============================================================================
# 2. EXTRACT MONTHLY NOWCAST (by_country) + BUILD common structures
# ==============================================================================

stopifnot(!is.null(Tensor_MF_TPRF_out$by_country))

y_now_by_country <- lapply(countries_eval, function(cc) {
  out_cc <- Tensor_MF_TPRF_out$by_country[[cc]]
  if (is.null(out_cc$y_nowcast)) stop("Missing y_nowcast for country: ", cc)
  as.numeric(out_cc$y_nowcast)
})
names(y_now_by_country) <- countries_eval

T_m_full <- length(y_now_by_country[[1]])
if (any(vapply(y_now_by_country, length, 1L) != T_m_full)) {
  stop("Monthly nowcast lengths differ across countries.")
}

# monthly dates aligned with nowcast
dates_m_full <- dates_m[seq_len(T_m_full)]

# in-sample months = 3*T_q_complete (quarters for which GDP observed)
n_in <- 3 * T_q_complete
if (T_m_full < n_in) stop("Not enough months: T_m_full < 3*T_q_complete.")

# period flag exactly like vector
period_vec <- rep("in-sample", T_m_full)
if (T_m_full > n_in) period_vec[(n_in + 1):T_m_full] <- "real-time"
period_vec <- factor(period_vec, levels = c("in-sample", "real-time"))

# indexes M1/M2/M3 (same as vector)
M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)


# ==============================================================================
# 3. PERIOD MASKS (PRE / COVID / POST) ON QUARTERS  [vector-style]
# ==============================================================================

start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

is_PRE   <- (dates_q >= start_eval  & dates_q <  covid_start)
is_COVID <- (dates_q >= covid_start & dates_q <= covid_end)
is_POST  <- (dates_q >  covid_end   & dates_q <= end_eval)
is_ALL   <- (dates_q >= start_eval  & dates_q <= end_eval)

period_labels <- c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")
period_masks  <- list(is_ALL, is_PRE, is_COVID, is_POST)


# ==============================================================================
# 4. PLOT: MONTHLY NOWCAST vs TRUE QUARTERLY GDP (FACET by country)
#    (same idea as vector; just repeated per country)
# ==============================================================================

df_true <- data.frame(
  date    = rep(dates_q, times = length(countries_eval)),
  country = rep(countries_eval, each = T_q_complete),
  value   = as.vector(y_true_q_all),
  series  = "Quarterly GDP"
)

df_now <- bind_rows(lapply(countries_eval, function(cc) {
  data.frame(
    date    = dates_m_full,
    country = cc,
    value   = y_now_by_country[[cc]],
    period  = period_vec
  )
}))

df_true$country <- factor(df_true$country, levels = countries_eval)
df_now$country  <- factor(df_now$country,  levels = countries_eval)

plot_nowcast <- ggplot() +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.15
  ) +
  geom_line(
    data = df_now,
    aes(x = date, y = value, color = period),
    linewidth = 0.9
  ) +
  geom_point(
    data = df_now %>% filter(period == "real-time"),
    aes(x = date, y = value),
    size = 1.6, color = "black"
  ) +
  geom_line(
    data = df_true,
    aes(x = date, y = value, color = "Quarterly GDP"),
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
    title = "Tensor MF-TPRF Monthly Nowcast (In-Sample + Real-Time) — all countries",
    x = "Date",
    y = "GDP growth"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(plot_nowcast)


# ==============================================================================
# 5. PERFORMANCE METRICS — RMSFE (M1, M2, M3) by country and period
#    (same helper as vector)
# ==============================================================================

rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

rmsfe_list <- list()

for (cc in countries_eval) {
  
  y_true_cc <- y_true_q_all[, cc]
  y_now_cc  <- y_now_by_country[[cc]]
  y_now_in  <- y_now_cc[1:n_in]
  
  for (p_id in seq_along(period_labels)) {
    lab  <- period_labels[p_id]
    mask <- period_masks[[p_id]]
    
    RMS_M1 <- rmsfe_period(mask, y_true_cc, y_now_in, M1_idx)
    RMS_M2 <- rmsfe_period(mask, y_true_cc, y_now_in, M2_idx)
    RMS_M3 <- rmsfe_period(mask, y_true_cc, y_now_in, M3_idx)
    
    rmsfe_list[[length(rmsfe_list) + 1L]] <- data.frame(
      country = cc, period = lab,
      M1 = RMS_M1, M2 = RMS_M2, M3 = RMS_M3
    )
  }
}

rmsfe_df <- bind_rows(rmsfe_list)

# print quick summary like vector
for (cc in countries_eval) {
  cat("\nTensor MF-TPRF RMSFE by period –", cc, "\n",
      "--------------------------------------------\n")
  
  tab_cc <- rmsfe_df %>% filter(country == cc)
  
  show_row <- function(lbl) {
    row <- tab_cc %>% filter(period == lbl)
    cat(lbl, ":  M1 =", round(row$M1,4), "  M2 =", round(row$M2,4), "  M3 =", round(row$M3,4), "\n")
  }
  
  show_row("Full sample")
  show_row("Pre-COVID")
  show_row("COVID period")
  show_row("Post-COVID")
}

# ==============================================================================
# 6. LaTeX table (same spirit as vector; one table all countries)
# ==============================================================================

rmsfe_pc <- rmsfe_df %>%
  filter(period != "Full sample") %>%
  mutate(period_short = dplyr::recode(
    period,
    "Pre-COVID"    = "PRE",
    "COVID period" = "COV",
    "Post-COVID"   = "POST"
  )) %>%
  pivot_wider(
    names_from  = period_short,
    values_from = c(M1, M2, M3),
    names_glue  = "{.value}_{period_short}"
  ) %>%
  arrange(factor(country, levels = countries_eval))

latex_tab <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Tensor MF-TPRF RMSFE by period (Pre, COVID, Post) -- all countries}\n",
  "\\label{tab:RMSFE_Tensor_MF_TPRF}\n",
  "\\scriptsize\n",
  "\\begin{tabular}{lccc ccc ccc}\n",
  "\\toprule\n",
  " & \\multicolumn{3}{c}{Pre-COVID} & \\multicolumn{3}{c}{COVID} & \\multicolumn{3}{c}{Post-COVID} \\\\\n",
  "Country & M1 & M2 & M3 & M1 & M2 & M3 & M1 & M2 & M3 \\\\\n",
  "\\midrule\n",
  paste(
    sprintf("%s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
            rmsfe_pc$country,
            rmsfe_pc$M1_PRE,  rmsfe_pc$M2_PRE,  rmsfe_pc$M3_PRE,
            rmsfe_pc$M1_COV,  rmsfe_pc$M2_COV,  rmsfe_pc$M3_COV,
            rmsfe_pc$M1_POST, rmsfe_pc$M2_POST, rmsfe_pc$M3_POST),
    collapse = "\n"
  ),
  "\n\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

cat("\n", latex_tab, "\n")


# ==============================================================================
# 7. SAVE GRAPH (same style as vector)
# ==============================================================================

file_graph_now <- file.path(
  path_graph,
  paste0(
    "T_MF_TPRF_FIT_tensor",
    "_cc-",      cc_lab,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    "_Lproxy-",  Lproxy,
    "_Lmidas-",  L_midas,
    "_pAR-",     p_ar,
    "_r1-",      r1,
    "_r2-",      r2,
    "_RobustF-", as.integer(isTRUE(params$Robust_F)),
    "_CovidM-",  as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-",  as.integer(isTRUE(params$covid_mask_q)),
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

cat("\n*** NOWCAST GRAPH SAVED TO ***\n", file_graph_now, "\n")







# ==============================================================================
#  TENSOR MF-TPRF — PSEUDO REAL-TIME RESULTS (VECTOR-STYLE LOAD + PLOTS + RMSFE)
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ==============================================================================
# 0) PATHS + LOAD PSEUDO-RT RESULTS (NO COPY/PASTE FILE NAME)
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
path_results <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph_tensor_RT")
if (!dir.exists(path_graph)) dir.create(path_graph, recursive = TRUE)

# ---- THESE MUST MATCH WHAT YOU USED WHEN SAVING ----
# (se preferisci, puoi derivarlo da res_rt$countries dopo il load, ma qui serve per costruire il filename)
countries_eval <- c("DE","FR","IT","ES")
cc_lab <- paste(countries_eval, collapse = "-")

# Qui assumo che questi oggetti esistano NEL TUO WORKSPACE del results script:
# params, N_m, N_q, roll_tensor (o almeno gli hyper_pre)
# Se non li hai in memoria, vedi nota in fondo per versione "auto-detect".
file_rt_tensor <- file.path(
  path_results,
  paste0(
    "T_MF_TPRF_RT_tensor",
    "_cc-", cc_lab,
    "_sel-", params$sel_method,
    "_Nm-", N_m,
    "_Nq-", N_q,
    "_Lproxy-", roll_tensor$hyper_pre$Lproxy,
    "_Lmidas-", roll_tensor$hyper_pre$L_midas,
    "_pAR-", roll_tensor$hyper_pre$p_AR,
    "_CvM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CvQ-", as.integer(isTRUE(params$covid_mask_q)),
    ".rds"
  )
)

if (!file.exists(file_rt_tensor)) stop("File pseudo-RT non trovato:\n", file_rt_tensor)
cat(">>> Loading pseudo-RT file:\n", file_rt_tensor, "\n\n")

res_rt <- readRDS(file_rt_tensor)

# ==============================================================================
# 1) EXTRACT META + SANITY CHECKS
# ==============================================================================

params_rt  <- res_rt$params
proxy_name <- if (!is.null(res_rt$proxy_name)) res_rt$proxy_name else "EA"

dates_m <- as.Date(res_rt$dates_m)
dates_q <- as.Date(res_rt$dates_q)

if (is.null(res_rt$pseudo_rt_all)) stop("Nel file manca pseudo_rt_all.")
df_rt_all <- res_rt$pseudo_rt_all

stopifnot(all(c("date","country","nowcast","month_in_quarter") %in% names(df_rt_all)))

df_rt_all <- df_rt_all %>%
  mutate(
    date = as.Date(date),
    type = factor(month_in_quarter, levels = c("M1","M2","M3"))
  ) %>%
  select(date, country, nowcast, type) %>%
  arrange(country, date, type)

countries_in_rt <- sort(unique(df_rt_all$country))
countries_eval  <- intersect(countries_eval, countries_in_rt)
if (length(countries_eval) == 0) stop("Nessun paese valido nei nowcast pseudo-RT.")

# True GDP
if (is.null(res_rt$Y_q_all)) stop("Nel file manca Y_q_all.")
Y_q_all <- res_rt$Y_q_all
if (is.null(colnames(Y_q_all))) stop("Y_q_all deve avere colnames (proxy + paesi).")
if (!(proxy_name %in% colnames(Y_q_all))) stop("proxy_name non in colnames(Y_q_all): ", proxy_name)

# allineamento Y_q_all vs dates_q
if (nrow(Y_q_all) != length(dates_q)) {
  stop("Mismatch Y_q_all vs dates_q: nrow(Y_q_all)=", nrow(Y_q_all),
       " length(dates_q)=", length(dates_q))
}

countries_y <- setdiff(colnames(Y_q_all), proxy_name)
countries_eval <- intersect(countries_eval, countries_y)
if (length(countries_eval) == 0) stop("Nessun paese in comune tra nowcast e Y_q_all.")

# ==============================================================================
# 2) BUILD TRUE GDP DATAFRAME (EVAL WINDOW)
# ==============================================================================

idx_eval_q <- which(dates_q >= params_rt$start_eval & dates_q <= params_rt$end_eval)
if (length(idx_eval_q) == 0) stop("Nessun trimestre nella evaluation window.")

dates_q_eval <- dates_q[idx_eval_q]
Y_q_eval     <- Y_q_all[idx_eval_q, , drop = FALSE]

df_true_all <- data.frame(
  date_q     = rep(dates_q_eval, times = length(countries_eval)),
  country    = rep(countries_eval, each = length(dates_q_eval)),
  GDP        = as.vector(as.matrix(Y_q_eval[, countries_eval, drop = FALSE])),
  quarter_id = rep(paste0(year(dates_q_eval), "Q", quarter(dates_q_eval)),
                   times = length(countries_eval)),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 3) PREP ROLLING NOWCAST: MAP MONTHS -> QUARTER_ID (EVAL WINDOW)
# ==============================================================================

df_rt_q <- df_rt_all %>%
  filter(country %in% countries_eval) %>%
  filter(date >= params_rt$start_eval & date <= params_rt$end_eval) %>%
  mutate(quarter_id = paste0(year(date), "Q", quarter(date))) %>%
  rename(date_m = date)

# ==============================================================================
# 4) JOIN TRUE GDP WITH NOWCASTS + PERIOD LABEL
# ==============================================================================

df_eval_all <- df_true_all %>%
  inner_join(df_rt_q, by = c("country","quarter_id")) %>%
  mutate(
    period = case_when(
      date_q <  params_rt$covid_start ~ "Pre-COVID",
      date_q >= params_rt$covid_start & date_q <= params_rt$covid_end ~ "COVID period",
      date_q >  params_rt$covid_end ~ "Post-COVID",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Pre-COVID","COVID period","Post-COVID"))
  )

if (nrow(df_eval_all) == 0) stop("Join vuoto: controlla date_q/dates_m e quarter_id mapping.")

# ==============================================================================
# 5) PLOT: ROLLING NOWCAST vs TRUE GDP (FACET PER PAESE)
# ==============================================================================

df_true_plot <- df_eval_all %>% select(country, date_q, GDP) %>% distinct()
df_rt_plot   <- df_eval_all %>% select(country, date_m, nowcast, type) %>% distinct() %>%
  rename(date = date_m)

df_true_plot$country <- factor(df_true_plot$country, levels = countries_eval)
df_rt_plot$country   <- factor(df_rt_plot$country,   levels = countries_eval)

plot_tensor_rt <- ggplot() +
  annotate(
    "rect",
    xmin = params_rt$covid_start, xmax = params_rt$covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.20
  ) +
  geom_line(
    data = df_true_plot,
    aes(x = date_q, y = GDP, color = "True GDP"),
    linewidth = 1.0
  ) +
  geom_line(
    data = df_rt_plot,
    aes(x = date, y = nowcast, color = type),
    linewidth = 0.8
  ) +
  geom_point(
    data = df_rt_plot,
    aes(x = date, y = nowcast, color = type),
    size = 1.8
  ) +
  scale_color_manual(
    values = c("True GDP" = "#D62728", "M1"="#1F77B4", "M2"="#2ECC71", "M3"="#F1C40F"),
    name = ""
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title    = "Tensor MF-TPRF – Rolling pseudo real-time nowcasts vs quarterly GDP",
    subtitle = "M1 (early), M2 (mid), M3 (end of quarter)",
    x = "Date", y = "GDP growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(plot_tensor_rt)

Lproxy_used <- res_rt$pseudo_rt_raw$hyper_pre$Lproxy
Lmidas_used <- res_rt$pseudo_rt_raw$hyper_pre$L_midas
pAR_used    <- res_rt$pseudo_rt_raw$hyper_pre$p_AR

file_graph_rt <- file.path(
  path_graph,
  paste0(
    "T_MF_TPRF_RT_Rolling_allCountries",
    "_Lproxy-", Lproxy_used,
    "_Lmidas-", Lmidas_used,
    "_pAR-",    pAR_used,
    "_",        format(params_rt$start_eval, "%Y-%m"),
    "_to_",     format(params_rt$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(filename = file_graph_rt, plot = plot_tensor_rt, width = 12, height = 8, dpi = 300)
cat("\n*** GRAPH SAVED TO ***\n", file_graph_rt, "\n")

# ==============================================================================
# 6) RMSFE BY COUNTRY, PERIOD, TYPE (M1/M2/M3) + FULL SAMPLE
# ==============================================================================

rmsfe_by_cp <- df_eval_all %>%
  filter(!is.na(period)) %>%
  group_by(country, period, type) %>%
  summarise(RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)), .groups = "drop")

rmsfe_full <- df_eval_all %>%
  group_by(country, type) %>%
  summarise(RMSFE = sqrt(mean((nowcast - GDP)^2, na.rm = TRUE)), .groups = "drop") %>%
  mutate(period = "Full sample")

rmsfe_all <- bind_rows(rmsfe_full, rmsfe_by_cp) %>%
  mutate(
    period  = factor(period, levels = c("Full sample","Pre-COVID","COVID period","Post-COVID")),
    type    = factor(type, levels = c("M1","M2","M3")),
    country = factor(country, levels = countries_eval)
  ) %>%
  arrange(country, period, type)

rmsfe_wide <- rmsfe_all %>%
  tidyr::pivot_wider(names_from = type, values_from = RMSFE) %>%
  arrange(country, period)

print(rmsfe_wide)

# ==============================================================================
# 7) LaTeX TABLE – RMSFE
# ==============================================================================

latex_tab_tensor_rt <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\scriptsize\n",
  "\\caption{Tensor MF-TPRF rolling RMSFE by country and period}\n",
  "\\label{tab:RMSFE_Tensor_MF_TPRF_Rolling}\n",
  "\\begin{tabular}{llccc}\n",
  "\\toprule\n",
  "Country & Period & M1 & M2 & M3 \\\\\n",
  "\\midrule\n",
  paste(
    sprintf("%s & %s & %.4f & %.4f & %.4f \\\\",
            rmsfe_wide$country,
            rmsfe_wide$period,
            rmsfe_wide$M1,
            rmsfe_wide$M2,
            rmsfe_wide$M3),
    collapse = "\n"
  ),
  "\n\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

cat(latex_tab_tensor_rt)

