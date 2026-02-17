# ==============================================================================
# 0. SETUP: PATHS, LIBRARIES, LOAD TENSOR T-MF-TPRF RESULTS
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
path_results <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs")
path_graph_tensor <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph_tensor")

if (!dir.exists(path_graph_tensor)) {
  dir.create(path_graph_tensor, recursive = TRUE)
}

# ===>>> metti qui esattamente il nome del file che hai salvato prima
file_tensor <- file.path(
  path_results,
"T-MF_TPRF_tensor_cc-DE-FR-IT-ES-EA_sel-none_Lproxy-1_Lmidas-2_Nm-34_Nq-11_pAR-1_robustF1_rtype-NW_CvM1_CvQ0_tgtcc-EA_r1-1_r2-1_2000-04_to_2025-10.rds"
)

res_tensor <- readRDS(file_tensor)

# Oggetti salvati
params             <- res_tensor$params
countries_eval     <- setdiff(res_tensor$countries, params$target_cc)
# ordine desiderato
country_order <- c("DE","FR","IT","ES")
# riordina countries_eval (e tieni solo quelli presenti)
countries_eval <- intersect(country_order, countries_eval)

tensor             <- res_tensor$tensor
Tensor_MF_TPRF_out <- res_tensor$Tensor_MF_TPRF
Lproxy             <- res_tensor$Lproxy
L_midas            <- res_tensor$L_midas

# Se hai salvato anche r:
if (!is.null(res_tensor$r)) {
  r <- res_tensor$r
} else {
  r <- c(NA_integer_, NA_integer_)
}

# EM diffs per la convergenza
EM_diffs <- res_tensor$EM_diffs


# ==============================================================================
# 1. TRUE QUARTERLY GDP FOR ALL COUNTRIES (FROM TENSOR)
# ==============================================================================

Y_tens      <- tensor$Y                        # [T_m x P x (N_m+N_q+1)]
gdp_col     <- tensor$target_col
dates_m_ts  <- as.Date(dimnames(Y_tens)[[1]])  # date mensili comuni
P           <- dim(Y_tens)[2]                  # numero paesi

# GDP tensoriale [T_m x P x 1]
gdp_tens <- Y_tens[, , gdp_col, drop = FALSE]

# Indici mesi che corrispondono a fine trimestre
mask_gdp <- !is.na(gdp_tens[, , 1])
idx_q    <- which(rowSums(mask_gdp) > 0)

# Date trimestrali comuni
dates_q_ts <- dates_m_ts[idx_q]

# Matrice T_q x P con GDP di ciascun paese
y_true_q_all <- gdp_tens[idx_q, , 1]          # [T_q x P]
T_q_complete <- nrow(y_true_q_all)

# Etichette paesi
# etichette + riordino coerente di y_true_q_all
stopifnot(ncol(y_true_q_all) == length(countries_eval))
colnames(y_true_q_all) <- countries_eval
y_true_q_all <- y_true_q_all[, countries_eval, drop = FALSE]



# ==============================================================================
# 2. MONTHLY NOWCAST (TENSOR MF-TPRF) + PERIODI  [CORRETTO: by_country]
# ==============================================================================

# 2.1 Estrai nowcast mensile per paese
stopifnot(!is.null(Tensor_MF_TPRF_out$by_country))

y_now_by_country <- lapply(countries_eval, function(cc) {
  out_cc <- Tensor_MF_TPRF_out$by_country[[cc]]
  if (is.null(out_cc) || is.null(out_cc$y_nowcast)) {
    stop("Manca y_nowcast per il paese: ", cc)
  }
  as.numeric(out_cc$y_nowcast)
})
names(y_now_by_country) <- countries_eval

# 2.2 Lunghezza mensile (assumo uguale per tutti; controllo)
T_m_full <- length(y_now_by_country[[1]])
if (any(vapply(y_now_by_country, length, 1L) != T_m_full)) {
  stop("Le lunghezze dei nowcast mensili non coincidono tra i paesi.")
}

if (T_m_full < 3 * T_q_complete) {
  stop("T_m_full ha meno mesi di 3 * T_q_complete.")
}

dates_m_full <- dates_m_ts[seq_len(T_m_full)]

# Numero di mesi in-sample (trimestri completi)
n_in <- 3 * T_q_complete

# Indici M1, M2, M3 dei trimestri in-sample
M1_idx <- seq(1, n_in, by = 3)
M2_idx <- seq(2, n_in, by = 3)
M3_idx <- seq(3, n_in, by = 3)

stopifnot(length(M1_idx) == T_q_complete)

# Periodi temporali sui TRIMESTRI
start_eval  <- params$start_eval
end_eval    <- params$end_eval
covid_start <- params$covid_start
covid_end   <- params$covid_end

is_PRE   <- (dates_q_ts >= start_eval  & dates_q_ts <  covid_start)
is_COVID <- (dates_q_ts >= covid_start & dates_q_ts <= covid_end)
is_POST  <- (dates_q_ts >  covid_end   & dates_q_ts <= end_eval)
is_ALL   <- (dates_q_ts >= start_eval  & dates_q_ts <= end_eval)

period_labels <- c("Full sample", "Pre-COVID", "COVID period", "Post-COVID")
period_masks  <- list(is_ALL, is_PRE, is_COVID, is_POST)


# ==============================================================================
# 3. RMSFE M1/M2/M3 PER PAESE E PERIODO  [CORRETTO: y_now_in per paese]
# ==============================================================================

rmsfe_period <- function(mask_q, y_true_q, y_now_in, M_idx) {
  if (!any(mask_q)) return(NA_real_)
  y_true_sub <- y_true_q[mask_q]
  y_now_sub  <- y_now_in[M_idx[mask_q]]
  sqrt(mean((y_true_sub - y_now_sub)^2, na.rm = TRUE))
}

rmsfe_list <- list()

for (j in seq_along(countries_eval)) {
  cc    <- countries_eval[j]
  y_cc  <- y_true_q_all[, j]
  y_now_cc <- y_now_by_country[[cc]]
  y_now_in <- y_now_cc[1:n_in]
  
  for (p_id in seq_along(period_labels)) {
    lab   <- period_labels[p_id]
    mask  <- period_masks[[p_id]]
    
    RMS_M1 <- rmsfe_period(mask, y_cc, y_now_in, M1_idx)
    RMS_M2 <- rmsfe_period(mask, y_cc, y_now_in, M2_idx)
    RMS_M3 <- rmsfe_period(mask, y_cc, y_now_in, M3_idx)
    
    rmsfe_list[[length(rmsfe_list) + 1L]] <- data.frame(
      country = cc,
      period  = lab,
      M1      = RMS_M1,
      M2      = RMS_M2,
      M3      = RMS_M3
    )
  }
}

rmsfe_df <- bind_rows(rmsfe_list)

# Matrice wide base: righe = (country, period), colonne = M1,M2,M3
rmsfe_wide <- rmsfe_df |>
  arrange(country, factor(period, levels = period_labels))

# ------------------------------------------------------------------------------
# 3.1 Tabella LaTeX multi-country: PRE / COVID / POST (M1/M2/M3)
# ------------------------------------------------------------------------------

rmsfe_pc_wide <- rmsfe_df |>
  filter(period != "Full sample") |>
  mutate(period_short = dplyr::recode(
    period,
    "Pre-COVID"    = "PRE",
    "COVID period" = "COV",
    "Post-COVID"   = "POST"
  )) |>
  select(country, period_short, M1, M2, M3) |>
  pivot_wider(
    names_from  = period_short,
    values_from = c(M1, M2, M3),
    names_glue  = "{.value}_{period_short}"
  ) |>
  arrange(country)

latex_tab_pc <- paste0(
  "\\begin{table}[htbp]\n",
  "\\centering\n",
  "\\caption{Tensor MF-TPRF RMSFE by period (Pre, COVID, Post) -- all countries}\n",
  "\\label{tab:RMSFE_TMF3PRF_periods}\n",
  "\\scriptsize\n",
  "\\begin{tabular}{l",
  "ccc",  # PRE
  "ccc",  # COV
  "ccc",  # POST
  "}\n",
  "\\toprule\n",
  " & \\multicolumn{3}{c}{Pre-COVID} & \\multicolumn{3}{c}{COVID period} & \\multicolumn{3}{c}{Post-COVID}\\\\\n",
  "Country & M1 & M2 & M3 & M1 & M2 & M3 & M1 & M2 & M3 \\\\\n",
  "\\midrule\n",
  paste(
    sprintf("%s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
            rmsfe_pc_wide$country,
            rmsfe_pc_wide$M1_PRE,  rmsfe_pc_wide$M2_PRE,  rmsfe_pc_wide$M3_PRE,
            rmsfe_pc_wide$M1_COV,  rmsfe_pc_wide$M2_COV,  rmsfe_pc_wide$M3_COV,
            rmsfe_pc_wide$M1_POST, rmsfe_pc_wide$M2_POST, rmsfe_pc_wide$M3_POST),
    collapse = "\n"
  ),
  "\n\\bottomrule\n",
  "\\end{tabular}\n",
  "\\end{table}\n"
)

cat(latex_tab_pc, "\n\n")


# ==============================================================================
# 4. HEATMAP RMSFE (COUNTRY x PERIOD, FACET PER M1/M2/M3)
# ==============================================================================

rmsfe_long <- rmsfe_df |>
  mutate(
    period  = factor(period, levels = period_labels),
    country = factor(country, levels = countries_eval)
  ) |>
  pivot_longer(
    cols      = c(M1, M2, M3),
    names_to  = "month_pos",
    values_to = "RMSFE"
  )

# 4.1 Heatmap FULL (Full + Pre + COVID + Post)
heat_all <- ggplot(rmsfe_long,
                   aes(x = period, y = country, fill = RMSFE)) +
  geom_tile(color = "grey30") +
  facet_wrap(~ month_pos, ncol = 3) +
  scale_fill_gradient(
    name  = "RMSFE",
    low   = "white",
    high  = "darkred",
    na.value = "grey90"
  ) +
  labs(
    title = "Tensor MF-TPRF – RMSFE by country, period and nowcast month",
    x = "Period",
    y = "Country"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(heat_all)

# 4.2 Heatmap solo Pre / COVID / Post
rmsfe_long_pc <- rmsfe_long |>
  filter(period != "Full sample")

heat_pc <- ggplot(rmsfe_long_pc,
                  aes(x = period, y = country, fill = RMSFE)) +
  geom_tile(color = "grey30") +
  facet_wrap(~ month_pos, ncol = 3) +
  scale_fill_gradient(
    name  = "RMSFE",
    low   = "white",
    high  = "darkred",
    na.value = "grey90"
  ) +
  labs(
    title = "Tensor MF-TPRF – RMSFE by country and period (Pre / COVID / Post)",
    x = "Period",
    y = "Country"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

print(heat_pc)

# Salvataggio heatmap
file_heat_all <- file.path(
  path_graph_tensor,
  paste0(
    "RMSFE_heatmap_all_",
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-",   params$p_AR,
    "_r1-",    r[1],
    "_r2-",    r[2],
    "_",       format(params$start_est, "%Y-%m"),
    "_to_",    format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

file_heat_pc <- file.path(
  path_graph_tensor,
  paste0(
    "RMSFE_heatmap_PreCovidPost_",
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-",   params$p_AR,
    "_r1-",    r[1],
    "_r2-",    r[2],
    "_",       format(params$start_est, "%Y-%m"),
    "_to_",    format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(file_heat_all, plot = heat_all, width = 10, height = 4.5, dpi = 300)
ggsave(file_heat_pc,  plot = heat_pc,  width = 10, height = 4.5, dpi = 300)

cat("*** HEATMAPS SAVED TO ***\n",
    file_heat_all, "\n",
    file_heat_pc,  "\n\n")


# ==============================================================================
# 5. PLOT NOWCAST vs TRUE GDP – FACET PER PAESE  [CORRETTO: nowcast per paese]
# ==============================================================================


# True GDP trimestrale
df_true <- data.frame(
  date    = rep(dates_q_ts, times = length(countries_eval)),
  country = rep(countries_eval, each = T_q_complete),
  value   = as.vector(y_true_q_all[, seq_along(countries_eval)]),
  series  = "Quarterly GDP"
)

# Periodo in-sample vs real-time (mensile)
period_vec <- rep("in-sample", T_m_full)
if (T_m_full > n_in) {
  period_vec[(n_in + 1):T_m_full] <- "real-time"
}
period_vec <- factor(period_vec, levels = c("in-sample", "real-time"))

# Nowcast mensile per paese
df_now <- bind_rows(
  lapply(countries_eval, function(cc) {
    data.frame(
      date    = dates_m_full,
      country = cc,
      value   = y_now_by_country[[cc]],
      period  = period_vec,
      series  = "Nowcast"
    )
  })
)

df_true$country <- factor(df_true$country, levels = countries_eval)
df_now$country  <- factor(df_now$country,  levels = countries_eval)

plot_tensor_nowcast <- ggplot() +
  annotate(
    "rect",
    xmin = covid_start, xmax = covid_end,
    ymin = -Inf, ymax = Inf,
    fill = "grey80", alpha = 0.15
  ) +
  geom_line(
    data = df_now,
    aes(x = date, y = value, color = period),
    linewidth = 0.8
  ) +
  geom_line(
    data = df_true,
    aes(x = date, y = value, color = series),
    linewidth = 1.0
  ) +
  scale_color_manual(
    values = c(
      "in-sample"     = "#1F77B4",
      "real-time"     = "#2ECC71",
      "Quarterly GDP" = "#D62728"
    ),
    breaks = c("in-sample", "real-time", "Quarterly GDP"),
    name = ""
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title = "Tensor MF-TPRF – Monthly nowcast vs Quarterly GDP",
    x = "Date",
    y = "GDP growth"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

print(plot_tensor_nowcast)

file_graph_nowcast <- file.path(
  path_graph_tensor,
  paste0(
    "T_MF_TPRF_Nowcast_allCountries_",
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-",   params$p_AR,
    "_",       format(params$start_est, "%Y-%m"),
    "_to_",    format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(
  filename = file_graph_nowcast,
  plot     = plot_tensor_nowcast,
  width    = 12,
  height   = 8,
  dpi      = 300
)

cat("*** NOWCAST PLOT SAVED TO ***\n", file_graph_nowcast, "\n\n")


# ==============================================================================
# 6. EM CONVERGENCE PLOT (DIFFS) + SAVE
# ==============================================================================

df_em <- data.frame(
  iter = seq_along(EM_diffs),
  diff = EM_diffs
)

plot_em <- ggplot(df_em, aes(x = iter, y = diff)) +
  geom_line() +
  geom_point(size = 1.8) +
  scale_y_log10() +
  labs(
    title = "EM convergence – Tensor model",
    x = "Iteration",
    y = "Max change in X (log scale)"
  ) +
  theme_minimal(base_size = 13)

print(plot_em)

file_graph_em <- file.path(
  path_graph_tensor,
  paste0(
    "EM_Convergence_Tensor_",
    "_Lproxy-", Lproxy,
    "_Lmidas-", L_midas,
    "_pAR-",   params$p_AR,
    "_r1-",    r[1],
    "_r2-",    r[2],
    "_",       format(params$start_est, "%Y-%m"),
    "_to_",    format(params$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(
  filename = file_graph_em,
  plot     = plot_em,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("*** EM CONVERGENCE PLOT SAVED TO ***\n", file_graph_em, "\n")



















library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ==============================================================================
# PATHS
# ==============================================================================
path_graph   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph_tensor_RT")

if (!dir.exists(path_graph)) dir.create(path_graph, recursive = TRUE)

# ==============================================================================
# 1) LOAD ONLY PSEUDO REAL-TIME RESULTS  (COPIA-INCOLLA QUI IL NOME FILE)
# ==============================================================================
file_rt_name <- "T-MF_TPRF_pseudoRT_tensor_cc-DE-FR-IT-ES_sel-none_Lproxy-1_Lmidas-1_Nm-34_Nq-11_pAR-1_CvM-TRUE_CvQ-FALSE_2017-01_to_2025-10.rds"
file_rt_tensor <- file.path(path_results, file_rt_name)

if (!file.exists(file_rt_tensor)) stop("File PSEUDO-RT non trovato:\n", file_rt_tensor)
cat(">>> Tensor pseudo-RT file:\n", file_rt_tensor, "\n\n")

res_rt <- readRDS(file_rt_tensor)

# Meta
params_rt  <- res_rt$params
proxy_name <- res_rt$proxy_name %||% "EA"   # se non salvato, default "EA"
dates_m    <- as.Date(res_rt$dates_m)
dates_q    <- as.Date(res_rt$dates_q)

# Nowcast rolling (PER PAESE)
stopifnot(!is.null(res_rt$pseudo_rt_all))
stopifnot(all(c("date","country","nowcast","month_in_quarter") %in% names(res_rt$pseudo_rt_all)))

df_rt_all <- res_rt$pseudo_rt_all %>%
  mutate(
    date = as.Date(date),
    type = factor(month_in_quarter, levels = c("M1","M2","M3"))
  ) %>%
  select(date, country, nowcast, type) %>%
  arrange(country, date, type)

countries <- sort(unique(df_rt_all$country))
cat(">>> Countries in pseudoRT:\n", paste(countries, collapse = ", "), "\n\n")

# ------------------------------------------------------------
# FIX ORDER COUNTRIES (DE, FR, IT, ES)
# ------------------------------------------------------------
country_order <- c("DE","FR","IT","ES")

# tieni solo quelli presenti nel dataset (senza introdurre NA)
countries <- intersect(country_order, countries)

# ==============================================================================
# 2) TRUE QUARTERLY GDP PER COUNTRY (DA Y_q_all SALVATO NEL PSEUDO-RT)
# ==============================================================================
if (is.null(res_rt$Y_q_all)) stop("Nel file pseudo-RT manca Y_q_all (serve per TRUE GDP per paese).")

Y_q_all <- res_rt$Y_q_all
if (is.null(colnames(Y_q_all))) stop("Y_q_all deve avere colnames (proxy + paesi).")
if (!(proxy_name %in% colnames(Y_q_all))) stop("proxy_name=", proxy_name, " non è in colnames(Y_q_all).")

# Paesi = colonne di Y_q_all tranne la proxy (interseco con quelli nei nowcast)
countries_y <- setdiff(colnames(Y_q_all), proxy_name)
countries   <- intersect(countries, countries_y)
if (length(countries) == 0) stop("Nessun paese in comune tra pseudo_rt_all e Y_q_all.")

# Check allineamento righe Y_q_all con dates_q
if (nrow(Y_q_all) != length(dates_q)) {
  stop("Y_q_all e dates_q non allineati: nrow(Y_q_all) = ", nrow(Y_q_all),
       " ma length(dates_q) = ", length(dates_q))
}

# finestra evaluation coerente con params_rt (su trimestri)
idx_eval_q <- which(dates_q >= params_rt$start_eval & dates_q <= params_rt$end_eval)
if (length(idx_eval_q) == 0) stop("Nessun trimestre nella evaluation window (start_eval/end_eval).")

dates_q_eval <- dates_q[idx_eval_q]
Y_q_eval     <- Y_q_all[idx_eval_q, , drop = FALSE]

df_true_all <- data.frame(
  date_q     = rep(dates_q_eval, times = length(countries)),
  country    = rep(countries, each = length(dates_q_eval)),
  GDP        = as.vector(as.matrix(Y_q_eval[, countries, drop = FALSE])),
  quarter_id = rep(paste0(year(dates_q_eval), "Q", quarter(dates_q_eval)),
                   times = length(countries)),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 3) PREP ROLLING NOWCAST PER JOIN SU QUARTER (PER PAESE!)
# ==============================================================================
df_rt_q <- df_rt_all %>%
  filter(country %in% countries) %>%
  mutate(quarter_id = paste0(year(date), "Q", quarter(date))) %>%
  filter(date >= params_rt$start_eval & date <= params_rt$end_eval)

# ==============================================================================
# 4) JOIN (country, quarter_id) + PERIOD LABEL
# ==============================================================================
df_eval_all <- df_true_all %>%
  inner_join(df_rt_q %>% rename(date_m = date), by = c("country","quarter_id")) %>%
  mutate(
    period = dplyr::case_when(
      date_q <  params_rt$covid_start ~ "Pre-COVID",
      date_q >= params_rt$covid_start & date_q <= params_rt$covid_end ~ "COVID period",
      date_q >  params_rt$covid_end ~ "Post-COVID",
      TRUE ~ NA_character_
    ),
    period = factor(period, levels = c("Pre-COVID","COVID period","Post-COVID"))
  )

# ==============================================================================
# 5) PLOT: ROLLING NOWCASTS vs TRUE GDP (FACET PER PAESE)
# ==============================================================================
df_true_plot <- df_eval_all %>% select(country, date_q, GDP) %>% distinct()
df_rt_plot   <- df_eval_all %>% select(country, date_m, nowcast, type) %>% distinct() %>%
  rename(date = date_m)

df_true_plot$country <- factor(df_true_plot$country, levels = countries)
df_rt_plot$country   <- factor(df_rt_plot$country,   levels = countries)
df_eval_all$country <- factor(df_eval_all$country, levels = countries)

plot_tensor_rt <- ggplot() +
  annotate("rect",
           xmin = params_rt$covid_start, xmax = params_rt$covid_end,
           ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.20) +
  geom_line(data = df_true_plot,
            aes(x = date_q, y = GDP, color = "True GDP"),
            linewidth = 1.0) +
  geom_line(data = df_rt_plot,
            aes(x = date, y = nowcast, color = type),
            linewidth = 0.8) +
  geom_point(data = df_rt_plot,
             aes(x = date, y = nowcast, color = type),
             size = 1.8) +
  scale_color_manual(
    values = c("True GDP" = "#D62728", "M1"="#1F77B4", "M2"="#2ECC71", "M3"="#F1C40F"),
    name = "Series"
  ) +
  facet_wrap(~ country, ncol = 2, scales = "free_y") +
  labs(
    title    = "Tensor MF-TPRF – Rolling Pseudo Real-Time Nowcasts vs Quarterly GDP",
    subtitle = "M1: early • M2: mid-quarter • M3: end-quarter",
    x = "Date", y = "GDP growth"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

print(plot_tensor_rt)

file_graph_rt <- file.path(
  path_graph,
  paste0(
    "T_MF_TPRF_RollingPseudoRT_allCountries_",
    "_Lproxy-", res_rt$pseudo_rt_raw$Lproxy_fix,
    "_Lmidas-", res_rt$pseudo_rt_raw$L_midas_fix,
    "_",        format(params_rt$start_eval, "%Y-%m"),
    "_to_",     format(params_rt$end_eval, "%Y-%m"),
    ".png"
  )
)

ggsave(filename = file_graph_rt, plot = plot_tensor_rt, width = 12, height = 8, dpi = 300)
cat("\n*** ROLLING TENSOR NOWCAST GRAPH SAVED TO ***\n", file_graph_rt, "\n")

# ==============================================================================
# 6) RMSFE BY COUNTRY, PERIOD, TYPE (M1/M2/M3)
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
    country = factor(country, levels = countries)
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
  "\\label{tab:RMSFE_Tensor_MF3PRF_Rolling_All}\n",
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



















# OLD 

# ==============================================================================
# MF-3PRF: 
#===============================================================================
#' PROXY
#' LOADINGS 
#' FATTORI
#' TARGET PREDICTION
#' add the high friquency predictors
# ==============================================================================


# X_hf <- X_em_tens_centered   
# X_lf <- X_em_agg_tens_centered 
# y_q  <- y_EA_q_centered
# Robust_F   = TRUE
# alpha      = 0.05
# robust_type = "NW"
# nw_lag      = 1


Tensor_MF_TPRF <- function(X_lf, X_hf, y_q,
                           Lproxy      = 1,
                           L_midas     = 1,
                           p_AR        = 1,
                           r           = c(1, 1),   # c(r1, r2)
                           Robust_F    = FALSE, 
                           alpha       = 0.10,
                           robust_type = c("White", "NW"),
                           nw_lag      = 1) {
  #--------------------------------------------------
  # Preliminari e dimensioni
  #--------------------------------------------------
  y_q  <- as.numeric(y_q)          # target trimestrale
  T_q  <- length(y_q)
  T_m  <- dim(X_hf)[1]             # numero di mesi
  
  if (L_midas < 1) stop("L_midas deve essere >= 1")
  if (L_midas > T_q) stop("L_midas non può superare T_q")
  
  #--------------------------------------------------
  # STEP 1–2 (Tensor TPRF a bassa frequenza):
  # autoproxy + A,B,F_low (già calcolati dalla funzione precedente)
  #--------------------------------------------------
  auto_out <- build_autoproxy_Tensor_3prf(
    X_lf   = X_lf,     # T_q x p1 x p2
    y_q    = y_q,      # T_q
    Lproxy = Lproxy,
    r      = r
  )
  
  Z_q   <- auto_out$Z      # T_q x Lproxy (proxy: y_q, e1, ..., e_{Lproxy-1})
  A     <- auto_out$A      # p1 x r1
  B     <- auto_out$B      # p2 x r2
  # F_low <- auto_out$F_hat  # T_q x r1 x r2 (fattori trimestrali, opzionale diagnostica)
  
  r1 <- ncol(A)
  r2 <- ncol(B)
  K  <- r1 * r2
  
  #--------------------------------------------------
  # STEP 2bis: fattori MENSILI dal tensore X_hf usando A,B
  #           F_hf[t,,] = (A'A)^{-1}A' X_hf[t,,] B(B'B)^{-1}
  #--------------------------------------------------
  A_star <- solve(t(A) %*% A) %*% t(A)  # r1 x p1
  B_star <- solve(t(B) %*% B) %*% t(B)  # r2 x p2
  
  dates_m <- dimnames(X_hf)[[1]]
  if (is.null(dates_m)) {
    dates_m <- as.character(seq_len(T_m))
  }
  
  F_hf <- array(
    NA_real_,
    dim = c(T_m, r1, r2),
    dimnames = list(
      dates_m,
      paste0("f1_", seq_len(r1)),
      paste0("f2_", seq_len(r2))
    )
  )
  
  for (t in seq_len(T_m)) {
    Xt <- X_hf[t, , ]                    # p1 x p2
    F_hf[t, , ] <- A_star %*% Xt %*% t(B_star)   # r1 x r2
  }
  
  # Vectorizzazione: F_vec[t,] = vec(F_hf[t,,])
  F_vec <- matrix(NA_real_, nrow = T_m, ncol = K)
  for (t in seq_len(T_m)) {
    F_vec[t, ] <- as.vector(F_hf[t, , ])
  }
  
  #--------------------------------------------------
  # STEP 3.1 – Fattori trimestrali (come nel MF-TPRF vettoriale)
  #--------------------------------------------------
  if (T_m < 3 * T_q) {
    stop("Ci sono meno mesi di quelli necessari per coprire tutti i trimestri di y_q.")
  }
  
  # trimestri 1..T_q (completi) → usiamo i primi 3*T_q mesi
  F1 <- F_vec[seq(1, 3 * T_q, by = 3), , drop = FALSE]  # mese 1 del trimestre
  F2 <- F_vec[seq(2, 3 * T_q, by = 3), , drop = FALSE]  # mese 2
  F3 <- F_vec[seq(3, 3 * T_q, by = 3), , drop = FALSE]  # mese 3
  
  # mesi extra → trimestre T_q+1 (ragged edge)
  rem <- T_m - 3 * T_q               # 0,1,2 (o 3) mesi del trimestre successivo
  F_next1 <- if (rem >= 1) F_vec[3 * T_q + 1, , drop = FALSE] else NULL
  F_next2 <- if (rem >= 2) F_vec[3 * T_q + 2, , drop = FALSE] else NULL
  F_next3 <- if (rem >= 3) F_vec[3 * T_q + 3, , drop = FALSE] else NULL
  
  #--------------------------------------------------
  # STEP 3.2 – U-MIDAS trimestrale con AR(p_AR) in y
  #--------------------------------------------------
  if (p_AR < 0) stop("p_AR deve essere >= 0")
  
  start_tau <- max(L_midas, p_AR) + 1
  T_eff     <- T_q - max(L_midas, p_AR)
  
  # variabile dipendente: y_τ per τ = start_tau,...,T_q
  y_dep <- y_q[start_tau:T_q]   # lunghezza T_eff
  
  # ---- blocco AR in y ----
  Y_lag <- NULL
  if (p_AR > 0) {
    for (j in 1:p_AR) {
      Y_lag <- cbind(
        Y_lag,
        y_q[(start_tau - j):(T_q - j)]
      )
    }
    colnames(Y_lag) <- paste0("y_lag", 1:p_AR)
  }
  
  # ---- blocco fattori U-MIDAS (F1, F2, F3 con lag 0..L_midas-1) ----
  Xreg_F <- NULL
  for (ell_id in 1:L_midas) {
    ell   <- ell_id - 1                     # ell = 0,...,L_midas-1
    idx   <- (start_tau - ell):(T_q - ell)  # lunghezza T_eff
    
    Xreg_F <- cbind(
      Xreg_F,
      F1[idx, , drop = FALSE],
      F2[idx, , drop = FALSE],
      F3[idx, , drop = FALSE]
    )
  }
  
  # Combino AR + fattori
  if (p_AR > 0) {
    Xreg <- cbind(Y_lag, Xreg_F)
  } else {
    Xreg <- Xreg_F
  }
  
  # OLS con intercetta: y = β0 + Y_lag ρ + Xreg_F β + errore
  X_tilde_3 <- cbind(1, Xreg)
  XtX_3     <- t(X_tilde_3) %*% X_tilde_3
  XtY_3     <- t(X_tilde_3) %*% y_dep
  
  beta_hat  <- solve(XtX_3, XtY_3)
  
  beta0     <- as.numeric(beta_hat[1])
  
  # primi p_AR coefficienti dopo l'intercetta = AR in y
  if (p_AR > 0) {
    rho_hat  <- as.numeric(beta_hat[2:(1 + p_AR)])   # ρ_1,...,ρ_{p_AR}
    beta_vec <- as.numeric(beta_hat[(2 + p_AR):length(beta_hat)])
  } else {
    rho_hat  <- numeric(0)
    beta_vec <- as.numeric(beta_hat[-1])
  }
  
  # beta_mat: righe = ℓ_id = 1..L_midas (ell = 0..L_midas-1)
  # colonne per ogni lag: [β_{ℓ,1} (K), β_{ℓ,2} (K), β_{ℓ,3} (K)]
  beta_mat <- matrix(beta_vec,
                     nrow = L_midas,
                     ncol = 3 * K,
                     byrow = TRUE)
  
  #--------------------------------------------------------
  # STEP 4 – Nowcasting mensile con AR(p_AR) in y
  #--------------------------------------------------------
  y_nowcast <- rep(NA_real_, T_m)
  
  # punto di partenza coerente con la stima di STEP 3.2
  start_tau <- max(L_midas, p_AR) + 1
  
  # 4.a) Trimestri COMPLETI (backtest pseudo real time)
  if (start_tau <= T_q) {
    for (tau in start_tau:T_q) {
      
      # mesi del trimestre τ
      month_idx <- ((tau - 1) * 3 + 1):(tau * 3)
      
      # -------- parte AR in y per il trimestre τ --------
      if (p_AR > 0) {
        # y_{τ-1}, ..., y_{τ-p_AR}
        y_lags_tau <- sapply(1:p_AR, function(j) y_q[tau - j])
        AR_part_tau <- sum(rho_hat * y_lags_tau)
      } else {
        AR_part_tau <- 0
      }
      
      # -------- parte fattoriale mese per mese --------
      for (m in 1:3) {   # m = 1,2,3 (mese nel trimestre τ)
        
        contrib <- 0
        
        # somma sui lag ℓ_id = 1..L_midas
        for (ell_id in 1:L_midas) {
          ell   <- ell_id - 1
          lag_q <- tau - ell
          if (lag_q < 1) next
          
          # per ogni mese mm = 1,2,3
          for (mm in 1:3) {
            
            # Ragged edge "interno" al trimestre:
            # - per ell=0 (trimestre τ stesso), includo solo i mesi mm <= m
            # - per i lag (ell>=1), includo sempre mm=1,2,3
            if (ell_id == 1 && mm > m) next
            
            F_qm <- switch(
              mm,
              `1` = F1[lag_q, ],
              `2` = F2[lag_q, ],
              `3` = F3[lag_q, ]
            )
            
            start_col  <- (mm - 1) * K + 1
            end_col    <- mm * K
            beta_block <- beta_mat[ell_id, start_col:end_col]
            
            contrib <- contrib + sum(F_qm * beta_block)
          }
        }
        
        # nowcast per il mese m del trimestre τ:
        y_nowcast[month_idx[m]] <- beta0 + AR_part_tau + contrib
      }
    }
  }
  
  # 4.b) Trimestre CORRENTE T_q+1 (se ho mesi extra)
  if (rem > 0) {
    
    tau_curr   <- T_q + 1
    month_curr <- 3 * T_q + seq_len(rem)   # indici mesi disponibili del trimestre T_q+1
    
    # -------- parte AR in y per il trimestre corrente --------
    if (p_AR > 0) {
      # y_{T_q+1-1}, ..., y_{T_q+1-p_AR} = y_{T_q}, y_{T_q-1}, ...
      y_lags_curr <- sapply(1:p_AR, function(j) y_q[tau_curr - j])
      AR_part_curr <- sum(rho_hat * y_lags_curr)
    } else {
      AR_part_curr <- 0
    }
    
    # -------- parte fattoriale con F_next1/F_next2/F_next3 --------
    for (m in 1:rem) {   # m = 1,2 (o 3) mesi osservati nel trimestre corrente
      
      contrib <- 0
      
      for (ell_id in 1:L_midas) {
        ell   <- ell_id - 1
        lag_q <- tau_curr - ell
        
        # Per i lag (ell>=1) uso i trimestri COMPLETI in 1..T_q
        if (ell_id >= 2 && (lag_q < 1 || lag_q > T_q)) next
        
        for (mm in 1:3) {
          
          if (ell_id == 1) {
            # Trimestre corrente T_q+1: uso solo mesi <= m
            if (mm > m) next
            
            F_qm <- switch(
              mm,
              `1` = F_next1,
              `2` = F_next2,
              `3` = F_next3
            )
            if (is.null(F_qm)) next   # mese non ancora osservato
            
          } else {
            # Trimestri laggati: uso sempre i 3 mesi F1,F2,F3
            F_qm <- switch(
              mm,
              `1` = F1[lag_q, ],
              `2` = F2[lag_q, ],
              `3` = F3[lag_q, ]
            )
          }
          
          start_col  <- (mm - 1) * K + 1
          end_col    <- mm * K
          beta_block <- beta_mat[ell_id, start_col:end_col]
          
          contrib <- contrib + sum(F_qm * beta_block)
        }
      }
      
      # nowcast per il mese m del trimestre T_q+1
      y_nowcast[month_curr[m]] <- beta0 + AR_part_curr + contrib
    }
  }
  
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  return(list(
    Z_q       = Z_q,      # tutte le proxy
    A         = A,
    B         = B,
    F_hf      = F_hf,     # fattori mensili tensoriali
    F1        = F1,
    F2        = F2,
    F3        = F3,
    beta0     = beta0,
    rho_hat   = rho_hat,
    beta_mat  = beta_mat,
    y_nowcast = y_nowcast
  ))
}






pseudo_realtime_Tensor_MF_TPRF <- function(
    X_tens_full,    # [T_m x P x (N_m + N_q)] – predittori STANDARDIZZATI (tensore)
    y_q,            # vettore trimestrale target EA (lunghezza T_q)
    params,
    dates_m,        # date mensili (lunghezza T_m, coerenti con X_tens_full)
    dates_q,        # date trimestrali (lunghezza T_q, 3° mese del trimestre)
    Unb,            # matrice [P x (N_m + N_q)] – ritardi di pubblicazione (mesi)
    agg_m,          # matrice [P x N_m] – schema aggregazione mensili → trimestri
    agg_q,          # matrice [P x N_q] – schema aggregazione trimestrali
    Lproxy_fix,     # Lproxy scelto FUORI (estimation sample)
    L_midas_fix,    # L_midas scelto FUORI (choose_UMIDAS_lag)
    tensor_info     # lista con meta-info: N_m, N_q, countries, vars
) {
  # -------------------------------------------------------------------
  # 0. META: dimensioni, attributi, controlli
  # -------------------------------------------------------------------
  N_m     <- tensor$n_M                     # da oggetto globale 'tensor'
  nQ_tot  <- tensor$n_Q
  
  # serie trimestrali (tutte) e conteggio N_q escludendo il target
  Q_series_all <- tensor$vars[(N_m + 1):(N_m + nQ_tot)]
  q_cols       <- which(!grepl(params$target, Q_series_all, ignore.case = TRUE))
  N_q          <- length(q_cols)
  
  countries <- tensor_info$countries
  var_names <- tensor_info$vars
  
  T_m <- dim(X_tens_full)[1]
  P   <- dim(X_tens_full)[2]
  K   <- dim(X_tens_full)[3]
  
  stopifnot(N_m + N_q == K)
  stopifnot(all(dim(Unb) == c(P, K)))
  
  # Attacco meta-info come attributi (servono per flatten/unflatten)
  attr(X_tens_full, "N_m")       <- N_m
  attr(X_tens_full, "N_q")       <- N_q
  attr(X_tens_full, "countries") <- countries
  attr(X_tens_full, "var_names") <- var_names
  
  # Vettori di aggregazione globali (ordine coerente con flatten)
  agg_m_global <- as.vector(t(agg_m))  # lunghezza = P * N_m
  agg_q_global <- as.vector(t(agg_q))  # lunghezza = P * N_q
  
  meta <- list(
    N_m       = N_m,
    N_q       = N_q,
    countries = countries,
    var_names = var_names
  )
  
  # -------------------------------------------------------------------
  # 1. Evaluation window
  # -------------------------------------------------------------------
  t_start <- which(dates_m == params$start_eval)
  t_end   <- which(dates_m == params$end_eval)
  
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("start_eval o end_eval non trovate in 'dates_m'.")
  }
  
  # contenitori nowcast per M1, M2, M3 (chiavi = date)
  now_M1 <- list()
  now_M2 <- list()
  now_M3 <- list()
  
  # -------------------------------------------------------------------
  # 2. LOOP sui mesi in real-time
  # -------------------------------------------------------------------
  for (tt in seq(t_start, t_end)) {
    
    date_t <- dates_m[tt]
    cat("\n>>> REAL-TIME at", as.character(date_t),
        "| Lproxy =", Lproxy_fix,
        "| L_midas =", L_midas_fix, "\n")
    
    # ------------------------------------------------
    # 2.1 Unbalancedness tensoriale fino a tt
    # ------------------------------------------------
    X_cut_tens <- unbalancedness_tensor(
      X_full    = X_tens_full,
      Unb       = Unb,
      current_t = tt
    )  # [tt x P x K]
    
    # se non c'è nessuna informazione utile → skip
    if (all(is.na(X_cut_tens))) next
    
    # maschera di osservazione: 1 = osservato, 0 = NA
    W_cut_tens <- ifelse(is.na(X_cut_tens), 0, 1)
    
    # ------------------------------------------------
    # 2.2 Standardizzazione (tensore)
    # ------------------------------------------------
    std_out   <- standardize_mat_with_na(X_cut_tens)
    X_cut_std <- std_out$X_scaled   # [tt x P x K]
    
    # ------------------------------------------------
    # 2.3 Flatten tensor → matrice per EM
    # ------------------------------------------------
    flat_cut <- flatten_tensor_to_matrix(
      X_tens = X_cut_std,
      W_tens = W_cut_tens,
      N_m    = N_m,
      N_q    = N_q,
      agg_q  = agg_q                 # matrice [P x N_q]
    )
    
    X_obs_mat  <- flat_cut$X_mat        # [T_m_cur x N_global]
    N_q_global <- flat_cut$N_q_global
    agg_q_glob <- flat_cut$agg_q_global # vettore globale trimestrali
    
    # ------------------------------------------------
    # 2.4 A_list per EM (vettoriale)
    # ------------------------------------------------
    A_cut <- A_list(
      X_na  = X_obs_mat,
      N_q   = N_q_global,
      agg_q = agg_q_glob
    )
    
    # ------------------------------------------------
    # 2.5 Init EM (Chen/Lam) sul tensore standardizzato
    # ------------------------------------------------
    init <- init_CL_ER(
      X_std   = X_cut_std,
      W       = W_cut_tens,
      kmax    = params$kmax,
      do_plot = FALSE
    )
    
    X_init_tens <- init$X_init          # [tt x P x K]
    r_vec       <- prod(init$r)         # rank vettoriale = r1*r2
    
    # ------------------------------------------------
    # 2.6 Flatten X_init (stesso ordine colonne di X_obs_mat)
    # ------------------------------------------------
    X_init_mat <- flatten_tensor_init(
      X_tens = X_init_tens,
      N_m    = N_m,
      N_q    = N_q
    )
    
    # ------------------------------------------------
    # 2.7 EM (completion) sul sotto-campione vettoriale
    # ------------------------------------------------
    EM_out <- EM_algorithm(
      X_init   = X_init_mat,
      X_obs    = X_obs_mat,
      A_list   = A_cut,
      r        = r_vec,
      max_iter = 50,
      tol      = 1e-4
    )
    
    X_em    <- EM_out$X_completed       # [T_m_cur x N_global]
    T_m_cur <- nrow(X_em)
    
    # attributi per poter rifare unflatten
    for (nm in c("N_m", "N_q", "countries", "var_names")) {
      attr(X_em, nm) <- attr(X_obs_mat, nm)
    }
    
    # ------------------------------------------------
    # 2.8 separo mensili / trimestrali e aggrego
    # ------------------------------------------------
    X_m_em <- X_em[, 1:flat_cut$N_m_global, drop = FALSE]
    X_q_em <- X_em[, (flat_cut$N_m_global + 1):flat_cut$N_global, drop = FALSE]
    
    X_mq_em <- agg_mq(X_m_em, agg_m_global)
    X_qq_em <- agg_qq(X_q_em, agg_q_global)
    
    X_em_agg <- cbind(X_mq_em, X_qq_em)
    
    # attributi anche per X_em_agg
    for (nm in c("N_m", "N_q", "countries", "var_names")) {
      attr(X_em_agg, nm) <- attr(X_em, nm)
    }
    
    # ------------------------------------------------
    # 2.9 Ritorno al tensore trimestrale
    # ------------------------------------------------
    X_em_tens      <- unflatten_matrix_to_tensor(X_em)      # [T_m_cur x P x (N_m+N_q)]
    X_em_agg_tens  <- unflatten_matrix_to_tensor(X_em_agg)  # [T_q_em  x P x (N_m+N_q)]
    T_q_em         <- dim(X_em_agg_tens)[1]
    
    # ------------------------------------------------
    # 2.10 Trimestri di PIL disponibili a date_t
    # ------------------------------------------------
    idx_pub <- which(dates_q < date_t)
    if (length(idx_pub) < 2) next  # troppo pochi trimestri osservati
    
    T_q_current <- tail(idx_pub, 1)
    T_q_use     <- min(T_q_current, T_q_em, length(y_q))
    
    if (T_q_use < 2) next
    
    y_q_cut <- y_q[1:T_q_use]
    
    X_lf_cut_tens <- X_em_agg_tens[1:T_q_use,  , , drop = FALSE]
    X_hf_cut_tens <- X_em_tens[1:T_m_cur,     , , drop = FALSE]
    
    # ------------------------------------------------
    # 2.11 Centering (per togliere intercetta in Tensor_MF_TPRF)
    # ------------------------------------------------
    lf_center     <- center_Y(X_lf_cut_tens)
    X_lf_centered <- lf_center$Y_centered
    
    hf_center     <- center_Y(X_hf_cut_tens)
    X_hf_centered <- hf_center$Y_centered
    
    y_q_centered  <- scale(y_q_cut, center = TRUE, scale = FALSE)
    
    # ------------------------------------------------
    # 2.12 TENSOR MF-TPRF (con Lproxy/L_midas fissi)
    # ------------------------------------------------
    T_MF_TPRF_RT <- Tensor_MF_TPRF(
      X_lf        = X_lf_centered,
      X_hf        = X_hf_centered,
      y_q         = as.numeric(y_q_centered),
      Lproxy      = Lproxy_fix,
      L_midas     = L_midas_fix,
      p_AR        = params$p_AR,
      r           = init$r,          # rank tensoriale corrente
      Robust_F    = params$Robust_F,
      alpha       = params$alpha,
      robust_type = params$robust_type,
      nw_lag      = params$nw_lag
    )
    
    y_rt_full <- T_MF_TPRF_RT$y_nowcast
    if (length(y_rt_full) == 0) next
    y_rt_last <- tail(y_rt_full, 1)
    
    # ------------------------------------------------
    # 2.13 Assegno il nowcast al mese M1/M2/M3 del trimestre target
    # ------------------------------------------------
    m_tr <- compute_m_tr(date_t, dates_q)  # 1, 2, 3 oppure NA
    if (is.na(m_tr)) next
    
    key <- as.character(date_t)
    if (m_tr == 1) now_M1[[key]] <- y_rt_last
    if (m_tr == 2) now_M2[[key]] <- y_rt_last
    if (m_tr == 3) now_M3[[key]] <- y_rt_last
  }
  
  # -------------------------------------------------------------------
  # 3. Output: Lproxy/L_midas usati + vettori M1/M2/M3 (nomi = date)
  # -------------------------------------------------------------------
  M1_vec <- unlist(now_M1)
  M2_vec <- unlist(now_M2)
  M3_vec <- unlist(now_M3)
  
  list(
    Lproxy_fix  = Lproxy_fix,
    L_midas_fix = L_midas_fix,
    M1          = M1_vec,
    M2          = M2_vec,
    M3          = M3_vec
  )
}


