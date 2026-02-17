# ==============================================================================
#             Dynamic Matrix Factor Models and the EM Algorithm: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#                               - Main Script - 
# ==============================================================================
# This script estimates a Dynamic Matrix Factor Model (DMFM) by the EM-algorithm 
# and performs a real-time nowcasting exercise for GDP growth in the Euro Area's 
# main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation and harmonization across countries
#   2. Construction of a 3D data tensor
#   3. Estimation of the number of latent factors and model lags
#   4. Model initialization and parameter estimation using the EM algorithm on 
#      the whole dataset
#   5. Pseudo real-time rolling nowcasting of GDP
#   6. Saving results for further use
# DMFM and DFM sections are modular and clearly separated for experimentation.
#
# ==============================================================================
# Author      : Davide Delfino
# Institution : Alma Mater Studiorum - University of Bologna
# Dataset     : Barigozzi & Lissona (2024), EA-MD-QD
#
# Disclaimer  : DMFM functions adapted from Barigozzi & Trapin (2024)
#
# ==============================================================================
# Notes       : Other EA countries such as BE, NL, AT, PT, IR, GR and different 
#               set of Monthly and Quarterly variables can be added
# ==============================================================================
# Script Type : Data Preparation / Estimation / Forecasting
# ==============================================================================

# ==============================================================================
# SET WORKING DIRECTORY
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "data")
path_func    <- file.path(path_main, "Matrix_FM/functions")
path_results <- file.path(path_main, "Matrix_FM/DMFM/results/outputs")
path_graph   <- file.path(path_main, "Matrix_FM/DMFM/results/graph")



# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

# Data Handling
library(tidyverse)
library(lubridate)
library(abind)

# Time Series Analysis
library(tseries)
library(zoo)
library(fBasics)
library(vars)

# Linear Algebra & Utilities
library(Matrix)
library(RSpectra)
library(MASS)
library(pracma)

# I/O
library(writexl)
library(readxl)
library(xtable)

# Plotting
library(ggplot2)
library(reshape2)
library(patchwork)


# ==============================================================================
# SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

source(file.path(path_func, "dmfm.prep.R"))
source(file.path(path_func, "dmfm.utils.R"))
source(file.path(path_func, "dmfm.initialization.R"))
source(file.path(path_func, "dmfm.estimation.R"))
source(file.path(path_func, "dmfm.kalman.R"))
source(file.path(path_func, "dmfm.results.R"))


# ==============================================================================
# STEP 1: DATA PREPARATION AND ALIGNMENT
# ==============================================================================

# Users's parameters
params <- list(
  start_est   = as.Date("2000-02-01"),
  start_eval  = as.Date("2017-01-01"),
  end_eval    = as.Date("2025-09-01"),
  covid_start = as.Date("2020-03-01"),
  covid_end   = as.Date("2021-07-01"),
  covid_mask  = TRUE,    # boolean, not string
  target      = "GDP",   # quarterly target variable
  sel_method  = "corr_threshold", # "none" | "corr_threshold" | "F-Test"
  
  # Variable selection parameters
  n_m        = 38,
  n_q        = 2,
  thr_m      = 0.10,
  thr_q      = 0.85,
  thr_F_test = 0.01
)


# Quarterly Treatment
R_mat <- matrix(c(
  2, -1, 0,  0,  0,
  3,  0, -1, 0,  0,
  2,  0,  0, -1, 0,
  1,  0,  0,  0, -1
), nrow = 4, byrow = TRUE)

q <- rep(0, 4)


kmax <- c(3, 7)  # Initial upper bounds for number of factors

# Country List
countries <- c("DE", "FR", "IT", "ES")

country <- "ES"
# ==============================================================================
# 4. PREPARE DATA FOR ALL COUNTRIES AND EXCTARCT THE ONE OF INTEREST
# ==============================================================================

all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
  covid_mask      = params$covid_mask
)


# ==============================================================================
# STEP 2: BUILD TENSOR DATASET (T × p1 × p2)
# ==============================================================================
# 2) Costruzione del tensore DMFM
tensor <- build_tensor_dmfm(all_countries, countries = countries, span = "union")

str(tensor$Y)  # [1] T_global  P1  P2
str(tensor$W)
tensor$vars      # nomi delle variabili comuni (base names)
tensor$dates[1]  # prima data
tail(tensor$dates, 1)  # ultima data



Y <- tensor$Y  # Data tensor
W <- tensor$W  # Mask (missing indicator)

# Standardize and demean
std <- standardize_Y(Y)
Y_std <- std$Y_scaled

# Check percentage of missing values
nan_percent_Y(Y_std)
gdp_idx <- which(dimnames(tensor$Y)[[3]] == params$target)



################################# !is.na GDP ###################################


ea_index <- which(dimnames(tensor$Y)[[2]] == country)
gdp_idx <- which(dimnames(tensor$Y)[[3]] == params$target)

# Estrai la serie del GDP (standardizzata)
gdp_ts <- tensor$Y[, ea_index, gdp_idx]


# Serie temporale con Date
gdp_df <- data.frame(
  Date = tensor$dates,              # oppure converti in trimestre con get_quarter()
  GDP  = gdp_ts
)

ggplot(gdp_df, aes(x = Date, y = GDP)) +
  geom_line(color = "#1f77b4", size = 1) +
  geom_point(data = subset(gdp_df, !is.na(GDP)), color = "black", size = 1.5) +
  labs(
    title = "Observed GDP Series from Tensor",
    subtitle = "Standardized GDP values for Italy (example)",
    x = "Date", y = "GDP (Standardized)"
  ) +
  theme_minimal(base_size = 13)


# ==============================================================================
# STEP 3: FACTOR ESTIMATION
# ==============================================================================

## === DMFM Pipeline === ##
k_hat <- mfm.f(Y_std, W, kmax)  # Estimate number of factors
f.lag_mar <- mar_model_selection_auto(Y_std, W, k_hat)  # Select lag order

## === DFM Alternative === ##
# k_hat <- dfm.f.vec(Y_std, W, kmax)
# f.lag_mar <- mar_model_selection_vec(Y_std, W, k_hat)


# ==============================================================================
# STEP 4: INITIALIZE MODEL INPUTS
# ==============================================================================

# k_hat <- c(1,3)

## === DMFM Initialization === ##
imp <- mfm.cl(Y_std, W, k_hat)
inputs <- dmfm.na.2.sv(imp$Y_imputed, k_hat, W, t = "dmfm")

## === DFM Alternative Initialization === ##
# imp <- mfm.cl.vec(Y_std, W, k_hat)
# inputs <- dmfm.na.2.sv.vec(imp$Y_imputed, k_hat, W, t = "dmfm")


# ==============================================================================
# STEP 5: FIT DMFM USING FULL DATASET
# ==============================================================================

out <- dmfm.na.em(
  . = inputs,
  Y = Y_std,
  k = k_hat,
  W = W,
  t = "dmfm"
)

# ============================================================================
# 12. SAVE DMFM RESULTS (FULL SAMPLE ESTIMATION)
# ============================================================================

# Cartella: /results/DMFM/IT/ (adatta se vuoi una struttura diversa)
path_results_dmfm_country <- file.path(path_results, "full_sample")
if (!dir.exists(path_results_dmfm_country)) {
  dir.create(path_results_dmfm_country, recursive = TRUE)
}

# Stringhe di servizio per il nome file
sel_method <- if (!is.null(params$sel_method)) params$sel_method else "none"
covid_mask_str <- if (isTRUE(params$covid_mask)) "TRUE" else "FALSE"

start_est_str <- format(params$start_est, "%Y-%m")
end_eval_str  <- format(params$end_eval, "%Y-%m")

# Oggetto da salvare: tutto quello che serve per rifare analisi/forecast
dmfm_full <- list(
  params  = params,
  
  # Tensor data (grezzi + mask + metadati)
  tensor  = list(
    Y     = Y,
    W     = W,
    dates = tensor$dates,
    vars  = tensor$vars
  ),
  
  # Standardizzazione (serve per tornare allo scale originale)
  standardization = std,      # quello restituito da standardize_Y()
  Y_std          = Y_std,
  
  # Selezione dimensione fattori e ordine MAR
  k_hat   = k_hat,
  f_lag   = f.lag_mar,
  
  # Imputazione iniziale e parametri di partenza EM
  imp      = imp,
  init_in  = inputs,
  
  # Risultati EM completi
  em_out   = out
)

# Nome file
file_dmfm_full <- file.path(
  path_results_dmfm_country,
  paste0(
    "DMFM_FullSample_",
    country,
    "_sel-",   sel_method,
    "_Covid-", covid_mask_str,
    "_",       start_est_str,
    "_to_",    end_eval_str,
    ".rds"
  )
)

saveRDS(dmfm_full, file_dmfm_full)

cat("\n*** SAVED DMFM FULL-SAMPLE RESULTS TO ***\n", file_dmfm_full, "\n")











# ==============================================================================
# STEP 6: ROLLING NOWCASTING (PSEUDO REAL-TIME)
# ==============================================================================


Type_m  <- all_countries$data[[country]]$TypeM
Type_q  <- all_countries$data[[country]]$TypeQ
Type    <- c(Type_m, Type_q)

Freq_m  <- all_countries$data[[country]]$freq_m
Freq_q  <- all_countries$data[[country]]$freq_q
Freq    <- c(Freq_m, Freq_q)

Unb_m   <- all_countries$data[[country]]$unb_m
Unb_q   <- all_countries$data[[country]]$unb_q
Unb     <- c(Unb_m, Unb_q)

Class_m <- all_countries$data[[country]]$ClassM
Class_q <- all_countries$data[[country]]$ClassQ
Class   <- c(Class_m, Class_q)

roll_nowcast <- rolling_nowcast_dmfm_by_release(
  Y_full = Y,
  W_full = W,
  k = k_hat,
  dates = tensor$dates,
  group = Unb,
  gdp_idx = gdp_idx,
  startEV = params$start_eval,
  endEV = params$end_eval,
  model_type = "dmfm",
  max_iter = 1000,
  eps = 1e-3
)


# ==============================================================================
# STEP 7: SAVE OUTPUT
# ==============================================================================

save(
  res, tensor, k_hat, W, inputs, std, out, roll_nowcast, 
  file = "results/DMFM/new_dmfm_rollnowcastLS_11_204_40MedvarCovidpluse-03.RData"
)

# (Optional) Load previously saved results
# load("Results/dmfm_rollnowcastLS_11_204_40var.RData")


# ==============================================================================
# 14. SAVE RESULTS – ONLY ROLLING NOWCAST (DFM)
# ==============================================================================

# ==============================================================================
# 14. SAVE RESULTS – ROLLING NOWCAST (DMFM)
# ==============================================================================

# Usiamo la stessa cartella del full sample: /results/DMFM/IT/
path_results_dmfm_country <- file.path(path_results, "rolling_nowcast")
if (!dir.exists(path_results_dmfm_country)) {
  dir.create(path_results_dmfm_country, recursive = TRUE)
}

# Stringhe di servizio (se già definite sopra puoi omettere questa parte)
sel_method    <- if (!is.null(params$sel_method)) params$sel_method else "none"
covid_mask_str <- if (isTRUE(params$covid_mask)) "TRUE" else "FALSE"

start_est_str <- format(params$start_est, "%Y-%m")
end_eval_str  <- format(params$end_eval, "%Y-%m")

# Contenuto da salvare: solo informazioni utili per analisi dei nowcast
dmfm_rolling <- list(
  country = country,
  params  = params,
  
  # Info sulle serie usate
  Freq = Freq,   # vettore "M"/"Q"
  Unb  = Unb,    # gruppi di unbalancedness
  
  # Info temporali
  dates_tensor = tensor$dates,
  gdp_idx      = gdp_idx,
  
  # Risultato della funzione rolling_nowcast_dmfm_by_release
  nowcast = roll_nowcast
)

# Nome file
file_dmfm_rolling <- file.path(
  path_results_dmfm_country,
  paste0(
    "DMFM_RollingNowcast_",
    country,
    "_sel-",  sel_method,
    "_Covid-", covid_mask_str,
    "_Nm-",   sum(Freq == "M"),
    "_Nq-",   sum(Freq == "Q"),
    "_",      start_est_str,
    "_to_",   end_eval_str,
    ".rds"
  )
)

saveRDS(dmfm_rolling, file_dmfm_rolling)

cat("\n*** SAVED DMFM ROLLING NOWCAST TO ***\n", file_dmfm_rolling, "\n")



# ==============================================================================
# END OF SCRIPT
# ==============================================================================




# ======================================================================
#  PLOTTING & EVALUATION: DMFM ROLLING NOWCAST (M1, M2, M3)
# ======================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ----------------------------------------------------------------------
# 0. Funzione per destandardizzare SOLO il GDP (scalare!)
# ----------------------------------------------------------------------
destandardize_gdp_scalar <- function(x_std, mean_Y, sd_Y, gdp_idx) {
  x_std * sd_Y[gdp_idx] + mean_Y[gdp_idx]
}

# ----------------------------------------------------------------------
# 1. Helper trimestri
# ----------------------------------------------------------------------
get_quarter <- function(date) {
  date <- as.Date(date)
  paste0(format(date, "%Y"), "Q", ceiling(as.numeric(format(date, "%m")) / 3))
}

quarter_to_date <- function(qstr) {
  year <- as.integer(substr(qstr, 1, 4))
  q    <- as.integer(substr(qstr, 6, 6))
  month_end <- q * 3
  as.Date(sprintf("%d-%02d-01", year, month_end))
}

# ----------------------------------------------------------------------
# 2. Impostazioni
# ----------------------------------------------------------------------
Y        <- tensor$Y
dates_m  <- tensor$dates
month_labels <- c("M1","M2","M3")
n_months <- 3

country_idx <- which(dimnames(Y)[[2]] == country)
if (length(country_idx) != 1) stop("Problema con country_idx.")

# ----------------------------------------------------------------------
# 3. I nomi nei nowcast diventano "YYYYQk"
# ----------------------------------------------------------------------
for (m in month_labels)
  names(roll_nowcast[[m]]) <- get_quarter(names(roll_nowcast[[m]]))

common_quarters <- sort(Reduce(intersect,
                               lapply(roll_nowcast[month_labels], names)))

if (length(common_quarters)==0) stop("No common quarters M1/M2/M3")

n_qtr <- length(common_quarters)

# ----------------------------------------------------------------------
# 4. Costruzione matrice nowcast standardizzati
# ----------------------------------------------------------------------
gdp_mat_std <- matrix(NA, n_qtr, 3,
                      dimnames=list(common_quarters, month_labels))

for (j in 1:3) {
  mj <- month_labels[j]
  for (i in 1:n_qtr) {
    q <- common_quarters[i]
    pred_vec <- roll_nowcast[[mj]][[q]]
    gdp_mat_std[i,j] <- pred_vec[country_idx]
  }
}

# ----------------------------------------------------------------------
# 5. True GDP standardizzato estratto dal tensore Y_std
# ----------------------------------------------------------------------
quarters_all <- get_quarter(dates_m)
gdp_true_std <- numeric(n_qtr)

for (i in 1:n_qtr) {
  q <- common_quarters[i]
  idx_q <- which(quarters_all == q)
  
  gdp_months_std <- Y[idx_q, country_idx, gdp_idx]
  non_na_idx <- which(!is.na(gdp_months_std))
  
  gdp_true_std[i] <- gdp_months_std[max(non_na_idx)]
}

names(gdp_true_std) <- common_quarters

# ----------------------------------------------------------------------
# 6.❗ DESTANDARDIZZAZIONE (SCALARE) CORRETTA
# ----------------------------------------------------------------------
gdp_true <- destandardize_gdp_scalar(
  gdp_true_std,
  mean_Y = std$mean,
  sd_Y   = std$sd,
  gdp_idx = gdp_idx
)

gdp_mat <- matrix(NA, n_qtr, 3,
                  dimnames=list(common_quarters, month_labels))

for (m in 1:3) {
  gdp_mat[,m] <- destandardize_gdp_scalar(
    gdp_mat_std[,m],
    mean_Y = std$mean,
    sd_Y   = std$sd,
    gdp_idx = gdp_idx
  )
}

# ----------------------------------------------------------------------
# 7. Data frame per plot
# ----------------------------------------------------------------------
quarter_dates <- quarter_to_date(common_quarters)

df_ts <- data.frame(
  Quarter     = factor(common_quarters, levels=common_quarters),
  QuarterDate = quarter_dates,
  TrueGDP     = gdp_true,
  M1 = gdp_mat[,1],
  M2 = gdp_mat[,2],
  M3 = gdp_mat[,3]
)

df_ts_long <- df_ts |>
  pivot_longer(cols=c("M1","M2","M3"),
               names_to="Model", values_to="Nowcast")

# ----------------------------------------------------------------------
# 8. PLOT – tutti insieme True vs M1 M2 M3
# ----------------------------------------------------------------------
g_all <- ggplot() +
  # TRUE GDP
  geom_line(
    data = df_ts,
    aes(x = Quarter, y = TrueGDP, color = "True GDP", group = 1),
    linewidth = 1.3
  ) +
  
  # NOWCAST M1/M2/M3
  geom_line(
    data = df_ts_long,
    aes(x = Quarter, y = Nowcast, color = Model, group = Model),
    linewidth = 1, linetype = "dashed"
  ) +
  
  scale_color_manual(
    values = c(
      "True GDP" = "black",
      "M1" = "#1f77b4",
      "M2" = "#ff7f0e",
      "M3" = "#2ca02c"
    ),
    breaks = c("True GDP", "M1", "M2", "M3")
  ) +
  labs(
    title = paste("DMFM – True GDP vs Nowcast M1/M2/M3 (", country, ")", sep = ""),
    x = "Quarter", y = "GDP (destandardized)",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g_all

# ----------------------------------------------------------------------
# 9. PLOT – facet M1/M2/M3 separati
# ----------------------------------------------------------------------
g_facet <- ggplot(df_ts_long, aes(x=Quarter, y=Nowcast)) +
  geom_line(color="#1f77b4", linewidth=1) +
  geom_line(aes(y=TrueGDP),
            data=df_ts, color="black", linewidth=1) +
  facet_wrap(~Model, ncol=1) +
  theme_minimal(base_size=14) +
  labs(
    title=paste("DMFM Nowcast per mese –", country),
    x="Quarter", y="GDP"
  ) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

g_facet

# ----------------------------------------------------------------------
# 10. Errori + RMSFE (destandardizzati!)
# ----------------------------------------------------------------------
df_errors <- df_ts |>
  select(Quarter, QuarterDate, TrueGDP, M1,M2,M3) |>
  pivot_longer(cols=c("M1","M2","M3"),
               names_to="Model", values_to="Nowcast") |>
  mutate(Error = Nowcast - TrueGDP)

compute_rmsfe <- function(x) sqrt(mean(x^2, na.rm=TRUE))

rmsfe_all <- df_errors |>
  group_by(Model) |>
  summarise(RMSFE = compute_rmsfe(Error))

print(rmsfe_all)

# BOXPLOT Errori
g_errors <- ggplot(df_errors, aes(x=Model, y=Error, fill=Model)) +
  geom_boxplot(alpha=0.8) +
  theme_minimal(base_size=14) +
  labs(title="Forecast Errors (DMFM)", y="Error", x="Model")

g_errors




#########################################



load("results/DMFM/dmfm_rollnowcastLS_11_204_40varCovidpluse-03.RData")


mod <- out$model


# 2. Ricostruisci idiosincratici
E_hat <- Y_std - mod$Y

e_DE_GDP <- na.omit(E_hat[, 3, 40])
acf(e_DE_GDP)
pacf(e_DE_GDP)

e_FR_GDP <- na.omit(E_hat[, 4, 40])

acf(e_FR_GDP)
pacf(e_FR_GDP)

# 5. Eventuale test di autocorrelazione
Box.test(e_DE_GDP, lag = 12, type = "Ljung-Box")
Box.test(e_FR_GDP, lag = 12, type = "Ljung-Box")
