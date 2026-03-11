# ==============================================================================
#                        Three Pass Regression Filter: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#
# ==============================================================================
# This script estimates a Three Pass regression filter (TPRF) for Mixed Frequency
# data and performs a pseudo real-time nowcasting exercise for GDP growth in the 
# Euro Area's main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation and harmonization across countries
#   2. NA's Imputation in predictors panel         --> Xiong & Pelger (2020)
#   3. Estimation of the number of latent factors  --> Ahn (2013)
#   4. Estimation of the number of Proxies         --> Kraemer & Sugiyama (2011)
#   5, Estimation of U-MIDAS and AR lags number    --> IC
#   6. Full sample estimation with MF-TPRF
#   7. Pseudo real-time rolling nowcast of GDP
#
# ==============================================================================
# Author          : Davide Delfino
# Institution     : Bocconi University
#
# ===============================================================================
# Main Reference  : Hepenstick & Marcellino (2018), MF-TPRF
# Dataset         : Barigozzi & Lissona (2024), EA-MD-QD
#
# ==============================================================================
# Notes : Other EA countries such as BE, NL, AT, PT, IR, GR and different set of
#         Monthly and Quarterly variables are available from the Dataset
#
# ==============================================================================
# Script Type     : Main script
# ==============================================================================

# ==============================================================================
# OUTPUT     
# - Fit     : better when COVID is not masked 
# - Nowcast : better when COVID is masked
# ==============================================================================





# ==============================================================================
# SET WORKING DIRECTORY & PATHS
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "data")
path_func    <- file.path(path_main, "Vector_FM/functions")
path_results <- file.path(path_main, "Vector_FM/MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "Vector_FM/MF-TPRF/results/graph")

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

## Data Handling
library(tidyverse)
library(lubridate)
library(abind)
library(dplyr)

## Time Series
library(tseries)
library(zoo)
library(MASS)
library(glmnet)
library(plsdof)

## test
library(sandwich)
library(lmtest)
library(car)

## I/O
library(readxl)

## Plotting
library(ggplot2)

## Resolve conflicts
library(conflicted)
conflict_prefer("select",  "dplyr", quiet = TRUE)
conflict_prefer("filter",  "dplyr", quiet = TRUE)

# ==============================================================================
# SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

# data preparation
source(file.path(path_func, "mf.tprf.utils.R"))
source(file.path(path_func, "mf.tprf.prep.R"))

# mf-tprf estimation
source(file.path(path_func, "mf.tprf.R"))

# imputation
source(file.path(path_func, "mf.tprf.imp.R"))

# factor selection
source(file.path(path_func, "mf.tprf.fs.R"))

# nowcasting
source(file.path(path_func, "mf.tprf.now.R"))


# ==============================================================================
# PARAMETERS & COUNTRY LIST & COUNTRY OF INTEREST SELECTION
# ==============================================================================

# Users's parameters
params <- list(
  start_est    = as.Date("2000-04-01"),   # growth rate 2001-01 NA
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-10-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,                    # Boolean: if empty == TRUE
  covid_mask_q = TRUE,                    # Boolean: if empty == TRUE
  target       = "GDP",                   # quarterly target variable
  
  # Variable selection parameters
  sel_method   = "none",                  # "none" | "corr_threshold" | "F-Test" | "LASSO"
  n_m          = 20,
  n_q          = 50,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,                        # 1=LASSO, 0=Ridge, (0,1)=Elastic Net
  
  # Factor and proxy selection parameters
  Kmax = 10,                               # max number of factors
  Zmax = 10,                               # max number of proxy
  
  # MF-TPRF parameters
  p_AR_max    = 4,                         # Max numbers of lags y in U-MIDAS
  Lmax        = 6,                         # Max numbers of lags in the Factors in U-MIDAS
  Robust_F    = FALSE,                      # Robust F test in the first step
  alpha       = 0.1,
  robust_type = "NW",                      # White--> HC | NW(Newey West)--> HAC
  nw_lag      = 1                          # Lag in the error autocorrelation in step 1
)

# Country List
countries <- c("DE", "FR", "IT", "ES")

# Selected country
country <- "FR"

# ==============================================================================
# 1 PREPARE DATA FOR ALL COUNTRIES AND EXCTARCT THE ONE OF INTEREST
# ==============================================================================
all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
  covid_mask_m    = params$covid_mask_m,
  covid_mask_q    = params$covid_mask_q
)


# Extract selected country's data
data     <- all_countries$data[[country]]$Data
dates_m  <- all_countries$data[[country]]$Dates
dates_q  <- all_countries$data[[country]]$DatesQ
series   <- all_countries$data[[country]]$Series


# ====================================
# 1.2 TARGET AND PREDICTORS EXTRACTION
# ====================================

# target vector extraction (GDP for "cc")
gdp_col     <- all_countries$data[[country]]$target_col
target_name <- paste0(params$target, "_", country)  # es: "GDP_FR"

# The quarterly growth in the first quarter is NA (NO EA before the 2000)
y   <- as.matrix(data[, gdp_col, drop = FALSE])       # monthly target
y_q <- as.matrix(y[!is.na(y)])                        # quarterly target

# predictors panel extraction (X for "cc")
X   <- as.matrix(data[, -gdp_col, drop = FALSE])      

# =======================
# 1.3 METADATA EXTRACTION
# =======================

# Panel dimensions
N_m     <- all_countries$data[[country]]$nM   # monthly indicators

N_q_tot      <- all_countries$data[[country]]$nQ 
q_series_all <- all_countries$data[[country]]$Series[(N_m + 1):(N_m + N_q_tot)]
q_cols       <- which(toupper(q_series_all) != toupper(target_name))
N_q          <- length(q_cols)                # quarterly indicators

N       <- N_m + N_q                          # Tot number of indicators


# Aggregation rule (stock & flows)
agg_m   <- all_countries$data[[country]]$agg_m
agg_q   <- all_countries$data[[country]]$agg_q[q_cols]
agg     <- c(agg_m, agg_q)

# Frequency (monthly or quarterly)
Freq_m  <- all_countries$data[[country]]$freq_m
Freq_q  <- all_countries$data[[country]]$freq_q[q_cols]
Freq    <- c(Freq_m, Freq_q)

# Indicators' releasing delay
Unb_m   <- all_countries$data[[country]]$unb_m
Unb_q   <- all_countries$data[[country]]$unb_q[q_cols]
Unb     <- c(Unb_m, Unb_q)

# Indicators' class (Real, Nominal, Financial, Survey)
Class_m <- all_countries$data[[country]]$ClassM
Class_q <- all_countries$data[[country]]$ClassQ[q_cols]
Class   <- c(Class_m, Class_q)

# ==============================================================================
# 2. NA's IMPUTATION          : XIONG & PELGER (2020)
# 3. LATENT FACTOR EXTRACTION : Ahn (2011)
# ==============================================================================

# Standardization
out_std <- standardize_with_na(X)
X_std   <- out_std$X_std

# Imputation 
imp_xp  <- init_XP_ER(X_std, params$Kmax)

X_xp    <- imp_xp$X_init     # Imputed Data
r_hat   <- imp_xp$r          # Number of latent factors in the imputed data

# ================================
# 3.1 AGGREGATION TO LOW FREQUENCY
# ================================

X_m_xp   <- X_xp[,1:N_m]
X_q_xp   <- X_xp[,(N_m+1):N]

X_mq_xp  <- agg_mq(X_m_xp, agg_m)
X_qq_xp  <- agg_qq(X_q_xp, agg_q)

X_xp_agg <- cbind(X_mq_xp, X_qq_xp)

# ==============================================================================
# 4. NUMBER OF PROXY : KREAMER & SUGIYAMA (2013) 
# ==============================================================================

pls.object <- select_L_autoproxy_3prf(X_xp_agg, y_q, Zmax = params$Zmax)
Lproxy     <- pls.object$L_opt

# ==============================================================================
# 5. NUMBER OF U-MIDAS LAGS
# ==============================================================================

lag_sel <- choose_UMIDAS_lag(
  X_lf        = X_xp_agg,
  X_hf        = X_xp,
  y_q         = y_q,
  Lmax        = params$Lmax,
  Lproxy      = Lproxy,
  p_AR_max    = params$p_AR_max,   
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
)

L_midas <- lag_sel$best_BIC$L
p_ar    <- lag_sel$best_BIC$p_AR

# ==============================================================================
# 6. MF-TPRF
# ==============================================================================

MF_TPRF_out <- MF_TPRF(
  X_lf        = X_xp_agg,
  X_hf        = X_xp,
  y_q         = y_q,
  Lproxy      = Lproxy,
  L_midas     = L_midas,
  p_AR        = p_ar,
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
)


# ==============================================================================
# 6.1 SAVE ALL RESULTS TO COUNTRY-SPECIFIC FOLDER
# ==============================================================================

## 14.1 Create folder for the country
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) {
  dir.create(path_country, recursive = TRUE)
}

## 14.2 Build filename
file_out <- file.path(
  path_country,
  paste0(
    "MF_TPRF_",      country,
    "_sel-",         params$sel_method,
    "_Nm-",          N_m,
    "_Nq-",          N_q,
    "_Lproxy-",      Lproxy,
    "_L_midas-",     L_midas,
    "_p-AR",         p_ar,
    "_Robust-F_",    params$Robust_F,
    "_Covid_m-",     params$covid_mask_m,
    "_Covid_q-",     params$covid_mask_q,
    ".rds"
  )
)

## 14.3 Save everything into one RDS
saveRDS(
  list(
    # Country info
    country = country,
    params  = params,
    
    # Imputation Results
    imp_xp     = imp_xp,
    X_xp       = X_xp,
    F_xp       = imp_xp$F_hat,
    lambda_xp  = imp_xp$lambda,
    X_mq_xp    = X_mq_xp,
    X_qq_xp    = X_qq_xp,
    X_xp_agg   = X_xp_agg,
    
    # MF-3PRF Results
    UMIDAS_lag = lag_sel,
    MF_TPRF    = MF_TPRF_out
  ),
  file_out
)

cat("\n*** SAVED MF-TPRF RESULTS TO ***\n", file_out, "\n")







# ==============================================================================
# 7. PSEUDO REAL-TIME FORECASTING EXERCISE
# ==============================================================================

# Pseudo real time nowcasting accounting for frequency mismatch and asyncronous
# release of macroeconomic indicators 

pseudo_realtime_raw <- pseudo_realtime_MF_TPRF_XP(
  X_full  = X,         # full mixed-frequency predictor matrix (with NA)
  y_q     = y_q,       # true quarterly target
  params  = params,    # all model parameters (include start_eval, end_eval, Lmax_midas)
  dates   = dates_m,   # monthly dates (for release schedule)
  dates_q = dates_q,   # quarterly publication dates
  Freq    = Freq,      # vector: monthly/quarterly frequency classification
  Unb     = Unb,       # unbalancedness info (lag di pubblicazione in mesi)
  agg_m   = agg_m,     # schema aggregazione mensile -> trimestrale
  agg_q   = agg_q      # schema aggregazione trimestrale
)

# piccola utility per rendere M1/M2/M3 in data frame puliti
list_to_df <- function(lst, tag) {
  if (length(lst) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  data.frame(
    date             = as.Date(names(lst)),
    nowcast          = as.numeric(unlist(lst)),
    month_in_quarter = tag,
    row.names        = NULL
  )
}


# ==============================================================================
# 16. SAVE REAL-TIME RESULTS (COUNTRY-SPECIFIC FOLDER)
# ==============================================================================

# Create directory for this country
path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) {
  dir.create(path_country, recursive = TRUE)
}

# Construct output filename
file_out <- file.path(
  path_country,
  paste0(
    "MF_TPRF_RT_",   country,
    "_sel-",         params$sel_method,
    "_Nm-",          N_m,
    "_Nq-",          N_q,
    "_Lproxy-",      pseudo_realtime_raw$hyper_pre$Lproxy,
    "_Lmidas-",      pseudo_realtime_raw$hyper_pre$L_midas,
    "_p-AR",         pseudo_realtime_raw$hyper_pre$p_AR,
    "_Robust-F_",    params$Robust_F,
    "_Covid_m-",     params$covid_mask_m,
    "_Covid_q-",     params$covid_mask_q,
    ".rds"
  )
)

# Remove existing file (fresh overwrite)
if (file.exists(file_out)) {
  file.remove(file_out)
}



df_M1 <- list_to_df(pseudo_realtime_raw$M1, "M1")
df_M2 <- list_to_df(pseudo_realtime_raw$M2, "M2")
df_M3 <- list_to_df(pseudo_realtime_raw$M3, "M3")

# tutti i nowcast insieme, ordinati per data e mese nel trimestre
df_pseudoRT_all <- rbind(df_M1, df_M2, df_M3)
if (nrow(df_pseudoRT_all) > 0L) {
  df_pseudoRT_all <- df_pseudoRT_all[order(df_pseudoRT_all$date,
                                           df_pseudoRT_all$month_in_quarter), ]
}

cat("file_out =\n", file_out, "\n")
nchar(file_out)
dir.exists(path_country)

# Save results in modo super riconoscibile
saveRDS(
  list(
    country          = country,
    params           = params,
    dates_m          = dates_m,
    dates_q          = dates_q,
    y_q              = y_q,
    
    # output grezzo della funzione pseudo_realtime_TPRF_EM (liste M1/M2/M3)
    pseudo_realtime_raw = pseudo_realtime_raw,
    
    # versioni tidy dei nowcast
    pseudo_realtime_M1  = df_M1,
    pseudo_realtime_M2  = df_M2,
    pseudo_realtime_M3  = df_M3,
    pseudo_realtime_all = df_pseudoRT_all
  ),
  file = file_out
)

cat("\n*** PSEUDO REAL-TIME RESULTS SAVED TO ***\n", file_out, "\n")

