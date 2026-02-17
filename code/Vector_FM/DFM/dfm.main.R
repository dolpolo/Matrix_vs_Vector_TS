# ==============================================================================
#             Dynamic Vector Factor Models and the EM Algorithm: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#                               - Main Script - 
# ==============================================================================
# This script estimates a Dynamic Vector Factor Model (DFM) by the EM-algorithm 
# and performs a pseudo real-time nowcasting exercise for GDP growth in the Euro
# Area's main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation and harmonization across countries
#   2. Estimation of the number of latent factors and model lags
#   4. Model initialization and parameter estimation using the EM algorithm on 
#      the whole dataset
#   5. Pseudo real-time rolling nowcasting of GDP
# The following DFM sections is clearly separated from the DMFM section reported
# in a related replication code.
#
# ==============================================================================
# Author      : Davide Delfino
# Institution : Alma Mater Studiorum - University of Bologna
# Dataset     : Barigozzi & Lissona (2024), EA-MD-QD
#
# Disclaimer  : DFM functions adapted from Bambura and Modugno (2014)
#
# ==============================================================================
# Notes       : Other EA countries such as BE, NL, AT, PT, IR, GR and different 
#               set of Monthly and Quarterly variables can be added
# ==============================================================================
# Script Type : Data Preparation / Estimation / Nowcasting
# ==============================================================================


# ==============================================================================
# 0. SET WORKING DIRECTORY & PATHS
# ==============================================================================

# Paths: data, functions, results, graphs
path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "data")
path_func    <- file.path(path_main, "Vector_FM/functions")
path_results <- file.path(path_main, "Vector_FM/DFM/results/outputs")
path_graph   <- file.path(path_main, "Vector_FM/DFM/results/graph")



# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
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


## I/O
library(readxl)

## Plotting
library(ggplot2)

## Resolve conflicts
library(conflicted)
conflict_prefer("select",  "dplyr", quiet = TRUE)
conflict_prefer("filter",  "dplyr", quiet = TRUE)



# ==============================================================================
# 2. SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

source(file.path(path_func, "dfm.prep.R"))
source(file.path(path_func, "dfm.fs.R"))
source(file.path(path_func, "dfm.kalman.R"))
source(file.path(path_func, "dfm.init.R"))
source(file.path(path_func, "dfm.em.R"))
source(file.path(path_func, "dfm.now.R"))


# ==============================================================================
# 3. PARAMETERS & COUNTRY LIST & COUNTRY OF INTEREST SELECTION
# ==============================================================================

# Users's parameters
params <- list(
  start_est    = as.Date("2000-04-01"),    # Accounting from missingenss in growth rates wrt 1999
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-10-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,                    # Boolean: if empty == TRUE
  covid_mask_q = FALSE,                    # Boolean: if empty == TRUE
  target       = "GDP",                    # quarterly target variable
  sel_method   = "none",         # "none" | "corr_threshold" | "t_test" | "LASSO"
  
  # Variable selection parameters
  n_m          = 38,
  n_q          = 10,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,                     # 1=LASSO, 0=Ridge, (0,1)=Elastic Net
  
  # Factor selection parameters
  Kmax = 10,                            # max number of factors
  pmax = 10,                            # max number of lags in var factors
  qmax = 10,                            # max number of lags in ar idio
  
  # DFM Parameters
  kappa = 1e-4,
  restr = "stock_flow"                  # "MM" or "stock_flow"
  
)

# Country List
countries <- c("DE", "FR", "IT", "ES")

# Selected country
country <- "FR"

# ==============================================================================
# 4. PREPARE DATA FOR ALL COUNTRIES AND EXCTARCT THE ONE OF INTEREST
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



# ==============================================================================
# 5. TARGET AND PREDICTORS IDENTIFICATION
# ==============================================================================

gdp_col <- all_countries$data[[country]]$target_col

y   <- as.matrix(data[, gdp_col, drop = FALSE])     # target
y_q <- y[!is.na(y)]                                 # quarterly target

# ==============================================================================
# 6. METADATA EXTRACTION
# ==============================================================================

NM      <- all_countries$data[[country]]$nM
NQ      <- all_countries$data[[country]]$nQ
N       <- NM + NQ

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

agg_m    <- all_countries$data[[country]]$agg_m
agg_q    <- all_countries$data[[country]]$agg_q
agg      <- c(agg_m, agg_q)



# ============================================================================
# 7. STANDARDIZATION (WITH NA)
# ============================================================================

out_std <- standardize_with_na(data)
X_std   <- out_std$X_std   # T × N matrix with NA


# ============================================================================
# 8. FACTOR SELECTION and VAR ORDER SELECTION
# ============================================================================

# All purpose estimator
cov_proxy_out <- all_purpose_covariance(X_std)

# Eigenvalue ratio
ER_proxy_out  <- select_num_factors_ER(cov_proxy_out$Sigma_tilde, Kmax = params$Kmax)
r             <- ER_proxy_out$r
eig           <- ER_proxy_out$eig

# Plot diagnostica selezione fattori
plot_factor_selection_diag(ER_proxy_out, Kmax = params$Kmax, r_sel = r)

# Stima fattori XP + plot
xp <- estimate_factors_XP(X_std, r = r)
plot_factors_XP(xp$F_hat, dates = dates_m, ncol = 1)

L_hat  <- xp$Lambda     # N x r
F_hat  <- xp$F_hat      # T x r
C_hat  <- xp$C_hat      # T x N  (componente comune, scala std)

# ------------------------------------------------------------------
# VAR order on factors
# ------------------------------------------------------------------
factor_var_ic <- select_p_var_ic(F_hat, pmax = params$pmax, include_const = TRUE)
p <- factor_var_ic$p_BIC

# ------------------------------------------------------------------
# AR order on idiosyncratic components
# ------------------------------------------------------------------
# 1) Residui idio
Ehat <- X_std - C_hat  # mantiene NA dove X_std è NA

# 2) Selezione ordine AR globale (consiglio: include_const = FALSE su dati standardizzati)
idio_ar_ic <- select_p_ar_ic_global(Ehat[, 1:NM, drop = FALSE],   # solo mensili (consigliato)
                                    pmax = params$qmax,
                                    include_const = FALSE,
                                    min_T = 30)

q <- idio_ar_ic$p_BIC

p <- 1
q <- 1
# ============================================================================
# 9. INITIALIZATION FOR THE DFM (Banbura–Modugno)
# ============================================================================

Init <- InitialCond(
  X_std    = X_std,
  r        = r,
  p_factor = p,
  q_idio   = q,
  NM       = NM,
  NQ       = NQ,
  restr    = params$restr,
  agg      = agg,              # length N or length NQ (agg_q)
  kappa    = params$kappa
)

# ============================================================================
# 11. EM ESTIMATION (DFM — BANBURA & MODUGNO)
# ============================================================================

res <- DFM_EM(
  X        = X_std,
  Init     = Init,
  max_iter = 200,
  tol      = 1e-3
)

# ============================================================================
# 12. DE-STANDARDIZE (SAVE MEAN/SD + OPTIONAL FITTED ON ORIGINAL SCALE)
# ============================================================================

# helpers
destd_vec <- function(x_std, mu, sd) as.numeric(mu + sd * x_std)

destd_mat <- function(X_std, mu, sd) {
  X_std <- as.matrix(X_std)
  stopifnot(length(mu) == ncol(X_std), length(sd) == ncol(X_std))
  sweep(sweep(X_std, 2, sd, `*`), 2, mu, `+`)
}

# store mean/sd used for standardization (needed to bring anything back to original scale)
mu_full <- out_std$mean
sd_full <- out_std$sd

# OPTIONAL but very useful: fitted values / smoothed signal on std + original scale
# (uses the same window as EM/Kalman routines if you pass n_lags_Q accordingly)
n_lags_Q_used <- if (params$restr == "MM") 5L else 3L   # if your framework fixes these; otherwise store from Init$meta$L

fit_obj <- dfm_kalman_fit(X_in = X_std, res = res, n_lags_Q = n_lags_Q_used)
X_fit_std  <- fit_obj$X_fit                     # T_eff x N (std scale)
X_fit_orig <- destd_mat(X_fit_std, mu_full, sd_full)

# GDP fitted / smoothed series (orig scale) if you want it ready
gdp_fit_std  <- X_fit_std[,  gdp_col]
gdp_fit_orig <- X_fit_orig[, gdp_col]

# ============================================================================
# 13. SAVE DFM RESULTS (FULL SAMPLE ESTIMATION) — INCLUDING ORIGINAL SCALE INFO
# ============================================================================

path_country <- file.path(path_results, country)
if (!dir.exists(path_country)) dir.create(path_country, recursive = TRUE)

covid_tag <- paste0("M", as.integer(params$covid_mask_m), "_Q", as.integer(params$covid_mask_q))

file_out_dfm <- file.path(
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

saveRDS(list(
  country = country,
  params  = params,
  
  # Input data
  X_raw   = data,      # original scale
  X_std   = X_std,     # standardized
  dates_m = dates_m,
  dates_q = dates_q,
  
  # Standardization maps (CRITICAL to go back to original scale)
  std_map = list(mean = mu_full, sd = sd_full),
  
  # Factor extraction (note: factors are latent; no "original scale")
  r     = r,
  F_hat = F_hat,
  L_hat = L_hat,
  
  # Initialization
  Init  = Init,
  
  # EM results
  A      = res$A,
  C      = res$C,
  Q      = res$Q,
  R      = res$R,
  Z0     = res$Z0,
  V0     = res$V0,
  loglik = res$loglik,
  crit   = res$crit,
  
  # Fitted values (std + original scale) + indices
  fit = list(
    n_lags_Q = n_lags_Q_used,
    t_index  = fit_obj$t_index,      # indices in the original monthly timeline
    X_fit_std  = X_fit_std,
    X_fit_orig = X_fit_orig,
    gdp_fit_std  = gdp_fit_std,
    gdp_fit_orig = gdp_fit_orig
  ),
  
  # Hyperselection diagnostics
  p = p,
  q = q,
  factor_var_ic = factor_var_ic,
  idio_ar_ic    = idio_ar_ic,
  ER_proxy_out  = ER_proxy_out
  
), file_out_dfm)

cat("\n*** SAVED DFM RESULTS TO ***\n", file_out_dfm, "\n")














# ==============================================================================
# 13. ROLLING REAL-TIME NOWCAST (DFM)
#   NOTE: fixed typos + save BOTH standardized and original-scale nowcasts
#         + save the rolling standardization maps used at each tt (for audit)
# ==============================================================================

now_dfm_re <- pseudo_realtime_DFM_EM_reestimate(
  X_full   = data,
  NQ       = NQ,
  params   = params,
  dates_m  = dates_m,
  dates_q  = dates_q,
  Freq     = Freq,
  Unb      = Unb,
  gdp_col  = gdp_col,
  agg      = agg,
  
  max_iter_em = 200,
  tol_em      = 1e-3,
  pmax        = params$pmax,
  qmax        = params$qmax
)

# ==============================================================================
# 14. SAVE RESULTS – ROLLING NOWCAST (DFM)  (std + original + metadata)
# ==============================================================================

path_results_dfm_country <- file.path(path_results, country)
if (!dir.exists(path_results_dfm_country)) dir.create(path_results_dfm_country, recursive = TRUE)

start_est_str <- format(params$start_est, "%Y-%m")
end_eval_str  <- format(params$end_eval, "%Y-%m")

sel_method <- if (!is.null(params$sel_method)) params$sel_method else "none"
covid_tag  <- paste0("M", as.integer(params$covid_mask_m), "_Q", as.integer(params$covid_mask_q))

file_dfm_rolling <- file.path(
  path_results_dfm_country,
  paste0(
    "DFM_RollingNowcast_", country,
    "_sel-", sel_method,
    "_Covid-", covid_tag,
    "_Nm-", NM,
    "_Nq-", NQ,
    "_rfix-", now_dfm_re$r_fix,
    "_pfix-", now_dfm_re$p_fix,
    "_", start_est_str,
    "_to_", end_eval_str,
    ".rds"
  )
)

# -------------------------------------------------------------------
# CONTENT TO SAVE
#   - now_dfm_re should ALREADY include:
#       r_fix, p_fix,
#       M1_std/M2_std/M3_std, M1_orig/M2_orig/M3_orig
#   - we add: std_map_by_t (optional), inputs metadata, and a compact "panel"
# -------------------------------------------------------------------

dfm_rolling <- list(
  country = country,
  params  = params,
  
  # predictor metadata (useful for reproducing ragged-edge)
  Freq = Freq,
  Unb  = Unb,
  agg  = agg,
  
  # time info
  dates_m = dates_m,
  dates_q = dates_q,
  gdp_col = gdp_col,
  
  # rolling nowcasts (std + original)
  nowcast = now_dfm_re
)

saveRDS(dfm_rolling, file_dfm_rolling)
cat("\n*** SAVED DFM ROLLING NOWCAST TO ***\n", file_dfm_rolling, "\n")


# ==============================================================================
# OPTIONAL (HIGHLY RECOMMENDED):
# build a clean dataframe for plots/eval and save it too
# ==============================================================================

rolling_to_df <- function(now_obj) {
  # now_obj = now_dfm_re
  make_df <- function(x, label, scale) {
    if (length(x) == 0) return(NULL)
    data.frame(date = as.Date(names(x)), value = as.numeric(x),
               vintage = label, scale = scale, row.names = NULL)
  }
  rbind(
    make_df(now_obj$M1_std,  "M1", "std"),
    make_df(now_obj$M2_std,  "M2", "std"),
    make_df(now_obj$M3_std,  "M3", "std"),
    make_df(now_obj$M1_orig, "M1", "orig"),
    make_df(now_obj$M2_orig, "M2", "orig"),
    make_df(now_obj$M3_orig, "M3", "orig")
  )
}

df_nowcasts_long <- rolling_to_df(now_dfm_re)

file_dfm_rolling_df <- sub("\\.rds$", "_longdf.rds", file_dfm_rolling)
saveRDS(df_nowcasts_long, file_dfm_rolling_df)
cat("\n*** SAVED LONG DF TO ***\n", file_dfm_rolling_df, "\n")
