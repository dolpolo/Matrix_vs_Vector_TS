# ==============================================================================
#                    Matrix Three Pass Regression Filter: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#
# ==============================================================================
# This script estimates a Matrix Three Pass regression filter (Matrix TPRF) for 
# Mixed Frequency data and performs a pseudo real-time nowcasting exercise for 
# GDP growth in the  Euro Area's main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation, tensor construction, and harmonization across countries
#   2. NA's Imputation in predictors panel         --> Cen & Lam (2025)
#   3. Estimation of the number of latent factors  --> Yu et al. (2020)
#   4. Estimation of the number of Proxies         --> Kraemer & Sugiyama (2011)
#   5, Estimation of U-MIDAS and AR lags number    --> IC
#   6. Full sample estimation with Matrix MF-TPRF
#   7. Pseudo real-time rolling nowcast of GDP
#
# ==============================================================================
# Author          : Davide Delfino
# Institution     : Bocconi University
#
# ===============================================================================
# Main Reference  : Barigozzi Delfino Marcellino (2026), Matrix MF-TPRF
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
# - Fit
# - Nowcast 
# ==============================================================================


# ==============================================================================
# SET WORKING DIRECTORY & PATHS
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "data")
path_func    <- file.path(path_main, "Matrix_FM/functions")
path_results <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph")

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
library(fBasics)
library(vars)
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
source(file.path(path_func, "matrix.mf.tprf.utils.R"))
source(file.path(path_func, "matrix.mf.tprf.prep.R"))

# Matrix mf-tprf estimation
source(file.path(path_func, "matrix.mf.tprf.R"))

# NANs imputation
source(file.path(path_func, "matrix.mf.tprf.imp.R"))

# factor selection
source(file.path(path_func, "matrix.mf.tprf.fs.R"))

# nowcasting
source(file.path(path_func, "matrix.mf.tprf.now.R"))

# ==============================================================================
# PARAMETERS & COUNTRY LIST
# ==============================================================================

# Users's parameters
params <- list(
  start_est    = as.Date("2000-04-01"),   # growth rate 2001-01 NA
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-10-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = FALSE,                    # Boolean: if empty == TRUE
  covid_mask_q = FALSE,                    # Boolean: if empty == TRUE
  target       = "GDP",                   # quarterly target variable
  target_cc    = "EA",                    # Aggregate from where to extract the proxy
  
  # Variable selection parameters
  sel_method   = "none",                   # "none" | "corr_threshold" | "F-Test" | "LASSO"
  n_m          = 30,
  n_q          = 20,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,                        # 1=LASSO, 0=Ridge, (0,1)=Elastic Net
  
  
  # Factor and proxy selection parameters
  Kmax     = c(1,3),                           # max number of factors
  Kmax_vec = 10,
  Zmax     = 10,                               # max number of proxy
  
  # MF-TPRF parameters
  p_AR_max    = 5,                         # Max numbers of lags y in U-MIDAS
  Lmax        = 5,                         # Max numbers of lags in the Factors in U-MIDAS
  Robust_F    = FALSE,                     # Robust F test in the first step
  alpha       = 0.1,
  robust_type = "NW",                      # White--> HC | NW(Newey West)--> HAC
  nw_lag      = 1                          # Lag in the error autocorrelation in step 1
)

# Country List
countries <- c("DE", "FR", "IT", "ES", "EA")

# Selected country
country <- "FR"

# ==============================================================================
# 1 PREPARE DATA FOR ALL COUNTRIES
# ==============================================================================
all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
  covid_mask_m    = params$covid_mask_m,
  covid_mask_q    = params$covid_mask_q
)

# ==================
# 1.1 PREPARE TENSOR
# ==================

build_tensor <- function(prep,
                         params,
                         var_scope = c("union", "intersection", "threshold"),
                         min_presence = NULL) {
  
  var_scope <- match.arg(var_scope)
  
  countries <- names(prep$data)
  data_list <- prep$data
  
  # ------------------------------------------------------------
  # 1. ASSE TEMPORALE GLOBALE da params (mensile)
  # ------------------------------------------------------------
  all_dates <- seq(params$start_est, params$end_eval, by = "month")
  T_global  <- length(all_dates)
  
  # ------------------------------------------------------------
  # 2. LISTE DI VARIABILI BASE (M e Q) PER OGNI PAESE
  # ------------------------------------------------------------
  series_list <- lapply(data_list, function(obj) obj$Series)
  
  monthly_base_list <- lapply(data_list, function(obj) {
    Series_m <- obj$Series[1:obj$nM]  # prime nM sono mensili
    remove_country_code(Series_m, countries = countries)
  })
  
  quarterly_base_list <- lapply(data_list, function(obj) {
    Series_q <- obj$Series[(obj$nM + 1):(obj$nM + obj$nQ)]  # poi nQ trimestrali
    remove_country_code(Series_q, countries = countries)
  })
  
  # ------------------------------------------------------------
  # 3. UNION / INTERSECTION / THRESHOLD SU MENSILI E TRIMESTRALI
  # ------------------------------------------------------------
  if (var_scope == "intersection") {
    base_M <- Reduce(intersect, monthly_base_list)
    base_Q <- Reduce(intersect, quarterly_base_list)
    
  } else if (var_scope == "union") {
    base_M <- sort(unique(unlist(monthly_base_list)))
    base_Q <- sort(unique(unlist(quarterly_base_list)))
    
  } else if (var_scope == "threshold") {
    
    if (is.null(min_presence)) {
      min_presence <- ceiling(length(countries) / 2)
    }
    
    count_M <- table(unlist(lapply(monthly_base_list, unique)))
    count_Q <- table(unlist(lapply(quarterly_base_list, unique)))
    
    base_M <- sort(names(count_M)[count_M >= min_presence])
    base_Q <- sort(names(count_Q)[count_Q >= min_presence])
  }
  
  # Ordine finale: prima tutte le mensili, poi tutte le trimestrali
  vars_ordered <- c(base_M, base_Q)
  
  if (length(vars_ordered) == 0L) {
    stop("Nessuna variabile selezionata dopo union/intersection/threshold.")
  }
  
  P1 <- length(countries)
  P2 <- length(vars_ordered)
  
  # ------------------------------------------------------------
  # 4. INIZIALIZZA TENSORE Y (DATI) E W (MASK)
  # ------------------------------------------------------------
  Y <- array(NA_real_,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               vars_ordered
             ))
  
  W <- array(0L,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               vars_ordered
             ))
  
  # ------------------------------------------------------------
  # 4bis. INIZIALIZZA METADATA [Paesi x Variabili]
  # ------------------------------------------------------------
  agg_mat   <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  class_mat <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  type_mat  <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  freq_mat  <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  unb_mat   <- matrix(NA_real_,       nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  
  # ------------------------------------------------------------
  # 5. RIEMPI Y, W E METADATA PAESE PER PAESE
  # ------------------------------------------------------------
  for (i in seq_along(countries)) {
    cc   <- countries[i]
    obj  <- data_list[[cc]]
    
    Data_cc   <- obj$Data        # matrice (T_cc x n_var_cc)
    Dates_cc  <- obj$Dates
    Series_cc <- obj$Series      # (mensili + trimestrali)
    
    base_cc <- remove_country_code(Series_cc, countries = countries)
    
    # Mappa date paese -> asse globale
    idx_time <- match(Dates_cc, all_dates)
    
    # Matrice temporanea per questo paese
    Y_cc <- matrix(NA_real_,
                   nrow = T_global, ncol = P2,
                   dimnames = list(as.character(all_dates), vars_ordered))
    
    for (j in seq_along(vars_ordered)) {
      v_base <- vars_ordered[j]
      
      col_idx <- which(base_cc == v_base)
      if (length(col_idx) == 0L) next
      
      # se più colonne con stesso base name → prendo la prima
      col_idx <- col_idx[1]
      
      # -----------------------------
      # 5a. DATI + MASK
      # -----------------------------
      Y_cc[idx_time, j] <- Data_cc[, col_idx]
      
      # -----------------------------
      # 5b. METADATA (M vs Q)
      # -----------------------------
      if (col_idx <= obj$nM) {
        # MENSILE
        k <- col_idx
        
        agg_mat[i, j]   <- obj$agg_m[k]
        class_mat[i, j] <- obj$ClassM[k]
        type_mat[i, j]  <- obj$TypeM[k]
        freq_mat[i, j]  <- obj$freq_m[k]
        unb_mat[i, j]   <- obj$unb_m[k]
        
      } else {
        # TRIMESTRALE
        k <- col_idx - obj$nM
        
        agg_mat[i, j]   <- obj$agg_q[k]
        class_mat[i, j] <- obj$ClassQ[k]
        type_mat[i, j]  <- obj$TypeQ[k]
        freq_mat[i, j]  <- obj$freq_q[k]
        unb_mat[i, j]   <- obj$unb_q[k]
      }
    }
    
    Y[, i, ] <- Y_cc
    W[, i, ] <- ifelse(is.na(Y_cc), 0L, 1L)
  }
  
  # ------------------------------------------------------------
  # 6. CLASSIFICA VARIABILI COME "M" / "Q"
  # ------------------------------------------------------------
  var_type <- rep(NA_character_, P2)
  names(var_type) <- vars_ordered
  var_type[vars_ordered %in% base_M] <- "M"
  var_type[vars_ordered %in% base_Q] <- "Q"
  
  n_M <- sum(var_type == "M")
  n_Q <- sum(var_type == "Q")
  
  idx_M <- which(var_type == "M")
  idx_Q <- which(var_type == "Q")
  
  vars_M <- vars_ordered[idx_M]
  vars_Q <- vars_ordered[idx_Q]
  
  # ------------------------------------------------------------
  # 7. SPLIT METADATA IN MENSILI / TRIMESTRALI
  # ------------------------------------------------------------
  agg_M   <- if (n_M > 0)   agg_mat[, idx_M, drop = FALSE]   else NULL
  agg_Q   <- if (n_Q > 0)   agg_mat[, idx_Q, drop = FALSE]   else NULL
  
  ClassM  <- if (n_M > 0) class_mat[, idx_M, drop = FALSE]   else NULL
  ClassQ  <- if (n_Q > 0) class_mat[, idx_Q, drop = FALSE]   else NULL
  
  TypeM   <- if (n_M > 0) type_mat[, idx_M, drop = FALSE]    else NULL
  TypeQ   <- if (n_Q > 0) type_mat[, idx_Q, drop = FALSE]    else NULL
  
  freq_M  <- if (n_M > 0) freq_mat[, idx_M, drop = FALSE]    else NULL
  freq_Q  <- if (n_Q > 0) freq_mat[, idx_Q, drop = FALSE]    else NULL
  
  unb_M   <- if (n_M > 0)  unb_mat[, idx_M, drop = FALSE]    else NULL
  unb_Q   <- if (n_Q > 0)  unb_mat[, idx_Q, drop = FALSE]    else NULL
  
  # ------------------------------------------------------------
  # 8. TARGET (base name = params$target, es. "GDP")
  # ------------------------------------------------------------
  target_base <- tolower(params$target)
  target_col  <- which(tolower(vars_ordered) == target_base)
  
  # ------------------------------------------------------------
  # 9. OUTPUT
  # ------------------------------------------------------------
  out <- list(
    Y          = Y,               # [T x Paesi x Variabili]
    W          = W,               # mask [T x Paesi x Variabili]
    dates      = all_dates,
    countries  = countries,
    vars       = vars_ordered,
    var_type   = var_type,        # "M" / "Q"
    n_M        = n_M,
    n_Q        = n_Q,
    idx_M      = idx_M,
    idx_Q      = idx_Q,
    vars_M     = vars_M,
    vars_Q     = vars_Q,
    target_col = target_col,
    
    # metadata globali [Paese x Variabile]
    agg        = agg_mat,
    class      = class_mat,
    type       = type_mat,
    freq       = freq_mat,
    unb        = unb_mat,
    
    # metadata separati M / Q (stile prepare_country_data)
    agg_M      = agg_M,
    agg_Q      = agg_Q,
    ClassM     = ClassM,
    ClassQ     = ClassQ,
    TypeM      = TypeM,
    TypeQ      = TypeQ,
    freq_M     = freq_M,
    freq_Q     = freq_Q,
    unb_M      = unb_M,
    unb_Q      = unb_Q
  )
  
  return(out)
}

tensor <- build_tensor(all_countries, params, "intersection")
# tensor <- build_tensor(all_countries, params, var_scope = "threshold", min_presence = 3)
Y <- tensor$Y
W <- tensor$W

nan_percent_Y(Y)

# ====================================
# 1.2 TARGET AND PREDICTORS EXTRAXTION
# ====================================

# target vector extraction (GDP growth)
gdp_col <- tensor$target_col

# The quarterly growth in the first quarter is NA (NO EA before the 2000)
y       <- Y[ , ,gdp_col, drop = FALSE]                               # Monthly Target
y_q     <- y[which(rowSums(!is.na(y[ , , 1])) > 0), , , drop = TRUE]  # Quarterly Target
y_EA_q  <- all_countries$proxy$Y_q                                    # Proxy GDP EA
y_q_all <- cbind(y_EA_q,y_q)

# predictors panel extraction
X     <- Y[, ,-gdp_col, drop = FALSE]    # predictors
W_x   <- W[, ,-gdp_col, drop = FALSE]    # predictors' mask

# =======================
# 1.3 METADATA EXTRACTION
# =======================

# dates 
dates_m <- as.Date(dimnames(tensor$Y)[[1]])                     # monthly dates
dates_q <- dates_m[which(rowSums(!is.na(y[, , 1])) > 0)]   # quarterly dates

# Panel dimensions
N_m     <- tensor$n_M                                   # monthly indicators

N_q_tot      <- tensor$n_Q
q_series_all <- tensor$vars[(N_m + 1):(N_m + N_q_tot)]
q_cols       <- which(!grepl(params$target, q_series_all, ignore.case = TRUE))
N_q          <- length(q_cols)                          # Quarterly series

N <- N_m + N_q                                          # Tot number of indicators

# Aggregation rule (stock & flows)
agg_m   <- tensor$agg_M
agg_q   <- tensor$agg_Q[, q_cols, drop = FALSE]
agg     <- cbind(agg_m, agg_q)

# Frequency (monthly or quarterly)
Freq_m  <- tensor$freq_M
Freq_q  <- tensor$freq_Q[,q_cols, drop = FALSE]
Freq    <- cbind(Freq_m, Freq_q)

# Indicators' releasing delay
Unb_m   <- tensor$unb_M
Unb_q   <- tensor$unb_Q[,q_cols, drop = FALSE]
Unb     <- cbind(Unb_m, Unb_q)

# Indicators' class (Real, Nominal, Financial, Survey)
Class_m <- tensor$ClassM
Class_q <- tensor$ClassQ[,q_cols, drop = FALSE]
Class   <- cbind(Class_m, Class_q)

# ==============================================================================
# 2. NA's IMPUTATION          : Cen & Lam (2025)
# 3. LATENT FACTOR EXTRACTION : Yu et al  (2020)
# ==============================================================================

# Standardization
out_std <- standardize_mat_with_na(X)
X_std   <- out_std$X_scaled

# Imputation 
imp_cl <- init_CL_Yu(X_std, W_x, params$Kmax)

X_cl   <- imp_cl$Y_init
r_hat  <- imp_cl$r

# ================================
# 3.1 AGGREGATION TO LOW FREQUENCY
# ================================

X_cl_q <- aggregate_tensor_to_quarterly(
  X_tens = X_cl,
  agg    = agg, 
  N_m    = N_m,
  N_q    = N_q
)

# ==============================================================================
# 4. NUMBER OF PROXY : KREAMER & SUGIYAMA (2013) 
# ==============================================================================

# vectorization
X_cl_q_vec  <- tensor_to_vector(X_tens = X_cl_q, N_m = N_m, N_q = N_q)

# IC. Kraemer & Sugiyama
pls.object <- select_L_autoproxy_3prf(X_cl_q_vec, y_EA_q, Zmax = params$Zmax)
Lproxy     <- pls.object$L_opt

# ==============================================================================
# 5. NUMBER OF U-MIDAS LAGS
# ==============================================================================

lag_sel <- choose_UMIDAS_grid_tensor_MF(
  X_lf      = X_cl_q,
  X_hf      = X_cl,
  y_q       = y_EA_q,              # target EA country
  r         = r_hat,               # CHECK !!!!!!
  Lproxy    = Lproxy,
  Lmax      = params$Lmax,
  p_AR_max  = params$p_AR_max,
  use_autoproxy = TRUE,
  y_proxy   = y_EA_q,              # proxy usata per costruire Z (EA GDP)
  standardize_proxy = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  ils_maxit = 100,
  ils_tol   = 1e-8
)

L_midas <- lag_sel$best_BIC$L
p_ar    <- lag_sel$best_BIC$p_AR


# ==============================================================================
# 6. Matrix MF-TPRF
# ==============================================================================

Tensor_MF_TPRF_out <- Tensor_MF_TPRF(
  X_lf        = X_cl_q,
  X_hf        = X_cl,
  Y_q_all     = y_q_all,
  proxy_name  = "EA",
  Lproxy      = Lproxy,
  L_midas     = L_midas,
  p_AR        = p_ar,             
  r           = r_hat,
  standardize_proxy = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  ils_maxit = 100,
  ils_tol   = 1e-8
)

# ====================
# 6.1 SAVE ALL RESULTS 
# ====================

path_tensor <- path_results
if (!dir.exists(path_tensor)) dir.create(path_tensor, recursive = TRUE)

r1 <- r_hat[1]; r2 <- r_hat[2]
cc_lab <- paste(countries, collapse = "-")

file_out <- file.path(
  path_tensor,
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

saveRDS(
  list(
    # General
    params    = params,
    countries = countries,
    tensor    = tensor,
    
    # --- Imputation / factor-selection stage (Cen&Lam + Yu style)
    out_std = out_std,          # (optional) contains scaling info; keep if you need to destandardize
    imp_cl  = imp_cl,           # output of init_CL_Yu()
    X_std   = X_std,            # (optional, can be huge) remove if not needed
    X_cl    = X_cl,             # monthly imputed tensor used in MF-TPRF
    X_cl_q  = X_cl_q,           # quarterly aggregated tensor used in Pass 1
    
    # Selected hyperparams
    r_hat   = r_hat,
    Lproxy  = Lproxy,
    L_midas = L_midas,
    p_ar    = p_ar,
    
    # Matrix MF-TPRF
    Tensor_MF_TPRF = Tensor_MF_TPRF_out
  ),
  file_out
)

cat("\n*** SAVED TENSOR MF-TPRF RESULTS TO ***\n", file_out, "\n")



# ==============================================================================
# 7. PSEUDO REAL-TIME FORECASTING EXERCISE 
# ==============================================================================
roll_tensor <- pseudo_realtime_Tensor_MF_TPRF(
  X_full       = X,          # predictors tensor monthly (with NA)
  W_full       = W_x,        # mask tensor 1/0 (same dim as X)
  Y_q_all_full = y_q_all,    # quarterly: EA + countries
  proxy_name   = "EA",
  params       = params,
  dates_m      = dates_m,    # monthly dates
  dates_q      = dates_q,    # quarterly dates
  agg          = agg,        # aggregation rules
  N_m          = N_m,
  N_q          = N_q,
  do_recalib   = TRUE
)

list_to_df_country <- function(lst_by_country, tag) {
  if (is.null(lst_by_country) || length(lst_by_country) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      country          = character(),
      nowcast          = numeric(),
      month_in_quarter = character()
    ))
  }
  
  out <- lapply(names(lst_by_country), function(cc) {
    x <- lst_by_country[[cc]]
    if (is.null(x) || length(x) == 0L) return(NULL)
    
    data.frame(
      date             = as.Date(names(x)),
      country          = cc,
      nowcast          = as.numeric(unlist(x)),
      month_in_quarter = tag
    )
  })
  
  dplyr::bind_rows(out)
}

df_M1 <- list_to_df_country(roll_tensor$M1, "M1")
df_M2 <- list_to_df_country(roll_tensor$M2, "M2")
df_M3 <- list_to_df_country(roll_tensor$M3, "M3")

df_pseudoRT_all <- dplyr::bind_rows(df_M1, df_M2, df_M3) |>
  dplyr::mutate(month_in_quarter = factor(month_in_quarter, levels = c("M1","M2","M3"))) |>
  dplyr::arrange(country, date, month_in_quarter)

# folder
if (!dir.exists(path_results)) dir.create(path_results, recursive = TRUE)

countries_eval <- setdiff(colnames(y_q_all), "EA")
cc_lab <- paste(countries_eval, collapse = "-")

file_out_rt <- file.path(
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

if (file.exists(file_out_rt)) file.remove(file_out_rt)

saveRDS(
  list(
    params      = params,
    proxy_name  = "EA",
    countries   = countries_eval,
    dates_m     = dates_m,
    dates_q     = dates_q,
    Y_q_all     = y_q_all,
    
    pseudo_rt_raw = roll_tensor,
    
    pseudo_rt_M1  = df_M1,
    pseudo_rt_M2  = df_M2,
    pseudo_rt_M3  = df_M3,
    pseudo_rt_all = df_pseudoRT_all
  ),
  file = file_out_rt
)

file_out_rt

RT_res <- readRDS(file_out_rt)







################################################################################
############################### VECTORIZATION ##################################
################################################################################

# ==============================================================================
# SET WORKING DIRECTORY & PATHS
# ==============================================================================

path_func_vec    <- file.path(path_main, "Vector_FM/functions")
path_results_vec <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs_vec")
path_graph_vec   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph_vec")

# ==============================================================================
# SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================
# data preparation
source(file.path(path_func_vec, "mf.tprf.utils.R"))
source(file.path(path_func_vec, "mf.tprf.prep.R"))

# mf-tprf estimation
source(file.path(path_func_vec, "mf.tprf.R"))

# imputation
source(file.path(path_func_vec, "mf.tprf.imp.R"))

# factor selection
source(file.path(path_func_vec, "mf.tprf.fs.R"))

# nowcasting
source(file.path(path_func_vec, "mf.tprf.now.R"))

# ==============================================================================
# 1 TENSOR VECTORIZATION
# ==============================================================================

tensor <- build_tensor(all_countries, params, "intersection")
data   <- tensor_to_vector(X_tens = tensor$Y, N_m = tensor$n_M, N_q = tensor$n_Q)

# ====================================
# 1.2 TARGET AND PREDICTORS EXTRACTION
# ====================================

all_names   <- colnames(data)
target_name <- paste0(country, "_", params$target)

gdp_col <- which(all_names == target_name)

y     <- as.matrix(data[, gdp_col, drop = FALSE])
y_q   <- as.matrix(y[!is.na(y[, 1]), , drop = FALSE])
X_vec <- as.matrix(data[, -gdp_col, drop = FALSE])

pred_names <- colnames(X_vec)

# =======================
# 1.3 METADATA EXTRACTION
# =======================

dates_m <- as.Date(dimnames(tensor$Y)[[1]])
dates_q <- dates_m[!is.na(y[, 1])]

# extract base variable names from vectorized columns: "FR_IP" -> "IP"
get_varname <- function(x) sub("^[^_]+_", "", x)

pred_var <- get_varname(pred_names)

# remove target variable from predictors for all countries, if present
keep <- pred_var != params$target

X_vec      <- X_vec[, keep, drop = FALSE]
pred_names <- pred_names[keep]
pred_var   <- pred_var[keep]

# base monthly and quarterly variable names
vars_m <- tensor$vars[1:tensor$n_M]
vars_q <- tensor$vars[(tensor$n_M + 1):(tensor$n_M + tensor$n_Q)]

# identify monthly and quarterly predictors in the vectorized panel
is_monthly   <- pred_var %in% vars_m
is_quarterly <- pred_var %in% vars_q

# dimensions after vectorization
N_m <- sum(is_monthly)
N_q <- sum(is_quarterly)
N   <- N_m + N_q

# aggregation rules aligned with X_vec columns
agg_m_lookup <- setNames(as.vector(tensor$agg_M), vars_m)
agg_q_lookup <- setNames(as.vector(tensor$agg_Q), vars_q)

agg_m <- matrix(agg_m_lookup[pred_var[is_monthly]], nrow = 1)
agg_q <- matrix(agg_q_lookup[pred_var[is_quarterly]], nrow = 1)

# frequency and release-delay rules aligned with X_vec columns
freq_m_lookup <- setNames(as.vector(tensor$freq_M), vars_m)
freq_q_lookup <- setNames(as.vector(tensor$freq_Q), vars_q)

unb_m_lookup  <- setNames(as.vector(tensor$unb_M), vars_m)
unb_q_lookup  <- setNames(as.vector(tensor$unb_Q), vars_q)

Freq <- c(
  freq_m_lookup[pred_var[is_monthly]],
  freq_q_lookup[pred_var[is_quarterly]]
)

Unb <- c(
  unb_m_lookup[pred_var[is_monthly]],
  unb_q_lookup[pred_var[is_quarterly]]
)
# ==============================================================================
# 2. NA's IMPUTATION          : XIONG & PELGER (2020)
# 3. LATENT FACTOR EXTRACTION : Ahn (2011)
# ==============================================================================

# Standardization
out_std <- standardize_with_na(X_vec)
X_std   <- out_std$X_std

# Imputation 
imp_xp  <- init_XP_ER(X_std, params$Kmax_vec)

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
path_country <- file.path(path_results_vec, country)
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
  X_full  = X_vec,      # full vectorized mixed-frequency predictor matrix (with NA)
  y_q     = y_q,        # true quarterly target
  params  = params,
  dates   = dates_m,    # monthly dates
  dates_q = dates_q,    # quarterly dates
  Freq    = Freq,       # frequency classification aligned with X_vec columns
  Unb     = Unb,        # publication delays aligned with X_vec columns
  agg_m   = agg_m,      # aggregation rules for monthly predictors
  agg_q   = agg_q       # aggregation rules for quarterly predictors
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
path_country <- file.path(path_results_vec, country)
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


