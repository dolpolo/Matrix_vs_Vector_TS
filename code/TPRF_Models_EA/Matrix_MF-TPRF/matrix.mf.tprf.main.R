# ==============================================================================
# Matrix MF-TPRF Main Script
# Mixed-Frequency Nowcasting for Euro Area GDP
# ------------------------------------------------------------------------------
# This script:
#   1. Prepares country data and builds the panel tensor
#   2. Extracts target, predictors, and metadata
#   3. Standardizes and imputes missing values
#   4. Selects proxies, factors, and lag orders
#   5. Estimates the Matrix MF-TPRF on the full sample
#   6. Saves results in a comparable and reproducible format
#   7. Runs the pseudo real-time nowcasting exercise
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_mat")

path_results  <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")

dir.create(path_results, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(tidyverse)
library(lubridate)
library(abind)
library(zoo)
library(MASS)
library(tseries)
library(fBasics)
library(vars)
library(glmnet)
library(plsdof)
library(sandwich)
library(lmtest)
library(car)
library(readxl)
library(conflicted)

conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)

# ==============================================================================
# 2. SOURCE FUNCTIONS
# ==============================================================================

source(file.path(path_func, "matrix.mf.tprf.utils.R"))
source(file.path(path_func, "matrix.mf.tprf.prep.R"))
source(file.path(path_func, "matrix.mf.tprf.imp.R"))
source(file.path(path_func, "matrix.mf.tprf.fs.R"))
source(file.path(path_func, "matrix.mf.tprf.R"))
source(file.path(path_func, "matrix.mf.tprf.now.R"))

# ==============================================================================
# 3. USER PARAMETERS
# ==============================================================================

params <- list(
  # Sample
  start_est    = as.Date("2000-04-01"),
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2026-02-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,
  covid_mask_q = TRUE,
  
  # Target
  target    = "GDP",
  target_cc = "EA",
  
  # Variable selection
  sel_method  = "LASSO",      # "corr" | "corr_threshold" | "F-Test" | "LASSO"
  n_m         = 20,
  n_q         = 5,
  thr_m       = 0.10,
  thr_q       = 0.85,
  thr_F_test  = 0.01,
  alpha_lasso = 1,
  
  # Factor / proxy selection
  Kmax     = c(3, 6),
  Zmax     = 5,
  
  # U-MIDAS / AR
  p_AR_max = 5,
  Lmax     = 5,
  
  # Step-1 inference
  Robust_F    = FALSE,
  alpha       = 0.10,
  robust_type = "NW",
  nw_lag      = 1
)

# Tags used in all saved outputs
model_name <- "matrix"
Size       <- get_size_tag(params$n_m, params$n_q)
sel        <- params$sel_method

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")

# ==============================================================================
# 4. DATA PREPARATION
# ==============================================================================

all_countries <- prepare_all_countries(
  countries    = countries,
  params       = params,
  path_raw     = path_data_raw,
  path_adj     = path_data_adj,
  covid_mask_m = params$covid_mask_m,
  covid_mask_q = params$covid_mask_q
)

tensor <- build_tensor(all_countries, params, var_scope = "union")
Y <- tensor$Y
W <- tensor$W

na_pct            <- nan_percent_Y(Y)
df_selection      <- selection_to_df(all_countries, params)
df_selection_wide <- selection_to_wide(df_selection)
# ==============================================================================
# 5. TARGET, PREDICTORS, AND METADATA
# ==============================================================================

gdp_col <- tensor$target_col

# Target
y       <- Y[, , gdp_col, drop = FALSE]
y_q     <- y[rowSums(!is.na(y[, , 1])) > 0, , , drop = TRUE]
y_EA_q  <- all_countries$proxy$Y_q
y_q_all <- cbind(y_EA_q, y_q)

# Predictors
X   <- Y[, , -gdp_col, drop = FALSE]
W_x <- W[, , -gdp_col, drop = FALSE]

# Dates
dates_m <- as.Date(dimnames(Y)[[1]])
dates_q <- dates_m[rowSums(!is.na(y[, , 1])) > 0]

# Dimensions
N_m <- tensor$n_M

N_q_tot      <- tensor$n_Q
q_series_all <- tensor$vars[(N_m + 1):(N_m + N_q_tot)]
q_cols       <- which(tolower(q_series_all) != tolower(params$target))
N_q          <- length(q_cols)
N            <- N_m + N_q

# Metadata
agg   <- cbind(tensor$agg_M,   tensor$agg_Q[,  q_cols, drop = FALSE])
Freq  <- cbind(tensor$freq_M,  tensor$freq_Q[, q_cols, drop = FALSE])
Unb   <- cbind(tensor$unb_M,   tensor$unb_Q[,  q_cols, drop = FALSE])
Class <- cbind(tensor$ClassM,  tensor$ClassQ[, q_cols, drop = FALSE])

# ==============================================================================
# 6. STANDARDIZATION, IMPUTATION, AND INITIAL FACTOR SELECTION
# ==============================================================================

out_std <- standardize_mat_with_na(X)
X_std   <- out_std$X_scaled

imp_cl <- init_CL_Yu(X_std, W_x, params$Kmax)

X_cl  <- imp_cl$Y_init
r_hat <- imp_cl$r

X_cl_q <- aggregate_tensor_to_quarterly(
  X_tens = X_cl,
  agg    = agg,
  N_m    = N_m,
  N_q    = N_q
)

# ==============================================================================
# 7. PROXY SELECTION
# ==============================================================================

X_cl_q_vec <- tensor_to_vector(
  X_tens = X_cl_q,
  N_m    = N_m,
  N_q    = N_q
)

pls_object <- select_L_autoproxy_3prf(
  X_cl_q_vec,
  y_EA_q,
  Zmax = params$Zmax
)

Lproxy <- pls_object$L_opt

# ==============================================================================
# 8. U-MIDAS AND AR LAG SELECTION
# ==============================================================================

lag_sel <- choose_UMIDAS_grid_tensor_MF(
  X_lf                      = X_cl_q,
  X_hf                      = X_cl,
  y_q                       = y_EA_q,
  r                         = r_hat,
  Lproxy                    = Lproxy,
  Lmax                      = params$Lmax,
  p_AR_max                  = params$p_AR_max,
  use_autoproxy             = TRUE,
  y_proxy                   = y_EA_q,
  standardize_proxy         = TRUE,
  orthonormalize_each_iter  = TRUE,
  orthonormalize_final_Z    = TRUE,
  ils_maxit                 = 100,
  ils_tol                   = 1e-8
)

L_midas    <- lag_sel$best_BIC$L
p_ar       <- lag_sel$best_BIC$p_AR
r_selected <- lag_sel$r_selected
ils_iter   <- lag_sel$pass1_obj_path

# ==============================================================================
# 9. FULL-SAMPLE ESTIMATION
# ==============================================================================

fit_tensor <- Tensor_MF_TPRF(
  X_lf                     = X_cl_q,
  X_hf                     = X_cl,
  Y_q_all                  = y_q_all,
  proxy_name               = "EA",
  Lproxy                   = Lproxy,
  L_midas                  = L_midas,
  p_AR                     = p_ar,
  r                        = r_hat,
  standardize_proxy        = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  ils_maxit                = 100,
  ils_tol                  = 1e-8
)

# ==============================================================================
# 10. SAVE FULL-SAMPLE RESULTS
# ==============================================================================

r1 <- r_selected[1]
r2 <- r_selected[2]

file_fit <- build_result_filename(
  path_out         = path_results,
  model            = model_name,
  stage            = "fit",
  Size             = Size,
  sel              = sel,
  countries        = countries,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = Lproxy,
  L_midas          = L_midas,
  p_ar             = p_ar,
  r1               = r1,
  r2               = r2,
  robust_f         = as.integer(isTRUE(params$Robust_F)),
  covid_m          = as.integer(isTRUE(params$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params$covid_mask_q)),
  ext              = "rds",
  timestamp        = TRUE,
  include_details  = TRUE
)
saveRDS(
  list(
    model_id = "T_MF_TPRF",
    model    = model_name,
    stage    = "full_sample",
    Size     = Size,
    sel      = sel,
    
    params     = params,
    countries  = countries,
    na_pct     = na_pct,
    selections = df_selection,
    selections_wide = df_selection_wide,
    sel_raw    = all_countries$sel,
    
    dates_m = dates_m,
    dates_q = dates_q,
    
    tensor = tensor,
    
    target = list(
      gdp_col    = gdp_col,
      y_q_all    = y_q_all,
      proxy_name = "EA"
    ),
    
    metadata = list(
      N_m   = N_m,
      N_q   = N_q,
      N     = N,
      agg   = agg,
      Freq  = Freq,
      Unb   = Unb,
      Class = Class
    ),
    
    preprocessing = list(
      out_std = out_std,
      X_std   = X_std,
      imp_cl  = imp_cl,
      X_cl    = X_cl,
      X_cl_q  = X_cl_q
    ),
    
    hyper = list(
      r_hat      = r_hat,
      r_selected = r_selected,
      Lproxy     = Lproxy,
      L_midas    = L_midas,
      p_ar       = p_ar,
      ils_iter   = ils_iter
    ),
    
    fit = fit_tensor
  ),
  file = file_fit
)

cat("\nSaved full-sample results to:\n", file_fit, "\n")

# ==============================================================================
# 11. PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

roll_tensor <- pseudo_realtime_tensor_mf_tprf_fixed_hyper(
  X_full          = X,
  W_full          = W_x,
  Unb             = Unb,
  Y_q_all_full    = y_q_all,
  proxy_name      = "EA",
  params          = params,
  dates_m         = dates_m,
  dates_q         = dates_q,
  agg             = agg,
  N_m             = N_m,
  N_q             = N_q,
  user_hyper_pre  = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL, r_targeted = NULL),
  user_hyper_post = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL, r_targeted = NULL)
)

df_M1 <- list_to_df_country(roll_tensor$M1, "M1")
df_M2 <- list_to_df_country(roll_tensor$M2, "M2")
df_M3 <- list_to_df_country(roll_tensor$M3, "M3")

df_pseudoRT_all <- bind_rows(df_M1, df_M2, df_M3) |>
  mutate(month_in_quarter = factor(month_in_quarter, levels = c("M1", "M2", "M3"))) |>
  arrange(country, date, month_in_quarter)

# ==============================================================================
# 12. SAVE PSEUDO REAL-TIME RESULTS
# ==============================================================================

countries_eval <- setdiff(colnames(y_q_all), "EA")

file_rt <- build_result_filename(
  path_out         = path_results,
  model            = model_name,
  stage            = "rt",
  Size             = Size,
  sel              = sel,
  countries        = countries_eval,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = roll_tensor$hyper_pre$Lproxy,
  L_midas          = roll_tensor$hyper_pre$L_midas,
  p_ar             = roll_tensor$hyper_pre$p_AR,
  r1               = NA,
  r2               = NA,
  robust_f         = as.integer(isTRUE(params$Robust_F)),
  covid_m          = as.integer(isTRUE(params$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params$covid_mask_q)),
  ext              = "rds",
  timestamp        = TRUE,
  include_details  = TRUE
)

saveRDS(
  list(
    model_id   = "T_MF_TPRF",
    model      = model_name,
    stage      = "pseudo_realtime",
    Size       = Size,
    sel        = sel,
    na_pct     = na_pct,
    selections = df_selection,
    selections_wide = df_selection_wide,
    sel_raw    = all_countries$sel,
    
    params      = params,
    proxy_name  = "EA",
    countries   = countries_eval,
    
    dates_m = dates_m,
    dates_q = dates_q,
    Y_q_all = y_q_all,
    
    metadata = list(
      N_m = N_m,
      N_q = N_q,
      N   = N,
      agg = agg,
      Unb = Unb
    ),
    
    hyper = list(
      pre  = roll_tensor$hyper_pre,
      post = roll_tensor$hyper_post
    ),
    
    pseudo_rt_raw = roll_tensor,
    pseudo_rt_M1  = df_M1,
    pseudo_rt_M2  = df_M2,
    pseudo_rt_M3  = df_M3,
    pseudo_rt_all = df_pseudoRT_all
  ),
  file = file_rt
)

cat("\nSaved pseudo real-time results to:\n", file_rt, "\n")



# ==============================================================================
# 9B. APPENDIX: FIXED-RANK ESTIMATION
# ==============================================================================

fixed_r <- c(2, 3)

fit_tensor_fixed <- Tensor_MF_TPRF_fixed(
  X_lf                     = X_cl_q,
  X_hf                     = X_cl,
  Y_q_all                  = y_q_all,
  proxy_name               = "EA",
  Lproxy                   = Lproxy,
  L_midas                  = L_midas,
  p_AR                     = p_ar,
  fixed_r                  = fixed_r,
  standardize_proxy        = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  ils_maxit                = 100,
  ils_tol                  = 1e-8
)

r1_fixed <- fit_tensor_fixed$r_selected[1]
r2_fixed <- fit_tensor_fixed$r_selected[2]

file_fixed <- file.path(
  path_results,
  paste0(
    "fit_fixed_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    "_r1-", r1_fixed,
    "_r2-", r2_fixed,
    ".rds"
  )
)

dir.create(dirname(file_fixed), recursive = TRUE, showWarnings = FALSE)
saveRDS(fit_tensor_fixed, file = file_fixed)

cat("\nSaved fixed-rank results to:\n", file_fixed, "\n")

