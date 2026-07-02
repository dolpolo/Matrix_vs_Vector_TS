# ==============================================================================
# Matrix MF-TPRF: Mixed-Frequency Targeted Matrix Factor Nowcasting
# ==============================================================================
# Author: Davide Delfino
#
# Purpose:
#   Estimates the Matrix MF-TPRF model for nowcasting Euro Area and national
#   GDP growth using a mixed-frequency macroeconomic panel.
#
# Main steps:
#   1. Data preparation, country-specific variable selection, and tensor setup.
#   2. Standardization and missing-value imputation.
#   3. Construction of scalar or multivariate target proxy systems.
#   4. Targeted matrix factor extraction and rank selection.
#   5. U-MIDAS lag selection and full-sample estimation.
#   6. Expanding-window pseudo real-time nowcasting by monthly vintage.
#
# Proxy designs:
#   - "scalar": EA GDP proxy with sequential residual proxies.
#   - "multivariate": standardized national GDP target vector.
#
# Outputs:
#   - Full-sample estimation object saved as .rds.
#   - Pseudo real-time country and EA nowcasts saved as .rds.
#   - Optional fixed-rank appendix output for the scalar specification.
#
# Required source files:
#   - matrix.mf.tprf.utils.R
#   - matrix.mf.tprf.prep.R
#   - matrix.mf.tprf.imp.R
#   - matrix.mf.tprf.fs.R
#   - matrix.mf.tprf.R
#   - matrix.mf.tprf.now.R
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_mat")

path_results <- file.path(
  path_main,
  "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs"
)

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
# 2. FUNCTIONS
# ==============================================================================

source(file.path(path_func, "matrix.mf.tprf.utils.R"))
source(file.path(path_func, "matrix.mf.tprf.prep.R"))
source(file.path(path_func, "matrix.mf.tprf.imp.R"))
source(file.path(path_func, "matrix.mf.tprf.fs.R"))
source(file.path(path_func, "matrix.mf.tprf.R"))
source(file.path(path_func, "matrix.mf.tprf.now.R"))


# ==============================================================================
# 3. PARAMETERS
# ==============================================================================

params <- list(
  
  start_est    = as.Date("2000-04-01"),
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2026-02-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,
  covid_mask_q = TRUE,
  
  target    = "GDP",
  target_cc = "EA",
  
  proxy_mode = "scalar",   # "scalar" | "multivariate"
  
  sel_method  = "LASSO",   # "LASSO" | "corr" 
  n_m         = 20,        # 20  |  25  | 30
  n_q         = 5,         # 5   |  15  | 30
  thr_m       = 0.10,
  thr_q       = 0.85,
  thr_F_test  = 0.01,
  alpha_lasso = 1,
  
  Kmax = c(3, 6),         # Large - corr ---> Kmax = c(2, 3)  !!!ATTENTION!!!
  Zmax = 5,
  
  p_AR_max = 5,
  Lmax     = 5,
  
  Robust_F    = FALSE,
  alpha       = 0.10,
  robust_type = "NW",
  nw_lag      = 1
  
)

proxy_mode <- match.arg(
  params$proxy_mode,
  choices = c("scalar", "multivariate")
)

model_name <- paste0("matrix_", proxy_mode)
Size       <- get_size_tag(params$n_m, params$n_q)
sel        <- params$sel_method

countries <- c(
  "DE", "FR", "IT", "ES",
  "NL", "BE", "AT", "PT",
  params$target_cc
)


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

tensor <- build_tensor(
  prep      = all_countries,
  params    = params,
  var_scope = "union"
)

Y <- tensor$Y
W <- tensor$W

na_pct            <- nan_percent_Y(Y)
df_selection      <- selection_to_df(all_countries, params)
df_selection_wide <- selection_to_wide(df_selection)


# ==============================================================================
# 5. TARGETS, PREDICTORS, AND METADATA
# ==============================================================================

gdp_col <- tensor$target_col

y <- Y[, , gdp_col, drop = FALSE]

is_quarter <- rowSums(!is.na(y[, , 1])) > 0

y_q <- y[is_quarter, , 1, drop = TRUE]
y_q <- as.matrix(y_q)

if (is.null(colnames(y_q))) {
  stop("Country names are missing from the quarterly GDP matrix.")
}

if (params$target_cc %in% colnames(y_q)) {
  stop("The aggregate target must not be included in y_q.")
}

rownames(y_q) <- as.character(dimnames(Y)[[1]][is_quarter])

y_EA_q <- as.numeric(all_countries$proxy$Y_q)

if (length(y_EA_q) != nrow(y_q)) {
  stop("EA GDP and national GDP samples have different lengths.")
}

y_q_all <- cbind(y_EA_q, y_q)

colnames(y_q_all)[1] <- params$target_cc
rownames(y_q_all) <- rownames(y_q)

X   <- Y[, , -gdp_col, drop = FALSE]
W_x <- W[, , -gdp_col, drop = FALSE]

dates_m <- as.Date(dimnames(Y)[[1]])
dates_q <- dates_m[is_quarter]

N_m <- tensor$n_M

N_q_tot      <- tensor$n_Q
q_series_all <- tensor$vars[(N_m + 1):(N_m + N_q_tot)]

q_cols <- which(
  tolower(q_series_all) != tolower(params$target)
)

N_q <- length(q_cols)
N   <- N_m + N_q

agg <- cbind(
  tensor$agg_M,
  tensor$agg_Q[, q_cols, drop = FALSE]
)

Freq <- cbind(
  tensor$freq_M,
  tensor$freq_Q[, q_cols, drop = FALSE]
)

Unb <- cbind(
  tensor$unb_M,
  tensor$unb_Q[, q_cols, drop = FALSE]
)

Class <- cbind(
  tensor$ClassM,
  tensor$ClassQ[, q_cols, drop = FALSE]
)


# ==============================================================================
# 6. STANDARDIZATION, IMPUTATION, AND QUARTERLY AGGREGATION
# ==============================================================================

out_std <- standardize_mat_with_na(X)
X_std   <- out_std$X_scaled

imp_cl <- init_CL_Yu(
  X_std,
  W_x,
  params$Kmax
)

X_cl  <- imp_cl$Y_init
r_hat <- imp_cl$r

X_cl_q <- aggregate_tensor_to_quarterly(
  X_tens = X_cl,
  agg    = agg,
  N_m    = N_m,
  N_q    = N_q
)

if (nrow(X_cl_q) != nrow(y_q_all)) {
  stop("Quarterly predictor tensor and GDP targets have different time dimensions.")
}

country_order <- dimnames(X_cl_q)[[2]]

if (is.null(country_order) || length(country_order) == 0L) {
  stop("Country names are missing from the row dimension of X_cl_q.")
}

if (params$target_cc %in% country_order) {
  stop("The aggregate target must not be included in the row dimension of X_cl_q.")
}

if (!setequal(country_order, setdiff(colnames(y_q_all), params$target_cc))) {
  stop("Country rows of X_cl_q and national GDP columns of y_q_all are not aligned.")
}

y_q_all <- y_q_all[
  ,
  c(params$target_cc, country_order),
  drop = FALSE
]

rmax_targeted <- params$Kmax


# ==============================================================================
# 7. PROXY SETUP
# ==============================================================================

proxy_selection <- NULL
Z_q_grid <- NULL

X_lf_grid <- X_cl_q
X_hf_grid <- X_cl
y_q_grid  <- y_EA_q

target_proxy_window <- NULL

if (proxy_mode == "scalar") {
  
  X_cl_q_vec <- tensor_to_vector(
    X_tens = X_cl_q,
    N_m    = N_m,
    N_q    = N_q
  )
  
  proxy_selection <- select_L_autoproxy_3prf(
    X_lf = X_cl_q_vec,
    y_q  = y_EA_q,
    Zmax = params$Zmax
  )
  
  Lproxy <- proxy_selection$L_opt
  
  use_autoproxy_grid <- TRUE
  pass1_orthonormalize_grid <- TRUE
  pass1_center_grid <- TRUE
  
} else {
  
  target_proxy_window <- prepare_multivariate_proxy_window(
    X_lf = X_cl_q,
    X_hf = X_cl,
    Y_q_all = y_q_all,
    country_cols = country_order
  )
  
  Z_q_grid <- target_proxy_window$Z_q
  
  X_lf_grid <- target_proxy_window$X_lf
  X_hf_grid <- target_proxy_window$X_hf
  y_q_grid  <- y_EA_q[target_proxy_window$q_idx]
  
  Lproxy <- ncol(Z_q_grid)
  
  use_autoproxy_grid <- FALSE
  pass1_orthonormalize_grid <- FALSE
  pass1_center_grid <- FALSE
  
  cat(
    "\nMultivariate Pass 1 / grid window:\n",
    "  First common target date: ",
    as.character(dates_q[target_proxy_window$q_first]),
    "\n",
    "  Last common target date : ",
    as.character(dates_q[target_proxy_window$q_last]),
    "\n",
    "  Targeting quarters      : ",
    target_proxy_window$T_q_targeting,
    "\n\n"
  )
}

# ==============================================================================
# 8. U-MIDAS LAG SELECTION
# ==============================================================================

lag_sel <- choose_UMIDAS_grid_tensor_MF(
  X_lf = X_lf_grid,
  X_hf = X_hf_grid,
  y_q  = y_q_grid,
  
  rmax     = rmax_targeted,
  Lproxy   = Lproxy,
  Lmax     = params$Lmax,
  p_AR_max = params$p_AR_max,
  
  use_autoproxy = use_autoproxy_grid,
  y_proxy       = y_EA_q,
  Z_q           = Z_q_grid,
  
  standardize_proxy        = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  
  pass1_orthonormalize_Z = pass1_orthonormalize_grid,
  pass1_center_Z         = pass1_center_grid,
  
  ils_maxit = 100,
  ils_tol   = 1e-8
)

L_midas    <- lag_sel$best_BIC$L
p_ar       <- lag_sel$best_BIC$p_AR
r_selected <- lag_sel$r_selected
ils_iter   <- lag_sel$pass1_obj_path


# ==============================================================================
# 9. FULL-SAMPLE ESTIMATION
# ==============================================================================

fit_tensor <- Tensor_MF_TPRF(
  X_lf    = X_cl_q,
  X_hf    = X_cl,
  Y_q_all = y_q_all,
  
  proxy_name   = params$target_cc,
  proxy_mode   = proxy_mode,
  forecast_mode = "both",
  
  Lproxy  = Lproxy,
  L_midas = L_midas,
  p_AR    = p_ar,
  rmax    = rmax_targeted,
  
  standardize_proxy        = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  
  ils_maxit = 100,
  ils_tol   = 1e-8
)

r_selected_fit <- fit_tensor$r_selected


# ==============================================================================
# 10. SAVE FULL-SAMPLE RESULTS
# ==============================================================================

file_fit <- build_result_filename(
  path_out        = path_results,
  model           = model_name,
  stage           = "fit",
  Size            = Size,
  sel             = sel,
  countries       = countries,
  N_m             = N_m,
  N_q             = N_q,
  Lproxy          = fit_tensor$meta$Lproxy,
  L_midas         = L_midas,
  p_ar            = p_ar,
  r1              = r_selected_fit[1],
  r2              = r_selected_fit[2],
  robust_f        = as.integer(isTRUE(params$Robust_F)),
  covid_m         = as.integer(isTRUE(params$covid_mask_m)),
  covid_q         = as.integer(isTRUE(params$covid_mask_q)),
  ext             = "rds",
  timestamp       = FALSE,
  include_details = FALSE
)

saveRDS(
  list(
    model_id   = "T_MF_TPRF",
    model      = model_name,
    proxy_mode = proxy_mode,
    stage      = "full_sample",
    Size       = Size,
    sel        = sel,
    
    params    = params,
    countries = countries,
    dates_m   = dates_m,
    dates_q   = dates_q,
    
    target = list(
      gdp_col    = gdp_col,
      proxy_name = params$target_cc,
      y_q_all    = y_q_all
    ),
    
    proxy = list(
      mode      = proxy_mode,
      columns   = fit_tensor$proxy_columns,
      Lproxy    = fit_tensor$meta$Lproxy,
      selection = proxy_selection
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
    
    selections = list(
      table = df_selection,
      wide  = df_selection_wide,
      raw   = all_countries$sel
    ),
    
    preprocessing = list(
      X_std  = X_std,
      imp_cl = imp_cl,
      X_cl   = X_cl,
      X_cl_q = X_cl_q
    ),
    
    hyper = list(
      r_impute        = r_hat,
      rmax_targeted   = rmax_targeted,
      r_selected_grid = r_selected,
      r_selected_fit  = r_selected_fit,
      Lproxy          = fit_tensor$meta$Lproxy,
      L_midas         = L_midas,
      p_ar            = p_ar,
      ils_iter        = ils_iter
    ),
    
    na_pct = na_pct,
    fit    = fit_tensor
  ),
  file = file_fit
)

cat("\nSaved full-sample results to:\n", file_fit, "\n")


# ==============================================================================
# 11. REAL-TIME DATA
# ==============================================================================

selection_end_pre  <- params$start_eval %m-% months(1)
selection_end_post <- params$covid_end

all_countries_rt_pre <- prepare_all_countries(
  countries     = countries,
  params        = params,
  path_raw      = path_data_raw,
  path_adj      = path_data_adj,
  covid_mask_m  = params$covid_mask_m,
  covid_mask_q  = params$covid_mask_q,
  selection_end = selection_end_pre
)

all_countries_rt_post <- prepare_all_countries(
  countries     = countries,
  params        = params,
  path_raw      = path_data_raw,
  path_adj      = path_data_adj,
  covid_mask_m  = params$covid_mask_m,
  covid_mask_q  = params$covid_mask_q,
  selection_end = selection_end_post
)

tensor_rt_pre <- build_tensor(
  prep      = all_countries_rt_pre,
  params    = params,
  var_scope = "union"
)

tensor_rt_post <- build_tensor(
  prep      = all_countries_rt_post,
  params    = params,
  var_scope = "union"
)

regime_data_pre <- make_regime_data(
  tensor_obj    = tensor_rt_pre,
  label         = "pre_evaluation_selection",
  selection_end = selection_end_pre,
  params        = params
)

regime_data_post <- make_regime_data(
  tensor_obj    = tensor_rt_post,
  label         = "post_covid_selection",
  selection_end = selection_end_post,
  params        = params
)

df_selection_rt_pre <- selection_to_df(
  all_countries_rt_pre,
  params
)

df_selection_rt_post <- selection_to_df(
  all_countries_rt_post,
  params
)

df_selection_wide_rt_pre <- selection_to_wide(
  df_selection_rt_pre
)

df_selection_wide_rt_post <- selection_to_wide(
  df_selection_rt_post
)

na_pct_rt_pre  <- nan_percent_Y(regime_data_pre$X_full)
na_pct_rt_post <- nan_percent_Y(regime_data_post$X_full)


# ==============================================================================
# 12. PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

roll_tensor <- pseudo_realtime_tensor_mf_tprf_fixed_hyper(
  X_full       = regime_data_pre$X_full,
  W_full       = regime_data_pre$W_full,
  Unb          = regime_data_pre$Unb,
  Y_q_all_full = y_q_all,
  
  proxy_name   = params$target_cc,
  proxy_mode   = proxy_mode,
  forecast_mode = "both",
  
  params  = params,
  dates_m = dates_m,
  dates_q = dates_q,
  
  agg = regime_data_pre$agg,
  N_m = regime_data_pre$N_m,
  N_q = regime_data_pre$N_q,
  
  regime_data_pre  = regime_data_pre,
  regime_data_post = regime_data_post,
  post_start_date  = params$covid_end %m+% months(1),
  
  user_hyper_pre = list(
    r_impute   = NULL,
    Lproxy     = NULL,
    p_AR       = NULL,
    L_midas    = NULL,
    r_targeted = NULL
  ),
  
  user_hyper_post = list(
    r_impute   = NULL,
    Lproxy     = NULL,
    p_AR       = NULL,
    L_midas    = NULL,
    r_targeted = NULL
  )
)


# ==============================================================================
# 13. REAL-TIME OUTPUT TABLES
# ==============================================================================

df_M1 <- list_to_df_country(roll_tensor$M1, "M1")
df_M2 <- list_to_df_country(roll_tensor$M2, "M2")
df_M3 <- list_to_df_country(roll_tensor$M3, "M3")

df_pseudoRT_all <- bind_rows(df_M1, df_M2, df_M3) |>
  mutate(
    month_in_quarter = factor(
      month_in_quarter,
      levels = c("M1", "M2", "M3")
    ),
    proxy_mode = proxy_mode
  ) |>
  arrange(country, date, month_in_quarter)

df_A_M1 <- data.frame(
  date             = as.Date(names(roll_tensor$aggregate$M1)),
  country          = params$target_cc,
  nowcast          = as.numeric(roll_tensor$aggregate$M1),
  month_in_quarter = "M1"
)

df_A_M2 <- data.frame(
  date             = as.Date(names(roll_tensor$aggregate$M2)),
  country          = params$target_cc,
  nowcast          = as.numeric(roll_tensor$aggregate$M2),
  month_in_quarter = "M2"
)

df_A_M3 <- data.frame(
  date             = as.Date(names(roll_tensor$aggregate$M3)),
  country          = params$target_cc,
  nowcast          = as.numeric(roll_tensor$aggregate$M3),
  month_in_quarter = "M3"
)

df_pseudoRT_A <- bind_rows(df_A_M1, df_A_M2, df_A_M3) |>
  mutate(
    month_in_quarter = factor(
      month_in_quarter,
      levels = c("M1", "M2", "M3")
    ),
    proxy_mode = proxy_mode
  ) |>
  arrange(date, month_in_quarter)


# ==============================================================================
# 14. SAVE REAL-TIME RESULTS
# ==============================================================================

file_rt <- build_result_filename(
  path_out        = path_results,
  model           = model_name,
  stage           = "rt",
  Size            = Size,
  sel             = sel,
  countries       = country_order,
  N_m             = regime_data_pre$N_m,
  N_q             = regime_data_pre$N_q,
  Lproxy          = roll_tensor$hyper_pre$Lproxy,
  L_midas         = roll_tensor$hyper_pre$L_midas,
  p_ar            = roll_tensor$hyper_pre$p_AR,
  r1              = NA,
  r2              = NA,
  robust_f        = as.integer(isTRUE(params$Robust_F)),
  covid_m         = as.integer(isTRUE(params$covid_mask_m)),
  covid_q         = as.integer(isTRUE(params$covid_mask_q)),
  ext             = "rds",
  timestamp       = FALSE,
  include_details = FALSE
)

saveRDS(
  list(
    model_id   = "T_MF_TPRF",
    model      = model_name,
    proxy_mode = proxy_mode,
    stage      = "pseudo_realtime",
    Size       = Size,
    sel        = sel,
    
    params     = params,
    proxy_name = params$target_cc,
    countries  = country_order,
    dates_m    = dates_m,
    dates_q    = dates_q,
    Y_q_all    = y_q_all,
    
    selection_protocol = list(
      full_sample = list(
        selection_end = params$end_eval,
        table         = df_selection,
        wide          = df_selection_wide,
        raw           = all_countries$sel
      ),
      pre = list(
        selection_end = selection_end_pre,
        table         = df_selection_rt_pre,
        wide          = df_selection_wide_rt_pre,
        raw           = all_countries_rt_pre$sel
      ),
      post = list(
        selection_end = selection_end_post,
        table         = df_selection_rt_post,
        wide          = df_selection_wide_rt_post,
        raw           = all_countries_rt_post$sel
      )
    ),
    
    metadata = list(
      full_sample = list(
        N_m    = N_m,
        N_q    = N_q,
        N      = N,
        K      = dim(X)[3],
        agg    = agg,
        Unb    = Unb,
        vars   = dimnames(X)[[3]],
        na_pct = na_pct
      ),
      pre = list(
        N_m    = regime_data_pre$N_m,
        N_q    = regime_data_pre$N_q,
        N      = regime_data_pre$N,
        K      = dim(regime_data_pre$X_full)[3],
        agg    = regime_data_pre$agg,
        Unb    = regime_data_pre$Unb,
        vars   = regime_data_pre$vars,
        na_pct = na_pct_rt_pre
      ),
      post = list(
        N_m    = regime_data_post$N_m,
        N_q    = regime_data_post$N_q,
        N      = regime_data_post$N,
        K      = dim(regime_data_post$X_full)[3],
        agg    = regime_data_post$agg,
        Unb    = regime_data_post$Unb,
        vars   = regime_data_post$vars,
        na_pct = na_pct_rt_post
      )
    ),
    
    hyper = list(
      pre  = roll_tensor$hyper_pre,
      post = roll_tensor$hyper_post
    ),
    
    regime_info = roll_tensor$regime_info,
    
    pseudo_rt_raw       = roll_tensor,
    pseudo_rt_M1        = df_M1,
    pseudo_rt_M2        = df_M2,
    pseudo_rt_M3        = df_M3,
    pseudo_rt_all       = df_pseudoRT_all,
    pseudo_rt_aggregate = df_pseudoRT_A
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

