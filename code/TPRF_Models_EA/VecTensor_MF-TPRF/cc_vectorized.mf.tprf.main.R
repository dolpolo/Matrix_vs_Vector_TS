# ==============================================================================
# Vectorized-from-Tensor MF-TPRF Main Script (FOR A SINGLE COUNTRY)
# ------------------------------------------------------------------------------
# This script:
#   1. Prepares country data and builds the tensor
#   2. Vectorizes the tensor into a mixed-frequency panel
#   3. Runs MF-TPRF for all countries on tensor-vectorized data
#   4. Builds cross-country plots and RMSFE tables
#   5. Runs a detailed single-country estimation
#   6. Runs a pseudo real-time exercise for one country
#   7. Saves outputs in a consistent format
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")

path_func_mat <- file.path(path_main, "functions/functions_mat")
path_func_vec <- file.path(path_main, "functions/functions_vec")

path_results_vec <- file.path(path_main, "TPRF_Models_EA/VecTensor_MF-TPRF/results/outputs_vec")
path_graph_vec   <- file.path(path_main, "TPRF_Models_EA/VecTensor_MF-TPRF/results/graph_vec")

dir.create(path_results_vec, recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_vec,   recursive = TRUE, showWarnings = FALSE)

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
# 2. SOURCE MATRIX FUNCTIONS
# Build tensor before loading vector functions with overlapping names
# ==============================================================================

source(file.path(path_func_mat, "matrix.mf.tprf.utils.R"))
source(file.path(path_func_mat, "matrix.mf.tprf.prep.R"))

# ==============================================================================
# 3. LOCAL HELPERS
# ==============================================================================

build_run_tag <- function(params) {
  paste0(
    "eval-", format(params$start_eval, "%Y%m"), "_", format(params$end_eval, "%Y%m"),
    "_sel-", params$sel_method,
    "_Nm-", params$n_m,
    "_Nq-", params$n_q,
    "_KmaxVec-", params$Kmax_vec,
    "_Zmax-", params$Zmax,
    "_Lmax-", params$Lmax,
    "_pARmax-", params$p_AR_max,
    "_CovidM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-", as.integer(isTRUE(params$covid_mask_q))
  )
}

list_to_df_nowcast <- function(lst, tag) {
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
# 4. USER PARAMETERS
# ==============================================================================

params <- list(
  start_est    = as.Date("2000-04-01"),
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-09-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,
  covid_mask_q = TRUE,
  target       = "GDP",
  target_cc    = "EA",
  
  sel_method   = "corr",
  n_m          = 30,
  n_q          = 30,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,
  
  Kmax         = c(3, 4),
  Kmax_vec     = 15,
  Zmax         = 5,
  
  p_AR_max     = 5,
  Lmax         = 5,
  Robust_F     = FALSE,
  alpha        = 0.10,
  robust_type  = "NW",
  nw_lag       = 1
)

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
country   <- "IT"

tag_run <- build_run_tag(params)
model_name <- "vectensor_country"
sel        <- params$sel_method
Size       <- get_size_tag(params$n_m, params$n_q)

# ==============================================================================
# 5. PREPARE DATA AND BUILD TENSOR
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
data   <- tensor_to_vector(X_tens = tensor$Y, N_m = tensor$n_M, N_q = tensor$n_Q)

# Optional diagnostic
nan_percent_Y(tensor$Y)

# ==============================================================================
# 6. SOURCE VECTOR FUNCTIONS
# ==============================================================================

source(file.path(path_func_vec, "mf.tprf.utils.R"))
source(file.path(path_func_vec, "mf.tprf.prep.R"))
source(file.path(path_func_vec, "mf.tprf.imp.R"))
source(file.path(path_func_vec, "mf.tprf.fs.R"))
source(file.path(path_func_vec, "mf.tprf.R"))
source(file.path(path_func_vec, "mf.tprf.now.R"))
source(file.path(path_func_vec, "mf.tprf.all_cc_results.R"))


# ==============================================================================
# 7. SINGLE-COUNTRY INPUTS
# ==============================================================================

country_inputs <- prepare_country_inputs_from_tensor_vectorized(
  tensor  = tensor,
  data    = data,
  country = country,
  params  = params
)

# ==============================================================================
# 8. SINGLE-COUNTRY FULL-SAMPLE ESTIMATION
# ==============================================================================

est_single <- estimate_country_mf_tprf(
  country_inputs = country_inputs,
  params         = params
)

fit_vec_from_tensor <- est_single$fit

# ==============================================================================
# 9. SAVE SINGLE-COUNTRY FULL-SAMPLE RESULTS
# ==============================================================================

path_country <- file.path(path_results_vec, country)
dir.create(path_country, recursive = TRUE, showWarnings = FALSE)

file_fit <- build_result_filename(
  path_out         = path_country,
  model            = model_name,
  stage            = "fit",
  Size             = Size,
  sel              = sel,
  countries        = country,
  N_m              = country_inputs$N_m,
  N_q              = country_inputs$N_q,
  Lproxy           = est_single$hyper$Lproxy,
  L_midas          = est_single$hyper$L_midas,
  p_ar             = est_single$hyper$p_ar,
  r1               = if (!is.null(est_single$hyper$r_hat)) est_single$hyper$r_hat else NA,
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
    model_id = "MF_TPRF_VECFROMTENSOR",
    model    = model_name,
    stage    = "full_sample",
    Size     = Size,
    sel      = sel,
    country  = country,
    params   = params,
    dates_m  = country_inputs$dates_m,
    dates_q  = country_inputs$dates_q,
    y_q      = country_inputs$y_q,
    metadata = list(
      N_m          = country_inputs$N_m,
      N_q          = country_inputs$N_q,
      N            = country_inputs$N,
      agg_m        = country_inputs$agg_m,
      agg_q        = country_inputs$agg_q,
      agg          = country_inputs$agg,
      freq         = country_inputs$freq,
      unb          = country_inputs$unb,
      series_names = country_inputs$series_names
    ),
    preprocessing = est_single$preprocessing,
    hyper         = est_single$hyper,
    fit           = est_single$fit
  ),
  file = file_fit
)

cat("\nSaved full-sample results to:\n", file_fit, "\n")

# ==============================================================================
# 10. SINGLE-COUNTRY PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

rt_single <- run_country_mf_tprf_rt(
  country_inputs = country_inputs,
  params         = params
)

# ==============================================================================
# 11. SAVE SINGLE-COUNTRY PSEUDO REAL-TIME RESULTS
# ==============================================================================

hyper_pre_rt <- if (!is.null(rt_single$raw$hyper_pre)) {
  rt_single$raw$hyper_pre
} else {
  list(Lproxy = NA, L_midas = NA, p_AR = NA)
}

file_rt <- build_result_filename(
  path_out         = path_country,
  model            = model_name,
  stage            = "rt",
  Size             = Size,
  sel              = sel,
  countries        = country,
  N_m              = country_inputs$N_m,
  N_q              = country_inputs$N_q,
  Lproxy           = hyper_pre_rt$Lproxy,
  L_midas          = hyper_pre_rt$L_midas,
  p_ar             = hyper_pre_rt$p_AR,
  r1               = if (!is.null(est_single$hyper$r_hat)) est_single$hyper$r_hat else NA,
  r2               = NA,
  robust_f         = as.integer(isTRUE(params$Robust_F)),
  covid_m          = as.integer(isTRUE(params$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params$covid_mask_q)),
  ext              = "rds",
  timestamp        = TRUE,
  include_details  = TRUE
)

if (file.exists(file_rt)) file.remove(file_rt)

saveRDS(
  list(
    model_id = "MF_TPRF_VECFROMTENSOR",
    model    = model_name,
    stage    = "pseudo_realtime",
    Size     = Size,
    sel      = sel,
    country  = country,
    params   = params,
    dates_m  = country_inputs$dates_m,
    dates_q  = country_inputs$dates_q,
    y_q      = country_inputs$y_q,
    pseudo_realtime_raw = rt_single$raw,
    pseudo_realtime_M1  = rt_single$M1,
    pseudo_realtime_M2  = rt_single$M2,
    pseudo_realtime_M3  = rt_single$M3,
    pseudo_realtime_all = rt_single$all
  ),
  file = file_rt
)

cat("\nSaved pseudo real-time results to:\n", file_rt, "\n")
