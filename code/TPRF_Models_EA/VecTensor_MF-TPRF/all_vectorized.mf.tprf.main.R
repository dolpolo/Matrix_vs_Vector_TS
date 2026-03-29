# ==============================================================================
# Vectorized-from-Tensor MF-TPRF Main Script
# ------------------------------------------------------------------------------
# This script:
#   1. Prepares country data and builds the tensor
#   2. Vectorizes the tensor into a mixed-frequency panel
#   3. Runs MF-TPRF for all countries on tensor-vectorized data
#   4. Builds cross-country plots and RMSFE tables
#   5. Saves outputs in a consistent format
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

tag_run    <- build_run_tag(params)
model_name <- "vectensor"
Size       <- get_size_tag(params$n_m, params$n_q)
sel        <- params$sel_method

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
# 7. CROSS-COUNTRY PIPELINE
# ==============================================================================

results_all <- run_all_countries_mf_tprf_from_tensor(
  countries    = countries,
  tensor       = tensor,
  data         = data,
  params       = params,
  path_results = path_results_vec
)

summary_all <- lapply(names(results_all), function(cc) {
  res <- results_all[[cc]]
  
  summarize_mf_tprf_country(
    MF_TPRF_res = list(
      MF_TPRF = res$full_sample$fit
    ),
    RT_res = list(
      pseudo_realtime_all = res$pseudo_realtime$all
    ),
    country = cc,
    dates_m = res$inputs$dates_m,
    dates_q = res$inputs$dates_q,
    y_q     = res$inputs$y_q,
    params  = res$params
  )
})
names(summary_all) <- names(results_all)

cross_out <- build_cross_country_outputs(
  summary_all = summary_all,
  params      = params
)

# ==============================================================================
# 8. CROSS-COUNTRY OUTPUTS
# ==============================================================================

print(cross_out$plot_nowcast_facet)
print(cross_out$plot_rt_facet)

cat(cross_out$latex_tab_insample_all)
cat(cross_out$latex_tab_rt_all)

# ==============================================================================
# 9. SAVE CROSS-COUNTRY OUTPUTS
# ==============================================================================

path_graph_all <- file.path(path_graph_vec, "ALL_COUNTRIES")
dir.create(path_graph_all, recursive = TRUE, showWarnings = FALSE)

file_graph_full <- file.path(
  path_graph_all,
  paste0(
    "plot_full_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    "_", tag_run,
    ".png"
  )
)

file_graph_rt <- file.path(
  path_graph_all,
  paste0(
    "plot_rt_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    "_", tag_run,
    ".png"
  )
)

file_tex_insample <- file.path(
  path_graph_all,
  paste0(
    "tab_insample_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    ".tex"
  )
)

file_tex_rt <- file.path(
  path_graph_all,
  paste0(
    "tab_rt_",
    model_name,
    "_Size-", Size,
    "_sel-", sel,
    ".tex"
  )
)

ggsave(file_graph_full, cross_out$plot_nowcast_facet, width = 14, height = 8, dpi = 300)
ggsave(file_graph_rt,   cross_out$plot_rt_facet,      width = 14, height = 8, dpi = 300)

writeLines(cross_out$latex_tab_insample_all, con = file_tex_insample)
writeLines(cross_out$latex_tab_rt_all,       con = file_tex_rt)

file_summary_cross <- build_result_filename(
  path_out         = path_results_vec,
  model            = model_name,
  stage            = "summary",
  Size             = Size,
  sel              = sel,
  countries        = names(results_all),
  N_m              = params$n_m,
  N_q              = params$n_q,
  Lproxy           = NA,
  L_midas          = NA,
  p_ar             = NA,
  r1               = NA,
  r2               = NA,
  robust_f         = as.integer(isTRUE(params$Robust_F)),
  covid_m          = as.integer(isTRUE(params$covid_mask_m)),
  covid_q          = as.integer(isTRUE(params$covid_mask_q)),
  ext              = "rds",
  timestamp        = FALSE,
  include_details  = TRUE
)

if (file.exists(file_summary_cross)) file.remove(file_summary_cross)

saveRDS(
  list(
    model_id               = "MF_TPRF_VECFROMTENSOR",
    model                  = model_name,
    stage                  = "cross_country",
    Size                   = Size,
    sel                    = sel,
    params                 = params,
    countries              = names(results_all),
    tag_run                = tag_run,
    
    df_now_full_all        = cross_out$df_now_full_all,
    df_quarterly_all       = cross_out$df_quarterly_all,
    plot_nowcast_facet     = cross_out$plot_nowcast_facet,
    tab_insample_all       = cross_out$tab_insample_all,
    latex_tab_insample_all = cross_out$latex_tab_insample_all,
    
    df_rt_all              = cross_out$df_rt_all,
    df_yq_eval_all         = cross_out$df_yq_eval_all,
    plot_rt_facet          = cross_out$plot_rt_facet,
    tab_rt_all             = cross_out$tab_rt_all,
    latex_tab_rt_all       = cross_out$latex_tab_rt_all,
    
    file_graph_full        = file_graph_full,
    file_graph_rt          = file_graph_rt,
    file_tex_insample      = file_tex_insample,
    file_tex_rt            = file_tex_rt
  ),
  file = file_summary_cross
)

cat("\nSaved cross-country summary to:\n", file_summary_cross, "\n")
cat("\nCross-country outputs saved to:\n")
cat(file_graph_full, "\n")
cat(file_graph_rt, "\n")
cat(file_tex_insample, "\n")
cat(file_tex_rt, "\n")