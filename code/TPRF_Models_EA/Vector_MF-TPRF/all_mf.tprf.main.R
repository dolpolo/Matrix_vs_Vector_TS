# ==============================================================================
# Vector MF-TPRF Main Script
# Mixed-Frequency Nowcasting for Euro Area GDP
# ------------------------------------------------------------------------------
# This script:
#   1. Prepares country-specific datasets
#   2. Runs the MF-TPRF for all countries
#   3. Builds cross-country plots and RMSFE tables
#   4. Runs a detailed single-country estimation
#   5. Runs a pseudo real-time exercise for one country
#   6. Saves all outputs in a consistent format
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_vec")

path_results  <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs")
path_graph    <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/graph")

dir.create(path_results, recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph,   recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. PACKAGES
# ==============================================================================

library(tidyverse)
library(lubridate)
library(abind)
library(zoo)
library(MASS)
library(tseries)
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

source(file.path(path_func, "mf.tprf.utils.R"))
source(file.path(path_func, "mf.tprf.prep.R"))
source(file.path(path_func, "mf.tprf.imp.R"))
source(file.path(path_func, "mf.tprf.fs.R"))
source(file.path(path_func, "mf.tprf.R"))
source(file.path(path_func, "mf.tprf.now.R"))
source(file.path(path_func, "mf.tprf.all_cc_results.R"))

# ==============================================================================
# 3. LOCAL HELPERS
# ==============================================================================

build_run_tag <- function(params) {
  paste0(
    "eval-", format(params$start_eval, "%Y%m"), "_", format(params$end_eval, "%Y%m"),
    "_sel-", params$sel_method,
    "_Nm-", params$n_m,
    "_Nq-", params$n_q,
    "_Kmax-", params$Kmax,
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
  
  sel_method   = "corr",
  n_m          = 30,
  n_q          = 30,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,
  
  Kmax         = 10,
  Zmax         = 5,
  
  p_AR_max     = 5,
  Lmax         = 5,
  
  Robust_F     = FALSE,
  alpha        = 0.10,
  robust_type  = "NW",
  nw_lag       = 1
)

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")
country   <- "IT"

tag_run <- build_run_tag(params)

# Tags used in all saved outputs
model_name <- "vector"
Size       <- get_size_tag(params$n_m, params$n_q)
sel        <- params$sel_method

# ==============================================================================
# 5. PREPARE DATA FOR ALL COUNTRIES
# ==============================================================================

all_countries <- prepare_all_countries(
  countries    = countries,
  params       = params,
  path_raw     = path_data_raw,
  path_adj     = path_data_adj,
  covid_mask_m = params$covid_mask_m,
  covid_mask_q = params$covid_mask_q
)

# ==============================================================================
# 6. CROSS-COUNTRY PIPELINE
# ==============================================================================

results_all <- run_all_countries_mf_tprf(
  countries     = countries,
  all_countries = all_countries,
  params        = params,
  path_results  = path_results
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
# 7. CROSS-COUNTRY OUTPUTS
# ==============================================================================

print(cross_out$plot_nowcast_facet)
print(cross_out$plot_rt_facet)

cat(cross_out$latex_tab_insample_all)
cat(cross_out$latex_tab_rt_all)

missing_by_country <- data.frame(
  country   = countries,
  n_missing = sapply(countries, function(cc) sum(is.na(all_countries$data[[cc]]$Data))),
  n_total   = sapply(countries, function(cc) length(all_countries$data[[cc]]$Data))
) %>%
  mutate(pct_missing = 100 * n_missing / n_total) %>%
  arrange(desc(pct_missing))

print(missing_by_country)

# ==============================================================================
# 8. SAVE CROSS-COUNTRY OUTPUTS
# ==============================================================================

path_graph_all <- file.path(path_graph, "ALL_COUNTRIES")
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
  path_out         = path_results,
  model            = model_name,
  stage            = "summary",
  Size             = Size,
  sel              = sel,
  countries        = countries,
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
    model_id               = "MF_TPRF",
    model                  = model_name,
    stage                  = "cross_country",
    Size                   = Size,
    sel                    = sel,
    params                 = params,
    countries              = countries,
    tag_run                = tag_run,
    missing_by_country     = missing_by_country,
    
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
cat("\nSaved cross-country outputs to:\n")
cat(file_graph_full, "\n")
cat(file_graph_rt, "\n")
cat(file_tex_insample, "\n")
cat(file_tex_rt, "\n")
