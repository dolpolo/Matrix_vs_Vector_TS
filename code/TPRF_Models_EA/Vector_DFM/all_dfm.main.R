# ==============================================================================
# Vector DFM Main Script
# DFM-EM for mixed-frequency GDP nowcasting
# Versions:
#   idio_spec = "auto" -> q selected automatically
#   idio_spec = "q0"   -> q fixed to 0
# ==============================================================================

# ==============================================================================
# 0. MAIN PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_dfm")

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

source(file.path(path_func, "dfm.utils.R"))
source(file.path(path_func, "dfm.prep.R"))
source(file.path(path_func, "dfm.fs.R"))
source(file.path(path_func, "dfm.kalman.R"))
source(file.path(path_func, "dfm.init.R"))
source(file.path(path_func, "dfm.em.R"))
source(file.path(path_func, "dfm.now.R"))
source(file.path(path_func, "dfm.all_cc_results.R"))

# ==============================================================================
# 3. LOCAL HELPERS
# ==============================================================================

build_run_tag_dfm <- function(params) {
  paste0(
    "eval-", format(params$start_eval, "%Y%m"), "_", format(params$end_eval, "%Y%m"),
    "_sel-", params$sel_method,
    "_Nm-", params$n_m,
    "_Nq-", params$n_q,
    "_Kmax-", params$Kmax,
    "_pmax-", params$pmax,
    "_qmax-", params$qmax,
    "_Idio-", params$idio_spec,
    "_CovidM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-", as.integer(isTRUE(params$covid_mask_q))
  )
}

get_size_tag_dfm <- function(N_m, N_q) {
  N_tot <- N_m + N_q
  
  if (N_tot <= 25) {
    "small"
  } else if (N_tot <= 50) {
    "medium"
  } else {
    "large"
  }
}

build_result_filename_dfm <- function(path_out,
                                      model,
                                      stage,
                                      Size,
                                      sel,
                                      countries,
                                      N_m,
                                      N_q,
                                      r = NA,
                                      p = NA,
                                      q = NA,
                                      idio_spec = "auto",
                                      covid_m,
                                      covid_q,
                                      ext = "rds",
                                      timestamp = TRUE) {
  
  stamp <- if (timestamp) paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")) else ""
  
  file.path(
    path_out,
    paste0(
      toupper(model), "_",
      stage,
      "_Size-", Size,
      "_sel-", sel,
      "_Idio-", idio_spec,
      "_cc-", paste(countries, collapse = "-"),
      "_Nm-", N_m,
      "_Nq-", N_q,
      "_r-", r,
      "_p-", p,
      "_q-", q,
      "_CovidM-", covid_m,
      "_CovidQ-", covid_q,
      stamp,
      ".",
      ext
    )
  )
}

# ==============================================================================
# 4. USER PARAMETERS
# ==============================================================================

params <- list(
  start_est    = as.Date("2000-04-01"),
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2026-02-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,
  covid_mask_q = TRUE,
  
  target       = "GDP",
  
  sel_method   = "corr",
  n_m          = 20,
  n_q          = 5,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,
  
  Kmax         = 10,
  pmax         = 5,
  qmax         = 5,
  
  # "auto" -> q selected automatically
  # "q0"   -> q fixed to 0
  idio_spec    = "q0",          # "auto" | "q0" | 
  
  kappa        = 1e-4,
  restr        = "stock_flow",
  
  max_iter     = 200,
  tol          = 1e-2
)

if (!params$idio_spec %in% c("auto", "q0")) {
  stop("params$idio_spec must be either 'auto' or 'q0'.")
}

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")

model_name <- "vector_dfm"
sel        <- params$sel_method
Size       <- get_size_tag_dfm(params$n_m, params$n_q)
tag_run    <- build_run_tag_dfm(params)

# ==============================================================================
# 5. OUTPUT PATHS
# ==============================================================================
# Separate folders for q-auto and q0, so outputs are not overwritten.

path_results <- file.path(
  path_main,
  "TPRF_Models_EA/Vector_DFM/results/outputs",
  params$idio_spec
)

path_graph <- file.path(
  path_main,
  "TPRF_Models_EA/Vector_DFM/results/graph",
  params$idio_spec
)

path_graph_rt <- file.path(
  path_main,
  "TPRF_Models_EA/Vector_DFM/results/graph_rt",
  params$idio_spec
)

dir.create(path_results,  recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph,    recursive = TRUE, showWarnings = FALSE)
dir.create(path_graph_rt, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. PREPARE DATA
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
# 7. CROSS-COUNTRY PIPELINE
# ==============================================================================

results_all <- run_all_countries_dfm(
  countries     = countries,
  all_countries = all_countries,
  params        = params,
  path_results  = path_results
)

summary_all <- lapply(names(results_all), function(cc) {
  res <- results_all[[cc]]
  
  summarize_dfm_country(
    DFM_res = res$full_sample,
    RT_res  = res$pseudo_realtime,
    country = cc,
    dates_m = res$inputs$dates_m,
    dates_q = res$inputs$dates_q,
    y_q     = res$inputs$y_q,
    params  = res$params
  )
})

names(summary_all) <- names(results_all)

cross_out <- build_cross_country_outputs_dfm(
  summary_all = summary_all,
  params      = params,
  model_label = paste0("Vector DFM (", params$idio_spec, ")")
)

# ==============================================================================
# 8. REAL-TIME PLOT VARIANTS
# ==============================================================================

rt_plot_variants <- build_rt_plot_variants(
  df_rt_all      = cross_out$df_rt_all,
  df_yq_eval_all = cross_out$df_yq_eval_all,
  params         = params,
  model_label    = paste0("Vector DFM (", params$idio_spec, ")")
)

rt_plot_files <- save_rt_plot_variants(
  rt_plots      = rt_plot_variants,
  path_graph_rt = path_graph_rt,
  model_name    = paste0(model_name, "_", params$idio_spec),
  Size          = Size,
  sel           = sel
)

# ==============================================================================
# 9. PRINT OUTPUTS
# ==============================================================================

print(cross_out$plot_nowcast_facet)
print(cross_out$plot_rt_facet)

print(rt_plot_variants$plot_all)
print(rt_plot_variants$plot_big4)
print(rt_plot_variants$plot_other4)
print(rt_plot_variants$plot_post8)

cat("\n================ DFM IN-SAMPLE LATEX ================\n")
cat(cross_out$latex_tab_insample_all)

cat("\n\n================ DFM ROLLING LATEX ================\n")
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
# 10. SAVE GRAPHS AND TEX
# ==============================================================================

path_graph_all <- file.path(path_graph, "ALL_COUNTRIES")
dir.create(path_graph_all, recursive = TRUE, showWarnings = FALSE)

file_graph_full <- file.path(
  path_graph_all,
  paste0("plot_full_", model_name, "_", tag_run, ".png")
)

file_graph_rt <- file.path(
  path_graph_all,
  paste0("plot_rt_", model_name, "_", tag_run, ".png")
)

file_tex_insample <- file.path(
  path_graph_all,
  paste0("tab_insample_", model_name, "_", tag_run, ".tex")
)

file_tex_rt <- file.path(
  path_graph_all,
  paste0("tab_rt_", model_name, "_", tag_run, ".tex")
)

ggsave(
  filename = file_graph_full,
  plot     = cross_out$plot_nowcast_facet,
  width    = 14,
  height   = 8,
  dpi      = 300
)

ggsave(
  filename = file_graph_rt,
  plot     = cross_out$plot_rt_facet,
  width    = 14,
  height   = 8,
  dpi      = 300
)

writeLines(cross_out$latex_tab_insample_all, con = file_tex_insample)
writeLines(cross_out$latex_tab_rt_all,       con = file_tex_rt)

# ==============================================================================
# 11. SAVE SUMMARY OBJECT
# ==============================================================================

q_label <- params$idio_spec

file_summary_cross <- build_result_filename_dfm(
  path_out   = path_results,
  model      = model_name,
  stage      = "summary",
  Size       = Size,
  sel        = sel,
  countries  = countries,
  N_m        = params$n_m,
  N_q        = params$n_q,
  r          = NA,
  p          = NA,
  q          = q_label,
  idio_spec  = params$idio_spec,
  covid_m    = as.integer(isTRUE(params$covid_mask_m)),
  covid_q    = as.integer(isTRUE(params$covid_mask_q)),
  ext        = "rds",
  timestamp  = FALSE
)

if (file.exists(file_summary_cross)) {
  stop(
    "Summary file already exists and was NOT overwritten:\n",
    file_summary_cross,
    "\n\nChange params$idio_spec, change the output folder, or delete the old file manually if you really want to replace it."
  )
}

saveRDS(
  list(
    model_id               = "DFM_EM",
    model                  = model_name,
    stage                  = "cross_country",
    Size                   = Size,
    sel                    = sel,
    idio_spec              = params$idio_spec,
    params                 = params,
    countries              = countries,
    tag_run                = tag_run,
    missing_by_country     = missing_by_country,
    
    hyper_by_country       = lapply(summary_all, function(x) list(
      full = x$hyper_full,
      rt   = x$hyper_rt
    )),
    
    hyper_full_all         = cross_out$hyper_full_all,
    hyper_rt_pre_all       = cross_out$hyper_rt_pre_all,
    hyper_rt_post_all      = cross_out$hyper_rt_post_all,
    
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
    
    rt_plot_variants       = rt_plot_variants,
    rt_plot_files          = rt_plot_files,
    
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

cat("\nSaved additional real-time plot variants to:\n")
cat(rt_plot_files$file_graph_rt_all, "\n")
cat(rt_plot_files$file_graph_rt_big4, "\n")
cat(rt_plot_files$file_graph_rt_other4, "\n")
cat(rt_plot_files$file_graph_rt_post8, "\n")