# ==============================================================================
# Matrix MF-TPRF - Selection Table Across Sizes
# ==============================================================================

library(dplyr)
library(tidyr)

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
path_results <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_func    <- file.path(path_main, "functions/functions_mat")

source(file.path(path_func, "matrix.mf.tprf.utils.R"))

model_name <- "matrix"
sel_method <- "LASSO"

# --- load small ---
file_fit_small <- find_result_file(
  path  = path_results,
  model = model_name,
  stage = "fit",
  Size  = "small",
  sel   = sel_method
)
res_small <- readRDS(file_fit_small)

# --- load medium ---
file_fit_medium <- find_result_file(
  path  = path_results,
  model = model_name,
  stage = "fit",
  Size  = "medium",
  sel   = sel_method
)
res_medium <- readRDS(file_fit_medium)

# --- load large ---
file_fit_large <- find_result_file(
  path  = path_results,
  model = model_name,
  stage = "fit",
  Size  = "large",
  sel   = sel_method
)
res_large <- readRDS(file_fit_large)

# --- build final latex table ---
latex_final <- build_latex_selection_table(
  sel_small  = res_small$selections_wide,
  sel_medium = res_medium$selections_wide,
  sel_large  = res_large$selections_wide
)

cat(latex_final)
