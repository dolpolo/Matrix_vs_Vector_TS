# ==============================================================================
# Matrix MF-TPRF - Selection Table Across Sizes
# ==============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

# ==============================================================================
# 1. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"

path_results <- file.path(
  path_main,
  "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs"
)

path_final_results <- file.path(
  path_main,
  "TPRF_Models_EA/Final_Tab_Graph"
)

path_func <- file.path(path_main, "functions/functions_mat")

dir.create(path_final_results, recursive = TRUE, showWarnings = FALSE)

source(file.path(path_func, "matrix.mf.tprf.utils.R"))

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================

model_name <- "matrix"
sel_method <- "LASSO"   # "LASSO" oppure "corr"

country_cols <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT")

# ==============================================================================
# 3. LOAD FIT OBJECTS
# ==============================================================================

load_fit_by_size <- function(Size) {
  file_fit <- find_result_file(
    path  = path_results,
    model = model_name,
    stage = "fit",
    Size  = Size,
    sel   = sel_method
  )
  
  message("\nLoaded ", Size, " fit from:\n", file_fit)
  
  res <- readRDS(file_fit)
  
  if (is.null(res$selections_wide)) {
    stop("Object for Size = ", Size, " does not contain `selections_wide`.")
  }
  
  res
}

res_small  <- load_fit_by_size("small")
res_medium <- load_fit_by_size("medium")
res_large  <- load_fit_by_size("large")

# ==============================================================================
# 4. BUILD LATEX TABLE
# ==============================================================================

latex_final <- build_latex_selection_table(
  sel_small    = res_small$selections_wide,
  sel_medium   = res_medium$selections_wide,
  sel_large    = res_large$selections_wide,
  country_cols = country_cols,
  caption      = paste0(
    "Macroeconomic predictors: country-specific selection across information-set sizes, ",
    sel_method, " preselection"
  ),
  label        = paste0(
    "tab:variables_selected_country_size_",
    tolower(sel_method)
  )
)

cat(latex_final)

# ==============================================================================
# 5. SAVE OUTPUT
# ==============================================================================

file_latex <- file.path(
  path_final_results,
  paste0("selection_table_across_sizes_", sel_method, ".tex")
)

writeLines(latex_final, file_latex)

message("\nSelection table saved to:\n", file_latex)
