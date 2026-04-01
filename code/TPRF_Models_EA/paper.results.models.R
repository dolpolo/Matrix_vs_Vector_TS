# ==============================================================================
# Final Results Script
# Matrix MF-TPRF vs Vector MF-TPRF vs Vectorized-from-Tensor MF-TPRF
# ------------------------------------------------------------------------------
# This script:
#   1. Loads cross-country summary objects for a chosen Size and sel
#   2. Prints model-specific plots and LaTeX tables
#   3. Builds comparison tables across models
#   4. Builds relative RMSFE tables
#   5. Saves final comparison outputs
# ==============================================================================

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ==============================================================================
# 1. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"

path_final_results <- file.path(path_main, "TPRF_Models_EA/Final_Tab_Graph")
dir.create(path_final_results, recursive = TRUE, showWarnings = FALSE)

path_matrix_results <- file.path(path_main, "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs")
path_vector_results <- file.path(path_main, "TPRF_Models_EA/Vector_MF-TPRF/results/outputs")
path_vectensor_results <- file.path(path_main, "TPRF_Models_EA/VecTensor_MF-TPRF/results/outputs_vec")

path_func     <- file.path(path_main, "functions/functions_mat")
source(file.path(path_func, "matrix.mf.tprf.utils.R"))
# ==============================================================================
# 2. CHOOSE CONFIGURATION
# ==============================================================================

Size <- "large"   # "small" | "medium" | "large"
sel  <- "corr"    # "corr"  | "LASSO" 

# Model names used in cross-country summaries
model_matrix    <- "matrix"
model_vector    <- "vector"
model_vectensor <- "vectensor"

# ==============================================================================
# 3. LOAD SUMMARY OBJECTS
# ==============================================================================

file_matrix_summary <- find_result_file(
  path  = path_matrix_results,
  model = model_matrix,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

file_vector_summary <- find_result_file(
  path  = path_vector_results,
  model = model_vector,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

file_vectensor_summary <- find_result_file(
  path  = path_vectensor_results,
  model = model_vectensor,
  stage = "summary",
  Size  = Size,
  sel   = sel
)

summary_matrix    <- readRDS(file_matrix_summary)
summary_vector    <- readRDS(file_vector_summary)
summary_vectensor <- readRDS(file_vectensor_summary)

cat("\nLoaded matrix summary from:\n", file_matrix_summary, "\n")
cat("\nLoaded vector summary from:\n", file_vector_summary, "\n")
cat("\nLoaded vectensor summary from:\n", file_vectensor_summary, "\n")

# Safety overwrite from loaded objects if available
Size <- if (!is.null(summary_matrix$Size)) summary_matrix$Size else Size
sel  <- if (!is.null(summary_matrix$sel))  summary_matrix$sel  else sel

# ==============================================================================
# 4. MODEL-SPECIFIC TABLES
# ==============================================================================

df_matrix_insample <- normalize_summary_table(
  summary_matrix$tab_insample_all,
  "Matrix"
)

df_matrix_rt <- normalize_summary_table(
  summary_matrix$tab_rt_all,
  "Matrix"
)

df_vector_insample <- normalize_summary_table(
  summary_vector$tab_insample_all,
  "Vector"
)

df_vector_rt <- normalize_summary_table(
  summary_vector$tab_rt_all,
  "Vector"
)

df_vectensor_insample <- normalize_summary_table(
  summary_vectensor$tab_insample_all,
  "VecTensor"
)

df_vectensor_rt <- normalize_summary_table(
  summary_vectensor$tab_rt_all,
  "VecTensor"
)

# Keep only common countries across models
common_countries_insample <- Reduce(
  intersect,
  list(
    unique(df_matrix_insample$country),
    unique(df_vector_insample$country),
    unique(df_vectensor_insample$country)
  )
)

common_countries_rt <- Reduce(
  intersect,
  list(
    unique(df_matrix_rt$country),
    unique(df_vector_rt$country),
    unique(df_vectensor_rt$country)
  )
)

df_matrix_insample    <- df_matrix_insample    %>% filter(country %in% common_countries_insample)
df_vector_insample    <- df_vector_insample    %>% filter(country %in% common_countries_insample)
df_vectensor_insample <- df_vectensor_insample %>% filter(country %in% common_countries_insample)

df_matrix_rt    <- df_matrix_rt    %>% filter(country %in% common_countries_rt)
df_vector_rt    <- df_vector_rt    %>% filter(country %in% common_countries_rt)
df_vectensor_rt <- df_vectensor_rt %>% filter(country %in% common_countries_rt)

# ==============================================================================
# 5. SHOW MODEL-SPECIFIC PLOTS
# ==============================================================================

print(summary_matrix$plot_nowcast_facet)
print(summary_matrix$plot_rt_facet)

print(summary_vector$plot_nowcast_facet)
print(summary_vector$plot_rt_facet)

print(summary_vectensor$plot_nowcast_facet)
print(summary_vectensor$plot_rt_facet)

# ==============================================================================
# 6. MODEL-SPECIFIC LATEX TABLES
# ==============================================================================

cat("\n================ MATRIX IN-SAMPLE ================\n")
cat(summary_matrix$latex_tab_insample_all)

cat("\n================ MATRIX ROLLING ================\n")
cat(summary_matrix$latex_tab_rt_all)

cat("\n================ VECTOR IN-SAMPLE ================\n")
cat(summary_vector$latex_tab_insample_all)

cat("\n================ VECTOR ROLLING ================\n")
cat(summary_vector$latex_tab_rt_all)

cat("\n================ VECTENSOR IN-SAMPLE ================\n")
cat(summary_vectensor$latex_tab_insample_all)

cat("\n================ VECTENSOR ROLLING ================\n")
cat(summary_vectensor$latex_tab_rt_all)

# ==============================================================================
# 7. COMPARISON TABLES
# ==============================================================================

comp_insample <- build_comparison_table(
  df_matrix_insample,
  df_vector_insample,
  df_vectensor_insample
)

comp_rt <- build_comparison_table(
  df_matrix_rt,
  df_vector_rt,
  df_vectensor_rt
)

latex_comp_insample <- comparison_to_latex(
  comp_insample,
  caption = paste0("In-sample RMSFE comparison across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:comparison_insample_models_", Size, "_", sel)
)

latex_comp_rt <- comparison_to_latex(
  comp_rt,
  caption = paste0("Rolling pseudo-real-time RMSFE comparison across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:comparison_rt_models_", Size, "_", sel)
)

cat("\n================ COMPARISON IN-SAMPLE ================\n")
print(comp_insample)
cat("\n", latex_comp_insample, "\n")

cat("\n================ COMPARISON ROLLING ================\n")
print(comp_rt)
cat("\n", latex_comp_rt, "\n")

# ==============================================================================
# 8. RELATIVE RMSFE TABLES
# ==============================================================================

rel_insample <- build_relative_table(comp_insample)
rel_rt       <- build_relative_table(comp_rt)

latex_rel_insample <- relative_to_latex(
  rel_insample,
  caption = paste0("Relative in-sample RMSFE across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:relative_insample_models_", Size, "_", sel)
)

latex_rel_rt <- relative_to_latex(
  rel_rt,
  caption = paste0("Relative rolling pseudo-real-time RMSFE across models (", Size, ", sel = ", sel, ")"),
  label   = paste0("tab:relative_rt_models_", Size, "_", sel)
)

cat("\n================ RELATIVE IN-SAMPLE ================\n")
print(rel_insample)
cat("\n", latex_rel_insample, "\n")

cat("\n================ RELATIVE ROLLING ================\n")
print(rel_rt)
cat("\n", latex_rel_rt, "\n")


# ==============================================================================
# 9. FINAL LATEX TABLE
# ==============================================================================

latex_final_custom <- build_final_large_style_latex(
  summary_matrix    = summary_matrix,
  summary_vector    = summary_vector,
  summary_vectensor = summary_vectensor,
  Size              = Size,
  sel               = sel
)

cat(latex_final_custom)

# ==============================================================================
# 10. SAVE FINAL OUTPUTS
# ==============================================================================

suffix_out <- paste0("Size-", Size, "_sel-", sel)

writeLines(
  summary_matrix$latex_tab_insample_all,
  file.path(path_final_results, paste0("matrix_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_matrix$latex_tab_rt_all,
  file.path(path_final_results, paste0("matrix_rt_", suffix_out, ".tex"))
)

writeLines(
  summary_vector$latex_tab_insample_all,
  file.path(path_final_results, paste0("vector_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_vector$latex_tab_rt_all,
  file.path(path_final_results, paste0("vector_rt_", suffix_out, ".tex"))
)

writeLines(
  summary_vectensor$latex_tab_insample_all,
  file.path(path_final_results, paste0("vectensor_insample_", suffix_out, ".tex"))
)
writeLines(
  summary_vectensor$latex_tab_rt_all,
  file.path(path_final_results, paste0("vectensor_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_comp_insample,
  file.path(path_final_results, paste0("comparison_insample_", suffix_out, ".tex"))
)
writeLines(
  latex_comp_rt,
  file.path(path_final_results, paste0("comparison_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_rel_insample,
  file.path(path_final_results, paste0("relative_insample_", suffix_out, ".tex"))
)
writeLines(
  latex_rel_rt,
  file.path(path_final_results, paste0("relative_rt_", suffix_out, ".tex"))
)

writeLines(
  latex_final_custom,
  file.path(path_final_results, paste0("final_custom_table_", suffix_out, ".tex"))
)

saveRDS(
  list(
    Size                 = Size,
    sel                  = sel,
    file_matrix_summary  = file_matrix_summary,
    file_vector_summary  = file_vector_summary,
    file_vectensor_summary = file_vectensor_summary,
    
    summary_matrix       = summary_matrix,
    summary_vector       = summary_vector,
    summary_vectensor    = summary_vectensor,
    
    df_matrix_insample   = df_matrix_insample,
    df_matrix_rt         = df_matrix_rt,
    df_vector_insample   = df_vector_insample,
    df_vector_rt         = df_vector_rt,
    df_vectensor_insample = df_vectensor_insample,
    df_vectensor_rt       = df_vectensor_rt,
    
    comp_insample        = comp_insample,
    comp_rt              = comp_rt,
    rel_insample         = rel_insample,
    rel_rt               = rel_rt
  ),
  file.path(path_final_results, paste0("final_model_comparison_", suffix_out, ".rds"))
)

cat("\nAll final outputs saved to:\n", path_final_results, "\n")
