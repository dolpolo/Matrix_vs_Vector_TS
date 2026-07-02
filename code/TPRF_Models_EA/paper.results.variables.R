# ==============================================================================
# Matrix MF-TPRF - Dated Predictor-Selection Tables
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(lubridate)
  library(readxl)
  library(glmnet)
})

# ==============================================================================
# 1. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"

path_results <- file.path(
  path_main,
  "TPRF_Models_EA/Matrix_MF-TPRF/results/outputs"
)

# Keep these identical to the paths used in the estimation script.
path_data_raw <- file.path(
  path_main,
  "data/raw"
)

path_data_adj <- file.path(
  path_main,
  "data/data_TR2"
)

path_func <- file.path(
  path_main,
  "functions/functions_mat"
)

path_output <- file.path(
  path_main,
  "TPRF_Models_EA/Final_Tab_Graph/dated_selection_tables"
)

dir.create(
  path_output,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==============================================================================
# 2. SOURCE FUNCTIONS
# ==============================================================================

source(file.path(path_func, "matrix.mf.tprf.utils.R"))
source(file.path(path_func, "matrix.mf.tprf.prep.R"))
source(file.path(path_func, "matrix.mf.tprf.now.R"))

# ==============================================================================
# 3. CONFIGURATION
# ==============================================================================

config_selection <- list(
  model         = "matrix_scalar",
  sel_method    = "LASSO",
  sizes         = c("small", "medium", "large"),
  countries     = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  path_results  = path_results,
  path_data_raw = path_data_raw,
  path_data_adj = path_data_adj
)

method_tag <- tolower(config_selection$sel_method)

method_label <- if (toupper(config_selection$sel_method) == "LASSO") {
  "LASSO"
} else {
  "corr"
}

# ==============================================================================
# 4. REBUILD INITIAL AND POST-COVID SELECTION REGIMES
# ==============================================================================

selection_report <- build_dated_selection_report(
  config = config_selection
)

# ==============================================================================
# 5. TABLE 1: INITIAL DATED SELECTION
# ==============================================================================

note_initial <- paste0(
  "\\textit{Notes:} The table reports the country-specific predictor sets ",
  "selected under the initial dated ", method_label,
  " regime. The screening sample ends in ",
  format_month_year(selection_report$dates$selection_end_initial),
  ", and the resulting sets are used from ",
  format_month_year(selection_report$dates$initial_regime_start),
  " until the scheduled post-COVID update. For each country, \\textbf{S}, ",
  "\\textbf{M}, and \\textbf{L} denote selection in the small, medium, and ",
  "large information sets, respectively; combined entries indicate selection ",
  "in more than one information set. \\textbf{Cl.} denotes variable class, ",
  "\\textbf{Cat.} distinguishes hard and soft indicators, \\textbf{Tr.} ",
  "reports the transformation code, \\textbf{Fr.} the sampling frequency, ",
  "and \\textbf{Del.} the approximate release delay in days. GDP is excluded."
)

latex_initial <- build_latex_selection_table(
  sel_small    = selection_report$initial$small$wide,
  sel_medium   = selection_report$initial$medium$wide,
  sel_large    = selection_report$initial$large$wide,
  country_cols = config_selection$countries,
  caption      = paste0(
    "Country-specific ", method_label,
    " predictor selection: initial regime."
  ),
  label = paste0(
    "tab:selection_initial_",
    method_tag
  ),
  note = note_initial
)

# ==============================================================================
# 6. TABLE 2: SCHEDULED POST-COVID UPDATE
# ==============================================================================

note_post_covid <- paste0(
  "\\textit{Notes:} The table reports the country-specific predictor sets ",
  "selected under the scheduled post-COVID ", method_label,
  " update. The screening sample ends in ",
  format_month_year(selection_report$dates$selection_end_updated),
  ", and the updated sets are used from ",
  format_month_year(selection_report$dates$updated_regime_start),
  " onward. For each country, \\textbf{S}, \\textbf{M}, and \\textbf{L} ",
  "denote selection in the small, medium, and large information sets, ",
  "respectively; combined entries indicate selection in more than one ",
  "information set. \\textbf{Cl.} denotes variable class, \\textbf{Cat.} ",
  "distinguishes hard and soft indicators, \\textbf{Tr.} reports the ",
  "transformation code, \\textbf{Fr.} the sampling frequency, and ",
  "\\textbf{Del.} the approximate release delay in days. GDP is excluded."
)

latex_post_covid <- build_latex_selection_table(
  sel_small    = selection_report$updated$small$wide,
  sel_medium   = selection_report$updated$medium$wide,
  sel_large    = selection_report$updated$large$wide,
  country_cols = config_selection$countries,
  caption      = paste0(
    "Country-specific ", method_label,
    " predictor selection: scheduled post-COVID update."
  ),
  label = paste0(
    "tab:selection_postcovid_",
    method_tag
  ),
  note = note_post_covid
)

# ==============================================================================
# 7. TABLE 3: PREDICTOR-SET COMPOSITION
# ==============================================================================

latex_dimensions <- build_latex_predictor_dimensions_table(
  tbl = selection_report$predictor_dimensions,
  caption = "Country-union predictor-set dimensions.",
  label = paste0(
    "tab:selection_composition_",
    method_tag
  ),
  evaluation_end = selection_report$fit_results$small$params$end_eval
)

# ==============================================================================
# 7B. TABLE 4: REAL-TIME PREDICTOR AVAILABILITY
# ==============================================================================

realtime_availability_tbl <- build_missing_matrix_table(
  selection_report = selection_report,
  var_scope        = "union"
)

latex_realtime_availability <- build_latex_missing_matrix_table(
  tbl              = realtime_availability_tbl,
  selection_report = selection_report,
  caption = paste0(
    "Real-time predictor availability across information sets, dated ",
    method_label,
    " selection regimes, and nowcast vintages."
  ),
  label = paste0(
    "tab:realtime_predictor_availability_",
    method_tag
  )
)

# ==============================================================================
# 7C. TABLE: PREDICTOR AVAILABILITY AT CALIBRATION VINTAGES
# ==============================================================================

calibration_availability_tbl <- build_calibration_missing_matrix_table(
  selection_report = selection_report,
  var_scope        = "union"
)

latex_calibration_availability <- build_latex_calibration_missing_matrix_table(
  tbl = calibration_availability_tbl,
  caption = paste0(
    "Predictor availability at the dated ",
    method_label,
    " calibration vintages."
  ),
  label = paste0(
    "tab:calibration_predictor_availability_",
    method_tag
  )
)

# ==============================================================================
# 8. TABLE 5: SELECTION TURNOVER
# ==============================================================================

latex_turnover <- build_latex_selection_turnover_table(
  tbl = selection_report$selection_turnover,
  caption = "Country-specific predictor-selection turnover across dated regimes.",
  label = paste0(
    "tab:country_predictor_turnover_",
    method_tag
  )
)

# ==============================================================================
# 9. TABLE 6: COUNTRY-SPECIFIC ADDED AND REMOVED PREDICTORS
# ==============================================================================

latex_changes <- build_latex_selection_changes_by_status_table(
  selection_transition = selection_report$selection_transition,
  country_order        = config_selection$countries,
  caption = paste0(
    "Country-specific predictor changes in the scheduled post-COVID ",
    method_label, " update."
  ),
  label = paste0(
    "tab:selection_changes_",
    method_tag
  )
)

# ==============================================================================
# 10. SAVE EACH TABLE AS A SEPARATE LATEX FILE
# ==============================================================================

writeLines(
  latex_initial,
  file.path(
    path_output,
    paste0("01_selection_initial_", method_tag, ".tex")
  )
)

writeLines(
  latex_post_covid,
  file.path(
    path_output,
    paste0("02_selection_postcovid_", method_tag, ".tex")
  )
)

writeLines(
  latex_dimensions,
  file.path(
    path_output,
    paste0("03_union_tensor_dimensions_", method_tag, ".tex")
  )
)

writeLines(
  latex_realtime_availability,
  file.path(
    path_output,
    paste0("04_realtime_predictor_availability_", method_tag, ".tex")
  )
)

writeLines(
  latex_turnover,
  file.path(
    path_output,
    paste0("05_country_predictor_turnover_", method_tag, ".tex")
  )
)

writeLines(
  latex_changes,
  file.path(
    path_output,
    paste0("06_selection_changes_", method_tag, ".tex")
  )
)

writeLines(
  latex_calibration_availability,
  file.path(
    path_output,
    paste0(
      "07_calibration_predictor_availability_",
      method_tag,
      ".tex"
    )
  )
)
# ==============================================================================
# 11. OPTIONAL: PRINT ONE TABLE AT A TIME
# ==============================================================================

cat(latex_initial)
cat(latex_post_covid)
cat(latex_dimensions)
cat(latex_realtime_availability)
cat(latex_calibration_availability)
cat(latex_turnover)
cat(latex_changes)



# ==============================================================================
# 12. UNIT-TARGET PROXY ROBUSTNESS TABLES
# ==============================================================================

path_output_unit_target <- file.path(
  path_output,
  "unit_target_proxy_tables"
)

dir.create(
  path_output_unit_target,
  recursive = TRUE,
  showWarnings = FALSE
)

unit_target_lasso <- make_unit_target_proxy_tables_from_files(
  path_results  = path_results,
  sel           = "LASSO",
  country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  model_unit    = "matrix_multivariate",
  model_aggregate = "matrix_scalar",
  digits = 3,
  caption_rmsfe = paste0(
    "Unit-target Matrix MF--TPRF RMSFE: LASSO-based selection"
  ),
  caption_relative = paste0(
    "Relative RMSFE of unit-target to aggregate-target Matrix MF--TPRF: ",
    "LASSO-based selection"
  ),
  label_rmsfe = "tab:unit_target_rmsfe_lasso",
  label_relative = "tab:unit_target_relative_rmsfe_lasso"
)

unit_target_corr <- make_unit_target_proxy_tables_from_files(
  path_results  = path_results,
  sel           = "corr",
  country_order = c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT"),
  model_unit    = "matrix_multivariate",
  model_aggregate = "matrix_scalar",
  digits = 3,
  caption_rmsfe = paste0(
    "Unit-target Matrix MF--TPRF RMSFE: correlation screening"
  ),
  caption_relative = paste0(
    "Relative RMSFE of unit-target to aggregate-target Matrix MF--TPRF: ",
    "correlation screening"
  ),
  label_rmsfe = "tab:unit_target_rmsfe_corr",
  label_relative = "tab:unit_target_relative_rmsfe_corr"
)

# ==============================================================================
# 13. SAVE UNIT-TARGET PROXY TABLES
# ==============================================================================

writeLines(
  unit_target_lasso$latex_rmsfe,
  file.path(
    path_output_unit_target,
    "08_unit_target_rmsfe_lasso.tex"
  )
)

writeLines(
  unit_target_lasso$latex_relative,
  file.path(
    path_output_unit_target,
    "09_unit_target_relative_rmsfe_lasso.tex"
  )
)

writeLines(
  unit_target_corr$latex_rmsfe,
  file.path(
    path_output_unit_target,
    "10_unit_target_rmsfe_corr.tex"
  )
)

writeLines(
  unit_target_corr$latex_relative,
  file.path(
    path_output_unit_target,
    "11_unit_target_relative_rmsfe_corr.tex"
  )
)

# ==============================================================================
# 14. PRINT LATEX TABLES AND UNDERLYING DATA
# ==============================================================================

cat(unit_target_lasso$latex_rmsfe)
cat(unit_target_lasso$latex_relative)

cat(unit_target_corr$latex_rmsfe)
cat(unit_target_corr$latex_relative)

print(unit_target_lasso$rmsfe_table)
print(unit_target_lasso$relative_table)

print(unit_target_corr$rmsfe_table)
print(unit_target_corr$relative_table)
