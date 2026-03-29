# ==============================================================================
#             Three Pass Regression Filter : A Simulation Study
#                               - Main Script - 
# ==============================================================================
# Data:
#
# ==============================================================================
# Purpose: 
#
# ==============================================================================
# Author      : Davide Delfino
#
# Institution : Bocconi University - Milan
#
# ==============================================================================
# Please don't hesitate to write to davide.delfino@unibocconi.it in case you need
# any kind clarification
# ==============================================================================
# Script Type : Simulation Study
# ==============================================================================

# ==============================================================================
# SET WORKING DIRECTORY
# ==============================================================================
if (dir.exists("/content/drive/MyDrive/Matrix_MF-TPRF/code")) {
  path_main <- "/content/drive/MyDrive/Matrix_MF-TPRF/code"
} else {
  path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
}

setwd(path_main)
path_func_sim <- file.path(path_main, "functions/functions_sim")
path_func_vec <- file.path(path_main, "functions/functions_vec")
path_func_mat <- file.path(path_main, "functions/functions_mat")

path_results <- file.path(path_main, "MonteCarlo_simulation/results/outputs")

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

# Data Handling
library(tidyverse)
library(lubridate)
library(abind)
library(forcats)

# Time Series Analysis
library(tseries)
library(zoo)

# I/O
library(writexl)
library(readxl)
library(xtable)

# Plotting
library(ggplot2)

# conflict
library(conflicted)
conflict_prefer("select", "dplyr", quiet = TRUE)
conflict_prefer("filter", "dplyr", quiet = TRUE)


# ==============================================================================
# LOAD USERS' FUNCTIONS LIBRARIES
# ==============================================================================

### Simulation
source(file.path(path_func_sim, "matrix.dgp.R"))
source(file.path(path_func_sim, "mc.sim.now.R"))


### TPRF
source(file.path(path_func_vec, "tprf.R"))

### MF-TPRF
# data preparation
source(file.path(path_func_vec, "mf.tprf.utils.R"))
source(file.path(path_func_vec, "mf.tprf.prep.R"))

# mf-tprf estimation
source(file.path(path_func_vec, "mf.tprf.R"))

# imputation
source(file.path(path_func_vec, "mf.tprf.imp.R"))

# factor selection
source(file.path(path_func_vec, "mf.tprf.fs.R"))

# nowcasting
source(file.path(path_func_vec, "mf.tprf.now.R"))

### MATRIX MF-TPRF
# data preparation
source(file.path(path_func_mat, "matrix.mf.tprf.utils.R"))
source(file.path(path_func_mat, "matrix.mf.tprf.prep.R"))

# Matrix mf-tprf estimation
source(file.path(path_func_mat, "matrix.mf.tprf.R"))

# NANs imputation
source(file.path(path_func_mat, "matrix.mf.tprf.imp.R"))

# factor selection
source(file.path(path_func_mat, "matrix.mf.tprf.fs.R"))

# nowcasting
source(file.path(path_func_mat, "matrix.mf.tprf.now.R"))

# ==============================================================================
# Params
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

Tm_in   <- if (length(args) >= 1) as.numeric(args[1]) else 402
p1_in   <- if (length(args) >= 2) as.numeric(args[2]) else 10
p2_in   <- if (length(args) >= 3) as.numeric(args[3]) else 100 
rhoG_in <- if (length(args) >= 4) as.numeric(args[4]) else 0.3
rhoE_in <- if (length(args) >= 5) as.numeric(args[5]) else 0.9


dgp_params <- make_dgp_params(
  Tm = Tm_in,                  # 102, 201, 402
  p1 = p1_in,                  # 10, 20
  p2 = p2_in,                  # 25, 50, 100
  r1 = 1,
  r2 = 1,
  v1 = 1,
  v2 = 4,
  rhoG = rhoG_in,              # 0.3, 0.95
  rhoE = rhoE_in,              # 0.0, 0.9
  share_monthly = 0.8,
  share_m_flow = 0.5,
  share_q_flow = 0.5
)

est_nowcast_params <- make_nowcast_estimation_params(
  Lproxy_vec  = 1,
  L_midas_vec = 1,
  p_AR_vec    = 1,
  Lproxy_mat  = 1,
  L_midas_mat = 1,
  p_AR_mat    = 1,
  r_mat_true  = c(1, 1),
  proxy_name  = "EA",
  standardize_proxy = TRUE,
  orthonormalize_each_iter = TRUE,
  orthonormalize_final_Z   = TRUE,
  ils_maxit = 100,
  ils_tol   = 1e-8,
  bi_intercept = TRUE,
  L_tprf = 1,
  use_tprf_standard = TRUE
)

# inspect one simulated dataset
sim <- simulate_matrix_dgp(dgp_params)

# ==============================================================================
# MONTE CARLO SIMULATION
# ==============================================================================
mc_res <- run_mc_simulation(
  dgp_params = dgp_params,
  est_params = est_nowcast_params,
  draws = dgp_params$draws,
  verbose = TRUE
)

# ------------------------------------------------------------------------------
# 1. Main summaries
# ------------------------------------------------------------------------------
mc_summary <- summarize_mc_results(mc_res)
mc_rmsfe_wide <- summarize_mc_rmsfe_wide(mc_res)
mc_mat_compact_rmsfe <- summarize_mc_mat_compact_rmsfe(mc_res)

# TPRF month-by-month summaries
mc_tprf_by_vintage <- summarize_tprf_forecasts_by_vintage(mc_res)
mc_tprf_h_wide     <- summarize_tprf_forecasts_wide(mc_res, value = "avg_horizon")
mc_tprf_rmsfe_wide <- summarize_tprf_forecasts_wide(mc_res, value = "RMSFE")

# ------------------------------------------------------------------------------
# 2. Main LaTeX tables
# ------------------------------------------------------------------------------
tex_mc_mat_compact <- mc_mat_compact_rmsfe_to_latex(
  mc_obj = mc_res,
  digits = 3,
  caption = "Matrix MF-TPRF RMSFE and relative performance by nowcast vintage",
  label = "tab:mc_mat_compact_rmsfe"
)

cat(tex_mc_mat_compact)

# ------------------------------------------------------------------------------
# 4. Plots
# ------------------------------------------------------------------------------
g_vint  <- plot_mc_by_vintage(mc_res)
g_vint

# ------------------------------------------------------------------------------
# Create MonteCarlo folder
# ------------------------------------------------------------------------------
path_mc <- file.path(path_results, "MonteCarlo")
if (!dir.exists(path_mc)) {
  dir.create(path_mc, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Build base filename
# ------------------------------------------------------------------------------
file_base <- paste0(
  "MC_nowcast",
  "_T",  dgp_params$Tm,
  "_p1", dgp_params$p1,
  "_p2", dgp_params$p2,
  "_D",  dgp_params$draws,
  "_rF", dgp_params$rhoF,
  "_rG", dgp_params$rhoG,
  "_rE", dgp_params$rhoE,
  "_R2", dgp_params$R2_BI,
  "_V",  est_nowcast_params$Lproxy_vec, "-", est_nowcast_params$L_midas_vec, "-", est_nowcast_params$p_AR_vec,
  "_M",  est_nowcast_params$Lproxy_mat, "-", est_nowcast_params$L_midas_mat, "-", est_nowcast_params$p_AR_mat,
  "_Tprf", est_nowcast_params$L_tprf,
  "_r",  paste(est_nowcast_params$r_mat_true, collapse = "")
)

# ------------------------------------------------------------------------------
# Save main LaTeX tables
# ------------------------------------------------------------------------------
write_mc_mat_compact_rmsfe_latex(
  mc_obj = mc_res,
  file = file.path(path_mc, paste0(file_base, "_summary_mat_compact_rmsfe.tex")),
  digits = 3,
  caption = "Matrix MF-TPRF RMSFE and relative performance by nowcast vintage",
  label = "tab:mc_mat_compact_rmsfe"
)

# ------------------------------------------------------------------------------
# Save main plots
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(path_mc, paste0(file_base, "_by_vintage.png")),
  plot = g_vint,
  width = 10,
  height = 6
)