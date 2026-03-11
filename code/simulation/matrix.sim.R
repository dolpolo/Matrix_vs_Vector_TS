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
path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)
path_func_sim <- file.path(path_main, "simulation/functions")
path_func_vec <- file.path(path_main, "Vector_FM/functions")
path_func_mat <- file.path(path_main, "Matrix_FM/functions")

path_results <- file.path(path_main, "simulation/results/outputs")

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
dgp_params <- make_dgp_params(
  Tm = 100,
  burn = 400,
  draws = 500,
  seed = 1712,
  n_m_flow = 10,
  n_m_stock = 5,
  n_q_flow = 3,
  n_q_stock = 2
)

est_tprf_params <- list(
  # TPRF
  L_tprf       = 1,
  
  # MF-TPRF
  Lproxy_vec   = 1,
  L_midas_vec  = 1,
  p_AR_vec     = 1,
  Kmax_vec     = 3,
  
  # Matric MF-TPRF
  Lproxy_mat   = 1,
  L_midas_mat  = 1,
  p_AR_mat     = 1,
  Kmax_mat     = c(1, 3)
  
)


# ==============================================================================
# DGP simulator
# ==============================================================================
sim <- simulate_matrix_dgp(dgp_params)

# ==============================================================================
# TPRF - Forecast
# ==============================================================================

tprf_out <- run_tprf_wrapper(
  sim_obj = sim,
  L = est_tprf_params$L_tprf
)

tprf_out$results
# ==============================================================================
# MF-TPRF - Nowcast M1,M2,M3
# ==============================================================================

mftprf_M1 <- run_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M1",
  Lproxy  = est_tprf_params$Lproxy_vec,
  L_midas = est_tprf_params$L_midas_vec,
  p_AR    = est_tprf_params$p_AR_vec,
  Kmax    = est_tprf_params$Kmax_vec
)

mftprf_M1

mftprf_M2 <- run_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M2",
  Lproxy  = est_tprf_params$Lproxy_vec,
  L_midas = est_tprf_params$L_midas_vec,
  p_AR    = est_tprf_params$p_AR_vec,
  Kmax    = est_tprf_params$Kmax_vec
)

mftprf_M2

mftprf_M3 <- run_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M3",
  Lproxy  = est_tprf_params$Lproxy_vec,
  L_midas = est_tprf_params$L_midas_vec,
  p_AR    = est_tprf_params$p_AR_vec,
  Kmax    = est_tprf_params$Kmax_vec
)

mftprf_M3
rbind(mftprf_M1, mftprf_M2, mftprf_M3)


# ==============================================================================
# MATRIX MF-TPRF - Nowcast M1,M2,M3
# ==============================================================================

matrix_M1 <- run_matrix_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M1",
  Lproxy  = est_tprf_params$Lproxy_mat,
  L_midas = est_tprf_params$L_midas_mat,
  p_AR    = est_tprf_params$p_AR_mat,
  Kmax    = est_tprf_params$Kmax_mat
)

matrix_M1$results

matrix_M2 <- run_matrix_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M2",
  Lproxy  = est_tprf_params$Lproxy_mat,
  L_midas = est_tprf_params$L_midas_mat,
  p_AR    = est_tprf_params$p_AR_mat,
  Kmax    = est_tprf_params$Kmax_mat
)

matrix_M3 <- run_matrix_mftprf_wrapper(
  sim_obj = sim,
  vintage = "M3",
  Lproxy  = est_tprf_params$Lproxy_mat,
  L_midas = est_tprf_params$L_midas_mat,
  p_AR    = est_tprf_params$p_AR_mat,
  Kmax    = est_tprf_params$Kmax_mat
)

rbind(matrix_M1$results, matrix_M2$results, matrix_M3$results)

# ==============================================================================
# MONTECARLO SIMULATION
# ==============================================================================

mc_res <- run_mc_simulation(
  dgp_params = dgp_params,
  est_params = est_tprf_params,
  draws = dgp_params$draws,
  verbose = TRUE
)

mc_summary <- summarize_mc_results(mc_res)
mc_summary

head(mc_res$raw)
str(mc_res$raw)
table(mc_res$raw$model, mc_res$raw$vintage)

g1 <- plot_mc_boxplot(mc_res)
g1

g2 <- plot_mc_rmse(mc_res)
g2

g3 <- plot_mc_by_vintage(mc_res)
g3

# ------------------------------------------------------------------------------
# 1. Create MonteCarlo folder
# ------------------------------------------------------------------------------

path_mc <- file.path(path_results, "MonteCarlo")
if (!dir.exists(path_mc)) {
  dir.create(path_mc, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Build base filename
# ------------------------------------------------------------------------------

file_base <- paste0(
  "MC",
  "_Tm-",    dgp_params$Tm,
  "_burn-",  dgp_params$burn,
  "_p1-",    dgp_params$p1,
  "_Nm-",    dgp_params$Nm,
  "_Nq-",    dgp_params$Nq,
  "_draws-", dgp_params$draws,
  "_rhoF-",  dgp_params$rhoF,
  "_rhoG-",  dgp_params$rhoG,
  "_rhoE-",  dgp_params$rhoE,
  "_R2-",    dgp_params$R2_BI,
  "_Ltprf-", est_tprf_params$L_tprf,
  "_LproxyV-", est_tprf_params$Lproxy_mc,
  "_LmidasV-", est_tprf_params$L_midas_mc,
  "_pARV-",    est_tprf_params$p_AR_mc,
  "_LproxyM-", est_tprf_params$Lproxy_mat,
  "_LmidasM-", est_tprf_params$L_midas_mat,
  "_pARM-",    est_tprf_params$p_AR_mat
)

# ------------------------------------------------------------------------------
# 3. Save RDS
# ------------------------------------------------------------------------------
# file_base <- "MC_Tm-400_burn-100_p1-5_Nm-60_Nq-20_draws-500_rhoF-0.9_rhoG-0.3_rhoE-0_R2-0.5_Ltprf-1_LproxyV-_LmidasV-_pARV-_LproxyM-1_LmidasM-1_pARM-1"
file_rds <- file.path(path_mc, paste0(file_base, ".rds"))

saveRDS(
  list(
    dgp_params      = dgp_params,
    est_tprf_params = est_tprf_params,
    mc_res          = mc_res,
    mc_summary      = mc_summary
  ),
  file_rds
)

# ------------------------------------------------------------------------------
# 4. Save plots
# ------------------------------------------------------------------------------

ggsave(file.path(path_mc, paste0(file_base, "_boxplot_error.png")),
       g1, width = 10, height = 6)

ggsave(file.path(path_mc, paste0(file_base, "_rmse.png")),
       g2, width = 10, height = 6)

ggsave(file.path(path_mc, paste0(file_base, "_by_vintage.png")),
       g3, width = 10, height = 6)

cat("\n*** SAVED MONTE CARLO RESULTS TO ***\n")
cat(file_rds, "\n")
