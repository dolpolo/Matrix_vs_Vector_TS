# ==============================================================================
# Vector MF-TPRF Main Script
# Mixed-Frequency Nowcasting for Euro Area GDP
# ------------------------------------------------------------------------------
# Single-country version with standardized result saving
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_vec")
path_utils    <- file.path(path_main, "functions/functions_mat")

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
  
  sel_method   = "LASSO",
  n_m          = 20,
  n_q          = 5,
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

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
country   <- "EA"

tag_run    <- build_run_tag(params)
model_name <- "vector"
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
# 6. SINGLE-COUNTRY INPUTS
# ==============================================================================

data    <- all_countries$data[[country]]$Data
dates_m <- all_countries$data[[country]]$Dates
dates_q <- all_countries$data[[country]]$DatesQ

gdp_col     <- all_countries$data[[country]]$target_col
target_name <- paste0(params$target, "_", country)

y   <- as.matrix(data[, gdp_col, drop = FALSE])
y_q <- as.matrix(y[!is.na(y)])

X <- as.matrix(data[, -gdp_col, drop = FALSE])

N_m <- all_countries$data[[country]]$nM

N_q_tot      <- all_countries$data[[country]]$nQ
q_series_all <- all_countries$data[[country]]$Series[(N_m + 1):(N_m + N_q_tot)]
q_cols       <- which(toupper(q_series_all) != toupper(target_name))
N_q          <- length(q_cols)
N            <- N_m + N_q

Size <- get_size_tag(N_m, N_q)

agg_m <- all_countries$data[[country]]$agg_m
agg_q <- all_countries$data[[country]]$agg_q[q_cols]
agg   <- c(agg_m, agg_q)

Freq_m <- all_countries$data[[country]]$freq_m
Freq_q <- all_countries$data[[country]]$freq_q[q_cols]
Freq   <- c(Freq_m, Freq_q)

Unb_m <- all_countries$data[[country]]$unb_m
Unb_q <- all_countries$data[[country]]$unb_q[q_cols]
Unb   <- c(Unb_m, Unb_q)

Class_m <- all_countries$data[[country]]$ClassM
Class_q <- all_countries$data[[country]]$ClassQ[q_cols]
Class   <- c(Class_m, Class_q)

# ==============================================================================
# 7. STANDARDIZATION, IMPUTATION, AND INITIAL FACTOR SELECTION
# ==============================================================================

out_std <- standardize_with_na(X)
X_std   <- out_std$X_std

imp_xp <- init_XP_ER(X_std, params$Kmax)

X_xp  <- imp_xp$X_init
r_hat <- imp_xp$r

X_m_xp <- X_xp[, 1:N_m, drop = FALSE]
X_q_xp <- X_xp[, (N_m + 1):N, drop = FALSE]

X_mq_xp  <- agg_mq(X_m_xp, agg_m)
X_qq_xp  <- agg_qq(X_q_xp, agg_q)
X_xp_agg <- cbind(X_mq_xp, X_qq_xp)

# ==============================================================================
# 8. PROXY AND LAG SELECTION
# ==============================================================================

pls_object <- select_L_autoproxy_3prf(
  X_xp_agg,
  y_q,
  Zmax = params$Zmax
)

Lproxy <- pls_object$L_opt

lag_sel <- choose_UMIDAS_lag(
  X_lf        = X_xp_agg,
  X_hf        = X_xp,
  y_q         = y_q,
  Lmax        = params$Lmax,
  Lproxy      = Lproxy,
  p_AR_max    = params$p_AR_max,
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
)

L_midas <- lag_sel$best_BIC$L
p_ar    <- lag_sel$best_BIC$p_AR

# ==============================================================================
# 9. FULL-SAMPLE ESTIMATION
# ==============================================================================

fit_vector <- MF_TPRF(
  X_lf        = X_xp_agg,
  X_hf        = X_xp,
  y_q         = y_q,
  Lproxy      = Lproxy,
  L_midas     = L_midas,
  p_AR        = p_ar,
  Robust_F    = params$Robust_F,
  alpha       = params$alpha,
  robust_type = params$robust_type,
  nw_lag      = params$nw_lag
)

# ==============================================================================
# 10. SAVE SINGLE-COUNTRY FULL-SAMPLE RESULTS
# ==============================================================================

path_country <- file.path(path_results, country)
dir.create(path_country, recursive = TRUE, showWarnings = FALSE)

file_fit <- build_result_filename(
  path_out         = path_country,
  model            = model_name,
  stage            = "fit",
  Size             = Size,
  sel              = sel,
  countries        = country,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = Lproxy,
  L_midas          = L_midas,
  p_ar             = p_ar,
  r1               = r_hat,
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
    model_id = "MF_TPRF",
    model    = model_name,
    stage    = "full_sample",
    Size     = Size,
    sel      = sel,
    country  = country,
    params   = params,
    dates_m  = dates_m,
    dates_q  = dates_q,
    y_q      = y_q,
    metadata = list(
      N_m   = N_m,
      N_q   = N_q,
      N     = N,
      agg_m = agg_m,
      agg_q = agg_q,
      Freq  = Freq,
      Unb   = Unb,
      Class = Class
    ),
    preprocessing = list(
      out_std  = out_std,
      imp_xp   = imp_xp,
      X_xp     = X_xp,
      X_mq_xp  = X_mq_xp,
      X_qq_xp  = X_qq_xp,
      X_xp_agg = X_xp_agg
    ),
    hyper = list(
      r_hat   = r_hat,
      Lproxy  = Lproxy,
      L_midas = L_midas,
      p_ar    = p_ar,
      pls     = pls_object,
      lag_sel = lag_sel
    ),
    fit = fit_vector
  ),
  file = file_fit
)

cat("\nSaved full-sample results to:\n", file_fit, "\n")

# ==============================================================================
# 11. PSEUDO REAL-TIME NOWCASTING
# ==============================================================================

pseudo_realtime_raw <- pseudo_realtime_MF_TPRF_XP(
  X_full  = X,
  y_q     = y_q,
  params  = params,
  dates   = dates_m,
  dates_q = dates_q,
  Freq    = Freq,
  Unb     = Unb,
  agg_m   = agg_m,
  agg_q   = agg_q,
  do_post_covid_recalibration = TRUE,
  user_hyper_pre  = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL),
  user_hyper_post = list(r_impute = NULL, Lproxy = NULL, p_AR = NULL, L_midas = NULL),
  verbose = TRUE
)

df_M1 <- list_to_df_nowcast(pseudo_realtime_raw$M1, "M1")
df_M2 <- list_to_df_nowcast(pseudo_realtime_raw$M2, "M2")
df_M3 <- list_to_df_nowcast(pseudo_realtime_raw$M3, "M3")

df_pseudoRT_all <- bind_rows(df_M1, df_M2, df_M3)
if (nrow(df_pseudoRT_all) > 0L) {
  df_pseudoRT_all <- df_pseudoRT_all %>%
    arrange(date, month_in_quarter)
}

# ==============================================================================
# 12. SAVE SINGLE-COUNTRY PSEUDO REAL-TIME RESULTS
# ==============================================================================

file_rt <- build_result_filename(
  path_out         = path_country,
  model            = model_name,
  stage            = "rt",
  Size             = Size,
  sel              = sel,
  countries        = country,
  N_m              = N_m,
  N_q              = N_q,
  Lproxy           = pseudo_realtime_raw$hyper_pre$Lproxy,
  L_midas          = pseudo_realtime_raw$hyper_pre$L_midas,
  p_ar             = pseudo_realtime_raw$hyper_pre$p_AR,
  r1               = r_hat,
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
    model_id = "MF_TPRF",
    model    = model_name,
    stage    = "pseudo_realtime",
    Size     = Size,
    sel      = sel,
    country  = country,
    params   = params,
    dates_m  = dates_m,
    dates_q  = dates_q,
    y_q      = y_q,
    pseudo_realtime_raw = pseudo_realtime_raw,
    pseudo_realtime_M1  = df_M1,
    pseudo_realtime_M2  = df_M2,
    pseudo_realtime_M3  = df_M3,
    pseudo_realtime_all = df_pseudoRT_all
  ),
  file = file_rt
)

cat("\nSaved pseudo real-time results to:\n", file_rt, "\n")