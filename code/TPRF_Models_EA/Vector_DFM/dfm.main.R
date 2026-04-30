# ==============================================================================
# Vector DFM Main Script
# Dynamic Factor Model estimated by EM for mixed-frequency GDP nowcasting
# ------------------------------------------------------------------------------
# Full-sample estimation with standardized saving structure
# ==============================================================================

# ==============================================================================
# 0. PATHS
# ==============================================================================

path_main <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

path_data_raw <- file.path(path_main, "data/raw")
path_data_adj <- file.path(path_main, "data/data_TR2")
path_func     <- file.path(path_main, "functions/functions_dfm")
path_utils    <- file.path(path_main, "functions/functions_mat")

path_results  <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/outputs")
path_graph    <- file.path(path_main, "TPRF_Models_EA/Vector_DFM/results/graph")

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
# 2. SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

source(file.path(path_func, "dfm.utils.R"))
source(file.path(path_func, "dfm.prep.R"))
source(file.path(path_func, "dfm.fs.R"))
source(file.path(path_func, "dfm.kalman.R"))
source(file.path(path_func, "dfm.init.R"))
source(file.path(path_func, "dfm.em.R"))
source(file.path(path_func, "dfm.now.R"))

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
    "_CovidM-", as.integer(isTRUE(params$covid_mask_m)),
    "_CovidQ-", as.integer(isTRUE(params$covid_mask_q))
  )
}

get_size_tag_dfm <- function(N_m, N_q) {
  N_tot <- N_m + N_q
  
  if (N_tot <= 25) {
    return("small")
  } else if (N_tot <= 50) {
    return("medium")
  } else {
    return("large")
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
                                      r,
                                      p,
                                      q,
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
      "_cc-", countries,
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

destd_vec <- function(x_std, mu, sd) {
  as.numeric(mu + sd * x_std)
}

destd_mat <- function(X_std, mu, sd) {
  X_std <- as.matrix(X_std)
  stopifnot(length(mu) == ncol(X_std), length(sd) == ncol(X_std))
  sweep(sweep(X_std, 2, sd, `*`), 2, mu, `+`)
}

# ==============================================================================
# 4. USER PARAMETERS
# ==============================================================================
# Hyperparameters selected from data:
#   - r : number of static factors
#   - p : VAR order for factors
#   - q : AR order for idiosyncratic components
#
# Tuning / design parameters fixed ex ante:
#   - Kmax, pmax, qmax : search bounds for model selection
#   - kappa            : numerical regularization in initialization
#   - restr            : mixed-frequency aggregation restriction
#   - max_iter, tol    : EM stopping rules
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
  
  sel_method   = "LASSO",
  n_m          = 20,
  n_q          = 5,
  thr_m        = 0.10,
  thr_q        = 0.85,
  thr_F_test   = 0.01,
  alpha_lasso  = 1,
  
  # Selection bounds
  Kmax         = 10,      # upper bound for factor search
  pmax         = 4,      # upper bound for VAR lag search on factors
  qmax         = 4,      # upper bound for AR lag search on idiosyncratic components
  
  # DFM specification
  kappa        = 1e-2,
  restr        = "stock_flow",   # or "MM"
  
  # EM controls
  max_iter     = 200,
  tol          = 1e-2
)

countries <- c("DE", "FR", "IT", "ES", "NL", "BE", "AT", "PT", "EA")
country   <- "FR"

tag_run    <- build_run_tag_dfm(params)
model_name <- "vector_dfm"
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
series  <- all_countries$data[[country]]$Series

gdp_col     <- all_countries$data[[country]]$target_col
target_name <- all_countries$data[[country]]$target_name

y   <- as.matrix(data[, gdp_col, drop = FALSE])
y_q <- as.matrix(y[!is.na(y)])

NM <- all_countries$data[[country]]$nM
NQ <- all_countries$data[[country]]$nQ
N  <- NM + NQ

Type_m <- all_countries$data[[country]]$TypeM
Type_q <- all_countries$data[[country]]$TypeQ
Type   <- c(Type_m, Type_q)

Freq_m <- all_countries$data[[country]]$freq_m
Freq_q <- all_countries$data[[country]]$freq_q
Freq   <- c(Freq_m, Freq_q)

Unb_m <- all_countries$data[[country]]$unb_m
Unb_q <- all_countries$data[[country]]$unb_q
Unb   <- c(Unb_m, Unb_q)

Class_m <- all_countries$data[[country]]$ClassM
Class_q <- all_countries$data[[country]]$ClassQ
Class   <- c(Class_m, Class_q)

agg_m <- all_countries$data[[country]]$agg_m
agg_q <- all_countries$data[[country]]$agg_q
agg   <- c(agg_m, agg_q)

Size <- get_size_tag_dfm(NM, NQ)

# ==============================================================================
# 7. STANDARDIZATION
# ==============================================================================
# In the DFM, the whole observed panel is standardized, including the GDP series,
# because GDP is part of the joint state-space system.
# ==============================================================================

out_std <- standardize_with_na(data)
X_std   <- out_std$X_std

mu_full <- out_std$mean
sd_full <- out_std$sd

# ==============================================================================
# 8. HYPERPARAMETER SELECTION
# ==============================================================================
# Selected from data:
#   (i)   r : number of factors
#   (ii)  p : VAR order for the factor dynamics
#   (iii) q : AR order for the idiosyncratic components
# ==============================================================================

# ------------------------------------------------------------------
# 8.1 Covariance proxy for factor-number selection
# ------------------------------------------------------------------
cov_proxy_out <- all_purpose_covariance(X_std)

# ------------------------------------------------------------------
# 8.2 Number of factors r via Eigenvalue Ratio
# ------------------------------------------------------------------
ER_proxy_out <- select_num_factors_ER(
  Sigma       = cov_proxy_out$Sigma_tilde,
  Kmax        = params$Kmax
)

r   <- ER_proxy_out$r
eig <- ER_proxy_out$eig

# Optional diagnostics
plot_factor_selection_diag(ER_proxy_out, Kmax = params$Kmax, r_sel = r)

# ------------------------------------------------------------------
# 8.3 Preliminary factor extraction (XP)
# ------------------------------------------------------------------
xp <- estimate_factors_XP(X_std, r = r)

L_hat <- xp$Lambda
F_hat <- xp$F_hat
C_hat <- xp$C_hat

# Optional diagnostics
plot_factors_XP(F_hat, dates = dates_m, ncol = 1)

# ------------------------------------------------------------------
# 8.4 Factor VAR order p via BIC
# ------------------------------------------------------------------
factor_var_ic <- select_p_var_ic(
  F_hat          = F_hat,
  pmax           = params$pmax,
  include_const  = TRUE
)

p <- factor_var_ic$p_BIC

# ------------------------------------------------------------------
# 8.5 Idiosyncratic AR order q via BIC
# ------------------------------------------------------------------
Ehat <- X_std - C_hat

idio_ar_ic <- select_p_ar_ic_global(
  Ehat[, 1:NM, drop = FALSE],   # monthly block only, standard choice
  pmax          = params$qmax,
  include_const = FALSE,
  min_T         = 30
)

q <- idio_ar_ic$p_BIC

# ==============================================================================
# 9. MODEL INITIALIZATION
# ==============================================================================
# The initialization uses the data-driven choices (r, p, q).
# ==============================================================================

Init <- InitialCond(
  X_std    = X_std,
  r        = r,
  p_factor = p,
  q_idio   = q,
  NM       = NM,
  NQ       = NQ,
  restr    = params$restr,
  agg      = agg,
  kappa    = params$kappa
)

# ==============================================================================
# 10. FULL-SAMPLE EM ESTIMATION
# ==============================================================================

fit_dfm <- DFM_EM(
  X        = X_std,
  Init     = Init,
  max_iter = params$max_iter,
  tol      = params$tol
)

# ==============================================================================
# 11. FITTED VALUES ON STANDARDIZED AND ORIGINAL SCALE
# ==============================================================================

n_lags_Q_used <- if (params$restr == "MM") 5L else 3L

fit_obj <- dfm_kalman_fit(
  X_in      = X_std,
  res       = fit_dfm,
  n_lags_Q  = n_lags_Q_used
)

X_fit_std  <- fit_obj$X_fit
X_fit_orig <- destd_mat(X_fit_std, mu_full, sd_full)

gdp_fit_std  <- X_fit_std[,  gdp_col]
gdp_fit_orig <- X_fit_orig[, gdp_col]

# ==============================================================================
# 12. SAVE SINGLE-COUNTRY FULL-SAMPLE RESULTS
# ==============================================================================
# The saved object mirrors the MF-TPRF structure:
#   model_id, model, stage, Size, sel, country, params, dates, metadata,
#   preprocessing, hyper, fit
# ==============================================================================

path_country <- file.path(path_results, country)
dir.create(path_country, recursive = TRUE, showWarnings = FALSE)

file_fit <- build_result_filename_dfm(
  path_out = path_country,
  model    = model_name,
  stage    = "fit",
  Size     = Size,
  sel      = sel,
  countries= country,
  N_m      = NM,
  N_q      = NQ,
  r        = r,
  p        = p,
  q        = q,
  covid_m  = as.integer(isTRUE(params$covid_mask_m)),
  covid_q  = as.integer(isTRUE(params$covid_mask_q)),
  ext      = "rds",
  timestamp= TRUE
)

saveRDS(
  list(
    model_id = "DFM_EM",
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
      N_m   = NM,
      N_q   = NQ,
      N     = N,
      Type  = Type,
      agg_m = agg_m,
      agg_q = agg_q,
      agg   = agg,
      Freq  = Freq,
      Unb   = Unb,
      Class = Class,
      series = series,
      gdp_col = gdp_col,
      target_name = target_name
    ),
    
    preprocessing = list(
      X_raw         = data,
      out_std       = out_std,
      X_std         = X_std,
      cov_proxy_out = cov_proxy_out,
      xp            = xp,
      Ehat          = Ehat
    ),
    
    hyper = list(
      r            = r,
      p            = p,
      q            = q,
      eig          = eig,
      ER_proxy_out = ER_proxy_out,
      factor_var_ic= factor_var_ic,
      idio_ar_ic   = idio_ar_ic,
      selected_from_data = c("r", "p", "q"),
      fixed_tuning = list(
        Kmax     = params$Kmax,
        pmax     = params$pmax,
        qmax     = params$qmax,
        kappa    = params$kappa,
        restr    = params$restr,
        max_iter = params$max_iter,
        tol      = params$tol
      )
    ),
    
    fit = list(
      Init          = Init,
      res           = fit_dfm,
      A             = fit_dfm$A,
      C             = fit_dfm$C,
      Q             = fit_dfm$Q,
      R             = fit_dfm$R,
      Z0            = fit_dfm$Z0,
      V0            = fit_dfm$V0,
      loglik        = fit_dfm$loglik,
      crit          = fit_dfm$crit,
      L_hat         = L_hat,
      F_hat         = F_hat,
      C_hat         = C_hat,
      std_map       = list(mean = mu_full, sd = sd_full),
      fit_obj       = fit_obj,
      n_lags_Q_used = n_lags_Q_used,
      X_fit_std     = X_fit_std,
      X_fit_orig    = X_fit_orig,
      gdp_fit_std   = gdp_fit_std,
      gdp_fit_orig  = gdp_fit_orig
    )
  ),
  file = file_fit
)

cat("\nSaved full-sample DFM results to:\n", file_fit, "\n")









# ==============================================================================
# 13. PSEUDO REAL-TIME NOWCASTING (DFM)
# ==============================================================================

pseudo_realtime_raw <- pseudo_realtime_DFM_EM_reestimate(
  X_full   = data,
  NQ       = NQ,
  params   = params,
  dates_m  = dates_m,
  dates_q  = dates_q,
  Freq     = Freq,
  Unb      = Unb,
  gdp_col  = gdp_col,
  agg      = agg,
  do_post_covid_recalibration = TRUE,
  user_hyper_pre  = list(r = NULL, p = NULL, q = NULL),
  user_hyper_post = list(r = NULL, p = NULL, q = NULL),
  max_iter_em = params$max_iter,
  tol_em      = params$tol,
  pmax        = params$pmax,
  qmax        = params$qmax,
  min_est_T   = 24,
  verbose     = TRUE
)

# ------------------------------------------------------------------------------
# Build compact data frames exactly as in MF-TPRF style
# Here we use ORIGINAL-SCALE nowcasts for the standard plotting/evaluation object
# ------------------------------------------------------------------------------

df_M1 <- list_to_df_nowcast(pseudo_realtime_raw$M1_orig, "M1")
df_M2 <- list_to_df_nowcast(pseudo_realtime_raw$M2_orig, "M2")
df_M3 <- list_to_df_nowcast(pseudo_realtime_raw$M3_orig, "M3")

df_pseudoRT_all <- bind_rows(df_M1, df_M2, df_M3)
if (nrow(df_pseudoRT_all) > 0L) {
  df_pseudoRT_all <- df_pseudoRT_all %>%
    arrange(date, month_in_quarter)
}

# ------------------------------------------------------------------------------
# Optional long data frame containing both standardized and original-scale values
# ------------------------------------------------------------------------------

rolling_to_df <- function(now_obj) {
  make_df <- function(x, label, scale) {
    if (length(x) == 0) return(NULL)
    data.frame(
      date    = as.Date(names(x)),
      value   = as.numeric(x),
      vintage = label,
      scale   = scale,
      row.names = NULL
    )
  }
  
  bind_rows(
    make_df(now_obj$M1_std,  "M1", "std"),
    make_df(now_obj$M2_std,  "M2", "std"),
    make_df(now_obj$M3_std,  "M3", "std"),
    make_df(now_obj$M1_orig, "M1", "orig"),
    make_df(now_obj$M2_orig, "M2", "orig"),
    make_df(now_obj$M3_orig, "M3", "orig")
  )
}

df_pseudoRT_long <- rolling_to_df(pseudo_realtime_raw)

# ==============================================================================
# 14. SAVE SINGLE-COUNTRY PSEUDO REAL-TIME RESULTS (DFM)
# ==============================================================================

path_country <- file.path(path_results, country)
dir.create(path_country, recursive = TRUE, showWarnings = FALSE)

file_rt <- build_result_filename_dfm(
  path_out = path_country,
  model    = model_name,
  stage    = "rt",
  Size     = Size,
  sel      = sel,
  countries= country,
  N_m      = NM,
  N_q      = NQ,
  r        = pseudo_realtime_raw$hyper_pre$r,
  p        = pseudo_realtime_raw$hyper_pre$p,
  q        = pseudo_realtime_raw$hyper_pre$q,
  covid_m  = as.integer(isTRUE(params$covid_mask_m)),
  covid_q  = as.integer(isTRUE(params$covid_mask_q)),
  ext      = "rds",
  timestamp= TRUE
)

saveRDS(
  list(
    model_id = "DFM_EM",
    model    = model_name,
    stage    = "pseudo_realtime",
    Size     = Size,
    sel      = sel,
    country  = country,
    params   = params,
    dates_m  = dates_m,
    dates_q  = dates_q,
    y_q      = y_q,
    
    metadata = list(
      N_m        = NM,
      N_q        = NQ,
      N          = NM + NQ,
      Freq       = Freq,
      Unb        = Unb,
      agg        = agg,
      gdp_col    = gdp_col,
      target_name= target_name
    ),
    
    pseudo_realtime_raw  = pseudo_realtime_raw,
    pseudo_realtime_M1   = df_M1,
    pseudo_realtime_M2   = df_M2,
    pseudo_realtime_M3   = df_M3,
    pseudo_realtime_all  = df_pseudoRT_all,
    pseudo_realtime_long = df_pseudoRT_long
  ),
  file = file_rt
)

cat("\nSaved pseudo real-time DFM results to:\n", file_rt, "\n")
