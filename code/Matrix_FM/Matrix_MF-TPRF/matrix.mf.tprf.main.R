# ==============================================================================
#             Matrix Mixed Frequency TPRF and the EM Algorithm: 
#         A Nowcasting Framework for Mixed Frequency Data in Euro Area
#                               - Main Script - 
# ==============================================================================
# This script estimates a Matrix MF-TPRF using the EM-algorithm to impute missing
# and performs a real-time nowcasting exercise for GDP growth in the Euro Area's 
# main economies:
# Germany (DE), France (FR), Italy (IT), and Spain (ES). 
# It includes:
#   1. Data preparation and harmonization across countries
#   2. Construction of a 3D data tensor
#   3. Estimation of the number of latent factors and model lags
#   4. Model initialization and parameter estimation using the EM algorithm on 
#      the whole dataset
#   5. Pseudo real-time rolling nowcasting of GDP
#   6. Saving results for further use
#
# ==============================================================================
# Author      : Davide Delfino
# Institution : Alma Mater Studiorum - University of Bologna
# Dataset     : Barigozzi & Lissona (2024), EA-MD-QD
#
# Disclaimer  : ---
#
# ==============================================================================
# Notes       : Other EA countries such as BE, NL, AT, PT, IR, GR and different 
#               set of Monthly and Quarterly variables can be added
# ==============================================================================
# Script Type : Data Preparation / Estimation / Nowcastinf / Matrix MF-TPRF
# ==============================================================================


# ==============================================================================
# 0. SET WORKING DIRECTORY & PATHS
# ==============================================================================

path_main    <- "C:/Users/david/Desktop/Paper/Matrix_vs_Vector_TS/code"
setwd(path_main)

# Paths: data, functions, results, graphs
path_data    <- file.path(path_main, "data")
path_func    <- file.path(path_main, "Matrix_FM/functions")
path_results <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/outputs")
path_graph   <- file.path(path_main, "Matrix_FM/Matrix_MF-TPRF/results/graph")



# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================

## Data Handling
library(tidyverse)
library(lubridate)
library(abind)
library(dplyr)

## Time Series
library(tseries)
library(zoo)
library(MASS)
library(fBasics)
library(vars)

## test
library(sandwich)
library(lmtest)
library(car)

## I/O
library(readxl)

## Plotting
library(ggplot2)

## Resolve conflicts
library(conflicted)
conflict_prefer("select",  "dplyr", quiet = TRUE)
conflict_prefer("filter",  "dplyr", quiet = TRUE)



# ==============================================================================
# 2. SOURCE USER-DEFINED FUNCTIONS
# ==============================================================================

source(file.path(path_func, "matrix.mf.tprf.prep.R"))
source(file.path(path_func, "matrix.mf.tprf.R"))
source(file.path(path_func, "matrix.mf.tprf.em.R"))
source(file.path(path_func, "matrix.mf.tprf.fs.R"))
source(file.path(path_func, "matrix.mf.tprf.now.R"))
source(file.path(path_func, "matrix.mf.tprf.utils.R"))



# ==============================================================================
# 3.1 PARAMETERS & COUNTRY LIST & COUNTRY OF INTEREST SELECTION
# ==============================================================================

# Users's parameters
params <- list(
  start_est    = as.Date("2000-04-01"),    # Accounting from missingenss in growth rates wrt 1999
  start_eval   = as.Date("2017-01-01"),
  end_eval     = as.Date("2025-10-01"),
  covid_start  = as.Date("2020-03-01"),
  covid_end    = as.Date("2021-07-01"),
  covid_mask_m = TRUE,                    # Boolean: if empty == TRUE
  covid_mask_q = FALSE,                    # Boolean: if empty == TRUE
  target       = "GDP",                    # quarterly target variable
  target_cc    = "EA",           # Aggregate from where to extract the proxy
  sel_method   = "none",         # "none" | "corr_threshold" | "F-Test" | "LASSO"
  
  # Variable selection parameters
  n_m        = 38,
  n_q        = 20,
  thr_m      = 0.10,
  thr_q      = 0.85,
  thr_F_test = 0.01,
  alpha_lasso  = 1,                        # 1=LASSO, 0=Ridge, (0,1)=Elastic Net
  
  
  # Factor and proxy selection parameters
  Kmax = 10,                            # max number of factors
  Zmax = 10,                            # max number of proxy
  
  
  # MF-TPRF parameters
  p_AR        = 1,                         # order AR(p) in U-MIDAS regression of y
  kmax        = c(1,5),                   # Max number of Factors c(p_1, p_2)
  Lmax        = 3,                         # Max numbers of lags in the Factors in U-MIDAS
  Robust_F    = TRUE,                      # Robust F test in the first step
  alpha       = 0.1,
  robust_type = "NW",                      # White--> HC | NW(Newey West)--> HAC
  nw_lag      = 1                          # Lag in the error autocorrelation in step 1
)

# Country List
countries <- c("DE", "FR", "IT", "ES", "EA")

# ==============================================================================
# 4. PREPARE DATA FOR ALL COUNTRIES
# ==============================================================================
all_countries <- prepare_all_countries(
  countries       = countries,
  params          = params,
  path            = path_data,
  covid_mask_m    = params$covid_mask_m,
  covid_mask_q    = params$covid_mask_q
)

# ==============================================================================
# 5. PREPARE TENSOR
# ==============================================================================

tensor <- build_tensor(all_countries, params, "intersection")

Y <- tensor$Y
W <- tensor$W

nan_percent_Y(Y)

# ==============================================================================
# 5. TARGET AND PREDICTORS IDENTIFICATION
# ==============================================================================

gdp_col <- tensor$target_col

y   <- Y[ , ,gdp_col, drop = FALSE]     # target
mask_gdp <- !is.na(y[ , , 1])
idx_q <- which(rowSums(mask_gdp) > 0)
y_q <- y[idx_q, , , drop = TRUE]       # [T_trimestri x Paesi]
y_q_centered <- scale(y_q,center = TRUE, scale = FALSE)


X     <- Y[, ,-gdp_col, drop = FALSE]    # predictors
W_x   <- W[, ,-gdp_col, drop = FALSE]    # predictors

# Proxy GDP EA
y_EA_q    <- all_countries$proxy$Y_q
dates_EA  <- as.Date(rownames(all_countries$proxy$Y_q))

# y_EA_q_centered <- scale(y_EA_q,center = TRUE, scale = FALSE)

y_q_all <- cbind(y_EA_q,y_q)
# y_q_centered_all <- cbind(y_EA_q_centered,y_q_centered)

# ==============================================================================
# 6. METADATA EXTRACTION
# ==============================================================================

## Dimensions
N_m     <- tensor$n_M
nQ_tot  <- tensor$n_Q

# Quarterly series (all)
Q_series_all <- tensor$vars[(N_m + 1):(N_m + nQ_tot)]

# Remove target from quarterly list
q_cols <- which(!grepl(params$target, Q_series_all, ignore.case = TRUE))
N_q    <- length(q_cols)

# Total number of predictors
N <- N_m + N_q

# PROBLEMA CON LA SELEZIONE DELLE VARIABILI
agg_m   <- tensor$agg_M
agg_q   <- tensor$agg_Q[, q_cols, drop = FALSE]
agg     <- cbind(agg_m, agg_q)

Type_m  <- tensor$TypeM
Type_q  <- tensor$TypeQ[,q_cols, drop = FALSE]
Type    <- cbind(Type_m, Type_q)

Freq_m  <- tensor$freq_M
Freq_q  <- tensor$freq_Q[,q_cols, drop = FALSE]
Freq    <- cbind(Freq_m, Freq_q)

Unb_m   <- tensor$unb_M
Unb_q   <- tensor$unb_Q[,q_cols, drop = FALSE]
Unb     <- cbind(Unb_m, Unb_q)

Class_m <- tensor$ClassM
Class_q <- tensor$ClassQ[,q_cols, drop = FALSE]
Class   <- cbind(Class_m, Class_q)

# --- VARS: costruisci la lista predittori coerente con X_tens_full ---
vars_m <- tensor$vars[tensor$idx_M]  # oppure tensor$vars_M se ce l'hai già
# trimestrali (base_Q) SENZA target
vars_q_all <- tensor$vars[(N_m + 1):(N_m + nQ_tot)]
vars_q_pred <- vars_q_all[q_cols]

vars_pred <- c(vars_m, vars_q_pred)

# tensor_info *coerente con i predittori*
tensor_info_pred <- list(
  n_M       = N_m,
  n_Q       = N_q,
  countries = tensor$countries,
  vars      = vars_pred
)


# ==============================================================================
# 7. STANDARDIZATION OF PREDICTORS
# ==============================================================================

out_std <- standardize_mat_with_na(X)
X_std   <- out_std$X_scaled


# ==============================================================================
# 9. EM Algorithm
# ==============================================================================

# 1) Tensor → matrice osservata
flat <- flatten_tensor_to_matrix(
  X_tens   = X_std,   # [T x Paesi x (N_m + N_q)] predittori STANDARDIZZATI
  W_tens   = W_x,
  N_m      = N_m,
  N_q      = N_q,
  agg_q    = agg_q
)

X_obs_mat  <- flat$X_mat   # con attributi N_m, N_q, countries, var_names
W_mat      <- flat$W_mat
agg_q_glob <- flat$agg_q_global

# 2) Inizializzazione Cen & Lam sul TENSORE
init <- init_CL_ER(
  X_std = X_std,   # stesso tensore [T x Paesi x (N_m + N_q)]
  W     = W_x,
  kmax  = params$kmax
)

X_init <- init$X_init   # tensore [T x Paesi x (N_m + N_q)]
r      <- init$r        # c(r1_hat, r2_hat)
r_vec  <- max(r)

# 3) Flatten dell'inizializzato: stesso ordine/nome colonne
X_init_mat <- flatten_tensor_init(
  X_tens   = X_init,
  N_m      = N_m,
  N_q      = N_q
)


################################################################################
# Find the number of predictors using xiong and pelger
X_init_mat <- flatten_tensor_init(
  X_tens   = X_std,
  N_m      = N_m,
  N_q      = N_q
)

# Initialization using Xiong and Pelger "all Purposes Estimator"
init <- init_XP_ER(X_init_mat)

X_init_mat <- init$X_init
r_vec      <- init$r

################################################################################

# 4) Matrici A_i globali
A_global <- A_list(
  X_na  = X_obs_mat,
  N_q   = flat$N_q_global,
  agg_q = agg_q_glob
)

# 5) EM vettoriale
em_out <- EM_algorithm(
  X_init   = X_init_mat,
  X_obs    = X_obs_mat,
  A_list   = A_global,
  r        = r_vec,
  max_iter = 100,
  tol      = 1e-4
)

X_em <- em_out$X_completed   # [T x N_global], con rownames/colnames ma senza attr

# 6) Copia attributi da X_obs_mat a X_em
for (nm in c("N_m", "N_q", "countries", "var_names")) {
  attr(X_em, nm) <- attr(X_obs_mat, nm)
}

# 7) Ritorno al tensore
X_em_tens <- unflatten_matrix_to_tensor(X_em)

X_m_em <- X_em[,1:flat$N_m_global]
X_q_em <- X_em[,(flat$N_m_global+1):flat$N_global]

# Vettore globale per le mensili (stesso ordine di X_m_em)
agg_m_global <- as.vector(t(agg_m))   # length = P1 * N_m

# Vettore globale per le trimestrali (stesso ordine di X_q_em)
agg_q_global <- as.vector(t(agg_q))   # length = P1 * N_q

X_mq_em <- agg_mq(X_m_em, agg_m_global)
X_qq_em <- agg_qq(X_q_em, agg_q_global)

X_em_agg <- cbind(X_mq_em, X_qq_em)

# copia attributi da X_em a X_em_agg
for (nm in c("N_m", "N_q", "countries", "var_names")) {
  attr(X_em_agg, nm) <- attr(X_em, nm)
}

# 7) Ritorno al tensore trimestrale
X_em_agg_tens <- unflatten_matrix_to_tensor(X_em_agg)

# ==============================================================================
# 8. FIND THE NUMBER OF PROXY AND U-MIDAS LAGS 
# ==============================================================================

# IC. Kraemer & Sugiyama
pls.object <- select_L_autoproxy_3prf(X_em_agg, y_EA_q, Zmax = params$Zmax)
Lproxy     <- pls.object$L_opt

# Number of U-MIDAS LAG
lag_sel <- choose_UMIDAS_lag_tensor_MF(X_em_agg_tens, X_em_tens, y_EA_q,
                                       r,
                                       Lmax   = params$Lmax,
                                       Lproxy = Lproxy,
                                       p_AR   = params$p_AR
)

L_midas <- lag_sel$lag_BIC


################################################################################

# ==============================================================================
# 9. Center EM output to don't include the constant in Tensor-MF-TPRF
# ==============================================================================

# X_em_agg_tens_centered_out <- center_Y(X_em_agg_tens)
# X_em_tens_centered_out     <- center_Y(X_em_tens)

# X_em_agg_tens_centered <- X_em_agg_tens_centered_out$Y_centered
# X_em_tens_centered <- X_em_tens_centered_out$Y_centered


# ==============================================================================
# 9. MF-TPRF
# ==============================================================================

Tensor_MF_TPRF_out <- Tensor_MF_TPRF(
  X_lf        = X_em_agg_tens,
  X_hf        = X_em_tens,
  Y_q_all     = y_q_all,
  proxy_name  = "EA",
  Lproxy      = Lproxy,
  L_midas     = L_midas,
  p_AR        = params$p_AR,
  r           = r
)


# ==============================================================================
# 14. SAVE ALL RESULTS (TENSOR VERSION)
# ==============================================================================

sanitize_label <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "-", x)   # sostituisce spazi, virgole, ecc. con '-'
  x <- gsub("-+", "-", x)              # evita multipli '-' consecutivi
  x <- gsub("(^-)|(-$)", "", x)        # toglie '-' iniziali/finali
  x
}

## 14.1 Cartella di output (già è globale, non per country)
path_tensor <- path_results
if (!dir.exists(path_tensor)) {
  dir.create(path_tensor, recursive = TRUE)
}

## 14.2 Helpers per label "pulite"
sel_method    <- sanitize_label(params$sel_method)
robust_F_lab  <- if (isTRUE(params$Robust_F)) "robustF1" else "robustF0"
robust_type   <- sanitize_label(params$robust_type)
covid_m_lab   <- if (isTRUE(params$covid_mask_m)) "CvM1" else "CvM0"
covid_q_lab   <- if (isTRUE(params$covid_mask_q)) "CvQ1" else "CvQ0"
target_cc_lab <- sanitize_label(params$target_cc)

# lista paesi compressa nel nome file
countries_lab <- sanitize_label(paste(countries, collapse = "-"))

# rank tensoriale
r1 <- r[1]
r2 <- r[2]

## 14.3 Build filename (NESSUN country singolo, ma "tensor")
file_out <- file.path(
  path_tensor,
  paste0(
    "T-MF_TPRF_tensor",
    "_cc-",      countries_lab,
    "_sel-",     sel_method,
    "_Lproxy-",  Lproxy,
    "_Lmidas-",  L_midas,
    "_Nm-",      N_m,
    "_Nq-",      N_q,
    "_pAR-",     params$p_AR,
    "_",         robust_F_lab,
    "_rtype-",   robust_type,
    "_",         covid_m_lab,
    "_",         covid_q_lab,
    "_tgtcc-",   target_cc_lab,
    "_r1-",      r1,
    "_r2-",      r2,
    "_",         format(params$start_est, "%Y-%m"),
    "_to_",      format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

## 14.4 Save everything into one RDS
saveRDS(
  list(
    # Info generali
    params      = params,
    countries   = countries,
    tensor      = tensor,
    
    # Risultati EM
    A_list_global = A_global,     # matrice A_list(...) che hai definito sopra
    X_init        = X_init,
    X_em          = X_em,
    X_em_tens     = X_em_tens,
    X_em_agg      = X_em_agg,
    X_em_agg_tens = X_em_agg_tens,
    F_EM          = em_out$F,
    Phi_EM        = em_out$Phi,
    EM_iter       = em_out$iterations,
    EM_diffs      = em_out$diffs,
    
    # Selezione fattori e lag
    r_hat       = r,
    Lproxy      = Lproxy,
    lag_sel     = L_midas,
    
    # Tensor MF-TPRF Results
    Tensor_MF_TPRF = Tensor_MF_TPRF_out
  ),
  file_out
)

cat("\n*** SAVED EM + TENSOR MF-TPRF RESULTS TO ***\n", file_out, "\n")



print(path_tensor)
print(file_out)
str(file_out)





# ==============================================================================
# 15. PSEUDO REAL-TIME FORECASTING EXERCISE (WITH RELEASE SCHEDULE)
# ==============================================================================

# Questo esercizio implementa il vero pseudo real-time MF-3PRF + EM,
# rispettando calendario dei rilasci, unbalancedness dati, revisione dei fattori
# e disponibilità delle variabili ad alta e bassa frequenza.
# X_full: matrice mensile flatten con T_m x N (predittori, SENZA GDP)
#         puoi usare X_obs_mat che hai costruito con flatten_tensor_to_matrix



# ==============================================================================
# PSEUDO REAL-TIME T-MF-TPRF: RUN + TIDY + SAVE (COUNTRY-SPECIFIC)
# ==============================================================================

library(dplyr)

# ------------------------------------------------------------
# 2) Call rolling / pseudo real-time
# ------------------------------------------------------------
roll_tensor <- pseudo_realtime_Tensor_MF_TPRF(
  X_tens_full = X,
  Y_q_all     = y_q_all,          # [T_q x (1+P)]  (EA + paesi)
  proxy_name  = "EA",
  params      = params,
  dates_m     = tensor$dates,
  dates_q     = dates_EA,
  Unb         = Unb,
  agg_m       = agg_m,
  agg_q       = agg_q,
  tensor_info = tensor_info_pred
)

# ------------------------------------------------------------
# 3) Utility: list(country -> named numeric date->value) -> df
# ------------------------------------------------------------
list_to_df_country <- function(lst_by_country, tag) {
  # lst_by_country: list(country -> named numeric) oppure list(country -> list(date->value))
  if (is.null(lst_by_country) || length(lst_by_country) == 0L) {
    return(data.frame(
      date             = as.Date(character()),
      country          = character(),
      nowcast          = numeric(),
      month_in_quarter = character(),
      row.names        = NULL
    ))
  }
  
  # se per qualche motivo tagli fuori eventuali NULL
  countries <- names(lst_by_country)
  if (is.null(countries)) {
    stop("lst_by_country non ha names(): mi aspetto list con nomi paese.")
  }
  
  out <- lapply(countries, function(cc) {
    lst <- lst_by_country[[cc]]
    if (is.null(lst) || length(lst) == 0L) return(NULL)
    
    vals  <- as.numeric(unlist(lst))
    dates <- as.Date(names(lst))
    
    if (any(is.na(dates))) {
      stop("Date non parsabili in lst_by_country[['", cc, "']]. ",
           "Mi aspetto names in formato 'YYYY-MM-DD'.")
    }
    
    data.frame(
      date             = dates,
      country          = cc,
      nowcast          = vals,
      month_in_quarter = tag,
      row.names        = NULL
    )
  })
  
  dplyr::bind_rows(out)
}

df_M1 <- list_to_df_country(roll_tensor$M1, "M1")
df_M2 <- list_to_df_country(roll_tensor$M2, "M2")
df_M3 <- list_to_df_country(roll_tensor$M3, "M3")

df_pseudoRT_all <- bind_rows(df_M1, df_M2, df_M3) |>
  mutate(
    month_in_quarter = factor(month_in_quarter, levels = c("M1", "M2", "M3"))
  ) |>
  arrange(country, date, month_in_quarter)

# ------------------------------------------------------------
# 4) Info dal tensore per naming / Nq coerente (GDP escluso)
# ------------------------------------------------------------
N_m_tensor   <- tensor$n_M
nQ_tot       <- tensor$n_Q
Q_series_all <- tensor$vars[(N_m_tensor + 1):(N_m_tensor + nQ_tot)]
q_cols       <- which(!grepl(params$target, Q_series_all, ignore.case = TRUE))
N_q_tensor   <- length(q_cols)

# Paesi valutati: togli "EA" se presente per qualsiasi motivo
countries_eval <- tensor$countries
countries_eval <- setdiff(countries_eval, "EA")

cc_str <- paste(countries_eval, collapse = "-")

# ------------------------------------------------------------
# 5) Output path + filename
# ------------------------------------------------------------
path_tensor_rt <- path_results
if (!dir.exists(path_tensor_rt)) dir.create(path_tensor_rt, recursive = TRUE)

file_out_rt <- file.path(
  path_tensor_rt,
  paste0(
    "T-MF_TPRF_pseudoRT_tensor_cc-", cc_str,
    "_sel-",    params$sel_method,
    "_Lproxy-", roll_tensor$Lproxy_fix,
    "_Lmidas-", roll_tensor$L_midas_fix,
    "_Nm-",     N_m_tensor,
    "_Nq-",     N_q_tensor,
    "_pAR-",    params$p_AR,
    "_CvM-",    params$covid_mask_m,
    "_CvQ-",    params$covid_mask_q,
    "_",        format(params$start_eval, "%Y-%m"),
    "_to_",     format(params$end_eval, "%Y-%m"),
    ".rds"
  )
)

cat("file_out_rt =\n", file_out_rt, "\n")

# ------------------------------------------------------------
# 6) Save (overwrite pulito)
# ------------------------------------------------------------
if (file.exists(file_out_rt)) file.remove(file_out_rt)

saveRDS(
  list(
    # Meta
    params     = params,
    countries  = countries_eval,
    proxy_name = "EA",
    dates_m    = tensor$dates,
    dates_q    = dates_EA,
    Y_q_all    = y_q_all,        # truth/proxy completo
    
    # raw output
    pseudo_rt_raw = roll_tensor,
    
    # tidy
    pseudo_rt_M1  = df_M1,
    pseudo_rt_M2  = df_M2,
    pseudo_rt_M3  = df_M3,
    pseudo_rt_all = df_pseudoRT_all
  ),
  file = file_out_rt
)

cat("\n*** TENSOR PSEUDO REAL-TIME RESULTS SAVED TO ***\n", file_out_rt, "\n")
