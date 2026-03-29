# ==============================================================================
# NOWCAST COMPARISON: BI vs Vec MF-TPRF vs Mat MF-TPRF
# ------------------------------------------------------------------------------
# Design choice:
#   We assume perfect completion of the missing entries induced by publication
#   delays. Hence, at each vintage M1, M2, M3, the feasible models are run on
#   the oracle-completed panels available up to that vintage.
#
# Consequences:
#   - No Cen-Lam or Xiong-Pelger completion step is used here
#   - No factor-rank selection step is used
#   - The BI benchmark is built from the true predictive factor f_true
#   - Standardization is done vintage-by-vintage on the completed HF panel
#   - A standard quarterly TPRF benchmark is also included
#   - The quarterly TPRF forecast horizon is adjusted over the current quarter:
#       M1 -> h = 2 if the previous quarter is not yet complete
#       M2 -> h = 1 once the previous quarter is complete
#       M3 -> h = 1
# ==============================================================================


# ==============================================================================
# ESTIMATION PARAMETERS
# Purpose:
#   Store only the parameters needed in the Monte Carlo nowcast comparison.
# Notes:
#   - No rank selection
#   - True matrix rank is imposed directly
#   - The standard quarterly TPRF benchmark is optional
# ==============================================================================
make_nowcast_estimation_params <- function(
    Lproxy_vec = 1,
    L_midas_vec = 1,
    p_AR_vec = 1,
    Lproxy_mat = 1,
    L_midas_mat = 1,
    p_AR_mat = 1,
    r_mat_true = c(1, 1),
    proxy_name = "EA",
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8,
    bi_intercept = TRUE,
    L_tprf = 1,
    use_tprf_standard = TRUE
) {
  stopifnot(length(r_mat_true) == 2)
  
  list(
    Lproxy_vec = Lproxy_vec,
    L_midas_vec = L_midas_vec,
    p_AR_vec = p_AR_vec,
    
    Lproxy_mat = Lproxy_mat,
    L_midas_mat = L_midas_mat,
    p_AR_mat = p_AR_mat,
    r_mat_true = r_mat_true,
    
    proxy_name = proxy_name,
    standardize_proxy = standardize_proxy,
    orthonormalize_each_iter = orthonormalize_each_iter,
    orthonormalize_final_Z = orthonormalize_final_Z,
    ils_maxit = ils_maxit,
    ils_tol = ils_tol,
    
    bi_intercept = bi_intercept,
    L_tprf = L_tprf,
    use_tprf_standard = use_tprf_standard
  )
}


# ==============================================================================
# HELPER: TARGET HORIZON FOR TPRF
# Purpose:
#   Determine the quarterly TPRF forecast horizon at each vintage.
# Notes:
#   - If the previous quarter is not yet complete under publication delays,
#     the standard quarterly TPRF must use h = 2
#   - Once the previous quarter becomes complete, it switches to h = 1
#   - With max monthly delay = 2:
#       M1 -> h = 2
#       M2 -> h = 1
#       M3 -> h = 1
# ==============================================================================
get_tprf_horizon_from_vintage <- function(sim_obj, vintage = c("M1", "M2", "M3")) {
  vintage <- match.arg(vintage)
  
  max_delay <- max(sim_obj$params$monthly_release_delay)
  
  vintage_month <- switch(
    vintage,
    M1 = 1L,
    M2 = 2L,
    M3 = 3L
  )
  
  if (vintage_month < max_delay) {
    h <- 2L
  } else {
    h <- 1L
  }
  
  h
}


# ==============================================================================
# HELPER: PREPARE QUARTERLY INPUTS FOR STANDARD TPRF
# Purpose:
#   Prepare the completed quarterly panel used by the standard quarterly TPRF.
# Notes:
#   - Uses only the fully aggregated quarterly panel
#   - The end of the usable quarterly training sample depends on the
#     vintage-specific horizon h
#   - If h = 2, training ends at q_star - 2
#   - If h = 1, training ends at q_star - 1
# ==============================================================================
prepare_tprf_quarterly_inputs <- function(sim_obj, vintage = c("M1", "M2", "M3")) {
  vintage <- match.arg(vintage)
  
  iset <- sim_obj$info_sets[[vintage]]
  
  X_lf_raw <- sim_obj$truth$X_quarterly_lf_raw
  y_q      <- sim_obj$truth$y_quarterly
  
  q_star <- iset$meta$q_star
  h      <- get_tprf_horizon_from_vintage(sim_obj, vintage)
  
  # Last quarter used in the TPRF training sample:
  # h = 1 -> use data up to q_star - 1
  # h = 2 -> use data up to q_star - 2
  q_train_end <- q_star - h
  
  if (q_train_end < 2) {
    stop("Not enough quarterly observations for TPRF standard at this vintage.")
  }
  
  X_train_raw <- X_lf_raw[1:q_train_end, , , drop = FALSE]
  
  # Vintage-safe standardization fitted on the available quarterly panel
  std_obj <- fit_standardization_tensor(X_train_raw)
  X_train_std <- apply_standardization_tensor(
    X_train_raw,
    center_mat = std_obj$center,
    scale_mat  = std_obj$scale
  )
  
  X_train_vec_std <- vectorize_panel(X_train_std)
  
  list(
    X_train_vec_std = X_train_vec_std,
    y_q             = y_q,
    q_star          = q_star,
    q_train_end     = q_train_end,
    h               = h,
    standardization = std_obj
  )
}


# ==============================================================================
# HELPER: ALIGN RESULT COLUMNS
# Purpose:
#   Harmonize output tables across models before row-binding.
# ==============================================================================
align_result_columns <- function(df) {
  if (is.null(df)) return(NULL)
  
  wanted <- c(
    "unit", "model", "vintage",
    "horizon_h", "q_train_end", "q_target",
    "y_hat", "y_true", "error"
  )
  
  missing <- setdiff(wanted, names(df))
  for (nm in missing) {
    df[[nm]] <- NA
  }
  
  df <- df[, wanted, drop = FALSE]
  df
}


# ==============================================================================
# HELPER: PREPARE COMPLETED INPUTS FOR A GIVEN VINTAGE
# Purpose:
#   Extract oracle-completed panels up to the chosen vintage and standardize
#   them in a vintage-safe way.
# Notes:
#   - HF completed panel differs across M1, M2, M3
#   - LF completed training panel is aligned with the completed quarters
#   - Standardization is fitted on the completed HF panel available at the
#     current vintage and then applied to both HF and LF panels
# ==============================================================================
prepare_completed_inputs_for_vintage <- function(sim_obj, vintage = c("M1", "M2", "M3")) {
  vintage <- match.arg(vintage)
  
  iset <- sim_obj$info_sets[[vintage]]
  
  # ---------------------------------------------------------------------------
  # 1. Raw oracle-completed panels up to the current vintage
  # ---------------------------------------------------------------------------
  X_hf_raw <- iset$oracle$X_hf_true_raw
  X_lf_raw <- iset$oracle$X_lf_true_raw
  
  stopifnot(length(dim(X_hf_raw)) == 3)
  stopifnot(length(dim(X_lf_raw)) == 3)
  
  # ---------------------------------------------------------------------------
  # 2. Vintage-safe standardization fitted on the completed HF panel
  # ---------------------------------------------------------------------------
  std_obj <- fit_standardization_tensor(X_hf_raw)
  
  X_hf_std <- apply_standardization_tensor(
    X_hf_raw,
    center_mat = std_obj$center,
    scale_mat  = std_obj$scale
  )
  
  X_lf_std <- apply_standardization_tensor(
    X_lf_raw,
    center_mat = std_obj$center,
    scale_mat  = std_obj$scale
  )
  
  # ---------------------------------------------------------------------------
  # 3. Vectorized versions for Vec MF-TPRF
  # ---------------------------------------------------------------------------
  X_hf_vec_std <- vectorize_panel(X_hf_std)
  X_lf_vec_std <- vectorize_panel(X_lf_std)
  
  list(
    matrix = list(
      X_hf_std = X_hf_std,
      X_lf_std = X_lf_std
    ),
    vector = list(
      X_hf_vec_std = X_hf_vec_std,
      X_lf_vec_std = X_lf_vec_std
    ),
    standardization = std_obj,
    oracle = iset$oracle,
    target = iset$target,
    meta = iset$meta
  )
}


# ==============================================================================
# HELPER: SAFE LAST NON-MISSING VALUE
# ==============================================================================
safe_last_non_na <- function(x) {
  x_ok <- x[!is.na(x)]
  if (length(x_ok) == 0) return(NA_real_)
  tail(x_ok, 1)
}


# ==============================================================================
# BI BENCHMARK: ALL UNITS, ONE VINTAGE
# Purpose:
#   Construct the infeasible best nowcast using the true predictive factor.
# Method:
#   For each unit i, estimate a quarterly regression of y_{i,q} on the oracle
#   BI regressors built from the true factor, and predict the target quarter
#   using the vintage-specific oracle regressor row.
# ==============================================================================
run_bi_wrapper <- function(
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    bi_intercept = TRUE
) {
  vintage <- match.arg(vintage)
  
  iset <- sim_obj$info_sets[[vintage]]
  
  X_bi_train  <- iset$oracle$X_bi_train
  x_bi_target <- iset$oracle$x_bi_target
  y_train     <- iset$target$y_train
  y_target    <- iset$target$y_target
  
  stopifnot(is.matrix(X_bi_train))
  stopifnot(is.numeric(x_bi_target))
  stopifnot(is.matrix(y_train))
  stopifnot(is.matrix(y_target))
  
  p1 <- ncol(y_train)
  out <- vector("list", p1)
  
  for (i in 1:p1) {
    y_i <- as.numeric(y_train[, i])
    
    keep <- complete.cases(X_bi_train) & !is.na(y_i)
    X_i  <- X_bi_train[keep, , drop = FALSE]
    y_i  <- y_i[keep]
    
    if (nrow(X_i) <= ncol(X_i)) {
      stop(sprintf("BI regression ill-posed for unit %d at vintage %s.", i, vintage))
    }
    
    if (bi_intercept) {
      X_i_fit      <- cbind(1, X_i)
      x_target_fit <- c(1, x_bi_target)
    } else {
      X_i_fit      <- X_i
      x_target_fit <- x_bi_target
    }
    
    beta_hat <- solve(crossprod(X_i_fit), crossprod(X_i_fit, y_i))
    y_hat_i  <- as.numeric(drop(x_target_fit %*% beta_hat))
    y_true_i <- as.numeric(y_target[1, i])
    
    out[[i]] <- data.frame(
      unit    = i,
      model   = "BI",
      vintage = vintage,
      y_hat   = y_hat_i,
      y_true  = y_true_i,
      error   = y_hat_i - y_true_i
    )
  }
  
  do.call(rbind, out)
}


# ==============================================================================
# STANDARD QUARTERLY TPRF: ONE UNIT, ONE VINTAGE, PREPARED INPUTS
# Purpose:
#   Run the standard quarterly TPRF for one unit at a given vintage using
#   quarterly inputs that have already been prepared.
# Notes:
#   - TPRF_h() is assumed to be already defined
#   - The forecast horizon h is chosen upstream and stored in prep
#   - The target quarter remains the current quarter q_star
# ==============================================================================
run_tprf_standard_one_unit_prepared <- function(
    prep,
    unit_id,
    vintage = c("M1", "M2", "M3"),
    L_tprf = 1
) {
  vintage <- match.arg(vintage)
  
  X_train <- prep$X_train_vec_std
  y_q     <- prep$y_q
  q_star  <- prep$q_star
  q_end   <- prep$q_train_end
  h       <- prep$h
  
  y_train <- as.numeric(y_q[1:q_end, unit_id])
  
  fit <- TPRF_h(
    X = X_train,
    y = y_train,
    L = L_tprf,
    h = h
  )
  
  y_hat  <- fit$y_forecast
  y_true <- as.numeric(y_q[q_star, unit_id])
  
  list(
    unit_id      = unit_id,
    vintage      = vintage,
    horizon_h    = h,
    q_train_end  = q_end,
    q_target     = q_star,
    y_hat        = y_hat,
    y_true       = y_true,
    error        = y_hat - y_true,
    fit          = fit
  )
}


# ==============================================================================
# STANDARD QUARTERLY TPRF: ALL UNITS, ONE VINTAGE, PREPARED INPUTS
# Purpose:
#   Run the standard quarterly TPRF for all units at a given vintage using
#   inputs prepared once upstream.
# Notes:
#   - Stores the horizon used at each vintage
#   - Stores the end of the quarterly training sample and the target quarter
# ==============================================================================
run_tprf_standard_wrapper_prepared <- function(
    prep,
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    L_tprf = 1
) {
  vintage <- match.arg(vintage)
  
  p1 <- sim_obj$dims$n_units
  out <- vector("list", p1)
  
  for (i in 1:p1) {
    tmp <- run_tprf_standard_one_unit_prepared(
      prep    = prep,
      unit_id = i,
      vintage = vintage,
      L_tprf  = L_tprf
    )
    
    out[[i]] <- data.frame(
      unit         = i,
      model        = "TPRF_STD",
      vintage      = vintage,
      horizon_h    = tmp$horizon_h,
      q_train_end  = tmp$q_train_end,
      q_target     = tmp$q_target,
      y_hat        = tmp$y_hat,
      y_true       = tmp$y_true,
      error        = tmp$error
    )
  }
  
  do.call(rbind, out)
}


# ==============================================================================
# STANDARD QUARTERLY TPRF: ONE UNIT, ONE VINTAGE
# Purpose:
#   Run the standard quarterly TPRF for one unit at a given vintage.
# Notes:
#   - TPRF_h() is assumed to be already defined
#   - The forecast horizon h is chosen automatically from the vintage
#   - The target quarter remains the current quarter q_star
# ==============================================================================
run_tprf_standard_one_unit <- function(
    sim_obj,
    unit_id,
    vintage = c("M1", "M2", "M3"),
    L_tprf = 1
) {
  vintage <- match.arg(vintage)
  
  prep <- prepare_tprf_quarterly_inputs(sim_obj, vintage = vintage)
  
  run_tprf_standard_one_unit_prepared(
    prep    = prep,
    unit_id = unit_id,
    vintage = vintage,
    L_tprf  = L_tprf
  )
}


# ==============================================================================
# STANDARD QUARTERLY TPRF: ALL UNITS, ONE VINTAGE
# Purpose:
#   Run the standard quarterly TPRF for all units at a given vintage.
# Notes:
#   - Stores the horizon used at each vintage
#   - Stores the end of the quarterly training sample and the target quarter
# ==============================================================================
run_tprf_standard_wrapper <- function(
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    L_tprf = 1
) {
  vintage <- match.arg(vintage)
  
  prep <- prepare_tprf_quarterly_inputs(sim_obj, vintage = vintage)
  
  run_tprf_standard_wrapper_prepared(
    prep    = prep,
    sim_obj = sim_obj,
    vintage = vintage,
    L_tprf  = L_tprf
  )
}


# ==============================================================================
# VECTOR MF-TPRF: ONE UNIT, ONE VINTAGE, PREPARED INPUTS
# Purpose:
#   Run Vec MF-TPRF for a single unit at a given nowcast vintage using
#   oracle-completed and vintage-standardized inputs that have already been
#   prepared upstream.
# ==============================================================================
run_vec_mftprf_one_unit_prepared <- function(
    prep,
    unit_id,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    Robust_F = FALSE,
    alpha = 0.10,
    robust_type = "NW",
    nw_lag = 1
) {
  vintage <- match.arg(vintage)
  
  X_hf_vec   <- prep$vector$X_hf_vec_std
  X_lf_vec   <- prep$vector$X_lf_vec_std
  y_q_train  <- as.numeric(prep$target$y_train[, unit_id])
  y_true     <- as.numeric(prep$target$y_target[, unit_id])
  
  stopifnot(is.matrix(X_hf_vec))
  stopifnot(is.matrix(X_lf_vec))
  stopifnot(unit_id >= 1, unit_id <= ncol(prep$target$y_train))
  
  fit <- MF_TPRF(
    X_lf        = X_lf_vec,
    X_hf        = X_hf_vec,
    y_q         = y_q_train,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    Robust_F    = Robust_F,
    alpha       = alpha,
    robust_type = robust_type,
    nw_lag      = nw_lag
  )
  
  y_hat <- safe_last_non_na(fit$y_nowcast)
  
  list(
    unit_id    = unit_id,
    vintage    = vintage,
    y_hat      = y_hat,
    y_true     = y_true,
    error      = y_hat - y_true,
    fit        = fit,
    X_hf_used  = X_hf_vec,
    X_lf_used  = X_lf_vec
  )
}


# ==============================================================================
# VECTOR MF-TPRF: ALL UNITS, ONE VINTAGE, PREPARED INPUTS
# Purpose:
#   Run Vec MF-TPRF for all units at a given nowcast vintage using inputs
#   prepared once upstream.
# ==============================================================================
run_vec_mftprf_wrapper_prepared <- function(
    prep,
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    Robust_F = FALSE,
    alpha = 0.10,
    robust_type = "NW",
    nw_lag = 1
) {
  vintage <- match.arg(vintage)
  
  p1 <- sim_obj$dims$n_units
  out <- vector("list", p1)
  
  for (i in 1:p1) {
    tmp <- run_vec_mftprf_one_unit_prepared(
      prep        = prep,
      unit_id     = i,
      vintage     = vintage,
      Lproxy      = Lproxy,
      L_midas     = L_midas,
      p_AR        = p_AR,
      Robust_F    = Robust_F,
      alpha       = alpha,
      robust_type = robust_type,
      nw_lag      = nw_lag
    )
    
    out[[i]] <- data.frame(
      unit    = i,
      model   = "VecMF_TPRF",
      vintage = vintage,
      y_hat   = tmp$y_hat,
      y_true  = tmp$y_true,
      error   = tmp$error
    )
  }
  
  do.call(rbind, out)
}


# ==============================================================================
# VECTOR MF-TPRF: ONE UNIT, ONE VINTAGE
# Purpose:
#   Run Vec MF-TPRF for a single unit at a given nowcast vintage using the
#   oracle-completed and vintage-standardized panels.
# ==============================================================================
run_vec_mftprf_one_unit <- function(
    sim_obj,
    unit_id,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    Robust_F = FALSE,
    alpha = 0.10,
    robust_type = "NW",
    nw_lag = 1
) {
  vintage <- match.arg(vintage)
  
  prep <- prepare_completed_inputs_for_vintage(sim_obj, vintage = vintage)
  
  run_vec_mftprf_one_unit_prepared(
    prep        = prep,
    unit_id     = unit_id,
    vintage     = vintage,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    Robust_F    = Robust_F,
    alpha       = alpha,
    robust_type = robust_type,
    nw_lag      = nw_lag
  )
}


# ==============================================================================
# VECTOR MF-TPRF: ALL UNITS, ONE VINTAGE
# ==============================================================================
run_vec_mftprf_wrapper <- function(
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    Robust_F = FALSE,
    alpha = 0.10,
    robust_type = "NW",
    nw_lag = 1
) {
  vintage <- match.arg(vintage)
  
  prep <- prepare_completed_inputs_for_vintage(sim_obj, vintage = vintage)
  
  run_vec_mftprf_wrapper_prepared(
    prep        = prep,
    sim_obj     = sim_obj,
    vintage     = vintage,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    Robust_F    = Robust_F,
    alpha       = alpha,
    robust_type = robust_type,
    nw_lag      = nw_lag
  )
}


# ==============================================================================
# MATRIX MF-TPRF: ALL UNITS, ONE VINTAGE, PREPARED INPUTS
# Purpose:
#   Run Matrix MF-TPRF jointly for all units at a given vintage using
#   oracle-completed and vintage-standardized tensor panels prepared upstream.
# ==============================================================================
run_mat_mftprf_wrapper_prepared <- function(
    prep,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    r_mat_true = c(1, 1),
    proxy_name = "EA",
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8
) {
  vintage <- match.arg(vintage)
  stopifnot(length(r_mat_true) == 2)
  
  X_hf_tens <- prep$matrix$X_hf_std
  X_lf_tens <- prep$matrix$X_lf_std
  y_train   <- prep$target$y_train
  y_true    <- prep$target$y_target
  z_train   <- prep$target$z_train
  
  stopifnot(length(dim(X_hf_tens)) == 3)
  stopifnot(length(dim(X_lf_tens)) == 3)
  stopifnot(is.matrix(y_train))
  stopifnot(is.matrix(y_true))
  
  p1 <- ncol(y_train)
  
  Y_q_all <- cbind(z_train, y_train)
  colnames(Y_q_all) <- c(proxy_name, paste0("unit", 1:p1))
  
  fit <- Tensor_MF_TPRF(
    X_lf        = X_lf_tens,
    X_hf        = X_hf_tens,
    Y_q_all     = Y_q_all,
    proxy_name  = proxy_name,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    r           = r_mat_true,
    standardize_proxy = standardize_proxy,
    orthonormalize_each_iter = orthonormalize_each_iter,
    orthonormalize_final_Z   = orthonormalize_final_Z,
    ils_maxit = ils_maxit,
    ils_tol   = ils_tol
  )
  
  res_list <- vector("list", p1)
  
  for (i in 1:p1) {
    cc_name  <- paste0("unit", i)
    y_hat_i  <- safe_last_non_na(fit$by_country[[cc_name]]$y_nowcast)
    y_true_i <- as.numeric(y_true[1, i])
    
    res_list[[i]] <- data.frame(
      unit    = i,
      model   = "MatMF_TPRF",
      vintage = vintage,
      y_hat   = y_hat_i,
      y_true  = y_true_i,
      error   = y_hat_i - y_true_i
    )
  }
  
  results_df <- do.call(rbind, res_list)
  
  list(
    results   = results_df,
    fit       = fit,
    r_true    = r_mat_true,
    X_hf_used = X_hf_tens,
    X_lf_used = X_lf_tens,
    Y_q_all   = Y_q_all
  )
}


# ==============================================================================
# MATRIX MF-TPRF: ALL UNITS, ONE VINTAGE
# Purpose:
#   Run Matrix MF-TPRF jointly for all units at a given vintage using the
#   oracle-completed and vintage-standardized tensor panels.
# ==============================================================================
run_mat_mftprf_wrapper <- function(
    sim_obj,
    vintage = c("M1", "M2", "M3"),
    Lproxy,
    L_midas,
    p_AR,
    r_mat_true = c(1, 1),
    proxy_name = "EA",
    standardize_proxy = TRUE,
    orthonormalize_each_iter = TRUE,
    orthonormalize_final_Z = TRUE,
    ils_maxit = 100,
    ils_tol = 1e-8
) {
  vintage <- match.arg(vintage)
  
  prep <- prepare_completed_inputs_for_vintage(sim_obj, vintage = vintage)
  
  run_mat_mftprf_wrapper_prepared(
    prep                     = prep,
    vintage                  = vintage,
    Lproxy                   = Lproxy,
    L_midas                  = L_midas,
    p_AR                     = p_AR,
    r_mat_true               = r_mat_true,
    proxy_name               = proxy_name,
    standardize_proxy        = standardize_proxy,
    orthonormalize_each_iter = orthonormalize_each_iter,
    orthonormalize_final_Z   = orthonormalize_final_Z,
    ils_maxit                = ils_maxit,
    ils_tol                  = ils_tol
  )
}


# ==============================================================================
# ONE MONTE CARLO REPLICATION
# Purpose:
#   Simulate one DGP draw and compare BI, Vec MF-TPRF, Mat MF-TPRF, and the
#   standard quarterly TPRF benchmark at the nowcast vintages M1, M2, M3.
# Notes:
#   - Perfect completion is assumed
#   - No completion step is estimated
#   - No rank selection is performed
#   - The quarterly TPRF benchmark adapts its forecast horizon over the quarter
# ==============================================================================
run_one_mc_replication <- function(dgp_params, est_params, rep_id) {
  dgp_params$seed <- dgp_params$seed + rep_id - 1
  sim <- simulate_matrix_dgp(dgp_params)
  
  vintages <- c("M1", "M2", "M3")
  
  # ---------------------------------------------------------------------------
  # Prepare completed inputs once per vintage
  # ---------------------------------------------------------------------------
  completed_prep <- setNames(vector("list", length(vintages)), vintages)
  for (v in vintages) {
    completed_prep[[v]] <- prepare_completed_inputs_for_vintage(
      sim_obj = sim,
      vintage = v
    )
  }
  
  # ---------------------------------------------------------------------------
  # Prepare quarterly TPRF inputs once per vintage
  # ---------------------------------------------------------------------------
  if (isTRUE(est_params$use_tprf_standard)) {
    tprf_prep <- setNames(vector("list", length(vintages)), vintages)
    for (v in vintages) {
      tprf_prep[[v]] <- prepare_tprf_quarterly_inputs(
        sim_obj = sim,
        vintage = v
      )
    }
  } else {
    tprf_prep <- NULL
  }
  
  # ---------------------------------------------------------------------------
  # BI across M1, M2, M3
  # ---------------------------------------------------------------------------
  bi_M1 <- run_bi_wrapper(
    sim_obj       = sim,
    vintage       = "M1",
    bi_intercept  = est_params$bi_intercept
  )
  
  bi_M2 <- run_bi_wrapper(
    sim_obj       = sim,
    vintage       = "M2",
    bi_intercept  = est_params$bi_intercept
  )
  
  bi_M3 <- run_bi_wrapper(
    sim_obj       = sim,
    vintage       = "M3",
    bi_intercept  = est_params$bi_intercept
  )
  
  # ---------------------------------------------------------------------------
  # Standard quarterly TPRF across M1, M2, M3
  # Horizon-adjusted:
  # M1 -> h = 2 if the previous quarter is not yet complete
  # M2 -> h = 1 once the previous quarter becomes complete
  # M3 -> h = 1
  # ---------------------------------------------------------------------------
  if (isTRUE(est_params$use_tprf_standard)) {
    tprf_M1 <- run_tprf_standard_wrapper_prepared(
      prep    = tprf_prep[["M1"]],
      sim_obj = sim,
      vintage = "M1",
      L_tprf  = est_params$L_tprf
    )
    
    tprf_M2 <- run_tprf_standard_wrapper_prepared(
      prep    = tprf_prep[["M2"]],
      sim_obj = sim,
      vintage = "M2",
      L_tprf  = est_params$L_tprf
    )
    
    tprf_M3 <- run_tprf_standard_wrapper_prepared(
      prep    = tprf_prep[["M3"]],
      sim_obj = sim,
      vintage = "M3",
      L_tprf  = est_params$L_tprf
    )
  } else {
    tprf_M1 <- NULL
    tprf_M2 <- NULL
    tprf_M3 <- NULL
  }
  
  # ---------------------------------------------------------------------------
  # Vector MF-TPRF across M1, M2, M3
  # ---------------------------------------------------------------------------
  vec_M1 <- run_vec_mftprf_wrapper_prepared(
    prep        = completed_prep[["M1"]],
    sim_obj     = sim,
    vintage     = "M1",
    Lproxy      = est_params$Lproxy_vec,
    L_midas     = est_params$L_midas_vec,
    p_AR        = est_params$p_AR_vec
  )
  
  vec_M2 <- run_vec_mftprf_wrapper_prepared(
    prep        = completed_prep[["M2"]],
    sim_obj     = sim,
    vintage     = "M2",
    Lproxy      = est_params$Lproxy_vec,
    L_midas     = est_params$L_midas_vec,
    p_AR        = est_params$p_AR_vec
  )
  
  vec_M3 <- run_vec_mftprf_wrapper_prepared(
    prep        = completed_prep[["M3"]],
    sim_obj     = sim,
    vintage     = "M3",
    Lproxy      = est_params$Lproxy_vec,
    L_midas     = est_params$L_midas_vec,
    p_AR        = est_params$p_AR_vec
  )
  
  # ---------------------------------------------------------------------------
  # Matrix MF-TPRF across M1, M2, M3
  # ---------------------------------------------------------------------------
  mat_M1 <- run_mat_mftprf_wrapper_prepared(
    prep                     = completed_prep[["M1"]],
    vintage                  = "M1",
    Lproxy                   = est_params$Lproxy_mat,
    L_midas                  = est_params$L_midas_mat,
    p_AR                     = est_params$p_AR_mat,
    r_mat_true               = est_params$r_mat_true,
    proxy_name               = est_params$proxy_name,
    standardize_proxy        = est_params$standardize_proxy,
    orthonormalize_each_iter = est_params$orthonormalize_each_iter,
    orthonormalize_final_Z   = est_params$orthonormalize_final_Z,
    ils_maxit                = est_params$ils_maxit,
    ils_tol                  = est_params$ils_tol
  )$results
  
  mat_M2 <- run_mat_mftprf_wrapper_prepared(
    prep                     = completed_prep[["M2"]],
    vintage                  = "M2",
    Lproxy                   = est_params$Lproxy_mat,
    L_midas                  = est_params$L_midas_mat,
    p_AR                     = est_params$p_AR_mat,
    r_mat_true               = est_params$r_mat_true,
    proxy_name               = est_params$proxy_name,
    standardize_proxy        = est_params$standardize_proxy,
    orthonormalize_each_iter = est_params$orthonormalize_each_iter,
    orthonormalize_final_Z   = est_params$orthonormalize_final_Z,
    ils_maxit                = est_params$ils_maxit,
    ils_tol                  = est_params$ils_tol
  )$results
  
  mat_M3 <- run_mat_mftprf_wrapper_prepared(
    prep                     = completed_prep[["M3"]],
    vintage                  = "M3",
    Lproxy                   = est_params$Lproxy_mat,
    L_midas                  = est_params$L_midas_mat,
    p_AR                     = est_params$p_AR_mat,
    r_mat_true               = est_params$r_mat_true,
    proxy_name               = est_params$proxy_name,
    standardize_proxy        = est_params$standardize_proxy,
    orthonormalize_each_iter = est_params$orthonormalize_each_iter,
    orthonormalize_final_Z   = est_params$orthonormalize_final_Z,
    ils_maxit                = est_params$ils_maxit,
    ils_tol                  = est_params$ils_tol
  )$results
  
  bi_M1   <- align_result_columns(bi_M1)
  bi_M2   <- align_result_columns(bi_M2)
  bi_M3   <- align_result_columns(bi_M3)
  
  vec_M1  <- align_result_columns(vec_M1)
  vec_M2  <- align_result_columns(vec_M2)
  vec_M3  <- align_result_columns(vec_M3)
  
  mat_M1  <- align_result_columns(mat_M1)
  mat_M2  <- align_result_columns(mat_M2)
  mat_M3  <- align_result_columns(mat_M3)
  
  tprf_M1 <- align_result_columns(tprf_M1)
  tprf_M2 <- align_result_columns(tprf_M2)
  tprf_M3 <- align_result_columns(tprf_M3)
  
  out <- rbind(
    bi_M1, bi_M2, bi_M3,
    vec_M1, vec_M2, vec_M3,
    mat_M1, mat_M2, mat_M3,
    tprf_M1, tprf_M2, tprf_M3
  )
  
  out$rep <- rep_id
  out
}


# ==============================================================================
# FULL MONTE CARLO LOOP
# Purpose:
#   Repeat the comparison across many DGP replications.
# ==============================================================================
run_mc_simulation <- function(dgp_params, est_params, draws, verbose = TRUE) {
  results_list <- vector("list", draws)
  failed_reps  <- integer(0)
  error_msgs   <- character(0)
  
  start_time <- Sys.time()
  
  for (b in 1:draws) {
    if (verbose) {
      cat("\n==============================\n")
      cat("Monte Carlo replication:", b, "of", draws, "\n")
      cat("==============================\n")
    }
    
    res_b <- tryCatch(
      {
        run_one_mc_replication(
          dgp_params = dgp_params,
          est_params = est_params,
          rep_id     = b
        )
      },
      error = function(e) {
        if (verbose) {
          cat("Replication", b, "FAILED:\n")
          cat(conditionMessage(e), "\n")
        }
        failed_reps <<- c(failed_reps, b)
        error_msgs  <<- c(error_msgs, conditionMessage(e))
        NULL
      }
    )
    
    results_list[[b]] <- res_b
  }
  
  results_list <- Filter(Negate(is.null), results_list)
  
  if (length(results_list) == 0) {
    stop("All Monte Carlo replications failed.")
  }
  
  results_df <- do.call(rbind, results_list)
  
  end_time <- Sys.time()
  
  list(
    raw = results_df,
    failed_reps = failed_reps,
    error_msgs = error_msgs,
    n_success = length(results_list),
    n_failed = length(failed_reps),
    start_time = start_time,
    end_time = end_time,
    elapsed = end_time - start_time
  )
}


# ==============================================================================
# MONTE CARLO RESULTS: SUMMARY BY MODEL AND VINTAGE
# Purpose:
#   Compute mean forecast error, MAE, MSFE, and RMSFE.
# Notes:
#   - Extra TPRF columns such as horizon_h, q_train_end and q_target are
#     ignored automatically by the aggregation
# ==============================================================================
summarize_mc_results <- function(mc_obj) {
  df <- if (is.list(mc_obj) && !is.null(mc_obj$raw)) mc_obj$raw else mc_obj
  
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or a Monte Carlo object with a $raw data.frame.")
  }
  
  required_cols <- c("model", "vintage", "error")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df$abs_error <- abs(df$error)
  df$sq_error  <- df$error^2
  
  out <- aggregate(
    cbind(error, abs_error, sq_error) ~ model + vintage,
    data = df,
    FUN = mean,
    na.rm = TRUE
  )
  
  names(out)[names(out) == "error"]     <- "mean_bias"
  names(out)[names(out) == "abs_error"] <- "MAE"
  names(out)[names(out) == "sq_error"]  <- "MSFE"
  out$RMSFE <- sqrt(out$MSFE)
  
  model_order   <- c("BI", "VecMF_TPRF", "MatMF_TPRF", "TPRF_STD")
  vintage_order <- c("M1", "M2", "M3")
  
  out$model   <- factor(out$model, levels = model_order)
  out$vintage <- factor(out$vintage, levels = vintage_order)
  out <- out[order(out$model, out$vintage), , drop = FALSE]
  
  out$model   <- as.character(out$model)
  out$vintage <- as.character(out$vintage)
  rownames(out) <- NULL
  
  out[, c("model", "vintage", "mean_bias", "MAE", "MSFE", "RMSFE")]
}


# ==============================================================================
# TPRF FORECASTS BY VINTAGE
# Purpose:
#   Summarize the standard quarterly TPRF forecasts across M1, M2, M3.
# Notes:
#   - Useful to check both the forecast path over the quarter and the average
#     horizon effectively used by the TPRF benchmark
# ==============================================================================
summarize_tprf_forecasts_by_vintage <- function(mc_obj) {
  df <- if (is.list(mc_obj) && !is.null(mc_obj$raw)) mc_obj$raw else mc_obj
  
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or a Monte Carlo object with a $raw data.frame.")
  }
  
  df <- df[df$model == "TPRF_STD", , drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("No TPRF_STD results found.")
  }
  
  required_cols <- c("vintage", "y_hat", "y_true", "error", "horizon_h")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required TPRF columns: ", paste(missing_cols, collapse = ", "))
  }
  
  out <- aggregate(
    cbind(y_hat, y_true, error, horizon_h) ~ vintage,
    data = df,
    FUN = mean,
    na.rm = TRUE
  )
  
  out$abs_error <- aggregate(abs(error) ~ vintage, data = df, FUN = mean, na.rm = TRUE)[, 2]
  out$sq_error  <- aggregate((error^2) ~ vintage, data = df, FUN = mean, na.rm = TRUE)[, 2]
  out$RMSFE     <- sqrt(out$sq_error)
  
  vintage_order <- c("M1", "M2", "M3")
  out$vintage <- factor(out$vintage, levels = vintage_order)
  out <- out[order(out$vintage), , drop = FALSE]
  out$vintage <- as.character(out$vintage)
  rownames(out) <- NULL
  
  names(out)[names(out) == "horizon_h"] <- "avg_horizon"
  
  out[, c("vintage", "avg_horizon", "y_hat", "y_true", "error", "abs_error", "sq_error", "RMSFE")]
}


# ==============================================================================
# TPRF FORECASTS IN WIDE FORMAT
# Purpose:
#   Return one selected TPRF quantity in wide format across M1, M2, M3.
# ==============================================================================
summarize_tprf_forecasts_wide <- function(mc_obj, value = c("y_hat", "y_true", "error", "avg_horizon", "RMSFE")) {
  value <- match.arg(value)
  
  df <- summarize_tprf_forecasts_by_vintage(mc_obj)
  
  value_col <- switch(
    value,
    y_hat = "y_hat",
    y_true = "y_true",
    error = "error",
    avg_horizon = "avg_horizon",
    RMSFE = "RMSFE"
  )
  
  out <- data.frame(
    model = paste0("TPRF_", value),
    M1 = NA_real_,
    M2 = NA_real_,
    M3 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (v in df$vintage) {
    out[1, v] <- df[df$vintage == v, value_col]
  }
  
  out
}


# ==============================================================================
# WIDE TABLE FOR A SELECTED METRIC
# rows    = models
# columns = M1, M2, M3
# ==============================================================================
summarize_mc_metric_wide <- function(mc_obj, metric = "RMSFE") {
  sum_df <- summarize_mc_results(mc_obj)
  
  if (!metric %in% names(sum_df)) {
    stop("Requested metric not found in summary table.")
  }
  
  wide_df <- reshape(
    sum_df[, c("model", "vintage", metric)],
    idvar = "model",
    timevar = "vintage",
    direction = "wide"
  )
  
  names(wide_df) <- sub(paste0("^", metric, "\\."), "", names(wide_df))
  
  wanted_cols <- c("model", "M1", "M2", "M3")
  existing_cols <- wanted_cols[wanted_cols %in% names(wide_df)]
  wide_df <- wide_df[, existing_cols, drop = FALSE]
  
  model_order <- c("BI", "VecMF_TPRF", "MatMF_TPRF", "TPRF_STD")
  wide_df$model <- factor(wide_df$model, levels = model_order)
  wide_df <- wide_df[order(wide_df$model), , drop = FALSE]
  wide_df$model <- as.character(wide_df$model)
  rownames(wide_df) <- NULL
  
  wide_df
}

summarize_mc_rmsfe_wide <- function(mc_obj) {
  summarize_mc_metric_wide(mc_obj, metric = "RMSFE")
}


# ==============================================================================
# MATRIX-FOCUSED RMSFE TABLES
# ==============================================================================
summarize_mc_model_metric_wide <- function(
    mc_obj,
    model_name = "MatMF_TPRF",
    metric = "RMSFE"
) {
  sum_df <- summarize_mc_results(mc_obj)
  
  if (!metric %in% names(sum_df)) {
    stop("Requested metric not found in summary table.")
  }
  
  df <- sum_df[sum_df$model == model_name, c("model", "vintage", metric), drop = FALSE]
  
  if (nrow(df) == 0) {
    stop("Requested model not found in summary table.")
  }
  
  wide_df <- reshape(
    df,
    idvar = "model",
    timevar = "vintage",
    direction = "wide"
  )
  
  names(wide_df) <- sub(paste0("^", metric, "\\."), "", names(wide_df))
  wanted_cols <- c("model", "M1", "M2", "M3")
  existing_cols <- wanted_cols[wanted_cols %in% names(wide_df)]
  wide_df <- wide_df[, existing_cols, drop = FALSE]
  
  rownames(wide_df) <- NULL
  wide_df
}


summarize_mc_relative_metric_wide <- function(
    mc_obj,
    numerator_model = "MatMF_TPRF",
    denominator_model = "VecMF_TPRF",
    metric = "RMSFE",
    row_name = NULL
) {
  sum_df <- summarize_mc_results(mc_obj)
  
  if (!metric %in% names(sum_df)) {
    stop("Requested metric not found in summary table.")
  }
  
  num_df <- sum_df[sum_df$model == numerator_model, c("vintage", metric), drop = FALSE]
  den_df <- sum_df[sum_df$model == denominator_model, c("vintage", metric), drop = FALSE]
  
  if (nrow(num_df) == 0) stop("Numerator model not found in summary table.")
  if (nrow(den_df) == 0) stop("Denominator model not found in summary table.")
  
  merged <- merge(
    num_df,
    den_df,
    by = "vintage",
    suffixes = c("_num", "_den"),
    all = FALSE
  )
  
  merged$ratio <- merged[[paste0(metric, "_num")]] / merged[[paste0(metric, "_den")]]
  
  if (is.null(row_name)) {
    row_name <- paste0(numerator_model, "/", denominator_model)
  }
  
  out <- data.frame(
    model = row_name,
    M1 = NA_real_,
    M2 = NA_real_,
    M3 = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (v in merged$vintage) {
    out[1, v] <- merged$ratio[merged$vintage == v]
  }
  
  out
}


summarize_mc_mat_compact_rmsfe <- function(mc_obj) {
  mat_lvl <- summarize_mc_model_metric_wide(
    mc_obj = mc_obj,
    model_name = "MatMF_TPRF",
    metric = "RMSFE"
  )
  mat_lvl$model <- "MatMF_TPRF"
  
  mat_vec <- summarize_mc_relative_metric_wide(
    mc_obj = mc_obj,
    numerator_model = "MatMF_TPRF",
    denominator_model = "VecMF_TPRF",
    metric = "RMSFE",
    row_name = "Mat/Vec"
  )
  
  mat_bi <- summarize_mc_relative_metric_wide(
    mc_obj = mc_obj,
    numerator_model = "MatMF_TPRF",
    denominator_model = "BI",
    metric = "RMSFE",
    row_name = "Mat/BI"
  )
  
  mat_tprf <- summarize_mc_relative_metric_wide(
    mc_obj = mc_obj,
    numerator_model = "MatMF_TPRF",
    denominator_model = "TPRF_STD",
    metric = "RMSFE",
    row_name = "Mat/TPRF"
  )
  
  out <- rbind(mat_lvl, mat_vec, mat_bi, mat_tprf)
  rownames(out) <- NULL
  out
}


# ==============================================================================
# GENERIC LATEX EXPORTER
# ==============================================================================
mc_summary_to_latex <- function(
    summary_df,
    digits = 3,
    caption = "Monte Carlo results",
    label = "tab:mc_results"
) {
  stopifnot(is.data.frame(summary_df))
  
  df <- summary_df
  num_cols <- vapply(df, is.numeric, logical(1))
  df[num_cols] <- lapply(
    df[num_cols],
    function(x) sprintf(paste0("%.", digits, "f"), x)
  )
  
  header <- paste(names(df), collapse = " & ")
  body <- apply(df, 1, function(row) paste(row, collapse = " & "))
  body <- paste0(body, " \\\\")
  
  align <- paste0("l", paste(rep("r", ncol(df) - 1), collapse = ""))
  
  latex <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\hline",
    paste0(header, " \\\\"),
    "\\hline",
    body,
    "\\hline",
    "\\end{tabular}",
    "\\end{table}"
  )
  
  paste(latex, collapse = "\n")
}

write_mc_summary_latex <- function(
    summary_df,
    file,
    digits = 3,
    caption = "Monte Carlo results",
    label = "tab:mc_results"
) {
  tex <- mc_summary_to_latex(
    summary_df = summary_df,
    digits = digits,
    caption = caption,
    label = label
  )
  writeLines(tex, con = file)
  invisible(tex)
}


# ==============================================================================
# LATEX EXPORTERS YOU ACTUALLY NEED
# ==============================================================================
mc_rmsfe_wide_to_latex <- function(
    mc_obj,
    digits = 3,
    caption = "Monte Carlo RMSFE by model and nowcast vintage",
    label = "tab:mc_rmsfe_wide"
) {
  mc_summary_to_latex(
    summary_df = summarize_mc_rmsfe_wide(mc_obj),
    digits = digits,
    caption = caption,
    label = label
  )
}

write_mc_rmsfe_wide_latex <- function(
    mc_obj,
    file,
    digits = 3,
    caption = "Monte Carlo RMSFE by model and nowcast vintage",
    label = "tab:mc_rmsfe_wide"
) {
  tex <- mc_rmsfe_wide_to_latex(
    mc_obj = mc_obj,
    digits = digits,
    caption = caption,
    label = label
  )
  writeLines(tex, con = file)
  invisible(tex)
}


mc_mat_compact_rmsfe_to_latex <- function(
    mc_obj,
    digits = 3,
    caption = "Matrix MF-TPRF RMSFE and relative performance by nowcast vintage",
    label = "tab:mc_mat_compact_rmsfe"
) {
  mc_summary_to_latex(
    summary_df = summarize_mc_mat_compact_rmsfe(mc_obj),
    digits = digits,
    caption = caption,
    label = label
  )
}

write_mc_mat_compact_rmsfe_latex <- function(
    mc_obj,
    file,
    digits = 3,
    caption = "Matrix MF-TPRF RMSFE and relative performance by nowcast vintage",
    label = "tab:mc_mat_compact_rmsfe"
) {
  tex <- mc_mat_compact_rmsfe_to_latex(
    mc_obj = mc_obj,
    digits = digits,
    caption = caption,
    label = label
  )
  writeLines(tex, con = file)
  invisible(tex)
}


# ==============================================================================
# LATEX EXPORTERS FOR TPRF MONTH-BY-MONTH RESULTS
# ==============================================================================
mc_tprf_forecasts_to_latex <- function(
    mc_obj,
    digits = 3,
    caption = "Standard TPRF forecasts across nowcast vintages",
    label = "tab:tprf_forecasts"
) {
  mc_summary_to_latex(
    summary_df = summarize_tprf_forecasts_by_vintage(mc_obj),
    digits = digits,
    caption = caption,
    label = label
  )
}

write_mc_tprf_forecasts_latex <- function(
    mc_obj,
    file,
    digits = 3,
    caption = "Standard TPRF forecasts across nowcast vintages",
    label = "tab:tprf_forecasts"
) {
  tex <- mc_tprf_forecasts_to_latex(
    mc_obj = mc_obj,
    digits = digits,
    caption = caption,
    label = label
  )
  writeLines(tex, con = file)
  invisible(tex)
}


# ==============================================================================
# PLOTS (KEEP ONLY ESSENTIAL ONES)
# ==============================================================================
plot_mc_rmsfe <- function(mc_obj) {
  sum_df <- summarize_mc_results(mc_obj)
  sum_df$spec <- paste(sum_df$model, sum_df$vintage, sep = " - ")
  
  ggplot(sum_df, aes(x = spec, y = RMSFE)) +
    geom_col() +
    theme_minimal() +
    labs(
      title = "Monte Carlo RMSFE by model and vintage",
      x = "Model / vintage",
      y = "RMSFE"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_mc_by_vintage <- function(mc_obj) {
  sum_df <- summarize_mc_results(mc_obj)
  
  ggplot(sum_df, aes(x = vintage, y = RMSFE, group = model, color = model)) +
    geom_line() +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "RMSFE across nowcast vintages",
      x = "Vintage",
      y = "RMSFE"
    )
}