# ==============================================================================
# TPRF
# ==============================================================================

run_tprf_wrapper <- function(sim_obj, L, train_quarters = NULL) {
  X_full <- sim_obj$truth$X_quarterly_lf_vec
  Y_full <- sim_obj$truth$y_quarterly
  
  stopifnot(is.matrix(X_full))
  stopifnot(is.matrix(Y_full))
  stopifnot(nrow(X_full) == nrow(Y_full))
  
  Tq_total <- nrow(X_full)
  p1 <- ncol(Y_full)
  
  if (is.null(train_quarters)) {
    train_quarters <- Tq_total - 1
  }
  
  if (train_quarters < 2) {
    stop("Need at least 2 training quarters.")
  }
  if (train_quarters >= Tq_total) {
    stop("train_quarters must be strictly smaller than the total number of quarters.")
  }
  
  X_train <- X_full[1:train_quarters, , drop = FALSE]
  
  results <- vector("list", p1)
  
  for (i in 1:p1) {
    y_full_i  <- as.numeric(Y_full[, i])
    y_train_i <- y_full_i[1:train_quarters]
    y_true_i  <- y_full_i[train_quarters + 1]
    
    fit_i <- TPRF(
      X = X_train,
      y = y_train_i,
      L = L
    )
    
    results[[i]] <- data.frame(
      unit    = i,
      model   = "TPRF",
      vintage = "forecast",
      y_hat   = fit_i$y_forecast,
      y_true  = y_true_i,
      error   = fit_i$y_forecast - y_true_i
    )
  }
  
  out_df <- do.call(rbind, results)
  
  list(
    results = out_df,
    X_train = X_train,
    train_quarters = train_quarters,
    L = L
  )
}


# ==============================================================================
# MF-TPRF
# ==============================================================================

get_vector_panel_metadata <- function(sim_obj) {
  p1 <- sim_obj$dims$n_units
  p2 <- sim_obj$dims$n_vars
  
  var_type_base <- sim_obj$params$var_type
  stopifnot(length(var_type_base) == p2)
  
  var_type_vec <- rep(var_type_base, times = p1)
  
  is_monthly   <- grepl("^monthly_", var_type_vec)
  is_quarterly <- grepl("^quarterly_", var_type_vec)
  
  # aggregation codes required by agg_mq / agg_qq:
  # 1 = stock, 2 = flow
  agg_code <- ifelse(
    grepl("_flow$", var_type_vec),
    2,
    1
  )
  
  agg_m <- agg_code[is_monthly]
  agg_q <- agg_code[is_quarterly]
  
  list(
    var_type_vec = var_type_vec,
    is_monthly = is_monthly,
    is_quarterly = is_quarterly,
    agg_m = agg_m,
    agg_q = agg_q,
    N_m = sum(is_monthly),
    N_q = sum(is_quarterly)
  )
}


run_mftprf_one_unit <- function(sim_obj, unit_id, vintage = c("M1", "M2", "M3"),
                                Lproxy, L_midas, p_AR,
                                Kmax = 5,
                                Robust_F = FALSE,
                                alpha = 0.10,
                                robust_type = "NW",
                                nw_lag = 1) {
  
  vintage <- match.arg(vintage)
  
  # -----------------------------
  # 1. Extract vintage data
  # -----------------------------
  X_mf_vec <- sim_obj$info_sets[[vintage]]$vector$X_mf_vec
  Y_q_full <- sim_obj$truth$y_quarterly
  
  meta_vec <- get_vector_panel_metadata(sim_obj)
  N_m <- meta_vec$N_m
  N_q <- meta_vec$N_q
  agg_m <- meta_vec$agg_m
  agg_q <- meta_vec$agg_q
  
  Tq_total <- nrow(Y_q_full)
  T_q_train <- Tq_total - 1
  
  y_q_train <- as.numeric(Y_q_full[1:T_q_train, unit_id])
  y_true    <- as.numeric(Y_q_full[T_q_train + 1, unit_id])
  
  # -----------------------------
  # 2. Standardization + imputation
  # -----------------------------
  out_std <- standardize_with_na(X_mf_vec)
  X_std   <- out_std$X_std
  
  imp_xp <- init_XP_ER(X_std, Kmax)
  X_xp   <- imp_xp$X_init
  
  # -----------------------------
  # 3. Split monthly / quarterly
  # -----------------------------
  X_m_xp <- X_xp[, meta_vec$is_monthly, drop = FALSE]
  X_q_xp <- X_xp[, meta_vec$is_quarterly, drop = FALSE]
  
  # -----------------------------
  # 4. Aggregate to quarterly LF
  # -----------------------------
  X_mq_xp <- agg_mq(X_m_xp, agg_m)
  X_qq_xp <- agg_qq(X_q_xp, agg_q)
  X_xp_agg <- cbind(X_mq_xp, X_qq_xp)
  
  # Keep only the complete quarters aligned with y_q_train
  X_lf_train <- X_xp_agg[1:T_q_train, , drop = FALSE]
  
  # -----------------------------
  # 5. Run MF-TPRF
  # -----------------------------
  fit <- MF_TPRF(
    X_lf        = X_lf_train,
    X_hf        = X_xp,
    y_q         = y_q_train,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    Robust_F    = Robust_F,
    alpha       = alpha,
    robust_type = robust_type,
    nw_lag      = nw_lag
  )
  
  y_hat <- tail(na.omit(fit$y_nowcast), 1)
  
  list(
    unit_id = unit_id,
    vintage = vintage,
    y_hat = y_hat,
    y_true = y_true,
    error = y_hat - y_true,
    fit = fit
  )
}

run_mftprf_wrapper <- function(sim_obj, vintage = c("M1", "M2", "M3"),
                               Lproxy, L_midas, p_AR,
                               Kmax = 5,
                               Robust_F = FALSE,
                               alpha = 0.10,
                               robust_type = "NW",
                               nw_lag = 1) {
  
  vintage <- match.arg(vintage)
  p1 <- sim_obj$dims$n_units
  
  out <- vector("list", p1)
  
  for (i in 1:p1) {
    tmp <- run_mftprf_one_unit(
      sim_obj = sim_obj,
      unit_id = i,
      vintage = vintage,
      Lproxy = Lproxy,
      L_midas = L_midas,
      p_AR = p_AR,
      Kmax = Kmax,
      Robust_F = Robust_F,
      alpha = alpha,
      robust_type = robust_type,
      nw_lag = nw_lag
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
# MATRIX MF-TPRF
# ==============================================================================

get_tensor_panel_metadata <- function(sim_obj) {
  p1 <- sim_obj$dims$n_units
  p2 <- sim_obj$dims$n_vars
  Nm <- sim_obj$dims$n_monthly_vars
  Nq <- sim_obj$dims$n_quarterly_vars
  
  var_type_base <- sim_obj$params$var_type
  stopifnot(length(var_type_base) == p2)
  
  # 1 = stock, 2 = flow
  agg_code_base <- ifelse(
    grepl("_flow$", var_type_base),
    2,
    1
  )
  
  # aggregate_tensor_to_quarterly wants agg of dimension P1 x V
  agg <- matrix(
    rep(agg_code_base, each = p1),
    nrow = p1,
    ncol = p2,
    byrow = FALSE
  )
  
  list(
    agg = agg,
    N_m = Nm,
    N_q = Nq
  )
}

run_matrix_mftprf_wrapper <- function(sim_obj, vintage = c("M1", "M2", "M3"),
                                      Lproxy, L_midas, p_AR,
                                      Kmax = c(1, 5),
                                      r = NULL,
                                      proxy_name = "EA",
                                      standardize_proxy = TRUE,
                                      orthonormalize_each_iter = TRUE,
                                      orthonormalize_final_Z = TRUE,
                                      ils_maxit = 100,
                                      ils_tol = 1e-8) {
  
  vintage <- match.arg(vintage)
  
  # -----------------------------
  # 1. Extract vintage data
  # -----------------------------
  X_mf     <- sim_obj$info_sets[[vintage]]$matrix$X_mf
  obs_mask <- sim_obj$info_sets[[vintage]]$matrix$obs_mask
  Y_q_full <- sim_obj$truth$y_quarterly
  z_q_full <- sim_obj$truth$z_quarterly
  
  Tq_total <- nrow(Y_q_full)
  p1 <- ncol(Y_q_full)
  
  T_q_train <- Tq_total - 1
  
  y_train <- Y_q_full[1:T_q_train, , drop = FALSE]
  y_true  <- Y_q_full[T_q_train + 1, , drop = FALSE]
  z_train <- z_q_full[1:T_q_train]
  z_true  <- z_q_full[T_q_train + 1]
  
  # -----------------------------
  # 2. Metadata
  # -----------------------------
  meta_tens <- get_tensor_panel_metadata(sim_obj)
  agg <- meta_tens$agg
  N_m <- meta_tens$N_m
  N_q <- meta_tens$N_q
  
  # -----------------------------
  # 3. Standardization + imputation
  # -----------------------------
  out_std <- standardize_mat_with_na(X_mf)
  X_std   <- out_std$X_scaled
  
  imp_cl <- init_CL_Yu(X_std, obs_mask, Kmax)
  X_cl   <- imp_cl$Y_init
  r_hat  <- if (is.null(r)) imp_cl$r else r
  
  # -----------------------------
  # 4. Aggregate to low frequency
  # -----------------------------
  X_cl_q <- aggregate_tensor_to_quarterly(
    X_tens = X_cl,
    agg    = agg,
    N_m    = N_m,
    N_q    = N_q
  )
  
  # Keep only complete quarters aligned with y_train
  X_lf_train <- X_cl_q[1:T_q_train, , , drop = FALSE]
  
  # -----------------------------
  # 5. Build Y_q_all
  # -----------------------------
  Y_q_all <- cbind(z_train, y_train)
  colnames(Y_q_all) <- c(proxy_name, paste0("unit", 1:p1))
  
  # -----------------------------
  # 6. Run Matrix MF-TPRF
  # -----------------------------
  fit <- Tensor_MF_TPRF(
    X_lf        = X_lf_train,
    X_hf        = X_cl,
    Y_q_all     = Y_q_all,
    proxy_name  = proxy_name,
    Lproxy      = Lproxy,
    L_midas     = L_midas,
    p_AR        = p_AR,
    r           = r_hat,
    standardize_proxy = standardize_proxy,
    orthonormalize_each_iter = orthonormalize_each_iter,
    orthonormalize_final_Z   = orthonormalize_final_Z,
    ils_maxit = ils_maxit,
    ils_tol   = ils_tol
  )
  
  # -----------------------------
  # 7. Extract final nowcast
  # -----------------------------
  res_list <- vector("list", p1)
  
  for (i in 1:p1) {
    cc_name <- paste0("unit", i)
    y_hat_i <- tail(na.omit(fit$by_country[[cc_name]]$y_nowcast), 1)
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
    results = results_df,
    fit = fit,
    r_hat = r_hat,
    X_hf_imputed = X_cl,
    X_lf_train = X_lf_train,
    Y_q_all = Y_q_all,
    z_true = z_true
  )
}



# ==============================================================================
# MONTECARLO SIMULATION
# ==============================================================================

run_one_mc_replication <- function(dgp_params, est_params, rep_id) {
  dgp_params$seed <- dgp_params$seed + rep_id - 1
  sim <- simulate_matrix_dgp(dgp_params)
  
  tprf_res <- run_tprf_wrapper(
    sim_obj = sim,
    L = est_params$L_tprf
  )$results
  tprf_res$model <- "TPRF"
  tprf_res$vintage <- "forecast"
  
  vec_M1 <- run_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M1",
    Lproxy = est_params$Lproxy_vec,
    L_midas = est_params$L_midas_vec,
    p_AR = est_params$p_AR_vec,
    Kmax = est_params$Kmax_vec
  )
  vec_M1$model <- "VecMF_TPRF"
  
  vec_M2 <- run_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M2",
    Lproxy = est_params$Lproxy_vec,
    L_midas = est_params$L_midas_vec,
    p_AR = est_params$p_AR_vec,
    Kmax = est_params$Kmax_vec
  )
  vec_M2$model <- "VecMF_TPRF"
  
  vec_M3 <- run_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M3",
    Lproxy = est_params$Lproxy_vec,
    L_midas = est_params$L_midas_vec,
    p_AR = est_params$p_AR_vec,
    Kmax = est_params$Kmax_vec
  )
  vec_M3$model <- "VecMF_TPRF"
  
  mat_M1 <- run_matrix_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M1",
    Lproxy = est_params$Lproxy_mat,
    L_midas = est_params$L_midas_mat,
    p_AR = est_params$p_AR_mat,
    Kmax = est_params$Kmax_mat
  )$results
  mat_M1$model <- "MatMF_TPRF"
  
  mat_M2 <- run_matrix_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M2",
    Lproxy = est_params$Lproxy_mat,
    L_midas = est_params$L_midas_mat,
    p_AR = est_params$p_AR_mat,
    Kmax = est_params$Kmax_mat
  )$results
  mat_M2$model <- "MatMF_TPRF"
  
  mat_M3 <- run_matrix_mftprf_wrapper(
    sim_obj = sim,
    vintage = "M3",
    Lproxy = est_params$Lproxy_mat,
    L_midas = est_params$L_midas_mat,
    p_AR = est_params$p_AR_mat,
    Kmax = est_params$Kmax_mat
  )$results
  mat_M3$model <- "MatMF_TPRF"
  
  out <- rbind(
    tprf_res,
    vec_M1, vec_M2, vec_M3,
    mat_M1, mat_M2, mat_M3
  )
  
  out$rep <- rep_id
  out
}

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
          rep_id = b
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
# MONTECARLO RESULTS
# ==============================================================================

summarize_mc_results <- function(mc_obj) {
  df <- if (is.list(mc_obj) && !is.null(mc_obj$raw)) mc_obj$raw else mc_obj
  
  if (!"error" %in% names(df)) {
    stop("Input must contain an 'error' column.")
  }
  
  if (!"abs_error" %in% names(df)) {
    df$abs_error <- abs(df$error)
  }
  
  if (!"sq_error" %in% names(df)) {
    df$sq_error <- df$error^2
  }
  
  out <- aggregate(
    cbind(error, abs_error, sq_error) ~ model + vintage,
    data = df,
    FUN = mean
  )
  
  names(out)[names(out) == "error"]     <- "mean_bias"
  names(out)[names(out) == "abs_error"] <- "MAE"
  names(out)[names(out) == "sq_error"]  <- "MSFE"
  
  out$RMSE <- sqrt(out$MSFE)
  
  out <- out[, c("model", "vintage", "mean_bias", "MAE", "MSFE", "RMSE")]
  out[order(out$model, out$vintage), ]
}

summarize_mc_by_unit <- function(mc_obj) {
  df <- mc_obj$raw
  
  out <- aggregate(
    cbind(error, abs_error, sq_error) ~ unit + model + vintage,
    data = df,
    FUN = mean
  )
  
  names(out)[names(out) == "error"]     <- "mean_bias"
  names(out)[names(out) == "abs_error"] <- "MAE"
  names(out)[names(out) == "sq_error"]  <- "MSFE"
  
  out$RMSE <- sqrt(out$MSFE)
  out[order(out$unit, out$model, out$vintage), ]
}



library(ggplot2)

plot_mc_boxplot <- function(mc_obj) {
  df <- mc_obj$raw
  df$spec <- paste(df$model, df$vintage, sep = " - ")
  
  ggplot(df, aes(x = spec, y = error)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = "Monte Carlo forecast errors",
      x = "Model / vintage",
      y = "Forecast error"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_mc_abs_error_boxplot <- function(mc_obj) {
  df <- mc_obj$raw
  df$spec <- paste(df$model, df$vintage, sep = " - ")
  
  ggplot(df, aes(x = spec, y = abs_error)) +
    geom_boxplot() +
    theme_minimal() +
    labs(
      title = "Monte Carlo absolute forecast errors",
      x = "Model / vintage",
      y = "Absolute error"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_mc_rmse <- function(mc_obj) {
  sum_df <- summarize_mc_results(mc_obj)
  sum_df$spec <- paste(sum_df$model, sum_df$vintage, sep = " - ")
  
  ggplot(sum_df, aes(x = spec, y = RMSE)) +
    geom_col() +
    theme_minimal() +
    labs(
      title = "Monte Carlo RMSE by model and vintage",
      x = "Model / vintage",
      y = "RMSE"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_mc_by_vintage <- function(mc_obj) {
  sum_df <- summarize_mc_results(mc_obj)
  
  ggplot(sum_df, aes(x = vintage, y = RMSE, group = model, color = model)) +
    geom_line() +
    geom_point(size = 2) +
    theme_minimal() +
    labs(
      title = "RMSE across vintages",
      x = "Vintage",
      y = "RMSE"
    )
}