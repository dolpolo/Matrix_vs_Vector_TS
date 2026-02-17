# ==============================================================================
# STANDARDIZATION with NA
# ==============================================================================

standardize_mat_with_na <- function(X) {
  T <- dim(X)[1]
  p1 <- dim(X)[2]
  p2 <- dim(X)[3]
  
  X_std <- array(NA, dim = dim(X), dimnames = dimnames(X))
  mean_X <- matrix(0, p1, p2)
  sd_X <- matrix(1, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      X_ij <- X[, i, j]
      mu <- mean(X_ij, na.rm = TRUE)
      sigma <- sd(X_ij, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      X_std[, i, j] <- (X_ij - mu) / sigma
      mean_X[i, j] <- mu
      sd_X[i, j] <- sigma
    }
  }
  
  return(list(X_scaled = X_std, mean = mean_X, sd = sd_X))
}


center_Y <- function(Y) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  Y_centered <- array(NA, dim = dim(Y))
  dimnames(Y_centered) <- dimnames(Y) 
  mean_Y <- matrix(0, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu <- mean(y_ij, na.rm = TRUE)
      Y_centered[, i, j] <- y_ij - mu
      mean_Y[i, j] <- mu
    }
  }
  
  return(list(Y_centered = Y_centered, mean = mean_Y))
}

decenter_Y <- function(Y_centered, mean_Y) {
  T <- dim(Y_centered)[1]
  p1 <- dim(Y_centered)[2]
  p2 <- dim(Y_centered)[3]
  
  Y_original <- array(NA, dim = dim(Y_centered))
  dimnames(Y_original) <- dimnames(Y_centered)
  
  for (t in 1:T) {
    Y_original[t,,] <- Y_centered[t,,] + mean_Y
  }
  
  return(Y_original)
}



# ==============================================================================
# COVID PERIOD MASKING
# ==============================================================================

mask_covid <- function(X, dates, class, start, end) {
  
  real <- tolower(class) %in% c("r", "real")
  per  <- dates >= start & dates <= end
  
  # LOGICAL mask (important!)
  M <- outer(per, real, FUN = "&")
  
  X[M] <- NA
  
  list(data = X, mask = M)
}




# ==============================================================================
# MONTHLY AND QUARTERLY VARIABLES TYPES OF AGGREGATION
# ==============================================================================
# Xm <- Ym
# Type_m <- leg_type[1:ncol(Ym)]
# j <- 1
agg_mq <- function(Xm, agg) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  # Numero di trimestri completi
  T_q <- floor(T / 3)
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    for (tau in 1:T_q) {
      m3 <- 3 * tau
      m2 <- m3 - 1
      m1 <- m3 - 2
      
      if (agg[j] == 2) {
        # FLOW → SOMMA
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
        
      } else if (agg[j] == 1) {
        # STOCK → MEDIA
        Xq[tau, j] <- mean(c(Xm[m1, j], Xm[m2, j], Xm[m3, j]), na.rm = TRUE)
        
      } else {
        stop("Aggregation code must be 1 (stock) or 2 (flow).")
      }
    }
  }
  
  # nomi delle colonne
  colnames(Xq) <- colnames(Xm)
  
  # nomi delle righe: usa il 3° mese di ogni trimestre (se disponibili)
  if (!is.null(rownames(Xm))) {
    rn <- rownames(Xm)
    idx_q <- seq(3, by = 3, length.out = T_q)
    rownames(Xq) <- rn[idx_q]
  }
  
  return(Xq)
}


agg_qq <- function(Xm, agg) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  # numero di trimestri completi
  T_q <- floor(T / 3)
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    for (tau in 1:T_q) {
      m3 <- 3 * tau
      m2 <- m3 - 1
      m1 <- m3 - 2
      
      if (agg[j] == 2) {
        # FLOW → somma
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
        
      } else if (agg[j] == 1) {
        # STOCK → media
        Xq[tau, j] <- mean(c(Xm[m1, j], Xm[m2, j], Xm[m3, j]), na.rm = TRUE)
      } else {
        stop("Aggregation code must be 1 (stock) or 2 (flow).")
      }
    }
  }
  
  # nomi colonne
  colnames(Xq) <- colnames(Xm)
  
  # nomi righe
  if (!is.null(rownames(Xm))) {
    rn <- rownames(Xm)
    idx_q <- seq(3, by = 3, length.out = T_q)
    rownames(Xq) <- rn[idx_q]
  }
  
  return(Xq)
}




# ==============================================================================
# MONTHLY AND QUARTERLY VARIABLES TYPES OF AGGREGATION
# ==============================================================================

# X <- Y_mq
# X <- Yq_nt
# K <- params$n_m
# K <- params$n_q
# lambda <- lambda_m
# lambda <- lambda_q

LASSO_set <- function(y, X, K, alpha = 1, nlambda = 2000, lambda_min_ratio = 1e-5) {
  Z  <- as.matrix(X)
  yy <- as.numeric(y)
  
  ok <- complete.cases(cbind(yy, Z))
  Z  <- Z[ok, , drop = FALSE]
  yy <- yy[ok]
  
  fit_path <- glmnet::glmnet(
    x = Z, y = yy,
    alpha = alpha, standardize = TRUE,
    nlambda = nlambda, lambda.min.ratio = lambda_min_ratio
  )
  
  # Ridge: nessuno zero -> scegli top-K per |beta|
  if (alpha == 0) {
    id <- ceiling(length(fit_path$lambda) / 2)
    return(fit_path$lambda[id])
  }
  
  id <- which.min(abs(fit_path$df - K))
  fit_path$lambda[id]
}

LASSO_select <- function(y, X, lambda, K, alpha = 1) {
  Z  <- as.matrix(X)
  yy <- as.numeric(y)
  
  ok <- complete.cases(cbind(yy, Z))
  Z  <- Z[ok, , drop = FALSE]
  yy <- yy[ok]
  
  fit <- glmnet::glmnet(Z, yy, alpha = alpha, standardize = TRUE, lambda = lambda)
  b   <- as.matrix(coef(fit))[, 1]
  b   <- b[setdiff(names(b), "(Intercept)")]
  
  # Ridge: nessuno zero -> scegli top-K per |beta|
  if (alpha == 0) {
    ord <- order(abs(b), decreasing = TRUE)
    return(names(b)[ord][1:min(K, length(b))])
  }
  
  # LASSO / Elastic Net: prendi non-zero
  idx_nz <- which(b != 0)
  sel <- names(b)[idx_nz]
  
  # se >K, taglia ai K più grandi per |beta|
  if (length(sel) > K) {
    b_sel <- b[sel]
    sel <- names(sort(abs(b_sel), decreasing = TRUE))[1:K]
  }
  
  sel
}



# ==============================================================================
# VARIABLE SELECTION BASED ON GDP CORRELATION, F-Test, LASSO
# ==============================================================================
# path   <- path_data
# covid_mask = TRUE
# cc <- "EA"

# Selects variables most correlated with GDP for each country.
# Returns base (country-independent) names for monthly and quarterly datasets.

select_vars <- function(countries, params, path) {
  
  target <- params$target
  
  # dates limits
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  sel_m <- list()
  sel_q <- list()
  
  for (cc in countries) {
    
    # ---- Load data ----
    dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
    dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
    
    # ---- select the correct time span
    dm  <-  dm %>%
      filter(Time >= start_lim & Time <= end_lim)
    dq  <-  dq %>%
      filter(Time >= start_lim & Time <= end_lim)
    
    # ---- Remove date column ----
    Ym <- as.data.frame(dm[, -1])
    Yq <- as.data.frame(dq[, -1])
    
    # ---- Legend ----
    leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet = "info")
    TR   <- floor(leg$TR)
    freq <- leg$Frequency
    agg  <- leg$Aggregation
    
    # ---- Identify target quarterly variable ----
    target_id <- grep(paste0("^", target, "_"), colnames(Yq))[1]
    names_q   <- colnames(Yq)
    
    # ---- Match quarterly dataset with legend ----
    idx_leg_in_q <- match(tolower(names_q), tolower(leg$Name))
    
    # ---- Identify quarterly series (all Q variables) ----
    freq_q <- freq[idx_leg_in_q]
    need_shift <- which(freq_q == "Q")
    
    # --------------------------------------------------------
    # SHIFT TRIMESTRALI
    # --------------------------------------------------------
    
    Yq_shifted_full <- as.data.frame(
      apply(Yq, 2, function(col){
        Y <- col                      # valori trimestrali
        kronecker(Y, c(NA, NA, 1))         # monthly expansion: NA,NA,Y
      })
    )
    
    # estrai SOLO i mesi finali dei trimestri (3, 6, 9, ...)
    rows_quarter_end <- seq(3, nrow(Yq_shifted_full), by = 3)
    Yq_shifted <- Yq_shifted_full[rows_quarter_end, , drop = FALSE]
    
    # --------------------------------------------------------
    # MONTHLY → QUARTERLY aggregation
    # --------------------------------------------------------
    
    Y_mq <- agg_mq(Ym, agg)   # dimensione Tq
    
    Tq <- nrow(Y_mq)
    
    # ---- Target quarterly series (aligned) ----
    y <- as.numeric(Yq_shifted[[target_id]])[1:Tq]
    
    # ---- Selection Method ----
    method <- params$sel_method
    
    #──────────────────────────────────────────────
    # METHOD 1: NONE (pick top n correlations)
    #──────────────────────────────────────────────
    if (method == "none") {
      
      # -- Monthly --
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- setNames(as.numeric(cm), colnames(Y_mq))
      vars_m_base <- names(sort(cm, decreasing=TRUE))[1:params$n_m]
      
      # -- Quarterly (shifted) --
      cq <- abs(cor(Yq_shifted[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- setNames(as.numeric(cq), colnames(Yq_shifted)[-target_id])
      vars_q_base <- names(sort(cq, decreasing=TRUE))[1:params$n_q]
      
      #──────────────────────────────────────────────
      # METHOD 2: CORRELATION THRESHOLD
      #──────────────────────────────────────────────
    } else if (method == "corr_threshold") {
      
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- setNames(as.numeric(cm), colnames(Y_mq))
      vars_m_base <- names(cm[cm >= params$thr_m])
      
      cq <- abs(cor(Yq_shifted[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- setNames(as.numeric(cq), colnames(Yq_shifted)[-target_id])
      vars_q_base <- names(cq[cq >= params$thr_q])
      
      #──────────────────────────────────────────────
      # METHOD 3: F-TEST
      #──────────────────────────────────────────────
    } else if (method == "t_test") {
      
      pvals_m <- sapply(1:ncol(Y_mq), function(j) {
        f <- summary(lm(y ~ Y_mq[, j]))$fstatistic
        pf(f[1], f[2], f[3], lower.tail = FALSE)
      })
      names(pvals_m) <- colnames(Y_mq)
      vars_m_base <- names(pvals_m[pvals_m <= params$thr_F_test])
      
      Yq_nt <- Yq_shifted[1:Tq, -target_id, drop=FALSE]
      
      pvals_q <- sapply(1:ncol(Yq_nt), function(j) {
        f <- summary(lm(y ~ Yq_nt[, j]))$fstatistic
        pf(f[1], f[2], f[3], lower.tail = FALSE)
      })
      names(pvals_q) <- colnames(Yq_nt)
      vars_q_base <- names(pvals_q[pvals_q <= params$thr_F_test])
      
      #──────────────────────────────────────────────
      # METHOD 4: LASSO
      #──────────────────────────────────────────────
    } else if (method == "LASSO") {
      
      alpha <- params$alpha_lasso
      
      # quarterly non-target block
      Yq_nt <- Yq_shifted[1:Tq, -target_id, drop = FALSE]
      
      # Monthly block
      lambda_m <- LASSO_set(y = y, X = Y_mq, K = params$n_m, alpha = alpha)
      vars_m_base <- LASSO_select(y = y, X = Y_mq, lambda = lambda_m,
                                  K = params$n_m, alpha = alpha)
      
      # Quarterly block
      lambda_q <- LASSO_set(y = y, X = Yq_nt, K = params$n_q, alpha = alpha)
      vars_q_base <- LASSO_select(y = y, X = Yq_nt, lambda = lambda_q,
                                  K = params$n_q, alpha = alpha)
    }
    #──────────────────────────────────────────────
    # FINALIZE
    #──────────────────────────────────────────────
    
    vars_m_cc <- vars_m_base
    vars_q_cc <- vars_q_base
    
    target_cc <- colnames(Yq)[target_id]
    
    sel_m[[cc]] <- vars_m_cc
    sel_q[[cc]] <- c(target_cc, vars_q_cc)
  }
  
  return(list(
    m = sel_m,
    q = sel_q
  ))
}


# ==============================================================================
# PREPARE DATA FOR A SINGLE COUNTRY
# ==============================================================================
# sel_m <- vars_sel$m
# sel_q <- vars_sel$q
prepare_country_data <- function(cc, params, sel_m, sel_q, path,
                                 covid_mask    = TRUE,
                                 covid_mask_m  = covid_mask,  # mask mensili
                                 covid_mask_q  = covid_mask)  # mask trimestrali
{
  # ============================================================
  # 1. LOAD DATA (MONTHLY LT + QUARTERLY LT)
  # ============================================================
  dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
  dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
  leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet = "info")
  
  # ============================================================
  # 2. LEGEND PARSING
  # ============================================================
  var_col   <- intersect("Name", names(leg))[1]
  leg_names <- tolower(trimws(leg[[var_col]]))
  leg_base  <- sub(paste0("_", tolower(cc), "$"), "", leg_names)
  
  TR        <- floor(leg$TR)
  leg_class <- leg$Class
  leg_type  <- leg$Type
  leg_freq  <- leg$Frequency
  leg_unb   <- leg$M1
  freq      <- leg$Frequency
  agg       <- leg$Aggregation
  
  # ============================================================
  # 3. DATES LIMITS
  # ============================================================
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  # ---- select the correct time span
  dm <- dm %>%
    filter(Time >= start_lim & Time <= end_lim)
  dq <- dq %>%
    filter(Time >= start_lim & Time <= end_lim)
  
  # ============================================================
  # 4. EXTRACT MONTHLY DATA
  # ============================================================
  date_m <- as.Date(dm[[1]])
  Ym_raw <- as.matrix(dm[, -1, drop = FALSE])
  Xm     <- Ym_raw[, sel_m[[cc]], drop = FALSE]
  
  # metadata mensili
  base_m  <- sub(paste0("_", cc, "$"), "", sel_m[[cc]])
  class_m <- leg_class[match(tolower(base_m), leg_base)]
  type_m  <- leg_type[match(tolower(base_m), leg_base)]
  freq_m  <- leg_freq[match(tolower(base_m), leg_base)]
  unb_m   <- leg_unb[match(tolower(base_m), leg_base)]
  agg_m   <- agg[match(tolower(base_m), leg_base)]
  
  # COVID mask monthly (controllata da covid_mask_m)
  if (covid_mask_m) {
    res_m <- mask_covid(Xm, date_m, class_m,
                        params$covid_start, params$covid_end)
    Xm     <- res_m$data
    mask_m <- res_m$mask
  } else {
    mask_m <- matrix(FALSE, nrow(Xm), ncol(Xm))
  }
  
  # ============================================================
  # 5. EXTRACT QUARTERLY DATA
  # ============================================================
  date_q <- as.Date(dq[[1]])
  Xq_all <- as.matrix(dq[, -1, drop = FALSE])
  
  # colonne trimestrali effettivamente presenti
  valid_cols <- sel_q[[cc]][sel_q[[cc]] %in% colnames(Xq_all)]
  
  Xq <- Xq_all[, valid_cols, drop = FALSE]
  Xq <- matrix(as.numeric(Xq), nrow = nrow(Xq_all),
               dimnames = list(NULL, valid_cols))
  
  # metadata trimestrali
  base_q  <- sub(paste0("_", cc, "$"), "", valid_cols)
  class_q <- leg_class[match(tolower(base_q), leg_base)]
  type_q  <- leg_type[match(tolower(base_q), leg_base)]
  freq_q  <- leg_freq[match(tolower(base_q), leg_base)]
  unb_q   <- leg_unb[match(tolower(base_q), leg_base)]
  agg_q   <- agg[match(tolower(base_q), leg_base)]
  
  # --------------------------------------------------------
  # 5B. IDENTIFICA LA TARGET TRA LE TRIMESTRALI
  #     (senza assumere che sia la prima colonna)
  # --------------------------------------------------------
  target_pattern <- paste0("^", tolower(params$target), "_")
  target_idx_q   <- grep(target_pattern, tolower(colnames(Xq)))
  
  # vettore classi usato SOLO per la maschera (target esclusa)
  class_q_mask <- class_q
  if (length(target_idx_q) > 0) {
    class_q_mask[target_idx_q] <- "target"  # qualcosa ≠ "r"/"real"
  }
  
  # date del terzo mese del trimestre
  date_q_third <- date_q %m+% months(2)
  
  # ============================================================
  # 6. CREATE FULL MONTHLY TIME AXIS
  # ============================================================
  date_full <- seq.Date(
    from = min(date_m),
    to   = max(date_m),
    by   = "month"
  )
  
  # ============================================================
  # 7. ALIGN MONTHLY
  # ============================================================
  Xm_full <- matrix(NA_real_, length(date_full), ncol(Xm))
  rownames(Xm_full) <- as.character(date_full)
  colnames(Xm_full) <- colnames(Xm)
  Xm_full[as.character(date_m), ] <- Xm
  
  mask_m_full <- matrix(FALSE, length(date_full), ncol(Xm))
  rownames(mask_m_full) <- as.character(date_full)
  colnames(mask_m_full) <- colnames(Xm)
  mask_m_full[as.character(date_m), ] <- mask_m
  
  # ============================================================
  # 8. ALIGN QUARTERLY (portate al 3° mese del trimestre)
  # ============================================================
  Xq_full <- matrix(NA_real_, length(date_full), ncol(Xq))
  rownames(Xq_full) <- as.character(date_full)
  colnames(Xq_full) <- colnames(Xq)
  
  # inserisci ciascun trimestre nel 3° mese
  for (t in seq_along(date_q_third)) {
    key <- as.character(date_q_third[t])
    if (key %in% rownames(Xq_full)) {
      Xq_full[key, ] <- Xq[t, ]
    }
  }
  
  # ---- COVID mask trimestrali (target esclusa tramite class_q_mask) ----
  if (covid_mask_q) {
    res_q <- mask_covid(Xq_full, date_full, class_q_mask,
                        params$covid_start, params$covid_end)
    Xq_full     <- res_q$data
    mask_q_full <- res_q$mask
  } else {
    mask_q_full <- matrix(FALSE, nrow(Xq_full), ncol(Xq_full))
    rownames(mask_q_full) <- rownames(Xq_full)
    colnames(mask_q_full) <- colnames(Xq_full)
  }
  
  # ============================================================
  # 9. APPLY FINAL TIME WINDOW
  # ============================================================
  idx <- which(date_full >= start_lim & date_full <= end_lim)
  
  Xm_out     <- Xm_full[idx, , drop = FALSE]
  Xq_out     <- Xq_full[idx, , drop = FALSE]
  mask_m_out <- mask_m_full[idx, , drop = FALSE]
  mask_q_out <- mask_q_full[idx, , drop = FALSE]
  dates_out  <- date_full[idx]
  
  dates_m_out <- date_m[date_m >= start_lim & date_m <= end_lim]
  
  # Date dei trimestri (3° mese) dentro la finestra
  dates_q_all <- date_q_third[date_q_third >= start_lim & date_q_third <= end_lim]
  dates_q_out <- dates_q_all
  
  # ============================================================
  # 10. IDENTIFY TARGET COLUMN (su tutte le serie)
  # ============================================================
  all_series <- c(sel_m[[cc]], sel_q[[cc]])
  
  target_pattern_all <- paste0("^", tolower(params$target), "_")
  target_col         <- grep(target_pattern_all, tolower(all_series))
  target_name        <- all_series[target_col]
  
  # ============================================================
  # 11. FINAL OUTPUT
  # ============================================================
  list(
    Data        = cbind(Xm_out, Xq_out),
    Dates       = dates_out,
    DatesM      = dates_m_out,
    DatesQ      = dates_q_out,
    Series      = all_series,
    nM          = length(sel_m[[cc]]),
    nQ          = length(sel_q[[cc]]),
    agg_m       = agg_m,
    agg_q       = agg_q,
    ClassM      = class_m,
    ClassQ      = class_q,
    TypeM       = type_m,
    TypeQ       = type_q,
    MaskM       = mask_m_out,
    MaskQ       = mask_q_out,
    freq_m      = freq_m,
    freq_q      = freq_q,
    unb_m       = unb_m,
    unb_q       = unb_q,
    idx_q       = ncol(Xm_out) + 1,
    target_col  = target_col,
    target_name = target_name
  )
}


# ==============================================================================
# GDP PROXY FROM EA 
# ==============================================================================

get_proxy_EA_target <- function(params, path) {
  
  cc_proxy <- params$target_cc   # "EA"
  target   <- params$target      # "GDP"
  
  dq  <- readxl::read_excel(file.path(path, paste0(cc_proxy, "dataQ_LT.xlsx")))
  
  date_q <- as.Date(dq[[1]])
  Yq     <- as.data.frame(dq[, -1, drop = FALSE])
  
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  keep_time <- (date_q >= start_lim & date_q <= end_lim)
  date_q    <- date_q[keep_time]
  Yq        <- Yq[keep_time, , drop = FALSE]
  
  target_pattern <- paste0("^", target, "_")
  target_id      <- grep(target_pattern, colnames(Yq))[1]
  if (is.na(target_id)) {
    stop("Target ", target, " not found in ", cc_proxy, " quarterly file.")
  }
  
  target_name <- colnames(Yq)[target_id]
  y_raw       <- as.numeric(Yq[[target_id]])
  
  # Terzo mese del trimestre
  date_q_third <- date_q %m+% months(2)
  
  # Rimuovi NA
  keep_nonNA <- !is.na(y_raw)
  y_raw      <- y_raw[keep_nonNA]
  date_q_third <- date_q_third[keep_nonNA]
  
  # Finestra finale
  keep_win <- (date_q_third >= start_lim & date_q_third <= end_lim)
  y_raw    <- y_raw[keep_win]
  date_q_third <- date_q_third[keep_win]
  
  # ---- MATRICE NOMINATA ----
  Y_q_proxy <- matrix(
    y_raw,
    ncol = 1,
    dimnames = list(
      as.character(date_q_third),
      cc_proxy
    )
  )
  
  list(
    country     = cc_proxy,
    Y_q         = Y_q_proxy,
    target_name = target_name
  )
}


# ==============================================================================
# WRAPPER FUNCTION: COMPLETE MULTI-COUNTRY DATA PREPARATION
# ==============================================================================


prepare_all_countries <- function(countries, params, path,
                                  covid_mask = TRUE,
                                  covid_mask_m,
                                  covid_mask_q) {
  
  # Paesi effettivi (senza l’area aggregata usata come proxy, es. "EA")
  countries_eff <- setdiff(countries, params$target_cc)
  
  # 1. Variable selection country-by-country (solo DE, FR, IT, ES)
  vars_sel <- select_vars(countries_eff, params, path)
  
  # 2. Prepare each country (DE, FR, IT, ES)
  out <- lapply(countries_eff, function(cc) {
    prepare_country_data(
      cc        = cc,
      params    = params,
      sel_m     = vars_sel$m,
      sel_q     = vars_sel$q,
      path      = path,
      covid_mask_m = params$covid_mask_m,
      covid_mask_q = params$covid_mask_q
    )
  })
  names(out) <- countries_eff
  
  # 3. Estrai anche la proxy aggregata (es. GDP EA)
  proxy_out <- get_proxy_EA_target(params, path)
  
  # 4. Output complessivo
  list(
    data   = out,        # lista: DE, FR, IT, ES
    sel    = vars_sel,   # variabili selezionate per ciascun paese
    proxy  = proxy_out   # PIL aggregato EA (y_q + dates_q + nome colonna)
  )
}


# =====================================================================
# BUILD MATRIX-VARIATE DATA (TENSOR) FROM PREPARED COUNTRY DATA
# =====================================================================


# prep <- all_countries

# Rimuove suffissi tipo "_IT", "_FR", "_DE", ecc.
remove_country_code <- function(x, countries) {
  pattern <- paste0("_(", paste(countries, collapse = "|"), ")$")
  sub(pattern, "", x, ignore.case = TRUE)
}

build_tensor <- function(prep,
                         params,
                         var_scope = c("union", "intersection")) {
  
  var_scope <- match.arg(var_scope)
  
  countries <- names(prep$data)
  data_list <- prep$data
  
  # ------------------------------------------------------------
  # 1. ASSE TEMPORALE GLOBALE da params (mensile)
  # ------------------------------------------------------------
  all_dates <- seq(params$start_est, params$end_eval, by = "month")
  T_global  <- length(all_dates)
  
  # ------------------------------------------------------------
  # 2. LISTE DI VARIABILI BASE (M e Q) PER OGNI PAESE
  # ------------------------------------------------------------
  series_list <- lapply(data_list, function(obj) obj$Series)
  
  monthly_base_list <- lapply(data_list, function(obj) {
    Series_m <- obj$Series[1:obj$nM]  # prime nM sono mensili
    remove_country_code(Series_m, countries = countries)
  })
  
  quarterly_base_list <- lapply(data_list, function(obj) {
    Series_q <- obj$Series[(obj$nM + 1):(obj$nM + obj$nQ)]  # poi nQ trimestrali
    remove_country_code(Series_q, countries = countries)
  })
  
  # ------------------------------------------------------------
  # 3. UNION / INTERSECTION SU MENSILI E TRIMESTRALI
  # ------------------------------------------------------------
  if (var_scope == "intersection") {
    base_M <- Reduce(intersect, monthly_base_list)
    base_Q <- Reduce(intersect, quarterly_base_list)
  } else {  # "union"
    base_M <- sort(unique(unlist(monthly_base_list)))
    base_Q <- sort(unique(unlist(quarterly_base_list)))
  }
  
  # Ordine finale: prima tutte le mensili, poi tutte le trimestrali
  vars_ordered <- c(base_M, base_Q)
  
  if (length(vars_ordered) == 0L) {
    stop("Nessuna variabile selezionata dopo union/intersection.")
  }
  
  P1 <- length(countries)
  P2 <- length(vars_ordered)
  
  # ------------------------------------------------------------
  # 4. INIZIALIZZA TENSORE Y (DATI) E W (MASK)
  # ------------------------------------------------------------
  Y <- array(NA_real_,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               vars_ordered
             ))
  
  W <- array(0L,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               vars_ordered
             ))
  
  # ------------------------------------------------------------
  # 4bis. INIZIALIZZA METADATA [Paesi x Variabili]
  # ------------------------------------------------------------
  agg_mat   <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  class_mat <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  type_mat  <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  freq_mat  <- matrix(NA_character_, nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  unb_mat   <- matrix(NA_real_,       nrow = P1, ncol = P2,
                      dimnames = list(countries, vars_ordered))
  
  # ------------------------------------------------------------
  # 5. RIEMPI Y, W E METADATA PAESE PER PAESE
  # ------------------------------------------------------------
  for (i in seq_along(countries)) {
    cc   <- countries[i]
    obj  <- data_list[[cc]]
    
    Data_cc   <- obj$Data        # matrice (T_cc x n_var_cc)
    Dates_cc  <- obj$Dates
    Series_cc <- obj$Series      # (mensili + trimestrali)
    
    base_cc <- remove_country_code(Series_cc, countries = countries)
    
    # Mappa date paese -> asse globale
    idx_time <- match(Dates_cc, all_dates)
    
    # Matrice temporanea per questo paese
    Y_cc <- matrix(NA_real_,
                   nrow = T_global, ncol = P2,
                   dimnames = list(as.character(all_dates), vars_ordered))
    
    for (j in seq_along(vars_ordered)) {
      v_base <- vars_ordered[j]
      
      col_idx <- which(base_cc == v_base)
      if (length(col_idx) == 0L) next
      
      # se più colonne con stesso base name → prendo la prima
      col_idx <- col_idx[1]
      
      # -----------------------------
      # 5a. DATI + MASK
      # -----------------------------
      Y_cc[idx_time, j] <- Data_cc[, col_idx]
      
      # -----------------------------
      # 5b. METADATA (M vs Q)
      # -----------------------------
      if (col_idx <= obj$nM) {
        # MENSILE
        k <- col_idx
        
        agg_mat[i, j]   <- obj$agg_m[k]
        class_mat[i, j] <- obj$ClassM[k]
        type_mat[i, j]  <- obj$TypeM[k]
        freq_mat[i, j]  <- obj$freq_m[k]
        unb_mat[i, j]   <- obj$unb_m[k]
        
      } else {
        # TRIMESTRALE
        k <- col_idx - obj$nM
        
        agg_mat[i, j]   <- obj$agg_q[k]
        class_mat[i, j] <- obj$ClassQ[k]
        type_mat[i, j]  <- obj$TypeQ[k]
        freq_mat[i, j]  <- obj$freq_q[k]
        unb_mat[i, j]   <- obj$unb_q[k]
      }
    }
    
    Y[, i, ] <- Y_cc
    W[, i, ] <- ifelse(is.na(Y_cc), 0L, 1L)
  }
  
  # ------------------------------------------------------------
  # 6. CLASSIFICA VARIABILI COME "M" / "Q"
  # ------------------------------------------------------------
  var_type <- rep(NA_character_, P2)
  names(var_type) <- vars_ordered
  var_type[vars_ordered %in% base_M] <- "M"
  var_type[vars_ordered %in% base_Q] <- "Q"
  
  n_M <- sum(var_type == "M")
  n_Q <- sum(var_type == "Q")
  
  idx_M <- which(var_type == "M")
  idx_Q <- which(var_type == "Q")
  
  vars_M <- vars_ordered[idx_M]
  vars_Q <- vars_ordered[idx_Q]
  
  # ------------------------------------------------------------
  # 7. SPLIT METADATA IN MENSILI / TRIMESTRALI
  # ------------------------------------------------------------
  agg_M   <- if (n_M > 0)   agg_mat[, idx_M, drop = FALSE]   else NULL
  agg_Q   <- if (n_Q > 0)   agg_mat[, idx_Q, drop = FALSE]   else NULL
  
  ClassM  <- if (n_M > 0) class_mat[, idx_M, drop = FALSE]   else NULL
  ClassQ  <- if (n_Q > 0) class_mat[, idx_Q, drop = FALSE]   else NULL
  
  TypeM   <- if (n_M > 0) type_mat[, idx_M, drop = FALSE]    else NULL
  TypeQ   <- if (n_Q > 0) type_mat[, idx_Q, drop = FALSE]    else NULL
  
  freq_M  <- if (n_M > 0) freq_mat[, idx_M, drop = FALSE]    else NULL
  freq_Q  <- if (n_Q > 0) freq_mat[, idx_Q, drop = FALSE]    else NULL
  
  unb_M   <- if (n_M > 0)  unb_mat[, idx_M, drop = FALSE]    else NULL
  unb_Q   <- if (n_Q > 0)  unb_mat[, idx_Q, drop = FALSE]    else NULL
  
  # ------------------------------------------------------------
  # 8. TARGET (base name = params$target, es. "GDP")
  # ------------------------------------------------------------
  target_base <- tolower(params$target)
  target_col  <- which(tolower(vars_ordered) == target_base)
  
  # ------------------------------------------------------------
  # 9. OUTPUT
  # ------------------------------------------------------------
  out <- list(
    Y          = Y,               # [T x Paesi x Variabili]
    W          = W,               # mask [T x Paesi x Variabili]
    dates      = all_dates,
    countries  = countries,
    vars       = vars_ordered,
    var_type   = var_type,        # "M" / "Q"
    n_M        = n_M,
    n_Q        = n_Q,
    idx_M      = idx_M,
    idx_Q      = idx_Q,
    vars_M     = vars_M,
    vars_Q     = vars_Q,
    target_col = target_col,
    
    # metadata globali [Paese x Variabile]
    agg        = agg_mat,
    class      = class_mat,
    type       = type_mat,
    freq       = freq_mat,
    unb        = unb_mat,
    
    # metadata separati M / Q (stile prepare_country_data)
    agg_M      = agg_M,
    agg_Q      = agg_Q,
    ClassM     = ClassM,
    ClassQ     = ClassQ,
    TypeM      = TypeM,
    TypeQ      = TypeQ,
    freq_M     = freq_M,
    freq_Q     = freq_Q,
    unb_M      = unb_M,
    unb_Q      = unb_Q
  )
  
  return(out)
}
# ==============================================================================
# NOT AVAILABLE PERCENTAGE
# ==============================================================================

nan_percent_Y <- function(Y) {
  total_values <- length(Y)
  total_NAs <- sum(is.na(Y))
  percent <- total_NAs / total_values * 100
  return(percent)
}
