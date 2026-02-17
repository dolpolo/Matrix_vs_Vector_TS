# ==============================================================================
# STANDARDIZATION AND DESTANDARDIZATION
# ==============================================================================

standardize_Y <- function(Y) {
  T <- dim(Y)[1]
  p1 <- dim(Y)[2]
  p2 <- dim(Y)[3]
  
  Y_std <- array(NA, dim = dim(Y), dimnames = dimnames(Y))
  mean_Y <- matrix(0, p1, p2)
  sd_Y <- matrix(1, p1, p2)
  
  for (i in 1:p1) {
    for (j in 1:p2) {
      y_ij <- Y[, i, j]
      mu <- mean(y_ij, na.rm = TRUE)
      sigma <- sd(y_ij, na.rm = TRUE)
      if (is.na(sigma) || sigma == 0) sigma <- 1
      Y_std[, i, j] <- (y_ij - mu) / sigma
      mean_Y[i, j] <- mu
      sd_Y[i, j] <- sigma
    }
  }
  
  return(list(Y_scaled = Y_std, mean = mean_Y, sd = sd_Y))
}

inverse_standardize_Y <- function(Y_scaled, mean_Y, sd_Y) {
  T <- dim(Y_scaled)[1]
  p1 <- dim(Y_scaled)[2]
  p2 <- dim(Y_scaled)[3]
  
  Y_original <- array(NA, dim = dim(Y_scaled), dimnames = dimnames(Y_scaled))
  
  for (t in 1:T) {
    Y_original[t,,] <- Y_scaled[t,,] * sd_Y + mean_Y
  }
  
  return(Y_original)
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

inverse_standardize_nowcasts <- function(nowcast_list, mean_Y, sd_Y) {
  lapply(nowcast_list, function(mat) {
    mat * sd_Y + mean_Y
  })
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
agg_mq <- function(Xm, Type_m) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  # numero di trimestri COMPLETI
  T_q <- floor(T / 3)
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    
    type_j <- tolower(Type_m[j])
    
    for (tau in 1:T_q) {
      m3  <- 3*tau
      m2  <- m3 - 1
      m1  <- m3 - 2
      
      if (type_j == "flow") {
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
        
      } else if (type_j == "average") {
        Xq[tau, j] <- (Xm[m1, j] + Xm[m2, j] + Xm[m3, j]) / 3
        
      } else if (type_j == "stock") {
        Xq[tau, j] <- Xm[m3, j]  # ultimo mese del trimestre
      }
    }
  }
  
  colnames(Xq) <- colnames(Xm)
  return(Xq)
}

agg_qq <- function(Xm, Types) {
  
  Xm <- as.matrix(Xm)
  T  <- nrow(Xm)
  N  <- ncol(Xm)
  
  T_q <- floor(T / 3)    # ⇐ CORRETTO, UGUALE A agg_mq()
  
  Xq <- matrix(NA, T_q, N)
  
  for (j in 1:N) {
    
    type_j <- tolower(Types[j])
    
    for (tau in 1:T_q) {
      m3 <- 3 * tau
      m2 <- m3 - 1
      m1 <- m3 - 2
      
      if (type_j == "flow")
        Xq[tau, j] <- Xm[m1, j] + Xm[m2, j] + Xm[m3, j]
      else if (type_j %in% c("average","avg"))
        Xq[tau, j] <- mean(Xm[c(m1, m2, m3), j])
      else if (type_j == "stock")
        Xq[tau, j] <- Xm[m3, j]
      else
        stop(paste("Tipo non riconosciuto:", colnames(Xm)[j]))
    }
  }
  
  colnames(Xq) <- colnames(Xm)
  Xq
}

# ==============================================================================
# VARIABLE SELECTION BASED ON GDP CORRELATION
# ==============================================================================

# Selects variables most correlated with GDP for each country.
# Returns base (country-independent) names for monthly and quarterly datasets.

select_vars <- function(countries, params, path) {
  
  target <- params$target
  
  sel_m <- list()
  sel_q <- list()
  
  for (cc in countries) {
    
    # ---- Load ----
    dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
    dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
    leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet = "info")
    
    # ---- Legend ----
    leg_type <- leg$Type
    
    # ---- Remove dates ----
    Ym <- as.data.frame(dm[, -1])
    Yq <- as.data.frame(dq[, -1])
    
    # ---- Target ----
    target_id <- grep(paste0("^", target, "_"), colnames(Yq))[1]
    target_q  <- as.numeric(Yq[[target_id]])
    
    # ---- Monthly -> Quarterly ----
    Y_mq <- agg_mq(Ym, leg_type[1:ncol(Ym)])
    
    Tq <- nrow(Y_mq)
    y  <- target_q[1:Tq]
    
    method <- params$sel_method
    
    #────────────────────────────────────────────────────
    # METHOD: NONE
    #────────────────────────────────────────────────────
    if (method == "none") {
      
      # Monthly
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- as.numeric(cm); names(cm) <- colnames(Y_mq)
      vars_m_base <- names(sort(cm, decreasing=TRUE))[1:params$n_m]
      
      # Quarterly
      cq <- abs(cor(Yq[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- as.numeric(cq); names(cq) <- colnames(Yq)[-target_id]
      vars_q_base <- names(sort(cq, decreasing=TRUE))[1:params$n_q]
      
      #────────────────────────────────────────────────────
      # METHOD: CORR_THRESHOLD
      #────────────────────────────────────────────────────
    } else if (method == "corr_threshold") {
      
      cm <- abs(cor(Y_mq, y, use="pairwise.complete.obs"))
      cm <- as.numeric(cm); names(cm) <- colnames(Y_mq)
      vars_m_base <- names(cm[cm >= params$thr_m])
      
      cq <- abs(cor(Yq[1:Tq, -target_id, drop=FALSE], y, use="pairwise.complete.obs"))
      cq <- as.numeric(cq); names(cq) <- colnames(Yq)[-target_id]
      vars_q_base <- names(cq[cq >= params$thr_q])
      
      #────────────────────────────────────────────────────
      # METHOD: F-TEST
      #────────────────────────────────────────────────────
    } else if (method == "F-Test") {
      
      # Monthly
      pvals_m <- sapply(1:ncol(Y_mq), function(j) {
        f <- summary(lm(Y_mq[, j] ~ y))$fstatistic
        pf(f[1], f[2], f[3], lower.tail = FALSE)
      })
      names(pvals_m) <- colnames(Y_mq)
      vars_m_base <- names(pvals_m[pvals_m <= params$thr_F_test])
      
      # Quarterly
      Yq_nt <- Yq[1:Tq, -target_id, drop=FALSE]
      pvals_q <- sapply(1:ncol(Yq_nt), function(j) {
        f <- summary(lm(Yq_nt[, j] ~ y))$fstatistic
        pf(f[1], f[2], f[3], lower.tail=FALSE)
      })
      names(pvals_q) <- colnames(Yq_nt)
      vars_q_base <- names(pvals_q[pvals_q <= params$thr_F_test])
    }
    
    #────────────────────────────────────────────────────
    # REINTRODUCI I SUFFISSI DEL PAESE
    #────────────────────────────────────────────────────
    
    vars_m_cc <- paste0(vars_m_base)      # le mensili NON hanno suffisso nei file Excel
    vars_q_cc <- paste0(vars_q_base)      # le trimestrali hanno già suffisso nei file Excel
    
    # GDP deve essere incluso come variabile trimestrale
    target_cc <- colnames(Yq)[target_id]  # es: "GDP_IT"
    
    sel_m[[cc]] <- vars_m_cc
    sel_q[[cc]] <- c(target_cc, vars_q_cc)
  }
  
  list(
    m = sel_m,
    q = sel_q
  )
}


# ==============================================================================
# PREPARE DATA FOR A SINGLE COUNTRY
# ==============================================================================
# path   <- path_data
# covid_mask = TRUE
# cc <- "FR"
# sel_m <- vars_sel$m
# sel_q <- vars_sel$q
prepare_country_data <- function(cc, params, sel_m, sel_q, path,
                                 covid_mask = TRUE) {
  
  # ============================================================
  # 1. LOAD DATA
  # ============================================================
  dm  <- read_excel(file.path(path, paste0(cc, "dataM_LT.xlsx")))
  dq  <- read_excel(file.path(path, paste0(cc, "dataQ_LT.xlsx")))
  leg <- read_excel(file.path(path, paste0(cc, "data.xlsx")), sheet="info")
  
  # ============================================================
  # 2. LEGEND PARSING
  # ============================================================
  var_col <- intersect("Name", names(leg))[1]
  leg_names <- tolower(trimws(leg[[var_col]]))
  leg_base  <- sub(paste0("_", tolower(cc), "$"), "", leg_names)
  leg_class <- leg$Class
  leg_type  <- leg$Type
  leg_freq  <- leg$Frequency
  leg_unb   <- leg$M1
  
  # ============================================================
  # 3. EXTRACT DATE LIMITS (evaluation → end estimation)
  # ============================================================
  start_lim <- params$start_est
  end_lim   <- params$end_eval
  
  # ============================================================
  # 4. EXTRACT MONTHLY DATA
  # ============================================================
  date_m <- as.Date(dm[[1]])
  Xm <- as.matrix(dm[, -1, drop=FALSE][, sel_m[[cc]], drop=FALSE])
  
  # metadata
  base_m  <- sub(paste0("_", cc, "$"), "", sel_m[[cc]])
  class_m <- leg_class[match(tolower(base_m), leg_base)]
  type_m  <- leg_type[match(tolower(base_m), leg_base)]
  freq_m  <- leg_freq[match(tolower(base_m), leg_base)]
  unb_m   <- leg_unb[match(tolower(base_m), leg_base)]
  
  # COVID mask monthly
  if (covid_mask) {
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
  Xq <- as.matrix(dq[, -1, drop=FALSE][, sel_q[[cc]], drop=FALSE])
  
  base_q  <- sub(paste0("_", cc, "$"), "", sel_q[[cc]])
  class_q <- leg_class[match(tolower(base_q), leg_base)]
  type_q  <- leg_type[match(tolower(base_q), leg_base)]
  freq_q  <- leg_freq[match(tolower(base_q), leg_base)]
  unb_q   <- leg_unb[match(tolower(base_q), leg_base)]
  
  # COVID mask quarterly (skip GDP)
  if (covid_mask) {
    Xq_nogdp      <- Xq[, -1, drop=FALSE]
    class_q_nogdp <- class_q[-1]
    res_q <- mask_covid(Xq_nogdp, date_q, class_q_nogdp,
                        params$covid_start, params$covid_end)
    Xq[, -1] <- res_q$data
    mask_q   <- cbind(FALSE, res_q$mask)
  } else {
    mask_q <- matrix(FALSE, nrow(Xq), ncol(Xq))
  }
  
  # ============================================================
  # 6. CREATE A UNIFIED MONTHLY TIME AXIS
  # ============================================================
  date_full <- seq.Date(
    from = min(date_m, date_q),
    to   = max(date_m, date_q),
    by = "month"
  )
  
  
  # ============================================================
  # 7. ALIGN MONTHLY TO FULL DATE AXIS
  # ============================================================
  Xm_full <- matrix(NA, nrow = length(date_full), ncol = ncol(Xm))
  rownames(Xm_full) <- as.character(date_full)
  colnames(Xm_full) <- colnames(Xm)
  Xm_full[as.character(date_m), ] <- Xm
  
  mask_m_full <- matrix(FALSE, nrow = length(date_full), ncol = ncol(Xm))
  rownames(mask_m_full) <- as.character(date_full)
  colnames(mask_m_full) <- colnames(Xm)
  mask_m_full[as.character(date_m), ] <- mask_m
  
  # ============================================================
  # 8. ALIGN QUARTERLY TO FULL DATE AXIS
  # (each quarterly observation placed in its real month)
  # ============================================================
  Xq_full <- matrix(NA, nrow = length(date_full), ncol = ncol(Xq))
  rownames(Xq_full) <- as.character(date_full)
  colnames(Xq_full) <- colnames(Xq)
  Xq_full[as.character(date_q), ] <- Xq
  
  mask_q_full <- matrix(FALSE, nrow = length(date_full), ncol = ncol(Xq))
  rownames(mask_q_full) <- as.character(date_full)
  colnames(mask_q_full) <- colnames(Xq)
  mask_q_full[as.character(date_q), ] <- mask_q
  
  # ============================================================
  # 9. APPLY FINAL TIME WINDOW (start_eval → end_est)
  # ============================================================
  idx <- which(date_full >= start_lim & date_full <= end_lim)
  
  Xm_out     <- Xm_full[idx, , drop=FALSE]
  Xq_out     <- Xq_full[idx, , drop=FALSE]
  mask_m_out <- mask_m_full[idx, , drop=FALSE]
  mask_q_out <- mask_q_full[idx, , drop=FALSE]
  dates_out  <- date_full[idx]
  
  idx_m <- which(date_m >= start_lim & date_m <= end_lim)
  dates_m_out <- date_m[idx_m]
  
  idx_q2 <- which(date_q >= start_lim & date_q <= end_lim)
  dates_q_out <- date_q[idx_q2]
  
  # ============================================================
  # 9B. IDENTIFY TARGET COLUMN (GENERALIZED, NOT ONLY GDP)
  # ============================================================
  
  # Pattern del target, es. "GDP" -> cercherà "gdp_" o "_gdp" o "gdp"
  target_pattern <- paste0("^", tolower(params$target), "_|_", tolower(params$target), "$")
  
  all_series <- c(sel_m[[cc]], sel_q[[cc]])   # nomi completi delle serie usate
  all_series_lower <- tolower(all_series)
  
  # Trova l’indice del target
  target_pos <- grep(target_pattern, all_series_lower)
  
  if (length(target_pos) == 0) {
    stop(paste0("Target '", params$target, "' non trovato tra le serie."))
  }
  
  # Indice finale dentro la matrice Data = [Xm_out | Xq_out]
  target_col <- target_pos
  target_name <- all_series[target_col]
  
  
  # ============================================================
  # 10. FINAL OUTPUT
  # ============================================================
  out <- list(
    Data   = cbind(Xm_out, Xq_out),
    Dates  = dates_out,
    DatesM = dates_m_out,
    DatesQ = dates_q_out,
    Series = c(sel_m[[cc]], sel_q[[cc]]),
    nM     = length(sel_m[[cc]]),
    nQ     = length(sel_q[[cc]]),
    TypeM  = type_m,
    ClassM = class_m,
    ClassQ = class_q,
    TypeQ  = type_q,
    MaskM  = mask_m_out,
    MaskQ  = mask_q_out,
    freq_m = freq_m,
    freq_q = freq_q,
    unb_m  = unb_m,
    unb_q  = unb_q,
    idx_q  = ncol(Xm_out) + 1,
    target_col  = target_col,
    target_name = target_name
    
  )
  
  return(out)
}


# ==============================================================================
# WRAPPER FUNCTION: COMPLETE MULTI-COUNTRY DATA PREPARATION
# ==============================================================================


prepare_all_countries <- function(countries, params, path, covid_mask=TRUE) {
  
  # 1. variable selection country-by-country
  vars_sel <- select_vars(countries, params, path)
  
  # 2. prepare each country
  out <- lapply(countries, function(cc) {
    prepare_country_data(
      cc        = cc,
      params    = params,
      sel_m     = vars_sel$m,
      sel_q     = vars_sel$q,
      path      = path,
      covid_mask = covid_mask
    )
  })
  
  names(out) <- countries
  
  list(
    data   = out,
    sel    = vars_sel
  )
}

# ==============================================================================
# BUILD METRIX-VARIATE DATA
# ==============================================================================

# Rimuove suffissi tipo "_IT", "_FR", "_DE", ecc.
remove_country_code <- function(x, countries) {
  pattern <- paste0("_(", paste(countries, collapse = "|"), ")$")
  sub(pattern, "", x, ignore.case = TRUE)
}

# ==============================================================================
# BUILD MATRIX-VARIATE DATA (TENSOR) FROM PREPARED COUNTRY DATA
# ==============================================================================

build_tensor_dmfm <- function(prep, 
                              countries = names(prep$data),
                              span = c("union", "intersection")) {
  
  span <- match.arg(span)
  data_list <- prep$data
  
  # ------------------------------------------------------------
  # 1. ASSE TEMPORALE GLOBALE
  # ------------------------------------------------------------
  date_list <- lapply(data_list, function(obj) obj$Dates)
  
  if (span == "union") {
    all_dates <- sort(unique(do.call(c, date_list)))
  } else {
    all_dates <- Reduce(intersect, date_list)
    all_dates <- sort(all_dates)
  }
  
  T_global <- length(all_dates)
  
  # ------------------------------------------------------------
  # 2. VARIABILI COMUNI (BASE NAME) + CLASSIFICAZIONE M / Q
  # ------------------------------------------------------------
  # 2a) base name complessive (mensili + trimestrali)
  series_list <- lapply(data_list, function(obj) obj$Series)
  base_series_list <- lapply(series_list, function(v) {
    remove_country_code(v, countries = countries)
  })
  common_base <- Reduce(intersect, base_series_list)
  common_base <- sort(common_base)
  
  if (length(common_base) == 0) {
    stop("Nessuna variabile comune tra i paesi, dopo la selezione country-specific.")
  }
  
  # 2b) base name mensili e trimestrali per ogni paese
  monthly_base_list <- lapply(data_list, function(obj) {
    Series_m <- obj$Series[1:obj$nM]  # prime nM sono mensili
    remove_country_code(Series_m, countries = countries)
  })
  
  quarterly_base_list <- lapply(data_list, function(obj) {
    Series_q <- obj$Series[(obj$nM + 1):(obj$nM + obj$nQ)]  # poi nQ trimestrali
    remove_country_code(Series_q, countries = countries)
  })
  
  # Intersezione solo sulle mensili
  common_M <- Reduce(intersect, monthly_base_list)
  # Intersezione solo sulle trimestrali
  common_Q <- Reduce(intersect, quarterly_base_list)
  
  # Etichetta ogni variabile comune come M / Q / mixed
  var_type <- rep("mixed", length(common_base))
  names(var_type) <- common_base
  
  var_type[common_base %in% common_M] <- "M"
  var_type[common_base %in% common_Q] <- "Q"
  
  n_M      <- sum(var_type == "M")
  n_Q      <- sum(var_type == "Q")
  n_mixed  <- sum(var_type == "mixed")
  
  # ------------------------------------------------------------
  # 3. INIZIALIZZA TENSORE Y (DATI) E W (MASK)
  # ------------------------------------------------------------
  P1 <- length(countries)     # numero di paesi
  P2 <- length(common_base)   # numero di variabili comuni (M + Q + mixed)
  
  Y <- array(NA_real_,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               common_base
             ))
  
  W <- array(0L,
             dim = c(T_global, P1, P2),
             dimnames = list(
               as.character(all_dates),
               countries,
               common_base
             ))
  
  # ------------------------------------------------------------
  # 4. RIEMPI Y E W PAESE PER PAESE
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
    Y_cc <- matrix(NA_real_, nrow = T_global, ncol = P2,
                   dimnames = list(as.character(all_dates), common_base))
    
    for (j in seq_along(common_base)) {
      v_base <- common_base[j]
      
      col_idx <- which(base_cc == v_base)
      if (length(col_idx) == 0) next
      col_idx <- col_idx[1]
      
      Y_cc[idx_time, j] <- Data_cc[, col_idx]
    }
    
    Y[, i, ] <- Y_cc
    W[, i, ] <- ifelse(is.na(Y_cc), 0L, 1L)
  }
  
  # ------------------------------------------------------------
  # 5. OUTPUT
  # ------------------------------------------------------------
  out <- list(
    Y         = Y,               # [T x Paesi x Variabili]
    W         = W,               # [T x Paesi x Variabili]
    dates     = all_dates,
    countries = countries,
    vars      = common_base,
    var_type  = var_type,        # "M", "Q" o "mixed" per ciascuna variabile
    n_M       = n_M,
    n_Q       = n_Q,
    n_mixed   = n_mixed          # dovrebbe essere 0 se non hai casi strani
  )
  
  return(out)
}