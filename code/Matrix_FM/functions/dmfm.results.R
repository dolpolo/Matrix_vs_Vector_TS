# ==============================================================================
# NOWCAST
# ==============================================================================

rolling_nowcast_dmfm_by_release <- function(Y_full, W_full, k, dates, group, gdp_idx, 
                                            startEV, endEV,
                                            model_type = "dmfm", max_iter = 1000, eps = 1e-3) {
  T <- dim(Y_full)[1]
  p1 <- dim(Y_full)[2]
  p2 <- dim(Y_full)[3]
  
  # Calcola gli indici temporali corrispondenti alle date
  t_start <- which(dates == startEV)
  t_end <- which(dates == endEV)
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("startEV o endEV non trovate nel vettore dates.")
  }
  
  # Output lists per i tre mesi del trimestre
  nowcast_M1 <- list()
  nowcast_M2 <- list()
  nowcast_M3 <- list()
  
  for (t in seq(t_start, t_end)) {
    
    # t = 288
    # t = 289
    # t = 290
    # t = 291 
    # t = 292
    
    cat(sprintf(">>> Rolling nowcast at t = %d (%s)\n", t, as.character(dates[t])))
    
    # Subset e maschere release
    Y_t <- Y_full[1:t, , , drop = FALSE]
    W_t <- W_full[1:t, , , drop = FALSE]
    
    # Standardizza il dataset fino a t
    std <- standardize_Y(Y_t)
    Y_t_std <- std$Y_scaled
    
    # pre_prova1 <- Y_t_std[,1,]
    # pre_prova2 <- Y_t_std[,1,]
    # pre_prova3 <- Y_t_std[,1,]
    # pre_prova4 <- Y_t_std[,1,]
    # pre_prova5 <- Y_t_std[,1,]
    
    # Applica regole di rilascio
    for (j in 1:p2) {
      freq_j <- group[j]
      delay <- freq_j - 1
      release_time <- t - delay
      if (release_time >= 1) {
        if (release_time + 1 <= t) {
          Y_t_std[(release_time + 1):t, , j] <- NA
          W_t[(release_time + 1):t, , j] <- 0
        }
      } else {
        Y_t_std[, , j] <- NA
        W_t[, , j] <- 0
      }
    }
    
    # prova1 <- Y_t_std[ ,1,]
    # prova2 <- Y_t_std[ ,1,]
    # prova3 <- Y_t_std[ ,1,]
    # prova4 <- Y_t_std[ ,1,]
    # prova5 <- Y_t_std[ ,1,]
    
    # EM Estimation
    imp <- mfm.cl(Y_t_std, W_t, k)
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W_full, t = model_type)
    fit <- dmfm.na.em(inputs, Y = Y_t_std, k = k, W = W_t, t = model_type, max.iter = max_iter, eps = eps)
    model <- fit$model
    
    # Ricostruisci il nowcast
    f_t <- model$ff[t, , ]
    Y_hat <- model$R %*% f_t %*% t(model$C)
    Y_hat_true <- Y_hat * std$sd + std$mean
    gdp_forecast <- Y_hat_true[, gdp_idx]
    
    # Salva per il mese del trimestre corrispondente
    m_trimestre <- ((as.integer(format(dates[t], "%m")) - 1) %% 3) + 1
    date_str <- as.character(dates[t])
    if (m_trimestre == 1) {
      nowcast_M1[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 2) {
      nowcast_M2[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 3) {
      nowcast_M3[[date_str]] <- gdp_forecast
    }
  }
  
  return(list(M1 = nowcast_M1, M2 = nowcast_M2, M3 = nowcast_M3))
}


# ==============================================================================
# PREDICTION
# ==============================================================================


# dmfm.predict----
dmfm.predict <- function(., Y, k, W, t){
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  O <- lapply(split(W,seq(nrow(W))), function(x){which(vec(x)==1)})
  if (t=="fm"){
    kout  <- kalman.na(list(t=.$B%x%.$A,
                            l=.$C%x%.$R,
                            q=.$Q%x%.$P,
                            h=.$K%x%.$H),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  } else if (t=="fv"){
    kout  <- kalman.na(list(t=.$BA,
                            l=.$CR,
                            q=.$QP,
                            h=.$KH),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }else{
    kout  <- kalman.na(list(t=.$BA,
                            l=.$C%x%.$R,
                            q=.$QP,
                            h=.$K%x%.$H),
                       apply(Y,1,vec), prod(k),
                       list(f=rep(0, prod(k)),
                            P=diag(prod(k))),
                       O)
  }
  
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
  .$Y  <- array(0, c(n, p[1], p[2]))
  for (i in 1:n){
    .$ff[i,,] <- matrix(kout$ff[,i],k[1],k[2])
    .$fs[i,,] <- matrix(kout$fs[,i],k[1],k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    if (t=="fv"){
      .$Y[i,,] <- matrix(.$CR%*%kout$ff[,i], p[1], p[2])
    }else{
      .$Y[i,,] <- .$R%*%.$ff[i,,]%*%t(.$C)
    }
  }
  
  return(.)
  
}

# ==============================================================================
# NOWCAST CON INTERCETTA (versione corretta)
# ==============================================================================

rolling_nowcast_dmfm_by_release <- function(Y_full, W_full, k, dates, group, gdp_idx, 
                                            startEV, endEV,
                                            model_type = "dmfm", max_iter = 1000, eps = 1e-3) {
  T  <- dim(Y_full)[1]
  p1 <- dim(Y_full)[2]
  p2 <- dim(Y_full)[3]
  
  # Indici temporali
  t_start <- which(dates == startEV)
  t_end   <- which(dates == endEV)
  if (length(t_start) == 0 || length(t_end) == 0) {
    stop("startEV o endEV non trovate nel vettore dates.")
  }
  
  # Output
  nowcast_M1 <- list()
  nowcast_M2 <- list()
  nowcast_M3 <- list()
  
  for (t in seq(t_start, t_end)) {
    cat(sprintf(">>> Rolling nowcast at t = %d (%s)\n", t, as.character(dates[t])))
    
    # Subset e maschere release fino a t
    Y_t <- Y_full[1:t, , , drop = FALSE]
    W_t <- W_full[1:t, , , drop = FALSE]
    
    # Standardizzazione fino a t
    std     <- standardize_Y(Y_t)
    Y_t_std <- std$Y_scaled
    
    # Regole di rilascio (osservazioni future -> NA/0)
    for (j in 1:p2) {
      freq_j <- group[j]
      delay <- freq_j - 1
      release_time <- t - delay
      if (release_time >= 1) {
        if (release_time + 1 <= t) {
          Y_t_std[(release_time + 1):t, , j] <- NA
          W_t[(release_time + 1):t, , j] <- 0
        }
      } else {
        Y_t_std[, , j] <- NA
        W_t[, , j] <- 0
      }
    }
    
    # EM Estimation (usa W_t anche nell'init!)
    imp    <- mfm.cl(Y_t_std, W_t, k)
    inputs <- dmfm.na.2.sv(imp$Y_imputed, k, W_t, t = model_type)   # <-- fix qui
    fit    <- dmfm.na.em(inputs, Y = Y_t_std, k = k, W = W_t, t = model_type,
                         max.iter = max_iter, eps = eps)
    model  <- fit$model
    
    # ---- Ricostruzione nowcast (in scala standardizzata) ----
    # Usa fs[t,,] (allo step finale equivale a f.u_t), e aggiungi l'intercetta alpha_s[t,]
    f_t      <- model$fs[t, , ]                     # <-- fs, non ff
    alpha_t  <- model$alpha_s[t, ]                  # vettore p1
    Y_hat_sd <- (alpha_t %*% t(rep(1, p2))) + model$R %*% f_t %*% t(model$C)
    
    # Destandardizza
    Y_hat_true <- Y_hat_sd * std$sd + std$mean
    
    # Estrai GDP (colonna gdp_idx)
    gdp_forecast <- Y_hat_true[, gdp_idx]           # lunghezza p1
    
    # Salva per mese del trimestre
    m_trimestre <- ((as.integer(format(dates[t], "%m")) - 1) %% 3) + 1
    date_str <- as.character(dates[t])
    if (m_trimestre == 1) {
      nowcast_M1[[date_str]] <- gdp_forecast
    } else if (m_trimestre == 2) {
      nowcast_M2[[date_str]] <- gdp_forecast
    } else {
      nowcast_M3[[date_str]] <- gdp_forecast
    }
  }
  
  list(M1 = nowcast_M1, M2 = nowcast_M2, M3 = nowcast_M3)
}
