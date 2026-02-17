# ==============================================================================
# EM ALGORITHM
# ==============================================================================

# ---------------------------- Dolp Method -------------------------------------
# Modifiche apportate alla funzione dmfm.na.em del Professor Trapin:

# 1. Aggiunta funzione solve_safe() per gestire matrici singolari con rcond + ginv
# 2. Sostituzione solve(R.1) e solve(C.1) con solve_safe(...) --> (ginv(...) se tol = 1e-6 )
# 3. Regolarizzazione di H e K: aggiunto 1e-6 * diag(...) per evitare singularità numerica
# 4. Introduzione di llk.history per tracciare andamento della log-verosimiglianza
# 5. Aggiunto controllo: warning se llk.new < llk.old (EM non dovrebbe peggiorare)
# 6. Stampa dell'evoluzione della log-likelihood e del delta
# 7. Return esteso: output = list(model = ., history = llk.history)

solve_safe <- function(A, tol = 1e-6) {
  if (rcond(A) < tol) {
    return(ginv(A))
  } else {
    return(solve(A))
  }
}

#' . <- inputs
#' Y = Y_aug
#' k = k
#' W = W_aug
#' t = t

dmfm.na.em <- function(., Y, k, W, t, max.iter = 100, eps = 1e-03){
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  M <- matrix(0, k[1]*k[2], k[1]*k[2])
  for (i in 1:k[1]){
    for (j in 1:k[2]){
      i.k1 <- matrix(0, k[1], 1)
      i.k2 <- matrix(0, k[2], 1)
      i.k1[i] <- 1
      i.k2[j] <- 1
      M <- M + (i.k1 %*% t(i.k2)) %x% t(i.k1 %*% t(i.k2))
    }
  }
  
  I1 <- diag(k[1])
  I2 <- diag(k[2])
  
  D.w  <- array(0, c(n, p[1]*p[2], p[1]*p[2]))
  D.wt <- array(0, c(n, p[1]*p[2], p[1]*p[2]))
  for (i in 1:n){
    D.w[i,,] <- diag(as.vector(vec(W[i,,])))
    D.wt[i,,] <- diag(as.vector(vec(t(W[i,,]))))
  }
  
  criterion <- TRUE
  iter <- 0
  llk.history <- data.frame(iter = integer(0), llk_old = numeric(0), llk_new = numeric(0), delta = numeric(0))
  
  while (criterion && iter < max.iter) {
    cat(sprintf(">>> Iterazione EM: %d\n", iter))
    .. <- .
    
    ### MAXIMIZATION STEP ###
    
    Y[is.na(Y)] <- 0
    
    ## MAXIMIZATION STEP ##
    
    ## Loadings (R - C) ##
    R.1 <- matrix(0, p[1]*k[1], p[1]*k[1])
    R.2 <- matrix(0, p[1]*k[1], 1)
    C.1 <- matrix(0, p[2]*k[2], p[2]*k[2])
    C.2 <- matrix(0, p[2]*k[2], 1)
    # for (i in 1:n){
    #   Y[i,,][is.na(Y[i,,])] <- 0
    #   for (r in 1:p[1]){
    #     for (q in 1:p[1]){
    #       R.1 <- R.1 + star(t(..$C)%*%D.w[i,(0:(p[2]-1))*p[1]+r,(0:(p[2]-1))*p[1]+q]%*%solve(..$K)%*%..$C, vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%x%(E.basis(r,q,p[1],p[1])%*%solve(..$H))
    #     }
    #   }
    #   R.2 <- R.2 + vec((W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%..$C%*%t(matrix(.$fs[i,,],k[1],k[2])))
    #   
    #   for (r in 1:p[2]){
    #     for (q in 1:p[2]){
    #       C.1 <- C.1 + star(t(..$R)%*%D.wt[i,(0:(p[1]-1))*p[2]+r,(0:(p[1]-1))*p[2]+q]%*%solve(..$H)%*%..$R, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,])%*%t(M))%x%(E.basis(r,q,p[2],p[2])%*%solve(..$K))
    #     }
    #   }
    #   C.2 <- C.2 + vec(t(W[i,,]*(solve(..$H)%*%Y[i,,]%*%solve(..$K)))%*%..$R%*%matrix(.$fs[i,,],k[1],k[2]))
    # }
    # .$R <- matrix(solve(R.1)%*%R.2, p[1], k[1])
    # .$C <- matrix(solve(C.1)%*%C.2, p[2], k[2])
    
    Y[is.na(Y)] <- 0
    for (i in 1:n){
      for (r in 1:p[1]){
        for (q in 1:p[1]){
          R.1 <- R.1 + star(t(..$C) %*% D.w[i,(0:(p[2]-1))*p[1]+r,(0:(p[2]-1))*p[1]+q] %*% solve(..$K) %*% ..$C,
                            vec(.$fs[i,,]) %*% t(vec(.$fs[i,,])) + .$Ps[i,,]) %x%
            (E.basis(r,q,p[1],p[1]) %*% solve(..$H))
        }
      }
      R.2 <- R.2 + vec((W[i,,] * (solve(..$H) %*% Y[i,,] %*% solve(..$K))) %*% ..$C %*% t(matrix(.$fs[i,,],k[1],k[2])))
    }
    .$R <- matrix(solve_safe(R.1) %*% R.2, p[1], k[1])
    
    for (i in 1:n){
      for (r in 1:p[2]){
        for (q in 1:p[2]){
          C.1 <- C.1 + star(t(.$R) %*% D.wt[i,(0:(p[1]-1))*p[2]+r,(0:(p[1]-1))*p[2]+q] %*% solve(..$H) %*% .$R,
                            M %*% (vec(.$fs[i,,]) %*% t(vec(.$fs[i,,])) + .$Ps[i,,]) %*% t(M)) %x%
            (E.basis(r,q,p[2],p[2]) %*% solve(..$K))
        }
      }
      C.2 <- C.2 + vec(t(W[i,,] * (solve(..$H) %*% Y[i,,] %*% solve(..$K))) %*% .$R %*% matrix(.$fs[i,,],k[1],k[2]))
    }
    .$C <- matrix(solve_safe(C.1) %*% C.2, p[2], k[2])
    
    ## Measurement variance (H - K) ##
    H. <- matrix(0, p[1], p[1])
    K. <- matrix(0, p[2], p[2])
    # for (i in 1:n){
    #   H. <- H. + t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
    #   H. <- H. - t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
    #   H. <- H. - t((W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
    #   H. <- H. + t(star(solve(..$K), D.w[i,,]%*%(.$C%x%.$R)%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(.$C%x%.$R)%*%t(D.w[i,,])))
    #   H. <- H. + (..$H%*%matrix(1,p[1],p[2])%*%..$K)%*%solve(..$K)%*%t(1-W[i,,])
    #   
    #   K. <- K. + t(t(W[i,,]*Y[i,,])%*%solve(..$H)%*%(W[i,,]*Y[i,,]))
    #   K. <- K. - t(t(W[i,,]*Y[i,,])%*%solve(..$H)%*%(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
    #   K. <- K. - t(t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$H)%*%(W[i,,]*Y[i,,]))
    #   K. <- K. + t(star(solve(..$H), D.wt[i,,]%*%(.$R%x%.$C)%*%(M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))%*%t(.$R%x%.$C)%*%t(D.wt[i,,])))
    #   K. <- K. + t(1-W[i,,])%*%solve(..$H)%*%(..$H%*%matrix(1,p[1],p[2])%*%..$K)
    # }
    # .$H <- diag(diag(H.)/(p[2]*n), p[1], p[1])
    # .$K <- diag(diag(K.)/(p[1]*n), p[2], p[2])
    for (i in 1:n){
      H. <- H. + t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
      H. <- H. - t((W[i,,]*Y[i,,])%*%solve(..$K)%*%t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      H. <- H. - t((W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$K)%*%t(W[i,,]*Y[i,,]))
      H. <- H. + t(star(solve(..$K), D.w[i,,]%*%(.$C%x%.$R)%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(.$C%x%.$R)%*%t(D.w[i,,])))
      H. <- H. + (..$H%*%matrix(1,p[1],p[2])%*%..$K)%*%solve(..$K)%*%t(1-W[i,,])
    }
    .$H <- diag(diag(H.)/(p[2]*n), p[1], p[1])
    for (i in 1:n){
      K. <- K. + t(t(W[i,,]*Y[i,,])%*%solve(.$H)%*%(W[i,,]*Y[i,,]))
      K. <- K. - t(t(W[i,,]*Y[i,,])%*%solve(.$H)%*%(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      K. <- K. - t(t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(.$H)%*%(W[i,,]*Y[i,,]))
      K. <- K. + t(star(solve(.$H), D.wt[i,,]%*%(.$R%x%.$C)%*%(M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))%*%t(.$R%x%.$C)%*%t(D.wt[i,,])))
      K. <- K. + t(1-W[i,,])%*%solve(.$H)%*%(.$H%*%matrix(1,p[1],p[2])%*%..$K)
    }
    .$K <- diag(diag(K.)/(p[1]*n), p[2], p[2])
    
    # Regularization to H and K to avoid singularities
    .$H <- .$H + 1e-6 * diag(nrow(.$H))
    .$K <- .$K + 1e-6 * diag(nrow(.$K))
    
    ## VAR parameters (A - B) ##
    if (t=="fm"){
      A.1 <- matrix(0, k[1], k[1])
      A.2 <- matrix(0, k[1], k[1])
      B.1 <- matrix(0, k[2], k[2])
      B.2 <- matrix(0, k[2], k[2])
      # for (i in 2:n){
      #   A.1 <- A.1 + star(solve(..$Q)%*%..$B, vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])
      #   A.2 <- A.2 + star(t(..$B)%*%solve(..$Q)%*%..$B, vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])
      #   B.1 <- B.1 + star(solve(..$P)%*%..$A, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(M))
      #   B.2 <- B.2 + star(t(..$A)%*%solve(..$P)%*%..$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,]))+.$Ps[i-1,,])%*%t(M))
      # }
      # .$A <- A.1%*%solve(A.2)
      # .$B <- B.1%*%solve(B.2)
      for (i in 2:n){
        A.1 <- A.1 + star(solve(..$Q)%*%..$B, vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])
        A.2 <- A.2 + star(t(..$B)%*%solve(..$Q)%*%..$B, vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])
      }
      .$A <- A.1%*%solve(A.2)
      for (i in 2:n){
        B.1 <- B.1 + star(solve(..$P)%*%.$A, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(M))
        B.2 <- B.2 + star(t(.$A)%*%solve(..$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,]))+.$Ps[i-1,,])%*%t(M))
      }
      .$B <- B.1%*%solve(B.2)
    }else{
      BA.1 <- matrix(0, prod(k), prod(k))
      BA.2 <- matrix(0, prod(k), prod(k))
      for (i in 2:n){
        BA.1 <- BA.1 + vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]
        BA.2 <- BA.2 + vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]
      }
      .$BA <- BA.1%*%solve(BA.2)
    }
    
    ## State variance (P - Q) ##
    if (t=="fm"){
      P. <- matrix(0, k[1], k[1])
      Q. <- matrix(0, k[2], k[2])
      # for (i in 2:n){
      #   P.1 <- star(solve(..$Q), (vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,]) )
      #   P.2 <- star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]))%*%t(.$A)
      #   P.3 <- .$A%*%t(star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])))
      #   P.4 <- .$A%*%star(t(.$B)%*%solve(..$Q)%*%.$B, (vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]))%*%t(.$A)
      #   P.  <- P. + P.1 - P.2 - P.3 + P.4
      #   Q.1 <- star(solve(..$P), M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))
      #   Q.2 <- star(solve(..$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M)))%*%t(.$B)
      #   Q.3 <- .$B%*%t(star(solve(..$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M))))
      #   Q.4 <- .$B%*%star(t(.$A)%*%solve(..$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])%*%t(M))%*%t(.$B)
      #   Q.  <- Q. + Q.1 - Q.2 - Q.3 + Q.4
      # }
      # .$P <- P./((n-1)*k[2])
      # .$Q <- Q./((n-1)*k[1])
      for (i in 2:n){
        P.1 <- star(solve(..$Q), (vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,]) )
        P.2 <- star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]))%*%t(.$A)
        P.3 <- .$A%*%t(star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])))
        P.4 <- .$A%*%star(t(.$B)%*%solve(..$Q)%*%.$B, (vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]))%*%t(.$A)
        P.  <- P. + P.1 - P.2 - P.3 + P.4
      }
      .$P <- P./((n-1)*k[2])
      for (i in 2:n){
        Q.1 <- star(solve(.$P), M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))
        Q.2 <- star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M)))%*%t(.$B)
        Q.3 <- .$B%*%t(star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M))))
        Q.4 <- .$B%*%star(t(.$A)%*%solve(.$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])%*%t(M))%*%t(.$B)
        Q.  <- Q. + Q.1 - Q.2 - Q.3 + Q.4
      }
      .$Q <- Q./((n-1)*k[1])
    }else{
      QP <- matrix(0, prod(k), prod(k))
      for (i in 2:n){
        QP <- QP + vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,] - (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(.$BA)
      }
      .$QP <- QP/(n-1)
    }
    
    # --- Update Q or QP ---
    if (t != "fm") {
      .$QP <- .$QP + 1e-6 * diag(nrow(.$QP))
    }
    ### EXPECTATION STEP ###
    O <- lapply(split(W, seq(nrow(W))), function(x){which(vec(x) == 1)})
    
    if (t == "fm") {
      kout <- kalman.na(list(
        t = .$B %x% .$A,
        l = .$C %x% .$R,
        q = .$Q %x% .$P,
        h = .$K %x% .$H
      ), apply(Y,1,vec), prod(k), list(f = vec(.$fs[1,,]), P = diag(prod(k))), O)
    } else {
      kout <- kalman.na(list(
        t = .$BA,
        l = .$C %x% .$R,
        q = .$QP,
        h = .$K %x% .$H
      ), apply(Y,1,vec), prod(k), list(f = vec(.$fs[1,,]), P = diag(prod(k))), O)
    }
    
    .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
    .$Pf <- .$Ps <- .$Cs <- array(0, c(n, prod(k), prod(k)))
    for (i in 1:n){
      .$ff[i,,] <- matrix(kout$ff[,i], k[1], k[2])
      .$fs[i,,] <- matrix(kout$fs[,i], k[1], k[2])
      .$Pf[i,,] <- kout$Pf[,,i]
      .$Ps[i,,] <- kout$Ps[,,i]
      .$Cs[i,,] <- kout$Cs[,,i]
    }
    
    ### CONVERGENCE CHECK ###
    llk.old <- dmfm.em.llk(.., Y, t)
    llk.new <- dmfm.em.llk(., Y, t)
    delta <- 2 * abs(llk.new - llk.old) / abs(llk.new + llk.old)
    
    llk.history <- rbind(llk.history, data.frame(iter = iter, llk_old = llk.old, llk_new = llk.new, delta = delta))
    
    cat(sprintf("Log-likelihood old: %.4f, new: %.4f, delta: %.6f\n", llk.old, llk.new, delta))
    
    if (llk.new < llk.old) {
      warning(sprintf("Likelihood decreased at iteration %d!", iter))
    }
    
    if (delta < eps) {
      criterion <- FALSE
      .$Y <- array(0, c(n, p[1], p[2]))
      for (i in 1:n){
        .$Y[i,,] <- .$R %*% .$fs[i,,] %*% t(.$C)
      }
    }
    
    iter <- iter + 1
  }
  
  return(list(model = ., history = llk.history))
}

# ==============================================================================
# CONVERGENCE CRITERIUM
# ==============================================================================

# dmfm.em.llk ----
# Safe log-likelihood for DMFM EM algorithm
# Prediction error gaussian log-likelihood used to control convergence
# usa as.matrix per sicurezza su .$Pf[i,,]

dmfm.em.llk <- function(., Y, t) {
  
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  
  if (t == "fv") {
    # Full vectorized version (non la usi quasi mai)
    V <- array(0, c(n, prod(p), prod(p)))
    m <- matrix(0, n, prod(p))
    llk <- rep(0, n)
    
    for (i in 1:n) {
      V[i,,] <- (.$C %x% .$R) %*% .$Pf[i,,] %*% t(.$C %x% .$R) + .$K %x% .$H
      m[i,]  <- (.$C %x% .$R) %*% vec(.$ff[i,,])
      llk[i] <- log(det(V[i,,])) + t(vec(Y[i,,]) - m[i,]) %*% solve(V[i,,]) %*% (vec(Y[i,,]) - m[i,])
    }
    
  } else {
    # Factorized version (quella che ti interessa)
    V <- array(0, c(n, prod(p), prod(p)))
    m <- array(0, c(n, p))
    d <- rep(0, n)
    llk <- rep(0, n)
    
    for (i in 1:n) {
      # Assicuro che .$Pf[i,,] sia una matrice numerica
      Pf_i <- as.matrix(.$Pf[i,,])
      
      # Costruisco la matrice di varianza predizione
      V[i,,] <- (.$C %x% .$R) %*% Pf_i %*% t(.$C %x% .$R) + .$K %x% .$H
      
      # Costruisco il valore medio predetto
      m[i,,] <- .$R %*% .$ff[i,,] %*% t(.$C)
      
      # Parte log-determinante per decomposizione
      d[i] <- log(det(solve(Pf_i) + (t(.$C) %*% solve(.$K) %*% .$C) %x% (t(.$R) %*% solve(.$H) %*% .$R))) +
        log(det(Pf_i)) + p[1] * log(det(.$K)) + p[2] * log(det(.$H))
      
      # Log-likelihood della previsione
      llk[i] <- d[i] + t(vec(Y[i,,] - m[i,,])) %*% solve(V[i,,]) %*% vec(Y[i,,] - m[i,,])
    }
  }
  
  out <- -0.5 * sum(llk)
  
  return(out)
}


# ==============================================================================
# EM ALGORITHM FOR DFM 
# ==============================================================================

dmfm.na.em.vec <- function(., Y, k, W, t = "dmfm", max.iter = 1000, eps = 1e-4) {
  n <- dim(Y)[1]
  p <- dim(Y)[-1]  # p[1] = 1 for DFM case
  
  M <- matrix(0, k[1] * k[2], k[1] * k[2])
  for (i in 1:k[1]) {
    for (j in 1:k[2]) {
      i.k1 <- matrix(0, k[1], 1)
      i.k2 <- matrix(0, k[2], 1)
      i.k1[i] <- 1
      i.k2[j] <- 1
      M <- M + (i.k1 %*% t(i.k2)) %x% t(i.k1 %*% t(i.k2))
    }
  }
  
  # Kalman filter initialization
  O <- lapply(1:dim(W)[1], function(t) which(W[t, 1, ] == 1))
  kout <- kalman.na.vec(list(
    t = .$BA,
    l = .$C,
    q = .$QP,
    h = .$K
  ),
  apply(Y, 1, vec), k[2],
  list(f = rep(0, k[2]), P = diag(k[2])), O)
  
  # Prepare storage
  .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
  .$Pf <- .$Ps <- .$Cs <- array(0, c(n, k[2], k[2]))
  .$Y  <- array(0, c(n, p[1], p[2]))
  
  for (i in 1:n) {
    .$ff[i,,] <- matrix(kout$ff[, i], k[1], k[2])
    .$fs[i,,] <- matrix(kout$fs[, i], k[1], k[2])
    .$Pf[i,,] <- kout$Pf[,,i]
    .$Ps[i,,] <- kout$Ps[,,i]
    .$Cs[i,,] <- kout$Cs[,,i]
    .$Y[i,,] <- matrix(.$fs[i,,] %*% t(.$C), 1, p[2])  # prediction
  }
  
  # Likelihood monitoring loop
  llk.history <- data.frame(iter = integer(0), llk_old = numeric(0), llk_new = numeric(0), delta = numeric(0))
  criterion <- TRUE
  iter <- 0
  
  while (criterion && iter < max.iter) {
    cat(sprintf("EM iteration: %d\n", iter))
    .. <- .
    
    # Update measurement matrix C
    C.1 <- matrix(0, p[2] * k[2], p[2] * k[2])
    C.2 <- matrix(0, p[2] * k[2], 1)
    for (i in 1:n) {
      fs_vec <- vec(.$fs[i,,])
      C.1 <- C.1 + kronecker(fs_vec %*% t(fs_vec) + .$Ps[i,,], diag(1))
      C.2 <- C.2 + kronecker(fs_vec, diag(1)) %*% vec(Y[i,,])
    }
    .$C <- matrix(solve_safe(C.1) %*% C.2, p[2], k[2])
    
    # Update measurement noise K
    K.sum <- matrix(0, p[2], p[2])
    for (i in 1:n) {
      Y.pred <- .$fs[i,,] %*% t(.$C)
      resid <- Y[i,,] - matrix(Y.pred, nrow = 1)
      K.sum <- K.sum + t(resid) %*% resid
    }
    .$K <- diag(diag(K.sum / n)) + 1e-6 * diag(p[2])
    
    # Update transition matrix BA and QP
    BA.num <- matrix(0, k[2], k[2])
    BA.den <- matrix(0, k[2], k[2])
    for (i in 2:n) {
      BA.num <- BA.num + t(.$fs[i - 1,,]) %*% .$fs[i,,]
      BA.den <- BA.den + t(.$fs[i - 1,,]) %*% .$fs[i - 1,,] + .$Ps[i - 1,,]
    }
    .$BA <- t(BA.num %*% solve_safe(BA.den))
    
    QP.sum <- matrix(0, k[2], k[2])
    for (i in 2:n) {
      pred <- .$fs[i - 1,,] %*% t(.$BA)
      resid <- .$fs[i,,] - pred
      QP.sum <- QP.sum + t(resid) %*% resid + .$Ps[i,,]
    }
    .$QP <- QP.sum / (n - 1) + 1e-6 * diag(k[2])
    
    # Expectation step
    kout <- kalman.na(list(
      t = .$BA,
      l = .$C,
      q = .$QP,
      h = .$K
    ),
    apply(Y, 1, vec), k[2],
    list(f = vec(.$fs[1,,]), P = diag(k[2])), O)
    
    for (i in 1:n) {
      .$ff[i,,] <- matrix(kout$ff[, i], k[1], k[2])
      .$fs[i,,] <- matrix(kout$fs[, i], k[1], k[2])
      .$Pf[i,,] <- kout$Pf[,,i]
      .$Ps[i,,] <- kout$Ps[,,i]
      .$Cs[i,,] <- kout$Cs[,,i]
    }
    
    # Convergence check
    llk.old <- dmfm.em.llk.vec(.., Y)
    llk.new <- dmfm.em.llk.vec(., Y)
    delta <- 2 * abs(llk.new - llk.old) / abs(llk.new + llk.old)
    
    llk.history <- rbind(llk.history, data.frame(iter = iter, llk_old = llk.old, llk_new = llk.new, delta = delta))
    cat(sprintf("LogLik: old = %.4f, new = %.4f, delta = %.6f\n", llk.old, llk.new, delta))
    
    if (delta < eps) {
      criterion <- FALSE
      for (i in 1:n) {
        .$Y[i,,] <- matrix(.$fs[i,,] %*% t(.$C), 1, p[2])
      }
    }
    
    iter <- iter + 1
  }
  
  return(list(model = ., history = llk.history))
}

# ==============================================================================
# CONVERGENCE CRITERIUM FOR DFM
# ==============================================================================

dmfm.em.llk.vec <- function(., Y) {
  n <- dim(Y)[1]
  p2 <- dim(Y)[3]  # p[2] nel caso vettoriale
  k <- dim(.$ff)[3]
  
  llk <- numeric(n)
  
  for (i in 1:n) {
    Pf_i <- as.matrix(.$Pf[i,,])
    
    # Predicted variance (dimensione: p2 x p2)
    V_i <- .$C %*% Pf_i %*% t(.$C) + .$K
    
    # Predicted mean (dimensione: 1 x p2)
    m_i <- .$ff[i,,] %*% t(.$C)
    
    # Observation
    y_i <- matrix(Y[i,,], nrow = 1)
    
    # Log-likelihood contribution
    residual <- y_i - m_i
    llk[i] <- log(det(V_i)) + residual %*% solve(V_i) %*% t(residual)
  }
  
  out <- -0.5 * sum(llk)
  return(out)
}


# ==============================================================================
# EM ALGORITHM WITH INTERCETPT
# ==============================================================================

# -------------------------------------------------------------------
# Utilities
# -------------------------------------------------------------------
solve_safe <- function(A, tol = 1e-6) {
  if (rcond(A) < tol) ginv(A) else solve(A)
}

# -------------------------------------------------------------------
# EM algorithm for DMFM with country-specific AR(1) intercepts
# -------------------------------------------------------------------
dmfm.na.em <- function(., Y, k, W, t, max.iter = 100, eps = 1e-03){
  # Dimensions
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  p1 <- p[1]; p2 <- p[2]; r <- prod(k)
  
  # Precompute M (as in original)
  M <- matrix(0, r, r)
  for (i in 1:k[1]){
    for (j in 1:k[2]){
      i.k1 <- matrix(0, k[1], 1); i.k1[i] <- 1
      i.k2 <- matrix(0, k[2], 1); i.k2[j] <- 1
      M <- M + (i.k1 %*% t(i.k2)) %x% t(i.k1 %*% t(i.k2))
    }
  }
  I1 <- diag(k[1]); I2 <- diag(k[2])
  
  # Masks (as in original)
  D.w  <- array(0, c(n, p1*p2, p1*p2))
  D.wt <- array(0, c(n, p1*p2, p1*p2))
  for (i in 1:n){
    D.w[i,,]  <- diag(as.vector(vec(W[i,,])))
    D.wt[i,,] <- diag(as.vector(vec(t(W[i,,]))))
  }
  
  criterion <- TRUE
  iter <- 0
  llk.history <- data.frame(iter = integer(0), llk_old = numeric(0),
                            llk_new = numeric(0), delta = numeric(0))
  
  while (criterion && iter < max.iter) {
    cat(sprintf(">>> Iterazione EM: %d\n", iter))
    .. <- .   # copy previous iterate
    
    # -------------------------
    # MAXIMIZATION STEP
    # -------------------------
    # Replace NAs with zeros (as in original, masking W handles availability)
    Y0 <- Y
    Y0[is.na(Y0)] <- 0
    
    # Build Y_tilde = Y - alpha_prev * 1'  (cruciale per R,C,H,K!)
    if (is.null(..$alpha_s)) {
      alpha_prev <- matrix(0, n, p1)
    } else {
      alpha_prev <- ..$alpha_s
    }
    Y_tilde <- array(0, c(n, p1, p2))
    one_p2 <- rep(1, p2)
    for (i in 1:n) {
      # sottraggo l'intercetta alla Y già con NAs->0
      Y_tilde[i,,] <- Y0[i,,] - alpha_prev[i, ] %*% t(one_p2)
    }
    
    ## ---- Loadings (R, C) ----
    R.1 <- matrix(0, p1*k[1], p1*k[1]); R.2 <- matrix(0, p1*k[1], 1)
    C.1 <- matrix(0, p2*k[2], p2*k[2]); C.2 <- matrix(0, p2*k[2], 1)
    
    for (i in 1:n){
      for (r1 in 1:p1){
        for (q1 in 1:p1){
          R.1 <- R.1 + star(t(..$C) %*% D.w[i,(0:(p2-1))*p1+r1,(0:(p2-1))*p1+q1] %*%
                              solve(..$K) %*% ..$C,
                            vec(.$fs[i,,]) %*% t(vec(.$fs[i,,])) + .$Ps[i,,]) %x%
            (E.basis(r1,q1,p1,p1) %*% solve(..$H))
        }
      }
      R.2 <- R.2 + vec((W[i,,] * (solve(..$H) %*% Y_tilde[i,,] %*% solve(..$K))) %*%
                         ..$C %*% t(matrix(.$fs[i,,],k[1],k[2])))
    }
    .$R <- matrix(solve_safe(R.1) %*% R.2, p1, k[1])
    
    for (i in 1:n){
      for (r2 in 1:p2){
        for (q2 in 1:p2){
          C.1 <- C.1 + star(t(.$R) %*% D.wt[i,(0:(p1-1))*p2+r2,(0:(p1-1))*p2+q2] %*%
                              solve(..$H) %*% .$R,
                            M %*% (vec(.$fs[i,,]) %*% t(vec(.$fs[i,,])) + .$Ps[i,,]) %*% t(M)) %x%
            (E.basis(r2,q2,p2,p2) %*% solve(..$K))
        }
      }
      C.2 <- C.2 + vec(t(W[i,,] * (solve(..$H) %*% Y_tilde[i,,] %*% solve(..$K))) %*%
                         .$R %*% matrix(.$fs[i,,],k[1],k[2]))
    }
    .$C <- matrix(solve_safe(C.1) %*% C.2, p2, k[2])
    
    ## ---- Measurement variance (H, K) ----
    H. <- matrix(0, p1, p1); K. <- matrix(0, p2, p2)
    
    for (i in 1:n){
      # tutte le occorrenze con Y -> Y_tilde (importantissimo)
      H. <- H. + t((W[i,,]*Y_tilde[i,,])%*%solve(..$K)%*%t(W[i,,]*Y_tilde[i,,]))
      H. <- H. - t((W[i,,]*Y_tilde[i,,])%*%solve(..$K)%*%t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      H. <- H. - t((W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(..$K)%*%t(W[i,,]*Y_tilde[i,,]))
      H. <- H. + t(star(solve(..$K),
                        D.w[i,,]%*%(.$C%x%.$R)%*%
                          (vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%
                          t(.$C%x%.$R)%*%t(D.w[i,,])))
      H. <- H. + (..$H%*%matrix(1,p1,p2)%*%..$K)%*%solve(..$K)%*%t(1-W[i,,])
    }
    .$H <- diag(diag(H.)/(p2*n), p1, p1)
    
    for (i in 1:n){
      K. <- K. + t(t(W[i,,]*Y_tilde[i,,])%*%solve(.$H)%*%(W[i,,]*Y_tilde[i,,]))
      K. <- K. - t(t(W[i,,]*Y_tilde[i,,])%*%solve(.$H)%*%(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C))))
      K. <- K. - t(t(W[i,,]*(.$R%*%.$fs[i,,]%*%t(.$C)))%*%solve(.$H)%*%(W[i,,]*Y_tilde[i,,]))
      K. <- K. + t(star(solve(.$H),
                        D.wt[i,,]%*%(.$R%x%.$C)%*%
                          (M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))%*%
                          t(.$R%x%.$C)%*%t(D.wt[i,,])))
      K. <- K. + t(1-W[i,,])%*%solve(.$H)%*%(.$H%*%matrix(1,p1,p2)%*%..$K)
    }
    .$K <- diag(diag(K.)/(p1*n), p2, p2)
    
    # Regularization (numeric stability)
    .$H <- .$H + 1e-6 * diag(nrow(.$H))
    .$K <- .$K + 1e-6 * diag(nrow(.$K))
    
    ## ---- VAR parameters (A,B vs BA) ----
    if (t == "fm"){
      A.1 <- matrix(0, k[1], k[1]); A.2 <- matrix(0, k[1], k[1])
      B.1 <- matrix(0, k[2], k[2]); B.2 <- matrix(0, k[2], k[2])
      for (i in 2:n){
        A.1 <- A.1 + star(solve(..$Q)%*%..$B, vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])
        A.2 <- A.2 + star(t(..$B)%*%solve(..$Q)%*%..$B, vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])
      }
      .$A <- A.1%*%solve(A.2)
      for (i in 2:n){
        B.1 <- B.1 + star(solve(..$P)%*%.$A, M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(M))
        B.2 <- B.2 + star(t(.$A)%*%solve(..$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,]))+.$Ps[i-1,,])%*%t(M))
      }
      .$B <- B.1%*%solve(B.2)
    } else {
      BA.1 <- matrix(0, r, r); BA.2 <- matrix(0, r, r)
      for (i in 2:n){
        BA.1 <- BA.1 + vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]
        BA.2 <- BA.2 + vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]
      }
      .$BA <- BA.1 %*% solve(BA.2)
    }
    
    ## ---- State variance (P,Q) or QP ----
    if (t == "fm"){
      P. <- matrix(0, k[1], k[1]); Q. <- matrix(0, k[2], k[2])
      for (i in 2:n){
        P.1 <- star(solve(..$Q), (vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,]) )
        P.2 <- star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,]))%*%t(.$A)
        P.3 <- .$A%*%t(star(solve(..$Q)%*%.$B, (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])))
        P.4 <- .$A%*%star(t(.$B)%*%solve(..$Q)%*%.$B, (vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,]))%*%t(.$A)
        P.  <- P. + P.1 - P.2 - P.3 + P.4
      }
      .$P <- P./((n-1)*k[2])
      for (i in 2:n){
        Q.1 <- star(solve(.$P), M%*%(vec(.$fs[i,,])%*%t(vec(.$fs[i,,]))+.$Ps[i,,])%*%t(M))
        Q.2 <- star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M)))%*%t(.$B)
        Q.3 <- .$B%*%t(star(solve(.$P)%*%.$A, M%*%((vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,]))+.$Cs[i-1,,])%*%t(M))))
        Q.4 <- .$B%*%star(t(.$A)%*%solve(.$P)%*%.$A, M%*%(vec(.$fs[i-1,,])%*%t(vec(.$fs[i-1,,])) + .$Ps[i-1,,])%*%t(M))%*%t(.$B)
        Q.  <- Q. + Q.1 - Q.2 - Q.3 + Q.4
      }
      .$Q <- Q./((n-1)*k[1])
    } else {
      QP <- matrix(0, r, r)
      for (i in 2:n){
        QP <- QP + vec(.$fs[i,,])%*%t(vec(.$fs[i,,])) + .$Ps[i,,] -
          (vec(.$fs[i,,])%*%t(vec(.$fs[i-1,,])) + .$Cs[i-1,,])%*%t(.$BA)
      }
      .$QP <- QP/(n-1)
      .$QP <- .$QP + 1e-6 * diag(nrow(.$QP))  # regularization
    }
    
    # -------------------------
    # EXPECTATION STEP
    # -------------------------
    n <- dim(Y0)[1]  # assicura che n sia definito qui
    O <- lapply(seq_len(n), function(i) which(vec(W[i,,]) == 1))
    
    # New intercept dynamics params (init if missing)
    if (is.null(.$rho))     .$rho     <- 0.98
    if (is.null(.$Q_alpha)) .$Q_alpha <- diag(1e-4, p1)
    
    # Build extended state-space: x_t = [alpha_t; f_t]
    Z_alpha <- kronecker(rep(1, p2), diag(p1))   # (p1*p2) x p1
    Z_f     <- (.$C %x% .$R)                     # (p1*p2) x r
    Z_ext   <- cbind(Z_alpha, Z_f)               # (p1*p2) x (p1 + r)
    
    T_ext <- Matrix::bdiag(.$rho * diag(p1), .$BA)
    Q_ext <- Matrix::bdiag(.$Q_alpha, .$QP)
    H_ext <- .$K %x% .$H
    
    state_dim <- p1 + r
    P0_block  <- Matrix::bdiag(1e6 * diag(p1), diag(r))
    start_state <- list(f = rep(0, state_dim),
                        P = as.matrix(P0_block))     # coercizione a base matrix
    
    # Coerce to base matrices for kalman.na
    T_ext_m <- as.matrix(T_ext)
    Q_ext_m <- as.matrix(Q_ext)
    
    kout <- kalman.na(list(
      t = T_ext_m,
      l = Z_ext,
      q = Q_ext_m,
      h = H_ext
    ), apply(Y0, 1, vec), state_dim, start_state, O)
    
    # Split outputs (save both full and factor-only covariances)
    .$alpha_f <- .$alpha_s <- matrix(0, n, p1)
    .$ff <- .$fs <- array(0, c(n, k[1], k[2]))
    
    .$Pf_full <- kout$Pf
    .$Ps_full <- kout$Ps
    .$Cs_full <- kout$Cs
    
    .$Pf <- .$Ps <- .$Cs <- array(0, c(n, r, r))
    
    for (i in 1:n){
      .$alpha_f[i, ] <- kout$ff[1:p1, i]
      .$alpha_s[i, ] <- kout$fs[1:p1, i]
      .$ff[i,,]      <- matrix(kout$ff[-(1:p1), i], k[1], k[2])
      .$fs[i,,]      <- matrix(kout$fs[-(1:p1), i], k[1], k[2])
      
      .$Pf[i,,] <- kout$Pf[-(1:p1), -(1:p1), i]
      .$Ps[i,,] <- kout$Ps[-(1:p1), -(1:p1), i]
      .$Cs[i,,] <- kout$Cs[-(1:p1), -(1:p1), i]
    }
    
    # ---- Update rho and Q_alpha using smoothed moments (EM-consistent) ----
    idx_a <- 1:p1
    S11 <- matrix(0, p1, p1)  # sum E[a_{t-1} a_{t-1}']
    S10 <- matrix(0, p1, p1)  # sum E[a_t a_{t-1}']
    S00 <- matrix(0, p1, p1)  # sum E[a_t a_t']
    
    for (tix in 2:n){
      Ps_t   <- .$Ps_full[,,tix]
      Ps_tm1 <- .$Ps_full[,,tix-1]
      Cs_tm1 <- .$Cs_full[,,tix-1]
      
      E_at_at     <- Ps_t[idx_a, idx_a]    + tcrossprod(.$alpha_s[tix, ])
      E_atm1_atm1 <- Ps_tm1[idx_a, idx_a]  + tcrossprod(.$alpha_s[tix-1, ])
      E_at_atm1   <- Cs_tm1[idx_a, idx_a]  + .$alpha_s[tix, ] %o% .$alpha_s[tix-1, ]
      
      S00 <- S00 + E_at_at
      S11 <- S11 + E_atm1_atm1
      S10 <- S10 + E_at_atm1
    }
    
    rho_hat <- sum(diag(S10)) / max(1e-12, sum(diag(S11)))   # scalar rho
    .$rho <- as.numeric(rho_hat)
    
    Qalpha_mat <- S00 - .$rho * t(S10) - S10 * .$rho + (.$rho^2) * S11
    .$Q_alpha  <- diag(pmax(1e-10, diag(Qalpha_mat/(n-1))))  # keep diagonal, non-negative
    
    # -------------------------
    # CONVERGENCE CHECK
    # -------------------------
    llk.old <- dmfm.em.llk(.., Y0, t)
    llk.new <- dmfm.em.llk(.,  Y0, t)
    delta <- 2 * abs(llk.new - llk.old) / max(1e-12, abs(llk.new + llk.old))
    
    llk.history <- rbind(llk.history,
                         data.frame(iter = iter, llk_old = llk.old,
                                    llk_new = llk.new, delta = delta))
    cat(sprintf("Log-likelihood old: %.4f, new: %.4f, delta: %.6f\n",
                llk.old, llk.new, delta))
    
    if (llk.new < llk.old) {
      warning(sprintf("Likelihood decreased at iteration %d!", iter))
    }
    
    if (delta < eps) {
      criterion <- FALSE
      .$Y <- array(0, c(n, p1, p2))
      for (i in 1:n){
        .$Y[i,,] <- (.$alpha_s[i, ] %*% t(one_p2)) + .$R %*% .$fs[i,,] %*% t(.$C)
      }
    }
    
    iter <- iter + 1
  }
  
  return(list(model = ., history = llk.history))
}

# -------------------------------------------------------------------
# Log-likelihood used for EM convergence (with intercept)
# -------------------------------------------------------------------
dmfm.em.llk <- function(., Y, t) {
  n <- dim(Y)[1]
  p <- dim(Y)[-1]
  p1 <- p[1]; p2 <- p[2]
  
  V <- array(0, c(n, prod(p), prod(p)))
  m <- array(0, c(n, p1, p2))
  d <- rep(0, n)
  llk <- rep(0, n)
  
  for (i in 1:n) {
    Pf_i <- as.matrix(.$Pf[i,,])
    
    V[i,,] <- (.$C %x% .$R) %*% Pf_i %*% t(.$C %x% .$R) + .$K %x% .$H
    
    # mean prediction includes alpha
    m[i,,] <- (.$alpha_s[i, ] %*% t(rep(1, p2))) + .$R %*% .$ff[i,,] %*% t(.$C)
    
    d[i] <- log(det(solve(Pf_i) + (t(.$C) %*% solve(.$K) %*% .$C) %x%
                      (t(.$R) %*% solve(.$H) %*% .$R))) +
      log(det(Pf_i)) + p1 * log(det(.$K)) + p2 * log(det(.$H))
    
    llk[i] <- d[i] + t(vec(Y[i,,] - m[i,,])) %*% solve(V[i,,]) %*% vec(Y[i,,] - m[i,,])
  }
  
  out <- -0.5 * sum(llk)
  return(out)
}

