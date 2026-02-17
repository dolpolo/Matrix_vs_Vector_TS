tr <- function(x){
  sum(diag(x))
}

TR <- function(A,A.h){
  tr((t(A)%*%A.h)%*%solve(t(A.h)%*%A.h)%*%(t(A.h)%*%A))/tr(t(A)%*%A)
}

DM <- function(A,A.h){
  norm(A.h%*%solve(t(A.h)%*%A.h)%*%t(A.h)-A%*%solve(t(A)%*%A)%*%t(A), "2")
}

star <- function(X,Y){
  m <- dim(X)[1]
  n <- dim(X)[2]
  p <- dim(Y)[1]/m
  q <- dim(Y)[2]/n
  out <- matrix(0, p, q)
  for (i in 1:m){
    for (j in 1:n){
      out <- out + X[i,j]*Y[(i-1)*p+(1:p),(j-1)*q+(1:q)]
    }
  }
  return(out)
}

O.hat <- function(Y,i,j,h){
  n <- dim(Y)[1]
  O <- 0
  for (t in 1:(n-h))
    O <- O + Y[t,,i]%*%t(Y[t+h,,j])
  return(O/(n-h))
}

E.basis <- function(i,j,p,q){
  out <- matrix(0, p, q)
  out[i,j] <- 1
  return(out)
}

