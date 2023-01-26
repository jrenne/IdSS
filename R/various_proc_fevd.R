variance.decomp <- function(IRFs){
  # This functions computes the variance decomposition associated to a set
  # of IRFs organised as variables x shocks x horizons x IRF draws (we must
  # have that the number of variables is equal to the number of shocks)
  n <- dim(IRFs)[1]
  H <- dim(IRFs)[3]
  N <- dim(IRFs)[4]
  Variance <- array(0,c(n,n,H,N))
  variance <- array(0,c(n,n,H,N,n))
  Select <- array(0,c(n,n,n))
  for (k in 1:N){
    Variance[,,1,k] = IRFs[,,1,k]%*%t(IRFs[,,1,k])
    for (i in 1:n){
      Select[i,i,i] <- 1
      A <- Select[,,i]
      variance[,,1,k,i] = IRFs[,,1,k]%*%A%*%t(IRFs[,,1,k])
    }
    for (h in 2:H){
      Variance[,,h,k] = IRFs[,,h,k]%*%t(IRFs[,,h,k])+Variance[,,h-1,k]
      for (i in 1:n){
        A <- Select[,,i]
        variance[,,h,k,i] = IRFs[,,h,k]%*%A%*%t(IRFs[,,h,k])+variance[,,h-1,k,i]
      }
    }
  }
  Variance_rep  <- replicate(n,Variance, simplify="array")
  vardecomp = variance / Variance_rep
  return(
    list(
      variance = variance,
      # variance(i,j,h,k,l): contribution of shock l to the forecast error
      # covariance between variable i and j for IRF draw k at horizon h
      # (if i=j then it's the contribution to the variance of i)
      Variance = Variance,
      # Variance(i,j,h,k): forecast error covariance between variable i
      # and j for IRF draw k at horizon h (if i=j then it's the
      # variance of i)
      vardecomp = vardecomp
      # vardecomp(i,j,h,k,l) = variance(i,j,h,k,l)/Variance(i,j,h,k)
    ))
}
