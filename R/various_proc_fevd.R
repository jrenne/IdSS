variance.decomp <- function(IRFs){
  # This functions computes the variance decomposition associated to a set
  # of IRFs organised as variables x shocks x horizons x IRF draws (we must
  # have that the number of variables is equal to the number of shocks)
  # Hence, IRFs has 4 dimensions. If there is only one draw (e.g. estimated model),
  #   then IRFs has only 3 dimensions.

  n <- dim(IRFs)[1]
  H <- dim(IRFs)[3]

  if(length(dim(IRFs))<4){
    # In that case, add the fourth dimension (remains 1)
    N <- 1
    IRFs <- array(IRFs,c(n,n,H,N))
  }

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

  if(N==1){
    variance  <- array(variance,c(n,n,H,n))
    Variance  <- array(Variance,c(n,n,H))
    vardecomp <- array(vardecomp,c(n,n,H,n))
  }

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
