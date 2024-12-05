
svar.fevmax <- function(Y,p,
                        nb.shocks,
                        names.of.shocks,
                        H1,
                        H2,
                        variable,
                        norm,
                        nb.periods.IRF,
                        bootstrap.replications, # This is used in the parametric bootstrap only
                        confidence.interval = 0.90, # expressed in pp.
                        indic.plot = 1 # Plots are displayed if = 1.
){
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # FEV maximization

  names.of.variables <- colnames(Y)

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  E <-list()
  for(i in 1:nb.shocks){
    E[[i]] <- array(0,c(n,n))
    E[[i]][variable[i],variable[i]] <- 1
  }

  ######
  y0.star <- rep(0,n*p)
  stdv.IRFs <- list()
  # do Choleski (we need the FULL IRFs)
  cholesky.res <- svar.ordering.all(Y,p,nb.periods.IRF,n)
  # Store results
  IRFs <- cholesky.res$IRFs
  est.VAR <- cholesky.res$est.VAR
  Phi     <- Acoef(est.VAR)
  B.hat   <- cholesky.res$B.hat
  cst     <- Bcoef(est.VAR)[,p*n+1]
  resids  <- residuals(est.VAR)
  Omega   <- var(resids)
  # if no bootstrap then simply use point estimate
  simulated.IRFs  <- replicate(1,IRFs, simplify="array")
  simulated.B.hat <- replicate(1,B.hat, simplify="array")
  simulated.Phi   <- replicate(1,Phi, simplify="array")
  # if bootstrap then generate and store simulated IRFs, B.hat and Phi
  if(bootstrap.replications>1){
    bootstrap.res <- param.bootstrap(y,p,nb.periods.IRF,n,bootstrap.replications,
                                     posit.of.shock = 0)
    simulated.IRFs  <- bootstrap.res$simulated.IRFs
    simulated.B.hat <- bootstrap.res$simulated.B.hat
    simulated.Phi   <- bootstrap.res$simulated.Phi}
  ####

  # Initialize Q as identity matrix
  Q <- replicate(bootstrap.replications,diag(n), simplify="array")
  # This is where we will store the IRFs
  IRFs.final <- array(NaN,c(n,n,nb.periods.IRF,bootstrap.replications))
  # This loop identifies the relevant Q for each bootstrap replication
  for (l in 1:bootstrap.replications){

    for(i in 1:nb.shocks){
      # Compute S
      WWW <- array(0,c(n-i+1,n-i+1))
      for (h in H1[i]:H2[i]){
        V99 <- simulated.IRFs[,,h,l]%*%Q[,i:n,l]
        # Notice that we use only columns 2 to n of Q:
        #  the first column selects the TFP surprise shock,
        #  which is the first shock in the Cholesky
        #  decomposition where TFP is ordered first.
        JJ <- (H2[i]+1-h)*t(V99)%*%E[[i]]%*%V99
        WWW <- WWW+JJ}
      r <- eigen(WWW)
      # Take the eigenvector with the highest eigenvalue
      eigvec <- matrix(r$vectors[,1],n-i+1,1)
      # We might need to adjust the sign
      if (simulated.IRFs[variable[i],,norm[i],l]%*%Q[,i:n,l]%*%eigvec>0){
        Q[,i,l] <- Q[,i:n,l]%*%eigvec
      }else{
        Q[,i,l] <- - Q[,i:n,l]%*%eigvec}
      Q[,(i+1):n,l] <- Null(Q[,1:i,l]) # we ensure that columns 3 to n
      #   are orthogonal to the first 2
    }
    for (t in 1:nb.periods.IRF){
      IRFs.final[,,t,l] <- simulated.IRFs[,,t,l]%*%Q[,,l]
    }
  }

  # compute some key moments of the simulated IRFs
  stdv.IRFs <- apply(IRFs.final,c(1,2,3),sd)
  CI.lower.bounds <-
    apply(IRFs.final,c(1,2,3),function(x){quantile(x,(1-confidence.interval)/2)})
  CI.upper.bounds <-
    apply(IRFs.final,c(1,2,3),function(x){quantile(x,1-(1-confidence.interval)/2)})
  CI.median <-
    apply(IRFs.final,c(1,2,3),function(x){quantile(x,0.5)})

  if(indic.plot==1){
    # Plot graphs
    for (j in 1:nb.shocks){
      par(mfrow=c(2,ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
      for(i in 1:n){
        plot(CI.median[i,j,],type="l",lwd=2,xlab="",ylab="",
             ylim=c(min(CI.lower.bounds[i,j,]),
                    max(CI.upper.bounds[i,j,])),
             main=paste("Effect of ",names.of.shocks[j]," on ",
                        names.of.variables[i],sep=""))
        abline(h=0,col="grey")
        lines(CI.lower.bounds[i,j,],col="red",lty=2,lwd=2)
        lines(CI.upper.bounds[i,j,],col="red",lty=2,lwd=2)}
    }
  }

  return(
    list(
      #IRFs = IRFs,
      simulated.IRFs = IRFs.final,
      stdv.IRFs = stdv.IRFs,
      CI.lower.bounds = CI.lower.bounds,
      CI.upper.bounds = CI.upper.bounds
    ))
}
