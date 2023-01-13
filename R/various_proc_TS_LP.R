# =============================================================
# MacroEconometrics course
# Various procedures
# Jean-Paul Renne, 2018
# =============================================================

autocov <- function(X,n){
  if(class(X)[1]!="matrix"){
    T <- length(X)
    X.1 <- X[1:(T-n)] - mean(X)
    X.2 <- X[(n+1):T] - mean(X)
    return(1/T * sum(X.1 * X.2))
  }else{
    T <- dim(X)[1]
    k <- dim(X)[2]
    vec.1 <- matrix(1,1,k)
    mean.X <- apply(X,2,mean)
    X.1 <- X[1:(T-n),] - t(matrix(mean.X,k,T-n))
    X.2 <- X[(n+1):T,] - t(matrix(mean.X,k,T-n))
    return(
      matrix(1/T * apply((X.1 %x% vec.1) * (vec.1 %x% X.2),2,sum),k,k)
    )
  }
}

NW.Weights <- function(q){
  if(q==0){
    weights <- 0
  }else{
    weights <- 1 - (1:q)/(q+1)
  }
  return(weights)
}


NW.LongRunVariance <- function(X,q){
  gamma0 <- autocov(X,0)
  LRV <- gamma0
  if(q>0){
    weights <- NW.Weights(q)
    for(i in 1:q){
      LRV <- LRV + weights[i] * (autocov(X,i) + t(autocov(X,i)))
    }
  }
  return(LRV)
}

tsls <- function(Y,X,Z,q){
  # Regress y on x, using z as an instrument.
  # Y, X and Z have to be matrices
  # The covariance matrix of parameter estimates is corrected for HAC (Newey-West)

  if(class(Y)[1]!="matrix"){
    print("Error: Y should be a matrix.")
    return(0)
  }
  if(class(X)[1]!="matrix"){
    print("Error: X should be a matrix.")
    return(0)
  }
  if(class(Z)[1]!="matrix"){
    print("Error: Z should be a matrix.")
    return(0)
  }

  # Determine sample size:
  T <- length(Y)

  # Compute the IV estimator:
  PZ <- solve(t(Z)%*%Z) %*% t(Z)
  Y.hat <- PZ %*% Y
  X.hat <- PZ %*% X
  PX.hat <- solve(t(X.hat)%*%X.hat) %*% t(X.hat)
  b.iv <- PX.hat %*% Y.hat

  # Residual estimates:
  eps <- Y - X %*% b.iv

  # Approximation of the covariance matrix of b.iv:
  Q <- T * solve(t(X) %*% Z %*% solve(t(Z)%*%Z) %*% t(Z) %*% X) %*% t(X.hat)

  Z.times.eps <- Z * matrix(eps,dim(Z)[1],dim(Z)[2])
  S <- NW.LongRunVariance(Z.times.eps,q)

  covmat.b.iv <- 1/T * Q %*% S %*% t(Q)

  return(list(
    b.iv=b.iv,
    covmat.b.iv = covmat.b.iv
  )
  )
}

svar.iv <- function(Y,Z,p,
                    names.of.variables,
                    nb.periods.IRF = 20,
                    z.AR.order=3, # This is used in the parametric bootstrap only
                    nb.bootstrap.replications = 0, # This is used in the parametric bootstrap only
                    confidence.interval = 0.90, # expressed in pp.
                    indic.plot = 1 # Plots are displayed if = 1.
){
  # This function implements the IV-SVAR approach.
  # We are interested in the IRF of Y following a shock on a structural shock (eta_1), Z being an instrument for eta_1
  # Y is modelled as a VAR(p) model.
  # The IRF are scaled in such a way that the contemporaneous impact of eta_1 on Y_1 is 1.

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  baseline.res <- svar.iv.aux(Y,Z,p,
                              names.of.variables,
                              nb.periods.IRF)
  IRFs <- baseline.res$IRFs


  if(nb.bootstrap.replications>0){
    all.simulated.IRFs.res <- array(NaN,c(nb.periods.IRF,n,nb.bootstrap.replications))

    # Computation of std dev of IRF by parametric Gaussian bootstrap
    # (see Stock and Watson, appendix A.2)

    # step 0: fit an AR process for z:
    est.AR.z <- arima(c(Z),order = c(z.AR.order,0,0))

    # step 1: parameterize a VAR for [y',z]':
    est.VAR <- baseline.res$est.VAR
    Phi <- Acoef(est.VAR)
    resids <- residuals(est.VAR)

    c.IV <- Bcoef(est.VAR)[,p*n+1]
    c.IV <- c(c.IV,0) # add one term for the instrument
    Phi.IV <- list()
    for(k in 1:max(p,z.AR.order)){
      aux <- matrix(0,n+1,n+1)
      if(k<=p){
        aux[1:n,1:n] <- Phi[[k]]
      }
      if(k <= z.AR.order){
        aux[n+1,n+1] <- est.AR.z$coef[k]
      }
      Phi.IV[[k]] <- aux
    }
    RESIDS.eps.plus.z <- cbind(resids,est.AR.z$residuals[(1+p):T])
    Sigma.IV <- var(RESIDS.eps.plus.z)
    cor.Sigma.IV <- cor(RESIDS.eps.plus.z)
    B.IV <- t(chol(Sigma.IV))
    y0.star <- NULL
    for(k in max(p,z.AR.order):1){
      y0.star <- c(y0.star,Y[k,],0)
    }

    all.B.tilde.1 <- NULL
    for(i in 1:nb.bootstrap.replications){
      simulated.y.plus.z <- simul.VAR(c.IV,Phi.IV,B.IV,nb.sim=dim(Y)[1],y0.star)
      simulated.Y <- simulated.y.plus.z[,1:n]
      simulated.Z <- matrix(simulated.y.plus.z[,n+1],ncol=1)

      simulated.res <- svar.iv.aux(simulated.Y,simulated.Z,p,
                                   names.of.variables,
                                   nb.periods.IRF)

      all.simulated.IRFs.res[,,i] <- simulated.res$IRFs
      all.B.tilde.1 <- cbind(all.B.tilde.1,simulated.res$B.tilde.1)
    }

    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    for(i in 1:n){
      all.IRFs.i <- matrix(all.simulated.IRFs.res[,i,],ncol=nb.bootstrap.replications)

      stdv.IRF <- apply(all.IRFs.i,1,sd)
      all.stdv.IRFs <- cbind(all.stdv.IRFs,stdv.IRF)

      all.CI.lower.bounds <- cbind(all.CI.lower.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,(1-confidence.interval)/2)}))
      all.CI.upper.bounds <- cbind(all.CI.upper.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,1-(1-confidence.interval)/2)}))
    }

  }else{
    all.simulated.IRFs.res <- NULL
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    all.B.tilde.1 <- NULL
  }

  if(indic.plot==1){
    par(mfrow=c(2,
                ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
    for(i in 1:n){
      plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
           ylim=c(min(all.CI.lower.bounds[,i]),
                  max(all.CI.upper.bounds[,i])),
           main=paste("Effect of shock on ",considered.variables[i],sep=""))
      abline(h=0,col="grey")
      lines(all.CI.lower.bounds[,i],col="red",lty=2,lwd=2)
      lines(all.CI.upper.bounds[,i],col="red",lty=2,lwd=2)
    }
  }

  return(
    list(
      IRFs = IRFs,
      all.simulated.IRFs.res = all.simulated.IRFs.res,
      all.stdv.IRFs = all.stdv.IRFs,
      Sigma.IV = Sigma.IV,
      cor.Sigma.IV = cor.Sigma.IV,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds,
      all.B.tilde.1 = all.B.tilde.1
    ))

}




svar.iv.aux <- function(Y,Z,p,
                        names.of.variables,
                        nb.periods.IRF){
  # This function implements the IV-SVAR approach.
  # We are interested in the IRF of Y following a shock on a structural shock (eta_1), Z being an instrument for eta_1
  # Y is modelled as a VAR(p) model.
  # The IRF are scaled in such a way that the contemporaneous impact of eta_1 on Y_1 is 1.

  if(class(Y)[1]!="matrix"){
    print("Error: Y should be a matrix.")
    return(0)
  }
  if(class(Z)[1]!="matrix"){
    print("Error: Z should be a matrix.")
    return(0)
  }

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  # Estimate a VAR(p) model for Y:
  colnames(Y) <- names.of.variables
  est.VAR <- VAR(Y,p=p)

  # get estimated Phi matrices:
  Phi <- Acoef(est.VAR)

  # Compute the covariance matrix of VAR residuals:
  resids <- residuals(est.VAR)
  #Omega <- var(resids)

  # Run TSLS estimations (eps[i] on eps[1], using Z as an instrument):
  B.tilde.1 <- NULL
  #stdv.B.tilde.1 <- NULL
  for(i in 1:n){
    YY <- matrix(resids[,i],ncol=1)
    XX <- matrix(resids[,1],ncol=1)
    ZZ <- matrix(Z[(1+p):length(Z)],ncol=1)
    eq.iv <- tsls(YY,XX,ZZ,q=3)
    B.tilde.1[i] <- eq.iv$b.iv
    #stdv.B.tilde.1[i] <- sqrt(eq.iv$covmat.b.iv)
  }

  # Compute IRFs
  B.tilde <- matrix(0,n,n)
  B.tilde[,1] <- B.tilde.1

  y0.star <- rep(0,dim(y)[2]*p)

  nb.sim <- nb.periods.IRF
  IRFs <- simul.VAR(c=rep(0,dim(y)[2]),
                    Phi,
                    B.tilde,
                    nb.sim,
                    y0.star,
                    indic.IRF = 1,
                    u.shock = c(1,rep(0,n-1)))

  return(list(
    B.tilde.1 = B.tilde.1,
    IRFs = IRFs,
    est.VAR = est.VAR
  ))
}


svar.ordering <- function(Y,p,
                          posit.of.shock = 1,
                          nb.periods.IRF = 20,
                          nb.bootstrap.replications = 0, # This is used in the parametric bootstrap only
                          confidence.interval = 0.90, # expressed in pp.
                          indic.plot = 1 # Plots are displayed if = 1.
){
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # Christiano, Eichembaum and Evans methodology

  names.of.variables <- colnames(Y)

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  baseline.res <- svar.ordering.aux(Y,p,
                                    posit.of.shock,
                                    nb.periods.IRF)
  IRFs <- baseline.res$IRFs

  if(nb.bootstrap.replications>0){
    all.simulated.IRFs.res <- array(NaN,c(nb.periods.IRF,n,nb.bootstrap.replications))

    # Computation of std dev of IRF by parametric Gaussian bootstrap
    # (see Stock and Watson, appendix A.2)

    # step 1: parameterize a VAR for [y',z]':
    est.VAR <- baseline.res$est.VAR
    Phi <- Acoef(est.VAR)
    B.hat <- baseline.res$B.hat
    cst <- Bcoef(est.VAR)[,p*n+1]

    y0.star <- NULL
    for(k in p:1){
      y0.star <- c(y0.star,Y[k,])
    }

    all.B.tilde.1 <- NULL
    for(i in 1:nb.bootstrap.replications){
      simulated.Y <- simul.VAR(cst,Phi,B.hat,nb.sim=dim(Y)[1],y0.star)
      colnames(simulated.Y) <- colnames(Y)

      simulated.res <- svar.ordering.aux(simulated.Y,p,
                                         posit.of.shock,
                                         nb.periods.IRF)

      all.simulated.IRFs.res[,,i] <- simulated.res$IRFs
    }

    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    for(i in 1:n){
      all.IRFs.i <- matrix(all.simulated.IRFs.res[,i,],ncol=nb.bootstrap.replications)

      stdv.IRF <- apply(all.IRFs.i,1,sd)
      all.stdv.IRFs <- cbind(all.stdv.IRFs,stdv.IRF)

      all.CI.lower.bounds <- cbind(all.CI.lower.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,(1-confidence.interval)/2)}))
      all.CI.upper.bounds <- cbind(all.CI.upper.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,1-(1-confidence.interval)/2)}))
    }

  }else{
    all.simulated.IRFs.res <- NULL
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
  }

  if(indic.plot==1){
    par(mfrow=c(2,
                ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
    for(i in 1:n){
      plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
           ylim=c(min(all.CI.lower.bounds[,i]),
                  max(all.CI.upper.bounds[,i])),
           main=paste("Effect of shock on ",names.of.variables[i],sep=""))
      abline(h=0,col="grey")
      if(nb.bootstrap.replications>0){
        lines(all.CI.lower.bounds[,i],col="red",lty=2,lwd=2)
        lines(all.CI.upper.bounds[,i],col="red",lty=2,lwd=2)
      }
    }
  }

  return(
    list(
      IRFs = IRFs,
      all.simulated.IRFs.res = all.simulated.IRFs.res,
      all.stdv.IRFs = all.stdv.IRFs,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds
    ))
}



svar.ordering.aux <- function(Y,p,
                              posit.of.shock,
                              nb.periods.IRF){
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # Christiano, Eichembaum and Evans methodology

  considered.variables <- colnames(Y)
  n <- dim(Y)[2]

  # Select number of lags in VAR models:
  #VARselect(y)

  est.VAR <- VAR(Y,p=p)

  # get estimated Phi matrices:
  Phi <- Acoef(est.VAR)

  # Compute the covariance matrix of VAR residuals:
  resids <- residuals(est.VAR)
  Omega <- var(resids)

  # Short-Run Restrictions ~ Cholesky:
  B.hat <- t(chol(Omega))

  y0.star <- rep(0,n*p)

  u.shock <- rep(0,n)
  u.shock[posit.of.shock] <- 1

  par(mfrow=c(2,4))

  IRFs <- simul.VAR(c=rep(0,dim(Y)[2]),
                    Phi,
                    B.hat,
                    nb.periods.IRF,
                    y0.star,
                    indic.IRF = 1,
                    u.shock = u.shock)

  return(list(
    IRFs = IRFs,
    est.VAR = est.VAR,
    B.hat = B.hat
  ))
}



make.jorda.irf <- function(Y,
                           posit.of.shock = 1,
                           nb.periods.IRF = 20,
                           nb.lags.endog.var.4.control=0,
                           indic.plot = 1, # Plots are displayed if = 1.
                           confidence.interval = 0.90){
  # This function implements Jorda's approach.

  considered.variables <- colnames(Y)
  T <- dim(Y)[1]
  n <- dim(Y)[2]

  mat.of.results <- matrix(NaN,nb.periods.IRF+1,n)
  all.stdv.IRFs <- matrix(NaN,nb.periods.IRF+1,n)
  for(h in 0:nb.periods.IRF){
    print(paste("Jorda's approach, Currently working on horizon h=",toString(h),sep=""))
    for(i in 1:n){
      y.i.t.p.h <- Y[(h+1+nb.lags.endog.var.4.control):T,i]
      X <- Y[(1+nb.lags.endog.var.4.control):(T-h),1:posit.of.shock]
      if(nb.lags.endog.var.4.control>0){
        for(k in 1:nb.lags.endog.var.4.control){
          X <- cbind(X,Y[(1+nb.lags.endog.var.4.control-k):(T-h-k),])
        }
      }
      X <- cbind(1,X)
      b <- solve(t(X)%*%X) %*% t(X) %*% y.i.t.p.h
      mat.of.results[h+1,i] <- b[posit.of.shock+1]

      # Approximation of the covariance matrix of b (Newey-West):
      Q <- T * solve(t(X) %*% X)
      e <- y.i.t.p.h - X %*% b
      X.times.eps <- X * matrix(e,dim(X)[1],dim(X)[2])
      S <- NW.LongRunVariance(X.times.eps,h+1)
      covmat.b <- 1/T * Q %*% S %*% t(Q)
      all.stdv.IRFs[h+1,i] <- sqrt(covmat.b[posit.of.shock+1,posit.of.shock+1])

      # Compute variance of the shock:
      if(i==posit.of.shock){
        variance.of.shock <- c(var(e))
      }
    }
  }
  IRFs <- mat.of.results / sqrt(variance.of.shock)
  all.stdv.IRFs <- all.stdv.IRFs / sqrt(variance.of.shock)

  if(indic.plot==1){
    par(mfrow=c(2,
                ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
    for(i in 1:n){
      lower.bound <- IRFs[,i] + qnorm((1-confidence.interval)/2) * all.stdv.IRFs[,i]
      upper.bound <- IRFs[,i] + qnorm(1-(1-confidence.interval)/2) * all.stdv.IRFs[,i]
      plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
           ylim=c(min(lower.bound),max(upper.bound)),
           main=paste("Effect of shock on ",considered.variables[i],sep=""))
      abline(h=0,col="grey")
      lines(lower.bound,col="red",lty=2,lwd=2)
      lines(upper.bound,col="red",lty=2,lwd=2)
    }
  }

  return(
    list(
      IRFs = IRFs,
      all.stdv.IRFs = all.stdv.IRFs
    )
  )
}


make.LPIV.irf <- function(Y,Z,
                          posit.of.shock = 1,
                          nb.periods.IRF = 20,
                          nb.lags.Y.4.control=0,
                          nb.lags.Z.4.control=0,
                          indic.plot = 1, # Plots are displayed if = 1.
                          confidence.interval = 0.90){
  # This function implements the LP-IV approach.
  # By convention, the shock has a unit effect on the first component of Y
  # Z is the instrumental variable.

  considered.variables <- colnames(Y)
  T <- dim(Y)[1]
  n <- dim(Y)[2]

  max.nb.lags.control <- max(nb.lags.Y.4.control,nb.lags.Z.4.control)

  IRFs <- matrix(NaN,nb.periods.IRF+1,n)
  all.stdv.IRFs <- matrix(NaN,nb.periods.IRF+1,n)

  for(h in 0:nb.periods.IRF){
    print(paste("LP-IV approach, Currently working on horizon h=",toString(h),sep=""))

    for(i in 1:n){
      ZZ <- matrix(Z[(1+max.nb.lags.control):(T-h),],ncol=1)
      X  <- matrix(Y[(1+max.nb.lags.control):(T-h),1],ncol=1)

      y.i.t.p.h <- matrix(Y[(h+1+max.nb.lags.control):T,i],ncol=1)
      if(nb.lags.Y.4.control>0){
        for(k in 1:nb.lags.Y.4.control){
          X  <- cbind(X, Y[(1+max.nb.lags.control-k):(T-h-k),])
          ZZ <- cbind(ZZ,Y[(1+max.nb.lags.control-k):(T-h-k),])
        }
      }
      if(nb.lags.Z.4.control>0){
        for(k in 1:nb.lags.Z.4.control){
          X  <- cbind(X, Z[(1+max.nb.lags.control-k):(T-h-k),])
          ZZ <- cbind(ZZ,Z[(1+max.nb.lags.control-k):(T-h-k),])
        }
      }
      X  <- cbind(1,X) # add a constant in the regression
      ZZ <- cbind(ZZ,1) # add a constant in the regression
      res.iv.est <- tsls(y.i.t.p.h,X,ZZ,h+1)
      IRFs[h+1,i] <- res.iv.est$b.iv[2] # Note: res.iv.est$b.iv[1] corresponds to the constant
      all.stdv.IRFs[h+1,i] <- sqrt(res.iv.est$covmat.b.iv[2,2])
    }
  }

  if(indic.plot==1){
    par(mfrow=c(2,
                ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
    for(i in 1:n){
      lower.bound <- IRFs[,i] + qnorm((1-confidence.interval)/2) * all.stdv.IRFs[,i]
      upper.bound <- IRFs[,i] + qnorm(1-(1-confidence.interval)/2) * all.stdv.IRFs[,i]
      plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
           ylim=c(min(lower.bound),max(upper.bound)),
           main=paste("Effect of shock on ",considered.variables[i],sep=""))
      abline(h=0,col="grey")
      lines(lower.bound,col="red",lty=2,lwd=2)
      lines(upper.bound,col="red",lty=2,lwd=2)
    }
  }

  return(
    list(
      IRFs = IRFs,
      all.stdv.IRFs = all.stdv.IRFs
    )
  )
}
