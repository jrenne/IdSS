#################################3

# Various procedures useful for #
# inference in SVARs            #

#################################3

#############################################################################################3
# SVAR identified with Cholesky decomposition (extends svar.ordering)
# with various inference methods



#' Identification of SVAR - Estimation of IRF
#'
#' This functions identifies a SVAR model and computes IRF. The estimation of structural shocks is done by Cholesky decomposition. Confidence intervals are obtained by boostrapping the estimated VAR model.
#'Optionally, the function performs inference on the IRFs using parametric bootstrap, non-parametric
#' bootstrap, Monte Carlo, or bootstrap-after-bootstrap procedures, and plots
#' confidence intervals if requested.
#' @param Y 	Numeric matrix of size \code{T×n} containing the endogenous variables (one column per variable).
#' @param p Integer for the lag order.
#' @param posit.of.shock An integer giving the position (column index in \code{Y})
#' of the structural shock of interest. Default is 1.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param inference Integer. Method used for inference:
#'   \itemize{
#'     \item 0 = no inference (default)
#'     \item 1 = parametric bootstrap
#'     \item 2 = non-parametric bootstrap
#'     \item 3 = Monte Carlo
#'     \item 4 = bootstrap-after-bootstrap
#'   }
#' @param confidence.interval A numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{IRFs}}{Matrix of estimated impulse response functions.}
#'   \item{\code{simulated.IRFs}}{Array of simulated IRFs across bootstrap or
#'     Monte Carlo replications (if \code{inference > 0}).}
#'   \item{\code{all.stdv.IRFs}}{Matrix of standard deviations of the simulated IRFs
#'     (if \code{inference > 0}).}
#'   \item{\code{all.CI.lower.bounds}}{Matrix of lower bounds of the confidence intervals.}
#'   \item{\code{all.CI.upper.bounds}}{Matrix of upper bounds of the confidence intervals.}
#'   \item{\code{all.CI.median}}{Matrix of medians of the simulated IRFs.}
#' }
#' @examples
#' library(vars)
#' library(IdSS)
#' data("USmonthlyExample")
#' y <- as.matrix(USmonthlyExample[considered.variables])
#' ## Monte-Carlo
#' res.svar.ordering <-
#'  svar.ordering(y,p=3, posit.of.shock = 5,
#'                  nb.periods.IRF = 20,
#'                  inference = 3,
#'                  nb.draws = 200,
#'                  confidence.interval = 0.90,
#'                  indic.plot = 1)
#' @export
svar.ordering <- function(Y,p,
                          posit.of.shock = 1,
                          nb.periods.IRF = 20,
                          inference = 0, # 0 -> no inference, 1 -> parametric bootstrap, 2 <- non-parametric bootstrap, 3 <- monte carlo
                          nb.draws = 100, # This is used if inference >0
                          confidence.interval = 0.90, # expressed in pp.
                          indic.plot = 1 # Plots are displayed if = 1.
){

  names.of.variables <- colnames(Y)

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  baseline.res <- svar.ordering.aux(Y,p,
                                    posit.of.shock,
                                    nb.periods.IRF)
  IRFs <- baseline.res$IRFs

  if (inference >0){
    # If inference, simulate a set of IRFs using one of the following methods:

    # Parametric bootstrap
    if (inference == 1){
    simulated.res <- param.bootstrap(Y,p,
                                    nb.periods.IRF,
                                    n,
                                    nb.draws,
                                    posit.of.shock)
    # store the set of simulated IRFs
    all.simulated.IRFs.res <- simulated.res$simulated.IRFs
    }
    # Parametric bootstrap
    if (inference == 2){
      simulated.res <- nonparam.bootstrap(Y,p,
                                                nb.periods.IRF,
                                                n,
                                                nb.draws,
                                                posit.of.shock)
      # store the set of simulated IRFs
      all.simulated.IRFs.res <- simulated.res$simulated.IRFs
    }
    # Monte Carlo simulations
    if (inference == 3){
      simulated.res <- montecarlo(Y,p,
                                   nb.periods.IRF,
                                   n,
                                   nb.draws,
                                   posit.of.shock)
      # store the set of simulated IRFs
      all.simulated.IRFs.res <- simulated.res$simulated.IRFs
    }
    if (inference == 4){
      simulated.res <- bootstrap.after.bootstrap(Y,p,
                                                 nb.periods.IRF,
                                                 n,
                                                 nb.draws,
                                                 posit.of.shock)
      # store the set of simulated IRFs
      all.simulated.IRFs.res <- simulated.res$simulated.IRFs
    }

    # compute some key moments of the simulated IRFs
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    all.CI.median <- NULL
    for(i in 1:n){
      all.IRFs.i <- matrix(all.simulated.IRFs.res[,i,],ncol=nb.draws)

      stdv.IRF <- apply(all.IRFs.i,1,sd)
      all.stdv.IRFs <- cbind(all.stdv.IRFs,stdv.IRF)

      all.CI.lower.bounds <- cbind(all.CI.lower.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,(1-confidence.interval)/2)}))
      all.CI.upper.bounds <- cbind(all.CI.upper.bounds,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,1-(1-confidence.interval)/2)}))
      all.CI.median <- cbind(all.CI.median,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,0.5)}))
    }

  }else{
  #if no inference, simply use the estimated IRFs

    all.simulated.IRFs.res <- NULL
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- IRFs
    all.CI.upper.bounds <- IRFs
    all.CI.median <- IRFs
  }

  # here we plot the IRFs (if asked)
  if(indic.plot==1){
    par(mfrow=c(2,ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
    for(i in 1:n){
      plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
           ylim=c(min(all.CI.lower.bounds[,i]),
                  max(all.CI.upper.bounds[,i])),
           main=paste("Effect of shock on ",names.of.variables[i],sep=""))
      abline(h=0,col="grey")
      if(inference>0){
        lines(all.CI.lower.bounds[,i],col="red",lty=2,lwd=2)
        lines(all.CI.upper.bounds[,i],col="red",lty=2,lwd=2)
      }
    }
  }

  return(
    list(
      IRFs = IRFs,
      simulated.IRFs = all.simulated.IRFs.res,
      all.stdv.IRFs = all.stdv.IRFs,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds,
      all.CI.median = all.CI.median
    ))
}

#############################################################################################3
# Generates the matrix VAR regressors of size (T-p)*(1+n*p)

#' @export
VAR.regressors <- function(Y,p){
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  X <- array(1,c(T-p,1))
  for (i in 1:p){
    X <- cbind(X,Y[(1+p-i):(T-i),])
  }
  return(X=X)
}

#############################################################################################3
# Transforms the Phis and constant terms into the vector phi

#' @export
VEC.coef <- function(Phi,cst,n,p){
  phi <- cst
  for (i in 1:p){
    phi <- cbind(phi,Phi[[i]])
  }
  phi <- t(phi)
  phi = matrix(phi,ncol=1)
  return(phi = phi)
}

#############################################################################################3
# Computes the (vector) mean and variance-covariance of Phi

#' @export
distrib.phi <- function(Y,p,Phi,cst,Omega,n){
  phi <- VEC.coef(Phi,cst,n,p)
  X <- VAR.regressors(Y,p)
  XX <- t(X)%*%X
  #Omega1 <- diag(diag(Omega))
  var.phi <- Omega%x%solve(XX)
  return(list(var.phi = var.phi,
              phi = phi))
}

#############################################################################################3
# Computes the (vector) mean and variance-covariance of Omega

#' @export
distrib.omega <- function(Omega,n,T){
  omega <- NULL
  for (i in 1:n){
    omega <- rbind(omega,matrix(Omega[i:n,i],ncol=1))
  }
  var.omega <- array(NaN,c(dim(omega)[1],dim(omega)[1]))
  ind.r <- 1 #keep track of row index
  ind.c <- 1 #keep track of column index
  for (i in 1:n){
    for (j in i:n){
      for (k in 1:n){
        for (l in k:n){
          # correlation between sigma_ij and sigma_kl
          var.omega[ind.r,ind.c] <- (Omega[i,k]*Omega[j,l]+Omega[i,l]*Omega[j,k])/T
          ind.c <- ind.c+1
        }
      }
      ind.r <- ind.r+1
      ind.c <- 1
    }
  }
  return(list(omega = omega,
              var.omega = var.omega))
}

#############################################################################################3
# Monte Carlo simulations

#' @export
montecarlo <- function(Y,p,
                       nb.periods.IRF,
                       n,
                       nb.draws,
                       posit.of.shock =0 # =0 to simulate for all shocks
                       ){
  # Estimate IRFs for all shocks or just one
  T <- dim(Y)[1]
  if (posit.of.shock !=0){
    baseline.res <- svar.ordering.aux(Y,p,
                                      posit.of.shock,
                                      nb.periods.IRF)
  }else{
    baseline.res <- svar.ordering.all(Y,p,
                                      nb.periods.IRF,
                                      n)
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[,p*n+1]
  resids <- residuals(est.VAR)
  Omega <- var(resids)

  simulated.IRFs  <- replicate(nb.draws, IRFs, simplify="array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify="array")
  simulated.Phi   <- replicate(nb.draws, Phi)

  # Get the asymptotic distribution of estimated coefficients
  distrib.phi.res <- distrib.phi(Y,p,Phi,cst,Omega,n)
  phi             <- distrib.phi.res$phi
  var.phi         <- distrib.phi.res$var.phi
  Phi.simul       <- list()
  distrib.omega.res <- distrib.omega(Omega,n,T)
  omega             <- distrib.omega.res$omega
  var.omega         <- distrib.omega.res$var.omega
  Omega.simul       <- array(0,c(n,n))

  i<-0
  xx<-0
  # Perform the simulations
  while(xx <nb.draws){
    i<-i+1
    # here draw B.hat
    omega.simul       <- mvrnorm(1, omega,var.omega)
    j <- 1
    for (k in 1:n){
      Omega.simul[k:n,k] <- omega.simul[j:(j+n-k)]
      Omega.simul[k,k:n] <- t(omega.simul[j:(j+n-k)])
      j <- j+n-k+1
    }
    r <- eigen(Omega.simul) # use this if some Omega.simul are not positive definite
    eig <- r$values
    if (length(eig[eig>0])!=n){
      break
    }else{
      xx<-xx+1
      # here draw Phi
      phi.simul       <- mvrnorm(1, phi, var.phi)
      B.hat.simul <- t(chol(Omega.simul))
      phi.simul.temp <- t(matrix(phi.simul,ncol=n,nrow=1+n*p))
      for (j in 1:p){ Phi.simul[[j]] <- phi.simul.temp[,(2+(j-1)*n):(1+j*n)]}
      # Simulate IRFs for all shocks or just one
      if (posit.of.shock!=0){
        IRFs.simul <- simul.VAR(c=rep(0,n),Phi.simul,
                                    B.hat.simul,
                                    nb.periods.IRF,
                                    y0.star=rep(0,n*p),
                                    indic.IRF = 1,
                                    diag(n)[,posit.of.shock])

        simulated.IRFs[,,xx] <- IRFs.simul
      }else{
      IRFs.simul <- simul.VAR.all(c=rep(0,n),n,Phi.simul,
                                  B.hat.simul,
                                  nb.periods.IRF,
                                  y0.star=rep(0,n*p),
                                  indic.IRF = 1,
                                  diag(n))

      simulated.IRFs[,,,xx] <- IRFs.simul
      }
      simulated.B.hat[,,xx] <- B.hat.simul
      simulated.Phi[,xx] <- Phi.simul
      #}
    }
  }
  return(list(simulated.IRFs = simulated.IRFs,
              simulated.B.hat = simulated.B.hat,
              simulated.Phi = simulated.Phi))
}

#############################################################################################3
# Gaussian parametric bootstrap simulations
# (see Stock and Watson, appendix A.2)

#' @export
param.bootstrap <- function(Y,p,
                            nb.periods.IRF,
                            n,
                            nb.draws,
                            posit.of.shock = 0 # =0 to simulate for all shocks
                            ){


  # Estimate IRFs for all shocks or just one
  if (posit.of.shock !=0){
    baseline.res <- svar.ordering.aux(Y,p,
                                      posit.of.shock,
                                      nb.periods.IRF)
  }else{
    baseline.res <- svar.ordering.all(Y,p,
                                      nb.periods.IRF,
                                      n)
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[,p*n+1]

  simulated.IRFs  <- replicate(nb.draws, IRFs, simplify="array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify="array")
  simulated.Phi   <- replicate(nb.draws, Phi)

  y0.star <- NULL
  for(k in p:1){
    y0.star <- c(y0.star,Y[k,])
  }

  # Perform the simulations
  for(i in 1:nb.draws){
    # Simulate T observations
    simulated.Y <- simul.VAR(cst,Phi,B.hat,nb.sim=dim(Y)[1],y0.star)
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock!=0){
    simulated.res <- svar.ordering.aux(simulated.Y,p,
                                       posit.of.shock,
                                       nb.periods.IRF)
    simulated.IRFs[,,i] <- simulated.res$IRFs
    }else{
    simulated.res <- svar.ordering.all(simulated.Y,p,
                                       nb.periods.IRF,
                                       n)
    simulated.IRFs[,,,i] <- simulated.res$IRFs
    }

    simulated.B.hat[,,i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi[,i] <- Acoef(simulated.est.VAR)
  }


  return(list(simulated.IRFs = simulated.IRFs,
              simulated.B.hat = simulated.B.hat,
              simulated.Phi = simulated.Phi))
}

#############################################################################################3
# Non-parametric bootstrap simulations

#' @export
nonparam.bootstrap <- function(Y,p,
                            nb.periods.IRF,
                            n,
                            nb.draws,
                            posit.of.shock
                            ){
  # Estimate IRFs for all shocks or just one
  if (posit.of.shock !=0){
    baseline.res <- svar.ordering.aux(Y,p,
                                      posit.of.shock,
                                      nb.periods.IRF)
  }else{
    baseline.res <- svar.ordering.all(Y,p,
                                      nb.periods.IRF,
                                      n)
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[,p*n+1]
  resids <- residuals(est.VAR)

  simulated.IRFs  <- replicate(nb.draws, IRFs, simplify="array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify="array")
  simulated.Phi   <- replicate(nb.draws, Phi)

  y0.star <- NULL
  for(k in p:1){
    y0.star <- c(y0.star,Y[k,])
  }

  # Perform the simulations
  for(i in 1:nb.draws){
    # draw set residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1],dim(Y)[1],replace=TRUE)
    simul.residuals <- resids[simul.index,]
    # Simulate T observations
    simulated.Y <- simul.VAR(cst,Phi,B = diag(n),nb.sim=dim(Y)[1],y0.star,indic.IRF = 0,eta = simul.residuals)
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock!=0){
      simulated.res <- svar.ordering.aux(simulated.Y,p,
                                         posit.of.shock,
                                         nb.periods.IRF)
      simulated.IRFs[,,i] <- simulated.res$IRFs
    }else{
      simulated.res <- svar.ordering.all(simulated.Y,p,
                                         nb.periods.IRF,
                                         n)
      simulated.IRFs[,,,i] <- simulated.res$IRFs
    }

    simulated.B.hat[,,i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi[,i] <- Acoef(simulated.est.VAR)
  }


  return(list(simulated.IRFs = simulated.IRFs,
              simulated.B.hat = simulated.B.hat,
              simulated.Phi = simulated.Phi))
}

#############################################################################################3
# Non-parametric bootstrap after bootstrap simulations

#' @export
bootstrap.after.bootstrap <- function(Y,p,
                               nb.periods.IRF,
                               n,
                               nb.draws,
                               posit.of.shock
){
  # Estimate IRFs for all shocks or just one
  if (posit.of.shock !=0){
    baseline.res <- svar.ordering.aux(Y,p,
                                      posit.of.shock,
                                      nb.periods.IRF)
  }else{
    baseline.res <- svar.ordering.all(Y,p,
                                      nb.periods.IRF,
                                      n)
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[,p*n+1]
  resids <- residuals(est.VAR)
  phi <- VEC.coef(Phi,c(n,0),n,p)


  simulated.IRFs  <- replicate(nb.draws, IRFs, simplify="array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify="array")
  simulated.Phi   <- replicate(nb.draws, Phi)
  simulated.phi   <- array(replicate(nb.draws, phi, simplify="array"),c(dim(phi)[1],nb.draws))

  y0.star <- NULL
  for(k in p:1){
    y0.star <- c(y0.star,Y[k,])
  }

  # FIRST BOOTSTRAP
  # Perform the simulations
  for(i in 1:nb.draws){
    # draw set of residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1],dim(Y)[1],replace=TRUE)
    simul.residuals <- resids[simul.index,]
    # Simulate T observations
    simulated.Y <- simul.VAR(cst,Phi,B = diag(n),nb.sim=dim(Y)[1],y0.star,indic.IRF = 0,eta = simul.residuals)
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock!=0){
      simulated.res <- svar.ordering.aux(simulated.Y,p,
                                         posit.of.shock,
                                         nb.periods.IRF)
      simulated.IRFs[,,i] <- simulated.res$IRFs
    }else{
      simulated.res <- svar.ordering.all(simulated.Y,p,
                                         nb.periods.IRF,
                                         n)
      simulated.IRFs[,,,i] <- simulated.res$IRFs
    }

    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi.temp <- Acoef(simulated.est.VAR)
    simulated.phi[,i] <- VEC.coef(simulated.Phi.temp,c(n,1),n,p)
  }

  # Compute the bias
  median.simulated.phi <- apply(simulated.phi,1,function(x){quantile(x,0.5)})
  bias <- median.simulated.phi-phi
  # Compute the unbiased Phi
  delta <-1
  for (l in 1:100){
    unbiased.phi <- phi-bias
    unbiased.Phi.temp <- t(matrix(unbiased.phi,ncol=n,nrow=1+n*p))
    unbiased.Phi <- list()
    for (j in 1:p){ unbiased.Phi[[j]] <- unbiased.Phi.temp[,(2+(j-1)*n):(1+j*n)]}
    unbiased.PHI <- make.PHI(unbiased.Phi)
    if (abs(eigen(unbiased.PHI)$values[1])<1){break}
    else{
      delta <- delta-0.01
      bias <- bias*delta
    }
  }
  # SECOND BOOTSTRAP
  # Perform the simulations
  for(i in 1:nb.draws){
    # draw set of residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1],dim(Y)[1],replace=TRUE)
    simul.residuals <- resids[simul.index,]
    # Simulate T observations
    simulated.Y <- simul.VAR(cst,unbiased.Phi,B = diag(n),nb.sim=dim(Y)[1],y0.star,indic.IRF = 0,eta = simul.residuals)
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock!=0){
      simulated.res <- svar.ordering.aux(simulated.Y,p,
                                         posit.of.shock,
                                         nb.periods.IRF)
    }else{
      simulated.res <- svar.ordering.all(simulated.Y,p,
                                         nb.periods.IRF,
                                         n)
    }

    simulated.B.hat[,,i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi.temp <- Acoef(simulated.est.VAR)
    simulated.phi <- VEC.coef(simulated.Phi.temp,c(n,1),n,p)
    # Compute the unbiased Phi
      unbiased.simulated.phi <- simulated.phi-bias
      unbiased.simulated.Phi.temp <- t(matrix(unbiased.simulated.phi,ncol=n,nrow=1+n*p))
      unbiased.simulated.Phi <- list()
      for (j in 1:p){ unbiased.simulated.Phi[[j]] <- unbiased.simulated.Phi.temp[,(2+(j-1)*n):(1+j*n)]}
      simulated.Phi[,i] <- unbiased.simulated.Phi

    if (posit.of.shock!=0){
      IRFs.simul <- simul.VAR(c=rep(0,n),simulated.Phi[,i],
                              simulated.B.hat[,,i],
                              nb.periods.IRF,
                              y0.star=rep(0,n*p),
                              indic.IRF = 1,
                              diag(n)[,posit.of.shock])

      simulated.IRFs[,,i] <- IRFs.simul
    }else{
      IRFs.simul <- simul.VAR.all(c=rep(0,n),n,simulated.Phi[,i],
                                  simulated.B.hat[,,i],
                                  nb.periods.IRF,
                                  y0.star=rep(0,n*p),
                                  indic.IRF = 1,
                                  diag(n))

      simulated.IRFs[,,,i] <- IRFs.simul
    }
    }

  return(list(simulated.IRFs = simulated.IRFs,
              simulated.B.hat = simulated.B.hat,
              simulated.Phi = simulated.Phi))
}

#############################################################################################3
# SVAR with Cholesky decomposition generating all IRFs (alternative to svar.ordering.aux)

#' @export
svar.ordering.all <- function(Y,p,
                              nb.periods.IRF,
                              n
){
  IRFs <- array(NaN,c(n,n,nb.periods.IRF))
  # Call svar.ordering.aux for all shocks
  for (j in 1:n){
    baseline.res <- svar.ordering.aux(Y,p,
                                      j,
                                      nb.periods.IRF)
    IRFs.j    <- baseline.res$IRFs
    dim(IRFs.j) <- c(nb.periods.IRF,n,1)
    IRFs[,j,] <- aperm(IRFs.j, c(2,3,1))
  }

  est.VAR <- baseline.res$est.VAR
  Phi     <- Acoef(est.VAR)
  B.hat   <- baseline.res$B.hat
  cst     <- Bcoef(est.VAR)[,p*n+1]
  resids  <- residuals(est.VAR)
  Omega   <- var(resids)

  return(list(
    IRFs = IRFs,
    est.VAR = est.VAR,
    Phi = Phi, # IRFs.draw*Q
    B.hat = B.hat,
    cst = cst,
    resids = resids,
    Omega = Omega # Q
  ))
}

#############################################################################################3
# Simulate IRFs for all shocks (one set of IRFs by shock)

#' @export
simul.VAR.all <- function(c=rep(0,n),n,Phi,
                          B.hat,
                          nb.periods.IRF,
                          y0.star,
                          indic.IRF = 1,
                          U.shocks = diag(n)
){
  IRFs <- array(NaN,c(n,n,nb.periods.IRF))
  # Call simul.VAR for all shocks
  for (j in 1:n){
    IRFs.j <- simul.VAR(c=rep(0,n),
                        Phi,
                        B.hat,
                        nb.periods.IRF,
                        y0.star,
                        indic.IRF = 1,
                        U.shocks[,j])
    IRFs.j <- t(IRFs.j)
    dim(IRFs.j) <- c(n,1,nb.periods.IRF)
    IRFs[,j,] <- IRFs.j
  }
  return(IRFs = IRFs)
}


