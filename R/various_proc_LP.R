# =============================================================
# MacroEconometrics course
# Various procedures
# Jean-Paul Renne, 2018
# =============================================================

# this function is used in NW.LongRunVariance()

#' @export
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

# Function used only in NW.LongRunVariance(), so we included directly into it.
# NW.Weights <- function(q){
#   if(q==0){
#     weights <- 0
#   }else{
#     weights <- 1 - (1:q)/(q+1)
#   }
#   return(weights)
# }


#' Newey–West Long-Run Variance Estimator
#'
#'Computes the Newey–West estimator of the long-run variance
#'(or long-run covariance matrix) of a univariate or
#'multivariate time series.
#'
#' @param X A vector or a matrix representing the time series.
#' If \code{X} is a matrix of size \eqn{T \times k}, the function returns a
#' \eqn{k \times k} long-run covariance matrix.
#'
#' @param q A non-negative integer specifying the Newey–West truncation lag.
#' If \code{q = 0}, the function returns the contemporaneous covariance
#' \code{autocov(X, 0)}.
#'
#' @return A scalar (if \code{X} is univariate) or a matrix (if \code{X} is
#' multivariate) giving the Newey–West long-run variance or covariance estimate.
#'
#' @examples
#' # Univariate example
#' set.seed(123)
#' x <- rnorm(200)
#' NW.LongRunVariance(x, q = 4)
#'
#' # Multivariate example
#' X <- cbind(rnorm(200), rnorm(200))
#' NW.LongRunVariance(X, q = 5)
#'
#' @export
NW.LongRunVariance <- function(X,q){
  gamma0 <- autocov(X,0)
  LRV <- gamma0
  if(q>0){
    weights <-  1 - (1:q)/(q+1)
    for(i in 1:q){
      LRV <- LRV + weights[i] * (autocov(X,i) + t(autocov(X,i)))
    }
  }
  return(LRV)
}

#' Two-Stage Least Squares (2SLS) estimation with HAC covariance
#'
#' @description
#' Computes the two-stage least squares (2SLS) estimator of a linear regression
#' model using instrumental variables, with a heteroskedasticity- and
#' autocorrelation-consistent (HAC) covariance matrix based on a Newey–West
#' long-run variance estimator.
#'
#' The model is:
#' \deqn{Y = X \beta + \varepsilon,}
#' where \eqn{X} may be endogenous and \eqn{Z} is a matrix of instruments.
#'
#' @param Y Numeric matrix of dimension \eqn{T \times 1} containing the dependent
#' variable.
#'
#' @param X Numeric matrix of dimension \eqn{T \times k} containing the regressors
#' (potentially endogenous).
#'
#' @param Z Numeric matrix of dimension \eqn{T \times m} containing the
#' instrumental variables.
#'
#' @param q Non-negative integer specifying the lag truncation parameter used in
#' the Newey–West long-run variance estimator.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{b.iv}: A \eqn{k \times 1} vector of 2SLS coefficient estimates.
#'   \item \code{covmat.b.iv}: A \eqn{k \times k} HAC covariance matrix of the
#'   estimator.
#' }
#'
#' @examples
#' set.seed(123)
#' T <- 200
#'
#' Z <- matrix(rnorm(T * 2), T, 2)
#' X <- cbind(1, Z[,1] + 0.5 * rnorm(T))
#' Y <- X %*% c(1, 2) + rnorm(T)
#'
#' res <- tsls(Y = as.matrix(Y), X = X, Z = Z, q = 4)
#' res$b.iv
#' res$covmat.b.iv
#'
#'
#' @export
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
  PZ <- Z %*% solve(t(Z)%*%Z) %*% t(Z)
  Y.hat <- PZ %*% Y
  X.hat <- PZ %*% X
  PX.hat <- solve(t(X.hat)%*%X.hat) %*% t(X.hat)
  b.iv <- PX.hat %*% Y.hat

  # Residual estimates:
  eps <- Y - X %*% b.iv

  # Approximation of the covariance matrix of b.iv:
  # Q <- T * solve(t(X) %*% Z %*% solve(t(Z)%*%Z) %*% t(Z) %*% X) %*% t(X.hat)
  Q <- T * solve(t(X) %*% Z %*% solve(t(Z)%*%Z) %*% t(Z) %*% X) %*%
    t(solve(t(Z)%*%Z) %*% t(Z) %*% X)

  Z.times.eps <- Z * matrix(eps,dim(Z)[1],dim(Z)[2])
  S <- NW.LongRunVariance(Z.times.eps,q)

  covmat.b.iv <- 1/T * Q %*% S %*% t(Q)

  return(list(
    b.iv=b.iv,
    covmat.b.iv = covmat.b.iv
  )
  )
}


#' Implementation of the IV-SVAR approach
#'
#'This function implements the IV-SVAR approach.
#'  We are interested in the IRF of Y following a structural shock (\eqn{\eta_1}),
#'  Z being an instrument for \eqn{eta_1}
#'  Y is modelled as a VAR(p) model.
#'  The IRF are scaled in such a way that the contemporaneous impact of eta_1 on Y_1 is 1.
#' @param Y A \code{T x n} numeric matrix
#'   of endogenous variables, with rows as time and columns as series.
#' @param Z A T × k (k ≥ 1) numeric matrix with external instrument(s) for the
#'   targeted structural shock (default: the first shock). If k > 1, columns are used
#'   jointly as instruments. Must be time-aligned with Y.
#' @param p Integer, the VAR lag order for \code{Y}.
#' @param names.of.variables Character vector of length \eqn{n} with the name of the endogenous variables.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param z.AR.order Integer \eqn{\ge 0}. AR order for \code{Z}. Used in the parametric bootstrap only.
#'  Default is \code{3}.
#' @param nb.bootstrap.replications An integer giving the number of bootstrap
#'   replications used to compute confidence intervals. Applies only when a
#'   parametric bootstrap is used. Default is 0 (no bootstrap).
#' @param confidence.interval  numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#'
#' @return
#' \item{IRFs}{A matrix of IRFs (horizons in rows).}
#' \item{all.simulated.IRFs.res}{Array  of
#'   bootstrap IRFs; \code{NULL} if no bootstrap.}
#' \item{all.stdv.IRFs}{Matrix of bootstrap standard
#'   deviations; \code{NULL} if no bootstrap.}
#' \item{Sigma.IV}{Empirical covariance matrix of the stacked innovation vector
#'   \eqn{(\hat{\varepsilon}_Y^\top,\; \hat{u}_Z)^\top} used in the bootstrap;
#'   \code{NULL} if no bootstrap.}
#' \item{cor.Sigma.IV}{Correlation matrix corresponding to \code{Sigma.IV};
#'   \code{NULL} if no bootstrap.}
#' \item{all.CI.lower.bounds, all.CI.upper.bounds}{
#'   Matrices of bootstrap percentile bands; \code{NULL} if no bootstrap.}
#' \item{all.B.tilde.1}{Matrix collecting the normalized impact vectors from
#'   each bootstrap draw (one column per replication); \code{NULL} if no
#'   bootstrap.}
#'
#'
#' @examples
#' data("USmonthly_1990_to_2012")
#'
#' indic.shock.name <- which(names(USmonthly_1990_to_2012)%in%c("FF4_TC","ED2_TC"))
#' Z <- as.matrix(USmonthly_1990_to_2012[,indic.shock.name])
#' considered.variables <- c("GS1","LIP","LCPI","EBP")
#' Y <- as.matrix(USmonthly_1990_to_2012[,considered.variables])
#' n <- length(considered.variables)
#'
#' svar.iv(Y,Z,p = 4,names.of.variables=considered.variables,
#'        nb.periods.IRF = 20,
#'        z.AR.order=1,
#'        nb.bootstrap.replications = 100,
#'        confidence.interval = 0.90,
#'        indic.plot=1)
#' @export
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



#' @export
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
    ZZ <- matrix(Z[(1+p):dim(Z)[1],],ncol=dim(Z)[2])
    eq.iv <- tsls(YY,XX,ZZ,q=3)
    B.tilde.1[i] <- eq.iv$b.iv
    #stdv.B.tilde.1[i] <- sqrt(eq.iv$covmat.b.iv)
  }

  # Compute IRFs
  B.tilde <- matrix(0,n,n)
  B.tilde[,1] <- B.tilde.1

  y0.star <- rep(0,dim(Y)[2]*p)

  nb.sim <- nb.periods.IRF
  IRFs <- simul.VAR(c=rep(0,dim(Y)[2]),
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


# Estimation of IRF
#
# This functions computes impulse response functions.
# The estimation of structural shocks is done by Cholesky decomposition.
# Confidence intervals are obtained by boostrapping the estimated VAR model.
#
# @param Y Numeric matrix of size \eqn{T \times n} containing the endogenous
#   variables (one column per variable). Column names are used for plot titles.
# @param p Integer for the lag order.
# @param posit.of.shock An integer giving the position (column index in \code{Y})
#   of the structural shock of interest. Default is 1.
# @param nb.periods.IRF An integer specifying the number of periods for which
#   impulse response functions are computed. Default is 20.
# @param nb.bootstrap.replications An integer giving the number of bootstrap
#   replications used to compute confidence intervals. Applies only when a
#   parametric bootstrap is used. Default is 0 (no bootstrap).
# @param confidence.interval A numeric value between 0 and 1 indicating the
#   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#   Default is 0.90.
# @param indic.plot Plots are displayed if = 1.
# @return A list with the following components:
# \describe{
#   \item{IRFs}{Matrix of estimated impulse response functions.}
#   \item{all.simulated.IRFs.res}{Array of simulated IRFs across bootstrap replications.}
#   \item{all.stdv.IRFs}{Matrix of standard deviations of the simulated IRFs.}
#   \item{all.CI.lower.bounds}{Matrix of lower bounds of bootstrap confidence intervals.}
#   \item{all.CI.upper.bounds}{Matrix of upper bounds of bootstrap confidence intervals.}
# }
#
# @references
# Christiano, L. J., Eichenbaum, M., and Evans, C. (1996).
# The effects of monetary policy shocks: Evidence from the flow of funds.
# The Review of Economics and Statistics, 78(1):16–34.
# @examples
# data("USmonthlyExample")
# considered.variables <- c("LIP","UNEMP","LCPI","LPCOM","FFR","NBR","TTR","M1")
# y <- as.matrix(USmonthlyExample[considered.variables])
# res.svar.ordering <- svar.ordering(y,p=3,
#                                    posit.of.shock = 5,
#                                    nb.periods.IRF = 20,
#                                    nb.bootstrap.replications = 100,
#                                    confidence.interval = 0.90, # expressed in pp.
#                                    indic.plot = 1 # Plots are displayed if = 1.
# )
# @export

# svar.ordering <- function(Y,p,
#                           posit.of.shock = 1,
#                           nb.periods.IRF = 20,
#                           nb.bootstrap.replications = 0, # This is used in the parametric bootstrap only
#                           confidence.interval = 0.90, # expressed in pp.
#                           indic.plot = 1 # Plots are displayed if = 1.
# ){
#   # This functions computes IRF (potentially with Confidence Intervals) using the
#   # Christiano, Eichembaum and Evans methodology
#
#   names.of.variables <- colnames(Y)
#
#   T <- dim(Y)[1]
#   n <- dim(Y)[2]
#
#   baseline.res <- svar.ordering.aux(Y,p,
#                                     posit.of.shock,
#                                     nb.periods.IRF)
#   IRFs <- baseline.res$IRFs
#
#   if(nb.bootstrap.replications>0){
#     all.simulated.IRFs.res <- array(NaN,c(nb.periods.IRF,n,nb.bootstrap.replications))
#
#     # Computation of std dev of IRF by parametric Gaussian bootstrap
#     # (see Stock and Watson, appendix A.2)
#
#     # step 1: parameterize a VAR for [y',z]':
#     est.VAR <- baseline.res$est.VAR
#     Phi <- Acoef(est.VAR)
#     B.hat <- baseline.res$B.hat
#     cst <- Bcoef(est.VAR)[,p*n+1]
#
#     y0.star <- NULL
#     for(k in p:1){
#       y0.star <- c(y0.star,Y[k,])
#     }
#
#     all.B.tilde.1 <- NULL
#     for(i in 1:nb.bootstrap.replications){
#       simulated.Y <- simul.VAR(cst,Phi,B.hat,nb.sim=dim(Y)[1],y0.star)
#       colnames(simulated.Y) <- colnames(Y)
#
#       simulated.res <- svar.ordering.aux(simulated.Y,p,
#                                          posit.of.shock,
#                                          nb.periods.IRF)
#
#       all.simulated.IRFs.res[,,i] <- simulated.res$IRFs
#     }
#
#     all.stdv.IRFs <- NULL
#     all.CI.lower.bounds <- NULL
#     all.CI.upper.bounds <- NULL
#     for(i in 1:n){
#       all.IRFs.i <- matrix(all.simulated.IRFs.res[,i,],ncol=nb.bootstrap.replications)
#
#       stdv.IRF <- apply(all.IRFs.i,1,sd)
#       all.stdv.IRFs <- cbind(all.stdv.IRFs,stdv.IRF)
#
#       all.CI.lower.bounds <- cbind(all.CI.lower.bounds,
#                                    apply(all.IRFs.i,1,
#                                          function(x){quantile(x,(1-confidence.interval)/2)}))
#       all.CI.upper.bounds <- cbind(all.CI.upper.bounds,
#                                    apply(all.IRFs.i,1,
#                                          function(x){quantile(x,1-(1-confidence.interval)/2)}))
#     }
#
#   }else{
#     all.simulated.IRFs.res <- NULL
#     all.stdv.IRFs <- NULL
#     all.CI.lower.bounds <- NULL
#     all.CI.upper.bounds <- NULL
#   }
#
#   if(indic.plot==1){
#     par(mfrow=c(2,
#                 ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
#     for(i in 1:n){
#       plot(IRFs[,i],type="l",lwd=2,xlab="",ylab="",
#            ylim=c(min(all.CI.lower.bounds[,i]),
#                   max(all.CI.upper.bounds[,i])),
#            main=paste("Effect of shock on ",names.of.variables[i],sep=""))
#       abline(h=0,col="grey")
#       if(nb.bootstrap.replications>0){
#         lines(all.CI.lower.bounds[,i],col="red",lty=2,lwd=2)
#         lines(all.CI.upper.bounds[,i],col="red",lty=2,lwd=2)
#       }
#     }
#   }
#
#   return(
#     list(
#       IRFs = IRFs,
#       all.simulated.IRFs.res = all.simulated.IRFs.res,
#       all.stdv.IRFs = all.stdv.IRFs,
#       all.CI.lower.bounds = all.CI.lower.bounds,
#       all.CI.upper.bounds = all.CI.upper.bounds
#     ))
# }


#' @export
svar.ordering.aux <- function(Y,p,
                              posit.of.shock,
                              nb.periods.IRF){
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # Christiano, Eichembaum and Evans methodology

  considered.variables <- colnames(Y)
  n <- dim(Y)[2]

  # Select number of lags in VAR models:
  #VARselect(Y)

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


#' Local-Projection IRFs à la Jordà (2005)
#'
#' Computes impulse response functions (IRFs) using Jordà’s local-projection
#' method. For each horizon \eqn{h = 0, \dots, H}, the routine runs a separate
#' regression of \eqn{y_{i,t+h}} on the contemporaneous shock series (chosen by
#' \code{posit.of.shock}) and optional lags of the endogenous variables, then
#' stacks the coefficients across horizons. Pointwise standard errors are
#' computed via a Newey–West estimator and responses are normalized by the
#' innovation variance of the shock series.
#'
#'@param Y Numeric matrix of size \eqn{T \times n} containing the endogenous
#'   variables (one column per variable). Column names are used for plot titles.
#' @param posit.of.shock An integer giving the position (column index in \code{Y})
#' of the structural shock of interest. Default is 1.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param nb.lags.endog.var.4.control Integer; number of lags of all endogenous
#'   variables to include as controls in each local projection (default 0).
#' @param confidence.interval A numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#'
#' @return A list with:
#' \item{IRFs}{Numeric matrix \eqn{(H+1) \times n} of IRFs normalized by the
#'   shock’s innovation standard deviation. Row \code{h+1} is horizon \code{h}.}
#' \item{all.stdv.IRFs}{Numeric matrix \eqn{(H+1) \times n} of pointwise
#'   (Newey–West) standard errors corresponding to \code{IRFs}.}
#' @export
#'
#' @references
#' Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
#' Projections. \emph{American Economic Review}, 95(1), 161–182.
#' https://doi.org/10.1257/0002828053828474
#'
#'
#' @examples
#' data("USmonthlyExample")
#' considered.variables <- c("LIP","UNEMP","LCPI","LPCOM","FFR","NBR","TTR","M1")
#' y <- as.matrix(USmonthlyExample[considered.variables])
#' res.jorda <- make.jorda.irf(y,posit.of.shock = 5,
#'                             nb.periods.IRF = 12, nb.lags.endog.var.4.control=3,
#'                           indic.plot = 1, # Plots are displayed if = 1.
#'                           confidence.interval = 0.90)
#' @export
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
    message(paste("Jorda's approach, Currently working on horizon h=",toString(h),sep=""))
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

#' Local Projection IV impulse response functions (LP-IV)
#'
#' @description
#' Computes impulse response functions (IRFs) using the Local Projection
#' Instrumental Variables (LP-IV) approach. For each horizon \eqn{h}, the function
#' estimates a separate IV regression of future outcomes on a shock proxy,
#' instrumented by \code{Z}, optionally controlling for lags of the endogenous
#' variables and the instrument.
#'
#' By convention, the shock is normalized to have a unit effect on the first
#' column of \code{Y}.
#' @details
#' Let \eqn{Y_t} be an \eqn{n}-dimensional vector of endogenous variables and
#' \eqn{Z_t} an instrument. For each horizon \eqn{h = 0, \dots, H} and each variable
#' \eqn{i = 1, \dots, n}, the function estimates:
#' \deqn{
#' Y_{i,t+h} = \alpha_{i,h} + \beta_{i,h} Y_{1,t} + \Gamma_{i,h}' W_t + u_{i,t+h},
#' }
#' using two-stage least squares, where:
#' \itemize{
#'   \item \eqn{Y_{1,t}} is the shock variable (first column of \code{Y}),
#'   \item \eqn{Z_t} is the instrument for \eqn{Y_{1,t}},
#'   \item \eqn{W_t} includes optional lags of \code{Y} and \code{Z},
#'   \item a constant is always included.
#' }
#'
#' The coefficient \eqn{\beta_{i,h}} is interpreted as the impulse response of
#' variable \eqn{i} at horizon \eqn{h}. Standard errors are obtained from the
#' HAC-corrected covariance matrix returned by \code{\link{tsls}}, using a
#' Newey–West lag length of \code{h + 1}.
#'
#' @param Y Numeric matrix of dimension \eqn{T \times n} containing the endogenous
#' variables. The first column is interpreted as the shock variable whose impulse
#' responses are traced out.
#'
#' @param Z Numeric matrix of dimension \eqn{T \times m} containing the instrument(s)
#' for the shock variable.
#'
#' @param nb.periods.IRF Non-negative integer specifying the maximum horizon
#' \eqn{H} for which impulse responses are computed.
#'
#' @param nb.lags.Y.4.control Non-negative integer specifying the number of lags
#' of \code{Y} included as control variables in each local projection.
#'
#' ' @param nb.lags.Z.4.control Non-negative integer specifying the number of lags
#' of \code{Z} included as control variables in each local projection.
#'
#' @param indic.plot Logical or integer (1/0). If equal to 1, impulse responses
#' with confidence intervals are plotted.
#'
#' @param confidence.interval Numeric scalar in \eqn{(0,1)} specifying the
#' confidence level used to construct pointwise confidence intervals.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{IRFs}: A matrix of dimension \eqn{(H+1) \times n} containing the
#'   estimated impulse responses for each variable and horizon.
#'   \item \code{all.stdv.IRFs}: A matrix of dimension \eqn{(H+1) \times n}
#'   containing the corresponding standard errors.
#' }
#'
#' @examples
#' data("USmonthly_1990_to_2012")
#' indic.shock.name <- which(names(USmonthly_1990_to_2012)%in%c("FF4_TC","ED2_TC"))
#' Z <- as.matrix(USmonthly_1990_to_2012[,indic.shock.name])
#' considered.variables <- c("GS1","LIP","LCPI","EBP")
#' Y <- as.matrix(USmonthly_1990_to_2012[,considered.variables])
#' n <- length(considered.variables)
#' res.LP.IV <- make.LPIV.irf(Y,Z,
#'                            nb.periods.IRF = 20,
#'                            nb.lags.Y.4.control=4,
#'                            nb.lags.Z.4.control=4,
#'                            indic.plot = 1, # Plots are displayed if = 1.
#'                            confidence.interval = 0.90)
#'
#' @export
make.LPIV.irf <- function(Y,Z,
                          #posit.of.shock = 1,
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
    message(paste("LP-IV approach, Currently working on horizon h=",toString(h),sep=""))

    for(i in 1:n){
      ZZ <- matrix(Z[(1+max.nb.lags.control):(T-h),],ncol=dim(Z)[2])
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
