# autocov <- function(X,n){
#   # X is a Txk matrix
#   # returns an estimator of E(X_{t-n}X_t)
#   T <- dim(X)[1]
#   k <- dim(X)[2]
#   vec.1 <- matrix(1,1,k)
#   mean.X <- apply(X,2,mean)
#   X.1 <- X[1:(T-n),] - t(matrix(mean.X,k,T-n))
#   X.2 <- X[(n+1):T,] - t(matrix(mean.X,k,T-n))
#   return(
#     matrix(1/T * apply((X.1 %x% vec.1) * (vec.1 %x% X.2),2,sum),k,k)
#   )
# }
#
# make.F <- function(phi){
#   # Make F matrix for an AR process
#   p <- length(phi)
#   F <- matrix(0,p,p)
#   F[1,] <- phi
#   if(p>1){
#     F[2:p,1:(p-1)] <- diag(p-1)
#   }
#   return(F)
# }
#
# make.dyn.mult <- function(phi,max.dyn.mult){
#   # Compute dynamic mutlipliers for an AR process
#   vec.dyn.mult <- NULL
#   F <- make.F(phi)
#   p <- length(phi)
#   F.j <- diag(p)
#   vec.j <- 0:max.dyn.mult
#   for(j in vec.j){
#     vec.dyn.mult <- c(vec.dyn.mult,F.j[1,1])
#     F.j <- F.j %*% F
#   }
#   return(vec.dyn.mult)
# }
#
#
# sim.arma <- function(c,phi,theta,sigma,T,y.0,nb.sim=1,make.IRF=0,X=NaN,beta=NaN){
#   # nb.sim samples of length T are simulated, all
#   # of them begin with y_0 = y.0, y.0 is of dimension p x 1
#   # sigma is a standard deviation
#   # Warning:
#   #     ->>>> for an AR process, set theta=1.
#   #     ->>>> that is, to simulate an ARMA(p,q), theta has to be of length (q+1).
#   # By default, only one simulation is done (nb.sim=1).
#   #
#   # If make.IRF=0,
#   #     then it means that the user wants to compute IRFs with a maximum horizon T.
#   # X is a matrix of exogenous variables; beta is a vector with appropriate dimensions.
#   p <- length(phi)
#   q <- length(theta)
#   if(length(y.0)!=p){
#     print("y.0 should have the same length as phi.")
#     return(NULL)
#   }
#
#   if(!is.na(X[1])){
#     X.beta <- X %*% matrix(beta,ncol=1)
#   }else{
#     X.beta <- rep(0,T)
#   }
#
#   if(make.IRF==0){# In that case, the user wants to perform standard simulations
#     eps <- sigma*matrix(rnorm(nb.sim*T),T,nb.sim) # These are all the shocks that will be used
#     #eps[1:q,] <- 0
#     eps.t <- matrix(0,q,nb.sim) # eps.t is of dimension (q x nb.sim)
#     # eps.t will change at each iteration, the i^th column corresponds to simulation i.
#     # At date t, the i^th column of eps.t contains, for simulation i: (epsilon[t],epsilon[t-1],...,epsilon[t-q+1])
#   }else{# In that case, the user wants to compute IRFs
#     eps <- matrix(0,T,nb.sim)
#     eps[1,] <- sigma # This is the initial impulsion
#     eps.t <- matrix(0,q,nb.sim) # eps.t is of dimension (q x nb.sim)
#   }
#   F <- make.F(phi) # This is the usual F matrix
#   Y <- NULL
#   y.00 <- matrix(y.0,p,nb.sim)
#   y.t <- y.00
#   for(t in 1:T){
#     if(q>1){
#       eps.t <- rbind(eps[t,],matrix(eps.t[1:(q-1),],q-1,nb.sim))
#     }else{
#       eps.t <- matrix(eps[t,],nrow=1)
#     }
#     theta.eps.t <- matrix(theta,nrow=1) %*% eps.t
#     theta.eps.t <- rbind(theta.eps.t,matrix(0,p-1,nb.sim))
#     y.t <- c(c,rep(0,p-1)) + F %*% y.t + theta.eps.t + X.beta[t]
#     Y <- rbind(Y,y.t[1,])
#   }
#   return(Y)
# }
#
#
#
#
# compute.fevd <- function(Phi,B,H){
#   # Forecast-error variance decomposition
#   p <- length(Phi)
#   n <- dim(Phi[[1]])[1]
#   y0.star <- rep(0,n*p)
#   c <- rep(0,n)
#   variance.decomp = array(0,c(H,n,n))
#   share.variance.decomp = array(0,c(H,n,n))
#   for(i in 1:n){
#     u.shock <- rep(0,n)
#     u.shock[i] <- 1
#     aux <- simul.VAR(c,Phi,B,nb.sim=H,y0.star,indic.IRF=1,u.shock=u.shock)
#     variance.decomp[,,i] <- apply(aux^2,2,cumsum)
#   }
#   # Compute shares:
#   tot.var <- apply(variance.decomp,c(1,2),sum)
#   for(i in 1:n){
#     share.variance.decomp[,,i] <- variance.decomp[,,i]/tot.var
#   }
#   return(list(variance.decomp=variance.decomp,
#               share.variance.decomp=share.variance.decomp))
# }
#
#
# log.l.gaussian <- function(eps,mu,sigma2){
#   # Function that computes the density of vector eps, where the elements
#   # of vector eps are Gaussian i.i.d. variables N(mu,sigma^2):
#   vec.log.f <- - 1/2*log(2*pi*sigma2) - (eps - mu)^2/(2*sigma2)
#   return(vec.log.f)
# }
#
#
# log.lik.armax <- function(THETA,Y,p,q,X=NaN){
#   c <- THETA[1]
#   phi <- THETA[2:(p+1)]
#   if(q>0){
#     theta <- c(1,THETA[(1+p+1):(1+p+q)])
#   }else{
#     theta <- 1
#   }
#   sigma <- THETA[1+p+q+1]
#   if(!is.na(X[1])){
#     if(is.null(dim(X))){
#       r <- 0
#     }else{
#       r <- dim(X)[2] - 1
#     }
#     beta <- THETA[(1+p+q+1+1):(1+p+q+1+(r+1))]
#   }else{
#     beta <- NULL
#   }
#   res <- armax.log.L(Y,c,phi,theta,sigma,X,beta)
#   return(-sum(res$vec.log.l))
# }
#
# armax.log.L <- function(Y,c,phi,theta,sigma,X=NaN,beta=NaN){
#   T <- length(Y)
#   p <- length(phi)
#   q <- length(theta)-1
#   if(is.nan(X[1])){
#     X.beta <- rep(0,T)
#   }else{
#     X.beta <- X %*% matrix(beta,ncol=1)
#     if(length(X.beta)!=T){
#       print("X is not of the same size as Y")
#       return(0)
#     }
#   }
#   if(q>=1){
#     vec.eps <- rep(0,q+1)
#     all.eps <- rep(0,p)
#     for(t in (p+1):T){
#       vec.eps <- c(Y[t] - c - sum(phi*Y[(t-1):(t-p)]) - sum(theta[2:(q+1)]*vec.eps[1:q])  - X.beta[t],vec.eps[1:q])
#       all.eps <- c(all.eps,vec.eps[1])
#     }
#   }else{# q == 0
#     if(p==0){
#       all.eps <- Y - c - X.beta
#     }else{
#       Y_1 <- Y[p:(T-1)]
#       if(p>1){
#         for(i in 2:p){
#           Y_1 <- cbind(Y_1,Y[(p+1-i):(T-i)])
#         }
#       }
#       all.eps <- c(
#         rep(0,p),
#         Y[(p+1):T] - c - Y_1 %*% matrix(phi,ncol=1) - X.beta[(p+1):T]
#       )
#     }
#   }
#   vec.log.l <- log.l.gaussian(all.eps,0,sigma^2)
#   return(list(
#     all.eps = all.eps,
#     vec.log.l = vec.log.l
#   )) # we return -log-lik because the optim procedure of R minimises functions
# }
#
# estim.armax <- function(Y,p,q,X=NaN){
#   # first step: estimate an AR(p)
#   T <- length(Y)
#   if(is.na(X[1])){
#     beta.0 <- NULL
#
#     if(p>=1){
#       y <- Y[(p+1):T]
#       Y_1 <- Y[p:(T-1)]
#       if(p>1){
#         for(i in 2:p){
#           Y_1 <- cbind(Y_1,Y[(p+1-i):(T-i)])
#         }
#       }
#       eq <- lm(y~Y_1)
#       c.0 <- eq$coefficients[1]
#       phi.0 <- eq$coefficients[2:(p+1)]
#       theta.0 <- rep(0,q)
#       sigma.0 <- sd(eq$residuals)
#     }else{
#       y <- Y
#       c.0 <- mean(y)
#       phi.0 <- 0
#       theta.0 <- rep(0,q)
#       sigma.0 <- sd(y)
#     }
#   }else{# exogenous variables
#     if(is.null(dim(X))){
#       r <- 0
#     }else{
#       r <- dim(X)[2] - 1
#     }
#     if(p>=1){
#       y <- Y[(p+1):T]
#       Y_1 <- Y[p:(T-1)]
#       if(p>1){
#         for(i in 2:p){
#           Y_1 <- cbind(Y_1,Y[(p+1-i):(T-i)])
#         }
#       }
#       M <- matrix(X,T,r+1)
#       eq <- lm(y~Y_1+M[(p+1):T,])
#       c.0 <- eq$coefficients[1]
#       phi.0 <- eq$coefficients[2:(p+1)]
#       theta.0 <- rep(0,q)
#       sigma.0 <- sd(eq$residuals)
#       beta.0 <- eq$coefficients[(1+p+1):(1+p+r+1)]
#     }else{
#       y <- Y
#       eq <- lm(y~X)
#       c.0 <- eq$coefficients[1]
#       phi.0 <- NULL
#       theta.0 <- rep(0,q)
#       sigma.0 <- sd(eq$residuals)
#       beta.0 <- eq$coefficients[(1+1):(1+r+1)]
#     }
#   }
#
#   THETA.0 <- c(c.0,phi.0,theta.0,sigma.0,beta.0)
#
#   MAXIT.NlMd <- 100*length(THETA.0)
#   MAXIT.BFGS <- 10
#
#   print("==================================================")
#   print("  ESTIMATING")
#   print("==================================================")
#
#   nb.iter <- 4
#   for(i in 1:nb.iter){
#     res.optim <- optim(par=THETA.0,
#                        fn=log.lik.armax,
#                        Y=Y,
#                        p=p,
#                        q=q,
#                        X=X,
#                        gr = NULL,
#                        method="Nelder-Mead",
#                        #method="CG",
#                        #method="BFGS",
#                        control=list(trace=FALSE,maxit=MAXIT.NlMd),hessian=FALSE)
#     THETA.0 <- res.optim$par
#
#     if(nb.iter==i){
#       hessian.TRUE <- TRUE
#     }else{
#       hessian.TRUE <- FALSE
#     }
#     res.optim <- optim(par=THETA.0,
#                        fn=log.lik.armax,
#                        Y=Y,
#                        p=p,
#                        q=q,
#                        X=X,
#                        gr = NULL,
#                        #method="Nelder-Mead",
#                        #method="CG",
#                        method="BFGS",
#                        control=list(trace=FALSE,maxit=MAXIT.BFGS),hessian=hessian.TRUE)
#     THETA.0 <- res.optim$par
#   }
#
#   print("  END OF ESTIMATION")
#   print("==================================================")
#   print("")
#   print("  RESULTS:")
#   print("  -----------------------")
#
#   THETA <- THETA.0
#   I <- solve(res.optim$hessian)
#
#   st.dev <- sqrt(diag(I))
#   t.ratio <- THETA/st.dev
#
#   res.matrix <- as.matrix(cbind(THETA,st.dev,t.ratio)) # to be printed
#
#   vec.names <- c("c")
#
#   if(p>0){
#     for(i in 1:p){
#       vec.names <- c(vec.names,paste("phi   t-",toString(i),sep=""))
#     }
#   }
#   if(q>0){
#     for(i in 1:q){
#       vec.names <- c(vec.names,paste("theta t-",toString(i),sep=""))
#     }
#   }
#   vec.names <- c(vec.names,"sigma")
#   if(!is.na(X[1])){
#     for(i in 0:r){
#       vec.names <- c(vec.names,paste("beta  t-",toString(i),sep=""))
#     }
#   }
#
#   rownames(res.matrix) <- vec.names
#
#   print(res.matrix)
#
#   print("==================================================")
#
#   c <- THETA[1]
#   phi <- THETA[2:(p+1)]
#   if(q>0){
#     theta <- c(1,THETA[(1+p+1):(1+p+q)])
#   }else{
#     theta <- 1
#   }
#   sigma <- THETA[1+p+q+1]
#   if(!is.na(X[1])){
#     if(is.null(dim(X))){
#       r <- 0
#     }else{
#       r <- dim(X)[2] - 1
#     }
#     beta <- THETA[(1+p+q+1+1):(1+p+q+1+(r+1))]
#   }else{
#     beta <- NULL
#   }
#
#   return(list(
#     c=c,phi=phi,theta=theta,sigma=sigma,beta=beta,
#     max.logl = -res.optim$value,I=I,st.dev=st.dev,THETA=THETA
#   ))
# }
#
#
# simul.garch <- function(zeta,alpha,delta,Nb.simul){
#   # Compute unconditional variance:
#   kappa <- zeta*(1-sum(delta))
#   Eu2 <- kappa/(1 - sum(alpha) - sum(delta))
#   m <- length(alpha)
#   r <- length(delta)
#   sigma2 <- zeta + sum(alpha)*Eu2/(1 - sum(delta))
#   vec.u2 <- rep(Eu2,m)
#   vec.h <- rep(Eu2,r)
#   all.u <- NULL
#   all.h <- NULL
#   for(i in 1:Nb.simul){
#     u <- sqrt(vec.h[1])*rnorm(1)
#     if(m>1){
#       vec.u2 <- c(u^2,vec.u2[1:(m-1)])
#     }else{
#       vec.u2 <- u^2
#     }
#     h <- zeta*(1-sum(delta)) + sum(delta * vec.h) + sum(alpha * vec.u2)
#     if(r>1){
#       vec.h <- c(h,vec.h[1:(r-1)])
#     }else{
#       vec.h <- h
#     }
#     all.u <- c(all.u,u)
#     all.h <- c(all.h,h)
#   }
#   return(list(all.u=all.u,all.h=all.h,Eu2=Eu2,kappa=kappa))
# }
#
# compute.garch <- function(theta,x,m0,r0){
#   zeta <- abs(theta[1])
#   if(m0>0){
#     m <- m0
#     alpha <- abs(theta[2:(m0+1)])
#   }else{
#     m <- 1
#     alpha <- 0
#   }
#   if(r0>0){
#     r <- r0
#     delta <- abs(theta[(m0+2):(1+m0+r0)])
#   }else{
#     r <- 1
#     delta <- 0
#   }
#   p <- max(m,r)
#   T <- length(x)
#   kappa <- zeta*(1-sum(delta))
#   Eu2 <- kappa/(1 - sum(alpha) - sum(delta))
#   vec.u2 <- x[p:(p-m+1)]^2
#   vec.h <- rep(Eu2,r)
#
#   h <- NULL
#   for(i in (p+1):T){
#     last.h <- kappa + sum(alpha * vec.u2) + sum(delta * vec.h)
#     h <- c(h,last.h)
#     if(m>1){
#       vec.u2 <- c(x[i]^2,vec.u2[1:(m-1)])
#     }else{
#       vec.u2 <- x[i]^2
#     }
#     if(r>1){
#       vec.h <- c(last.h,vec.h[1:(r-1)])
#     }else{
#       vec.h <- last.h
#     }
#   }
#   aux <- sum(alpha) + sum(delta)
#   if(aux >= 1){
#     logf <- 10000000000000000000 * (aux - 1)
#   }else{
#     logf <- -sum(-1/2*log(2*pi*h) - x[(p+1):T]^2/(2*h))
#   }
#
#   names.param <- "zeta"
#   if(m0>0){
#     names.param <- c(names.param,paste("alpha",1:m0,sep=""))
#   }
#   if(r0>0){
#     names.param <- c(names.param,paste("delta",1:r0,sep=""))
#   }
#   return(list(logf=logf,h = c(rep(NaN,p),h),names.param=names.param))
# }



simul.VAR <- function(c,Phi,B,nb.sim,y0.star,indic.IRF=0,u.shock=0,eta=NaN){
  # This function simulates a VAR model, initial condition = y0.star
  # Phi is a list, each element of which is a Phi_i matrix. Hence it has p elements if we consider a VAR(p)
  # If eta is different from NaN then the structural shocks are those defined in this
  #      matrix of dimension (nb.sim x n)
  p <- length(Phi)
  n <- dim(Phi[[1]])[1]
  # Check that right dimension for y0.star:
  if((indic.IRF==0)&(length(y0.star)!=(n*p))){
    print("The dimension of y0.star should be (np x 1) where p is the number of lags in the VAR and n is the dimension of y")
    return(0)
  }
  if((indic.IRF!=0)&(length(u.shock)!=n)){
    print("If you want to compute IRFs, u.shock has to be of length n, where n is the number of dependent variables")
  }
  PHI <- make.PHI(Phi)
  c.star <- c(c,rep(0*c,p-1))
  B.star <- matrix(0,n*p,n)
  B.star[1:n,1:n] <- B
  y <- y0.star
  Y <- NULL
  if(class(eta)[1]=="matrix"){
    eps <- eta %*% t(B.star)
    nb.sim <- dim(eta)[1]
  }else{
    eps <- matrix(rnorm(nb.sim*n),nb.sim,n) %*% t(B.star)
  }
  for(t in 1:nb.sim){
    if(indic.IRF==0){
      y <- c.star + PHI %*% y + c(eps[t,])
    }else{
      if(t==1){
        y <- B.star %*% c(u.shock)
      }else{
        y <- PHI %*% y
      }
    }
    Y <- rbind(Y,c(y))
  }
  return(Y[,1:n])
}

simul.VARMA <- function(Model,nb.sim,Y0,eta0,indic.IRF=0){
  # Model is a list containing:
  # Mu (vector of constants),
  # Phi (array of autoregressive matrices),
  # Theta (array of MA matrices),
  # C (matrix of dimension n x n, this is the mixing matrix)
  # distri (which characterizes the distribution of the shocks)
  # Y0   contains the initial values of Y,   it has to be of dimension (p x n) * 1 (concatenation of Y_1,...,Y_p)
  # eta0 contains the initial values of eps, it has to be of dimension (q x n) * 1 (concatenation of eta_1,...,eta_q)
  # Notations:
  # n is the dimension of Y, p is the AR order, q is the MA order.
  n <- dim(Model$Phi)[1]
  p <- dim(Model$Phi)[3]
  if(is.matrix(Model$Theta)==1){
    q <- 1
  }else{
    q <- dim(Model$Theta)[3]
  }

  MU <- c(Model$Mu,rep(0,n*(p-1)))
  PHI <- make.PHI(Model$Phi)
  THETA <- cbind(diag(n),-matrix(Model$Theta,nrow=n))
  THETA <- rbind(THETA,matrix(0,n*(p-1),n*(q+1)))
  CC <- diag(q+1) %x% Model$C

  y <- Y0
  eta <- eta0

  if(indic.IRF==1){
    eta.simul <- matrix(0,nb.sim,n)
    eta <- 0*eta
    MU <- 0*MU
  }else{
    eta.simul <- simul.distri(Model$distri,nb.sim)
  }
  eta.simul[1,] <- eta0[1:n]

  Y <- NULL
  EPS <- NULL
  V <- NULL
  ETA <- NULL
  for(t in 1:nb.sim){
    eta <- c(eta.simul[t,],eta[1:(n*q)])
    eps <- CC %*% eta
    v <- THETA %*% eps
    y <- MU + PHI %*% y + v
    Y <- cbind(Y,y)
    ETA <- cbind(ETA,eta)
    EPS <- cbind(EPS,eps)
    V <- cbind(V,v)
  }
  return(list(Y=Y,EPS=EPS,ETA=ETA,V=V))
}

make.PHI <- function(Phi){
  if(class(Phi)=="list"){
    p <- length(Phi)
    n <- dim(Phi[[1]])[1]
  }else{
    p <- dim(Phi)[3]
    n <- dim(Phi)[1]
  }
  PHI <- matrix(0,n*p,n*p)
  if(p>1){
    PHI[(n+1):(n*p),1:((p-1)*n)] <- diag((p-1)*n)
  }
  for(i in 1:p){
    if(class(Phi)=="list"){
      PHI[1:n,((i-1)*n+1):(i*n)] <- Phi[[i]]
    }else{
      PHI[1:n,((i-1)*n+1):(i*n)] <- Phi[,,i]
    }
  }
  return(PHI)
}

simul.distri <- function(distri,nb.sim,basic.drawings=NaN){
  # Simulation of independent shocks
  eps <- NULL
  nb.var <- length(distri$type) # number of variables
  if(is.na(basic.drawings[1])){
    U <- matrix(runif(nb.var*nb.sim),nb.sim,nb.var)
  }else{
    if(dim(basic.drawings)[1]!=nb.sim){
      print("Dimension of basic.drawing not consistent with nb.sim")
      return(0)
    }else{
      U <- basic.drawings
    }
  }

  for(i in 1:nb.var){
    if(distri$type[i]=="gaussian"){
      eps.i <- qnorm(U[,i])
    }else if(distri$type[i]=="mixt.gaussian"){
      eps.i <- qmixt(alpha=U[,i],distri$mu[i],distri$sigma[i],distri$p[i])
    }else if(distri$type[i]=="student"){
      nu <- distri$df[i]
      eps.i <- qt(U[,i],df = nu)/sqrt(nu/(nu-2))
    }
    eps <- cbind(eps,eps.i)
  }
  return(eps)
}


