
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

