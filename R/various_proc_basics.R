
simul.VAR <- function(c,Phi,B,nb.sim,y0.star=NaN,indic.IRF=0,u.shock=0,eta=NaN){
  # This function simulates a VAR model, initial condition = y0.star
  # Phi is a list, each element of which is a Phi_i matrix. Hence it has p elements if we consider a VAR(p)
  # If eta is different from NaN then the structural shocks are those defined in this
  #      matrix of dimension (nb.sim x n)

  # Make sure Phi is a list, convert it if this is not the case:
  if(class(Phi)[1]!="list"){
    if(class(Phi)[1]=="matrix"){
      Phi <- list(Phi)
    }else{# In thisd case, Phi is an array
      Phi <- lapply(seq_len(dim(Phi)[3]), function(k) Phi[ , , k])
    }
  }

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

simul.VARMA <- function(Model,nb.sim,Y0=NaN,eta0,indic.IRF=0){
  # Model is a list containing:
  # c (vector of constants),
  # Phi (array of autoregressive matrices),
  # Theta (array of MA matrices),
  # C, or B (matrix of dimension n x n, this is the mixing matrix)
  # distri (which characterizes the distribution of the shocks),
  #    By default, "distri" is set to standard normal distributions.
  # Y0   contains the initial values of Y,   it has to be of dimension (p x n) * 1 (concatenation of Y_1,...,Y_p)
  # eta0 contains the initial values of eps, it has to be of dimension (q x n) * 1 (concatenation of eta_1,...,eta_q)
  # Notations:
  # n is the dimension of Y, p is the AR order, q is the MA order.
  # !!! Important Note !!! -----------------------------------------------------
  # In Model, if B is specified, and not C, then the VARMA specification is:
  # y_t = c + Phi1·y_{t-1} + ... + Phip·y_{t-p} +
  #       B·eta_t + Theta1·B·eta_{t-1} + ... + Thetaq·B·eta_{t-q}
  # otherwise, if C is specified, the VARMA specification is:
  # y_t = c + Phi1·y_{t-1} + ... + Phip·y_{t-p} +
  #       C·eta_t - Theta1·C·eta_{t-1} - ... - Thetaq·C·eta_{t-q}
  # (Notice the change in the signs in the MA part)
  # ----------------------------------------------------------------------------
  n <- dim(Model$Phi)[1]
  p <- dim(Model$Phi)[3]

  # The model is specified via B or C?:
  indic_B_and_not_C <- !is.null(Model$B)

  if(is.null(Model$Theta)){
    q <- 0
    THETA <- diag(n)
  }else{
    if (is.matrix(Model$Theta) == 1) {
      q <- 1
      THETA <- cbind(diag(n),
                     ((+1)*indic_B_and_not_C + (-1)*!indic_B_and_not_C)*
                       matrix(Model$Theta, nrow = n))
    } else {
      q <- dim(Model$Theta)[3]
      THETA <- cbind(diag(n),
                     ((+1)*indic_B_and_not_C + (-1)*!indic_B_and_not_C)*
                       matrix(Model$Theta, nrow = n))
    }
  }
  THETA <- rbind(THETA, matrix(0, n * (p - 1), n * (q + 1)))

  if(indic_B_and_not_C){
    CC <- diag(q + 1) %x% Model$B
  }else{
    CC <- diag(q + 1) %x% Model$C
  }

  # Set type of distributions if missing:
  if(is.null(Model$distri)){
    Model$distri <- list(type=rep("gaussian",n))
  }

  MU <- c(Model$c, rep(0, n * (p - 1)))
  PHI <- make.PHI(Model$Phi)

  if(is.na(Y0[1])){
    if(indic.IRF==1){
      Y0 <- rep(0,n*p)
    }else{
      Y0 <- solve(diag(n*p) - PHI) %*% MU
    }
  }

  y <- Y0
  eta <- eta0
  if (indic.IRF == 1) {
    eta.simul <- matrix(0, nb.sim, n)
    eta <- 0 * eta
    MU <- 0 * MU
  } else {
    eta.simul <- simul.distri(Model$distri, nb.sim)
  }
  eta.simul[1, ] <- eta0[1:n]
  Y <- NULL
  EPS <- NULL
  V <- NULL
  ETA <- NULL
  for (t in 1:nb.sim) {
    if(q > 0){
      eta <- c(eta.simul[t, ], eta[1:(n * q)])
    }else{
      eta <- eta.simul[t, ]
    }
    eps <- CC %*% eta
    v <- THETA %*% eps
    y <- MU + PHI %*% y + v
    Y <- cbind(Y, y)
    ETA <- cbind(ETA, eta)
    EPS <- cbind(EPS, eps)
    V <- cbind(V, v)
  }
  return(list(Y = Y, EPS = EPS, ETA = ETA, V = V))
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


make_variance_decompo <- function(Phi,B,maxHorizon,
                                  mfrow=NaN,
                                  names.var=NaN,
                                  names.shock=NaN){
  # This function produces a plot showing the variance decomposition
  # associated with a VAR model defined by autoregressive matrices Phi,
  # and an impact matrix B.
  # 'mfrow' defines the plot layout.
  # 'names.var' allows to provide variable names (1,...,n by default).
  # 'names.shock' allows to provide shock names (1,...,n by default).

  n <- dim(B)[1] # dimension of the state vector
  if(is.na(names.var[1])){
    names.var <- paste("Variable ", 1:n)
  }
  if(is.na(names.shock[1])){
    names.shock <- paste("Shock ", 1:n)
  }

  IRFs <- array(NaN,c(n,n,maxHorizon))
  for(i in 1:n){
    u.shock <- rep(0,n)
    u.shock[i] <- 1
    Y <- simul.VAR(c=matrix(0,n,1),Phi=Phi,B = B,nb.sim=maxHorizon,
                   indic.IRF = 1,u.shock = u.shock)
    IRFs[,i,] <- t(Y)
  }
  res <- variance.decomp(IRFs)$vardecomp

  if(is.na(mfrow[1])){
    par(mfrow=c(1,n)) # define plot margins
  }else{
    par(mfrow=mfrow) # define plot margins
  }
  par(mar = c(6, 4, 4, 2) + 0.1) # define plot margins

  for(variable in 1:n){
    res_variable <-
      res[variable,variable,,]
    res_variable <- apply(res_variable,1,cumsum)
    res_variable <- rbind(0,res_variable)
    plot(1:maxHorizon,rep(0,maxHorizon),ylim=c(0,1),col="white",las=1,
         xlab="Horizon",ylab="Variance share",main=names.var[variable])
    for(k in 1:n){
      polygon(c(1:maxHorizon,maxHorizon:1),
              c(res_variable[k,],rev(res_variable[k+1,])),col=k+1)
    }

    # Draw legend below the plot
    legend("bottom",
           inset = -0.3,              # how far below (negative pushes it out)
           horiz = TRUE,              # horizontal layout
           bty = "n",                 # no box
           pch = 15,                  # filled square symbols
           col = 2:(n+1),                 # colors
           legend = names.shock,
           xpd = TRUE,                # allow drawing outside plot area
           pt.cex = 1.5, cex = 0.9)
  }
  return(0)
}

