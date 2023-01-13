


# ===============================================
# ===============================================
# Residual test
# ===============================================
# ===============================================



BaiNg.normality.test <- function(Y,kurt=3,skew=0,q=4){
  X <- Y[!is.na(Y)]
  T <- length(X)

  # Compute first four central moments:
  mu.1 <- mean(X)
  XX <- (X - mu.1)
  mu.2 <- mean(XX^2)
  sigma <- sqrt(mu.2)
  mu.3 <- mean(XX^3)
  mu.4 <- mean(XX^4)

  # Skewness test:
  tau   <- mu.3 / sigma^3
  kappa <- mu.4 / sigma^4
  alpha <- matrix(c(1,-3*sigma^2,-3*sigma*tau/2),nrow=1)
  Z <- cbind(XX^3-mu.3,XX,XX^2-sigma^2)
  Gamma <- NW.LongRunVariance(Z,q)
  alpha.2 <- alpha[1:2]
  Gamma.22 <- Gamma[1:2,1:2]
  var.tau <- 1/T * matrix(alpha.2,nrow=1) %*% Gamma.22 %*% matrix(alpha.2,ncol=1) / sigma^6
  pi.3 <- (tau - skew) / sqrt(var.tau)
  pval.3 <- 2*(1-pnorm(abs(pi.3)))

  # Kurtosis test:
  Beta <- matrix(c(1,-4*mu.3,-2*sigma^2*kappa),nrow=1)
  W <- cbind(XX^4-mu.4,XX,XX^2-sigma^2)
  Omega <- NW.LongRunVariance(W,q)
  var.kappa <- 1/T * Beta %*% Omega %*% t(Beta) / sigma^8
  pi.4 <- (kappa - kurt)/sqrt(var.kappa)
  pval.4 <- 2*(1-pnorm(abs(pi.4)))

  # Joint test:
  pi.3.4 <- pi.3^2 + pi.4^2
  pval.3.4 <- 1-pchisq(pi.3.4,df=2)

  return(list(
    var.tau = var.tau,
    var.kappa = var.kappa,
    pi.3 = pi.3,
    pi.4 = pi.4,
    pval.3 = pval.3,
    pval.4 = pval.4,
    pi.3.4 = pi.3.4,
    pval.3.4 = pval.3.4
  ))
}

LV.normality.test <- function(Y){
  X <- Y[!is.na(Y)]
  T <- length(X)

  # Compute first four central moments:
  mu.1 <- mean(X)
  XX <- (X - mu.1)
  mu.2 <- mean(XX^2)
  mu.3 <- mean(XX^3)
  mu.4 <- mean(XX^4)

  gamma.0 <- autocov(matrix(X,ncol=1),n = 0)

  F.3 <- gamma.0^3
  F.4 <- gamma.0^4

  for(i in 1:(T-1)){
    gamma.j <- autocov(matrix(X,ncol=1),n = i)
    F.3 <- F.3 + 2*gamma.j^3
    F.4 <- F.4 + 2*gamma.j^4
  }

  G <- T*mu.3^2/(6 * F.3) + T*(mu.4 - 3*mu.2^2)^2/(24 * F.4)

  # Joint test:
  pval <- 1-pchisq(G,df=2)

  return(list(
    G = G,
    pval = pval
  ))
}

Mardia.multi.normality <- function(X){
  X.bar <- apply(X,2,mean)
  T <- dim(X)[1]
  d <- dim(X)[2]
  X.demeaned <- X - matrix(1,T,1) %*% matrix(X.bar,nrow=1)
  S <- 1/T * t(X.demeaned) %*% X.demeaned
  S_1 <- solve(S)

  Big.Matrix <- X.demeaned %*% S_1 %*% t(X.demeaned)
  # The entries on the diagonal are the D_i^2, the extra-diagonal entries are the D_ij's

  D_i.2 <- diag(Big.Matrix)
  k.d <- 1/T * sum(sqrt(D_i.2)^4)
  b.d <- 1/T^2 * sum(Big.Matrix^3)

  test.stat.b.d <- 1/6 * T * b.d
  test.stat.k.d <- (k.d - d*(d+2))/sqrt(8*d*(d+2)/T)

  df <- d*(d+1)*(d+2)/6

  pval.b.d <- 1-pchisq(test.stat.b.d,df)
  pval.k.d <- 1-pnorm(test.stat.k.d)

  return(list(
    test.stat.b.d = test.stat.b.d,
    test.stat.k.d = test.stat.k.d,
    pval.b.d = pval.b.d,
    pval.k.d = pval.k.d
  ))
}

Lutkepohl_Theilen.normality.test <- function(X){
  X.bar <- apply(X,2,mean)
  T <- dim(X)[1]
  d <- dim(X)[2]
  X.demeaned <- X - matrix(1,T,1) %*% matrix(X.bar,nrow=1)
  S <- 1/T * t(X.demeaned) %*% X.demeaned
  S_1 <- solve(S)

  L <- chol(S)
  L_1 <- solve(L)

  Z <- X.demeaned %*% t(L_1)

  Z3 <- Z^3
  aux <- apply(Z3,2,mean)
  b.1.L <- T/6 * sum(aux^2)

  Z4.exc <- Z^4 - 3
  aux <- apply(Z4.exc,2,mean)
  b.2.L <- T/24 * sum(aux^2)

  D.L <- b.1.L + b.2.L
  pval <- 1 - pchisq(D.L,2*d)

  return(list(
    D.L = D.L,
    pval = pval
  ))
}







# ===============================================
# ===============================================
# SSVARMA ML Estimation
# ===============================================
# ===============================================




ML.inversion.loglik <- function(param0,Model,Y){
  logl <- ML.inversion.loglik.aux(param0,Model,Y)
  return(-sum(Re(logl)))
}


ML.inversion.loglik.aux <- function(param0,Model,Y){

  save(param0,file="results/tempo.res")

  Model.aux <- Model

  n <- dim(Model$Phi)[1]
  p <- dim(Model$Phi)[3]

  if(length(param0)>2*n^2+3*n){
    # we also estimate Phi
    Model.aux$Phi[,,] <- param0[1:(n^2*p)]
    param <- param0[(n^2*p+1):length(param0)]
  }else{
    param <- param0
  }

  Theta.C <- make.Theta.C(param,n)
  Model.aux$C <- Theta.C$C
  Model.aux$Theta <- array(Theta.C$Theta,c(n,n,1))
  Model.aux$Mu <- 0 * Model$Mu
  if(length(Theta.C$param.4.distri)>0){
    Model.aux$distri$type  <- rep("mixt.gaussian",n)

    theta1 <- Theta.C$param.4.distri[1:n]
    theta2 <- Theta.C$param.4.distri[(n+1):(2*n)]
    theta3 <- Theta.C$param.4.distri[(2*n+1):(3*n)]

    VEC.PARAM <-  make.p.mu.sigma(theta1,theta2,theta3)
    p     <- VEC.PARAM$p
    mu    <- VEC.PARAM$mu
    sigma <- VEC.PARAM$sigma

    Model.aux$distri$sigma <- sigma
    Model.aux$distri$p     <- p
    Model.aux$distri$mu    <- mu
  }

  estim.shocks <- estim.struct.shocks(Y,Model.aux)
  EPS.est <- estim.shocks$EPS.est
  eigen.theta <- eigen(Theta.C$Theta)$values
  eigen.theta.larger.than.one <- abs(eigen.theta[abs(eigen.theta)>1])
  logl <- g(EPS.est,Model.aux) -
    log(prod(eigen.theta.larger.than.one))*(length(eigen.theta.larger.than.one)>0)
  logl[logl < -10^(100)] <- -10000
  return(logl)
}


grad.ML.inversion.loglik <- function(param0,Model,Y,eps){
  K <- length(param0)
  grad.Logl <- NULL
  logl.0 <- ML.inversion.loglik.aux(param0,Model,Y)
  for(i in 1:K){
    param.i <- param0
    param.i[i] <- param.i[i] + eps[i]
    logl.i <- ML.inversion.loglik.aux(param.i,Model,Y)
    grad.Logl <- cbind(grad.Logl,(logl.i - logl.0)/eps[i])
  }
  return(grad.Logl)
}


estim.MA.inversion <- function(param.ini,Model,V,
                               MAXIT.nlminb=30,MAXIT.NlMd=100,nb.loops=1,
                               indic.Blaschke=0,
                               MAXIT.nlminb.BM=30,MAXIT.NlMd.BM=100,nb.loops.BM=1,
                               indic.print=1,indic.compute.cov.mat=1){

  n <- dim(Model$Phi)[1]
  p.order <- dim(Model$Phi)[3]

  if(length(param.ini)==2*n^2){
    indic.estim.mixtures <- 0
    indic.estim.phi <- 0
  }else{
    if(length(param.ini)==2*n^2+3*n){
      indic.estim.mixtures <- 1
      indic.estim.phi <- 0
    }else{
      indic.estim.mixtures <- 1
      indic.estim.phi <- 1
    }
  }

  print("param.ini:")
  print(param.ini)

  res.estim <- estim.MA.inversion.aux(param.ini,Model,V,
                                      MAXIT.nlminb,MAXIT.NlMd,nb.loops,indic.print)
  all.loglik <- res.estim$max.log.lik

  all.res.estim <- NULL
  all.res.estim$res.estim.1 <- res.estim

  if(indic.Blaschke==1){
    all.MA <- compute.all.forms(res.estim$Theta,res.estim$C)
    all.MA$all.Theta <- Re(all.MA$all.Theta)
    all.MA$all.C <- Re(all.MA$all.C)
    nb.cases <- dim(all.MA$all.Theta)[3]
    for(i in 1:(nb.cases-1)){#the last element of all.MA (last layers) correspond to the original estimation
      if(indic.print==1){
        print("======================================================")
        print(paste("End of estimation ",toString(i)," out of ",toString(nb.cases),sep=""))
        print("======================================================")
      }
      if((indic.estim.mixtures==1)&(indic.estim.phi==0)){
        param.ini.i <- c(
          all.MA$all.Theta[,,i],
          all.MA$all.C[,,i],
          param.ini[(2*n^2+1):(length(param.ini))]
        )
      }
      if((indic.estim.mixtures==0)&(indic.estim.phi==0)){
        param.ini.i <- c(
          all.MA$all.Theta[,,i],
          all.MA$all.C[,,i]
        )
      }
      if((indic.estim.mixtures==1)&(indic.estim.phi==1)){
        param.ini.i <- c(
          param.ini[1:(n^2*p.order)],
          all.MA$all.Theta[,,i],
          all.MA$all.C[,,i],
          param.ini[(n^2*p.order + 2*n^2+1):(length(param.ini))]
        )
      }
      print("param.ini.i:")
      print(param.ini.i)
      save(param.ini.i,file="essai.res")
      res.estim.i <- estim.MA.inversion.aux(param.ini.i,Model,V,
                                            MAXIT.nlminb.BM,MAXIT.NlMd.BM,nb.loops.BM,indic.print)
      all.loglik <- c(all.loglik,res.estim.i$max.log.lik)

      eval(parse(text =
                   gsub(" ","",
                        paste("all.res.estim$res.estim.",toString(i+1)," <- res.estim.i"))))

      if(res.estim.i$max.log.lik>res.estim$max.log.lik){
        res.estim <- res.estim.i
      }
    }
  }

  res.estim$all.res.estim <- all.res.estim

  res.estim$all.loglik <- all.loglik

  if(indic.compute.cov.mat==1){
    # compute covariance matrix of estimates

    if(indic.print==1){print("Computing Covariance matrix of estimator")}

    res.optim <- optim(res.estim$param,ML.inversion.loglik,
                       Model = Model,
                       Y=V,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="CG",
                       #method="BFGS",
                       control=list(trace=FALSE,maxit=5),hessian=TRUE)

    J_n <- res.optim$hessian

    # Look for unused parameters:
    indic.no.move <- NULL
    for(i in 1:dim(J_n)[1]){
      if(sum(abs(J_n[i,]))==0){
        indic.no.move <- c(indic.no.move,i)
      }
    }
    if(length(indic.no.move)>0){
      for(i in indic.no.move){
        J_n[i,i] <- 1
      }
    }

    res.J.eig <- eigen(J_n)
    J_n_1 <- res.J.eig$vectors %*% diag(1/res.J.eig$values) %*% t(res.J.eig$vectors)

    dlogl <- grad.ML.inversion.loglik(res.estim$param,Model,Y=V,
                                      eps = pmax(abs(res.estim$param)*.001,0.00001))

    V_n <- t(dlogl) %*% dlogl

    indic.no.move <- NULL
    for(i in 1:dim(V_n)[1]){
      if(sum(abs(V_n[i,]))==0){
        indic.no.move <- c(indic.no.move,i)
      }
    }
    if(length(indic.no.move)>0){
      for(i in indic.no.move){
        V_n[i,i] <- 1
      }
    }

    MV <- J_n_1

    save(J_n,J_n_1,V_n,res.estim,file="results/J_n.sav")

    st.dev <- sqrt(diag(MV))
    t.stud <- res.estim$param/sqrt(diag(MV))
    param <- res.estim$param
    if(indic.print==1){
      print(
        cbind(param,st.dev,t.stud)
      )}

    res.estim$MV <- MV

  }

  return(res.estim)
}



estim.MA.inversion.aux <- function(param.ini,Model,V,
                                   MAXIT.nlminb,MAXIT.NlMd,nb.loops,indic.print){

  n <- dim(Model$Phi)[1]
  p.order <- dim(Model$Phi)[3]

  if(length(param.ini)==2*n^2){
    indic.estim.mixtures <- 0
    indic.estim.phi <- 0
  }else{
    if(length(param.ini)==2*n^2+3*n){
      indic.estim.mixtures <- 1
      indic.estim.phi <- 0
    }else{
      indic.estim.mixtures <- 1
      indic.estim.phi <- 1
    }
  }

  if(indic.print==1){print(paste("Model with n=",toString(n),sep=""))}

  logl.ini <- ML.inversion.loglik(param.ini,Model,V)
  if(indic.print==1){print(paste("Log-lik at param.ini: ",
                                 toString(logl.ini),sep=""))}

  for(ii in 1:nb.loops){
    # Use Optimx ==========================================
    if(indic.print==1){print("Use of nlminb")}

    if(MAXIT.nlminb>0){
      maxlik <- optimx(param.ini,ML.inversion.loglik,
                       Model = Model,
                       Y=V,
                       method = "nlminb",
                       control=list(trace=0, maxit=MAXIT.nlminb, kkt = FALSE))
      if(maxlik$value < logl.ini){
        param.est <- c(as.matrix(maxlik)[1:length(param.ini)])
        param.ini <- param.est
        logl.ini <- maxlik$value
      }
    }

    if(indic.print==1){print("Use of Nelder-Mead")}
    res.optim <- optim(param.ini,ML.inversion.loglik,
                       Model = Model,
                       Y=V,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="CG",
                       #method="BFGS",
                       control=list(trace=FALSE,maxit=MAXIT.NlMd))
    if(res.optim$value < logl.ini){
      param.est <- res.optim$par
      param.ini <- param.est
      logl.ini <- res.optim$value
    }else{
      param.est <- param.ini
    }

    if(indic.estim.phi==0){
      Phi.est <- Model$Phi
      res.est <- make.Theta.C(param.est,n)
      Theta.est <- res.est$Theta
      C.est <- res.est$C
    }else{
      Phi.est <- array(param.est[1:(n^2*p.order)],c(n,n,p.order))
      res.est <- make.Theta.C(param.est[(n^2*p.order+1):length(param.est)],n)
      Theta.est <- res.est$Theta
      C.est <- res.est$C
    }

    if(length(res.est$param.4.distri)>0){
      # This change in variable esnures that E(X)=0 and Var(X)=1 is feasible:
      theta1 <- res.est$param.4.distri[1:n]
      theta2 <- res.est$param.4.distri[(n+1):(2*n)]
      theta3 <- res.est$param.4.distri[(2*n+1):(3*n)]

      VEC.PARAM <-  make.p.mu.sigma(theta1,theta2,theta3)
      p     <- VEC.PARAM$p
      mu    <- VEC.PARAM$mu
      sigma <- VEC.PARAM$sigma

      distri.est <- list(
        type = rep("mixt.gaussian",n),
        sigma = sigma,
        p = p,
        mu = mu
      )
    }else{
      distri.est <- Model$distri
    }


    if(indic.print==1){
      print(paste("Phi at the end of iteration ",toString(ii),":",sep=""))
      print(Phi.est)
      print(paste("Theta at the end of iteration ",toString(ii),":",sep=""))
      print(Theta.est)
      print(paste("C at the end of iteration ",toString(ii),":",sep=""))
      print(C.est)
      print(paste("Distri at the end of iteration ",toString(ii),":",sep=""))
      print(distri.est)
      print(paste("Log-likelihood at the end of iteration ",toString(ii),":",sep=""))
      print(ML.inversion.loglik(param.ini,Model,V))}

  }

  max.log.lik <- - ML.inversion.loglik(param.ini,Model,V)

  return(list(
    param=param.est,
    Phi=Phi.est,
    Theta=Theta.est,
    C=C.est,
    distri=distri.est,
    max.log.lik = max.log.lik
  ))

}


compute.gamma <- function(param){
  n <- length(param)/3

  theta1 <- param[1:n]
  theta2 <- param[(n+1):(2*n)]
  theta3 <- param[(2*n+1):(3*n)]

  VEC.PARAM <-  make.p.mu.sigma(theta1,theta2,theta3)
  p     <- VEC.PARAM$p
  mu    <- VEC.PARAM$mu
  sigma <- VEC.PARAM$sigma

  return(c(mu,sigma,p)
  )
}


grad.compute.gamma <- function(param,eps=.0000001){
  res <- NULL
  gamma.0 <- compute.gamma(param)
  for(i in 1:length(param)){
    param.i <- param
    param.i[i] <- param.i[i] + eps
    gamma.i <- compute.gamma(param.i)
    res <- cbind(res,(gamma.i-gamma.0)/eps)
  }
  return(res)
}


compute.phi.theta.C.distri <- function(param,n){

  p <- c((length(param) - 2*n^2 - 3*n)/n^2)

  Phi <- array(param[1:(n^2*p)],c(n,n,p))

  param.remain <- param[(n^2*p+1):length(param)]
  Theta.C <- make.Theta.C(param.remain,n)

  Theta <- Theta.C$Theta
  C     <- Theta.C$C

  theta1 <- Theta.C$param.4.distri[1:n]
  theta2 <- Theta.C$param.4.distri[(n+1):(2*n)]
  theta3 <- Theta.C$param.4.distri[(2*n+1):(3*n)]

  VEC.PARAM <-  make.p.mu.sigma(theta1,theta2,theta3)
  p     <- VEC.PARAM$p
  mu    <- VEC.PARAM$mu
  sigma <- VEC.PARAM$sigma

  return(c(Phi,Theta,C,mu,sigma,p)
  )
}


grad.compute.phi.theta.C.distri <- function(param,n,eps=.0000001){
  res <- NULL
  gamma.0 <- compute.phi.theta.C.distri(param,n)
  for(i in 1:length(param)){
    param.i <- param
    param.i[i] <- param.i[i] + eps
    gamma.i <- compute.phi.theta.C.distri(param.i,n)
    res <- cbind(res,(gamma.i-gamma.0)/eps)
  }
  return(res)
}








# ===============================================
# ===============================================
# SSVARMA 2SLS-GMM Estimation
# ===============================================
# ===============================================


compute.theoretical.moments <- function(C0,C1,kappa.3,kappa.4,u,v){
  # This function computes theoretical moments of the form:
  # E[(u'Z_t + v'Z_{t-1})^q], q = 2, 3 and 4,
  # where Z_t = C0 eta_t + C1 eta_{t-1}, with E(eta)=0 and Var(eta)=Id.
  #
  # C0 and C1 are matrices of dimension n x n.
  # u and v are of dimension k x n, where is the number of combinations u'Z_t + v'Z_{t-1} to consider
  # kappa3 is a n-dimensional vector containing the third-order cumulants
  # kappa4 is a n-dimensional vector containing the fourth-order cumulants

  n <- dim(C0)[1]

  k <- dim(u)[1]

  matrix.kappa.3 <- matrix(1,k,1) %*% matrix(kappa.3,nrow=1)
  matrix.kappa.4 <- matrix(1,k,1) %*% matrix(kappa.4,nrow=1)

  uC0 <- u %*% C0
  uC1 <- u %*% C1
  vC0 <- v %*% C0
  vC1 <- v %*% C1

  moments.2 <- apply( uC0^2 + (uC1 + vC0)^2 + vC1^2,1,sum)
  moments.3 <- apply((uC0^3 + (uC1 + vC0)^3 + vC1^3)*matrix.kappa.3,1,sum)
  moments.4 <- apply((uC0^4 + (uC1 + vC0)^4 + vC1^4)*matrix.kappa.4,1,sum) + 3*moments.2^2

  return(list(
    moments.2 = moments.2,
    moments.3 = moments.3,
    moments.4 = moments.4,
    vec.h.theo = matrix(c(moments.2,moments.3,moments.4),nrow=1)
  ))
}


compute.empirical.moments <- function(Z0,Z1,u,v){
  # This function computes theoretical moments of the form:
  # E[(u'Z_t + v'Z_{t-1})^q], q = 2, 3 and 4,
  # where Z_t = C0 eta_t + C1 eta_{t-1}, with E(eta)=0 and Var(eta)=Id.
  #
  # C0 and C1 are matrices of dimension n x n.
  # u and v are of dimension k x n, where is the number of combinations u'Z_t + v'Z_{t-1} to consider
  # kappa3 is a n-dimensional vector containing the third-order cumulants
  # kappa4 is a n-dimensional vector containing the fourth-order cumulants

  n <- dim(Z0)[2]

  k <- dim(u)[1]

  uZ0_vZ1 <- Z0 %*% t(u) + Z1 %*% t(v)

  moments.2 <- apply(uZ0_vZ1^2,2,mean)
  moments.3 <- apply(uZ0_vZ1^3,2,mean)
  moments.4 <- apply(uZ0_vZ1^4,2,mean)

  return(list(
    moments.2 = moments.2,
    moments.3 = moments.3,
    moments.4 = moments.4,
    matrix.h_t = cbind(uZ0_vZ1^2,uZ0_vZ1^3,uZ0_vZ1^4)
  ))
}


compute.distance.moments <- function(param,Z,u,v,Omega,
                                     indic.3rd.4th.order.moments=1){
  # If indic.3rd.4th.order.moments==1 (default) then:
  #       the estimation is carried out on third AND fourth order moments.
  # If indic.3rd.4th.order.moments==3 (default) then:
  #       the estimation is carried out on third-order moments only.
  # If indic.3rd.4th.order.moments==4 (default) then:
  #       the estimation is carried out on fourth-order moments only.

  n <- dim(Z)[2]

  C0 <- matrix(param[1:n^2],n,n)
  THETA <- matrix(param[(n^2+1):(2*n^2)],n,n)
  C1 <- - THETA %*% C0

  if(indic.3rd.4th.order.moments==1){
    kappa.3 <- param[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- 1 + kappa.3^2 + abs(param[(2*n^2+n+1):(2*n^2+2*n)]) - 3
  }else if(indic.3rd.4th.order.moments==3){
    kappa.3 <- param[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- pmin(0,1 + kappa.3^2 - 3)
  }else if(indic.3rd.4th.order.moments==4){
    kappa.4 <- 1 + abs(param[(2*n^2+1):(2*n^2+n)])
    kappa.3 <- 0*kappa.4
  }

  T <- dim(Z)[1] - 1
  Z0 <- Z[2:(T+1),]
  Z1 <- Z[1:T,]

  empirical.moments <- compute.empirical.moments(Z0,Z1,u,v)
  theoretical.moments <- compute.theoretical.moments(C0,C1,kappa.3,kappa.4,u,v)

  if(indic.3rd.4th.order.moments==1){
    empirical.moments <-
      c(empirical.moments$moments.2,empirical.moments$moments.3,empirical.moments$moments.4)
    theoretical.moments <-
      c(theoretical.moments$moments.2,theoretical.moments$moments.3,theoretical.moments$moments.4)
  }else if(indic.3rd.4th.order.moments==3){
    empirical.moments <-
      c(empirical.moments$moments.2,empirical.moments$moments.3)
    theoretical.moments <-
      c(theoretical.moments$moments.2,theoretical.moments$moments.3)
  }else if(indic.3rd.4th.order.moments==4){
    empirical.moments <-
      c(empirical.moments$moments.2,empirical.moments$moments.4)
    theoretical.moments <-
      c(theoretical.moments$moments.2,theoretical.moments$moments.4)
  }

  vec.distance <- empirical.moments - theoretical.moments

  distance <- matrix(vec.distance,nrow=1) %*% Omega %*% matrix(vec.distance,ncol=1)

  return(1000*distance)
}


compute.g.function <- function(Y,Model,indic.constant=1){
  # This function computes function g defined in Appenddic E.2 of the paper.
  # Works only for VARMA(1,1) models.
  # If indic.constant==1, the constant (in the VAR mode) is estimated
  n <- dim(Y)[2]
  T <- dim(Y)[1]
  p <- dim(Model$Phi)[3]

  Big.PHI <- matrix(Model$Phi,n,p*n)
  Big.PHI <- rbind(
    Big.PHI,
    cbind(diag(n*(p-1)),matrix(0,n*(p-1),n))
  )
  Big.Mu <- c(Model$Mu,rep(0,n*(p-1)))
  m <- (solve(diag(n*p) - Big.PHI) %*% Big.Mu)[1:n]

  # ======================================
  # Block corresponding to E(y - m)=0:
  first.block <- apply(Y,2,function(x){mean(x,na.rm=TRUE)}) - m

  # ======================================
  # Block corresponding to zero covariances:
  matrix.m <- matrix(m,ncol=1) %*% matrix(1,1,T)

  Z <- t(Y) - matrix.m
  Y_i <- Y
  for(i in 1:p){
    Y_i <- rbind(
      rep(NaN,n),
      Y_i[1:(T-1),]
    )
    Z <- Z - Model$Phi[,,i] %*% (t(Y_i) - matrix.m)
  }

  Y_1 <- rbind(rep(NaN,n),matrix(Y[1:(T-1),],T-1,n))
  Y_i <- Y_1
  Y_lags <- NULL
  for(i in 1:p){
    Y_i <- rbind(rep(NaN,n),matrix(Y_i[1:(T-1),],T-1,n))
    Y_lags <- cbind(Y_lags,Y_i)
  }
  # Y_lags is of dimension T x (p*n)
  # Its first column block is Y_2
  second.block <- 1/T * Z[,(p+1+1):T] %*% Y_lags[(p+1+1):T,]

  if(indic.constant==1){
    mean.g_t   <- c(first.block,second.block)
    matrix.g_t <- cbind(
      t(t(Y) - matrix.m)[(p+1+1):T,],
      (matrix(1,1,n*p) %x% t(Z)[(p+1+1):T,]) * (Y_lags[(p+1+1):T,] %x% matrix(1,1,n))
    )
  }else{
    mean.g_t   <- c(second.block)
    matrix.g_t <- (matrix(1,1,n*p) %x% t(Z)[(p+1+1):T,]) *
      (Y_lags[(p+1+1):T,] %x% matrix(1,1,n))
  }
  return(list(
    mean.g_t = mean.g_t,
    matrix.g_t = matrix.g_t
  ))
}


compute.h.function <- function(Y,Model,C0,C1,kappa.3,kappa.4,u,v,
                               indic.3rd.4th.order.moments=1){
  # This function computes function h defined in Appendix E.2 of the paper.
  # Works only for VARMA(1,1) models.
  # If indic.3rd.4th.order.moments==1 (default) then:
  #       the estimation is carried out on third AND fourth order moments.
  # If indic.3rd.4th.order.moments==3 (default) then:
  #       the estimation is carried out on third-order moments only.
  # If indic.3rd.4th.order.moments==4 (default) then:
  #       the estimation is carried out on fourth-order moments only.

  n <- dim(Y)[2]
  T <- dim(Y)[1]
  p <- dim(Model$Phi)[3]

  nb.combinations <- dim(u)[1]

  Big.PHI <- matrix(Model$Phi,n,p*n)
  Big.PHI <- rbind(
    Big.PHI,
    cbind(diag(n*(p-1)),matrix(0,n*(p-1),n))
  )
  Big.Mu <- c(Model$Mu,rep(0,n*(p-1)))
  m <- (solve(diag(n*p) - Big.PHI) %*% Big.Mu)[1:n]
  matrix.m <- matrix(m,ncol=1) %*% matrix(1,1,T)

  Z0 <- t(Y) - matrix.m
  Y_i <- Y
  for(i in 1:p){
    Y_i <- rbind(
      rep(NaN,n),
      matrix(Y_i[1:(T-1),],T-1,n)
    )
    Z0 <- Z0 - Model$Phi[,,i] %*% (t(Y_i) - matrix.m)
  }
  Z0 <- t(Z0)
  Z1 <- matrix(Z0[(p+1):(T-1),],T-(p+1),n)
  # remove NaNs:
  Z0 <- matrix(Z0[(p+1+1):T,],T-(p+1),n)

  empir <- compute.empirical.moments(Z0,Z1,u,v)
  theor <- compute.theoretical.moments(C0,C1,kappa.3,kappa.4,u,v)

  if(indic.3rd.4th.order.moments==1){
    mean.h_t <- c(
      empir$moments.2 - theor$moments.2,
      empir$moments.3 - theor$moments.3,
      empir$moments.4 - theor$moments.4
    )
    TT <- dim(empir$matrix.h_t)[1]
    matrix.h_t <- empir$matrix.h_t - matrix(1,TT,1) %*% matrix(theor$vec.h.theo,nrow=1)

  }else if(indic.3rd.4th.order.moments==3){
    mean.h_t <- c(
      empir$moments.2 - theor$moments.2,
      empir$moments.3 - theor$moments.3
    )
    TT <- dim(empir$matrix.h_t)[1]
    matrix.h_t <- empir$matrix.h_t[,1:(2*nb.combinations)] -
      matrix(1,TT,1) %*% matrix(theor$vec.h.theo[1:(2*nb.combinations)],nrow=1)

  }else if(indic.3rd.4th.order.moments==4){
    mean.h_t <- c(
      empir$moments.2 - theor$moments.2,
      empir$moments.4 - theor$moments.4
    )
    TT <- dim(empir$matrix.h_t)[1]
    indic.columns <- c(1:nb.combinations,(2*nb.combinations+1):(3*nb.combinations))
    matrix.h_t <- empir$matrix.h_t[,indic.columns] -
      matrix(1,TT,1) %*% matrix(theor$vec.h.theo[indic.columns],nrow=1)
  }

  matrix.h_t.complete <- matrix(NaN,T,dim(matrix.h_t)[2])
  matrix.h_t.complete[(T - TT + 1):T,] <- matrix.h_t

  return(
    list(
      mean.h_t = mean.h_t,
      matrix.h_t = matrix.h_t,
      matrix.h_t.complete = matrix.h_t.complete
    )
  )
}


estim.VARMAp1.2SLS.GMM <- function(Y,u,v,
                                   nb.iterations.gmm=6,
                                   maxitNM=2000,
                                   indic.Blaschke=1,
                                   indic.print=1,
                                   lag.NW=3,
                                   indic.estim.phi=1,
                                   indic.constant = 0,
                                   indic.3rd.4th.order.moments=1,
                                   param.ini=NULL,p=1,addit.IV=0){
  # This function estimates a VARMA(p,1) process for Y
  # Y is of dimension n x T.
  # If indic.estim.phi==1, a VMA(1) model is estimated.
  # If indic.3rd.4th.order.moments==1 (default) then:
  #       the estimation is carried out on third AND fourth order moments.
  # If indic.3rd.4th.order.moments==3 (default) then:
  #       the estimation is carried out on third-order moments only.
  # If indic.3rd.4th.order.moments==4 (default) then:
  #       the estimation is carried out on fourth-order moments only.
  # If param.ini is not NULL, then this vector of param is used as starting values of
  #       the numerical optim problem

  if((indic.estim.phi==0)&(p>0)){
    print("Cannot have (indic.estim.phi==0)&(p>0)")
    print("If you want (indic.estim.phi==0), set p=0")
    return(0)
  }

  n <- dim(Y)[2]
  T <- dim(Y)[1]

  # ===================================================
  # First step of the estimation (IV approach):
  if(indic.estim.phi==1){
    res.estim.IV <- TSLS.estim(Y,p,1,indic.constant = indic.constant,addit.IV = addit.IV)
    Phi.IV <- res.estim.IV$Phi.est
    Mu.IV <- res.estim.IV$Mu.est
    Big.PHI <- matrix(Phi.IV,n,p*n)
    Big.PHI <- rbind(
      Big.PHI,
      cbind(diag(n*(p-1)),matrix(0,n*(p-1),n))
    )
    Big.Mu <- c(Mu.IV,rep(0,n*(p-1)))
    m.IV <- (solve(diag(n*p) - Big.PHI) %*% Big.Mu)[1:n]
  }else{
    Phi.IV <- array(0,c(n,n,1))
    Mu.IV <- rep(0,n)
    m.IV <- rep(0,n)
  }

  Model.IV <- list(
    Mu = Mu.IV,
    Phi = array(Phi.IV,c(n,n,max(p,1))) # only for VARMA(1,1)
  )

  # ===================================================
  # Second step of the estimation (GMM approach):

  # (i) compute Z0 and Z1:
  Y_1 <- rbind(rep(NaN,n),matrix(Y[1:(T-1),],T-1,n))

  matrix.m <- matrix(m.IV,ncol=1) %*% matrix(1,1,T)

  Z0 <- t(Y) - matrix.m
  Y_i <- Y
  if(p>0){
    for(i in 1:p){
      Y_i <- rbind(
        rep(NaN,n),
        matrix(Y_i[1:(T-1),],T-1,n)
      )
      Z0 <- Z0 - Phi.IV[,,i] %*% (t(Y_i) - matrix.m)
    }
    Z0 <- t(Z0)
  }else{
    Z0 <- Y - t(matrix.m)
  }
  Z1 <- matrix(Z0[(p+1):(T-1),],T-(p+1),n)
  # remove NaNs:
  Z0 <- matrix(Z0[(p+1+1):T,],T-(p+1),n)

  # (ii) Estimate a first model based on GMM (matching sample Cov and 1st-order autoCov matrix):
  Cov <- var(Z0)
  Z_01 <- cbind(Z0,Z1)
  aux <- var(Z_01)
  AutoCov <- aux[(n+1):(2*n),1:n]
  res.gmm <- gmm.estim.sigma.theta(matrix(Cov,n,n),matrix(AutoCov,n,n),MAXIT.nlminb=50)
  Theta.gmm <- res.gmm$Theta

  Sigma.gmm <- res.gmm$Sigma
  C.gmm <- t(chol(Sigma.gmm))

  kappa.3.ini <- rep(0.1,n)
  kappa.4.ini <- rep(kappa.4.tstud,n)

  if(is.null(param.ini)){
    if(indic.3rd.4th.order.moments==1){
      param.beta <- c(C.gmm,Theta.gmm,kappa.3.ini,kappa.4.ini + 2 - kappa.3.ini^2)
    }else if(indic.3rd.4th.order.moments==3){
      param.beta <- c(C.gmm,Theta.gmm,kappa.3.ini)
    }else if(indic.3rd.4th.order.moments==4){
      param.beta <- c(C.gmm,Theta.gmm,kappa.4.ini + 2)
    }
  }else{
    param.beta <- param.ini
  }

  # (iii) Initialization of the weighting matrix:
  res.compute.asympt <- compute.asymptotic.distri.2SLS.GMM(Y,
                                                           matrix.g_t  = res.estim.IV$matrix.g_t,
                                                           Q.Pi        = res.estim.IV$Q.Pi,
                                                           param.alpha = res.estim.IV$param.alpha.est,
                                                           param.beta  = param.beta,
                                                           u,v,lag.NW,indic.estim.phi,
                                                           indic.3rd.4th.order.moments,p=p)
  Omega <- res.compute.asympt$Omega.star

  # Overweight order-2 moments to begin with:
  nb.combi <- dim(u)[1]
  Omega[1:nb.combi,1:nb.combi] <- 10*Omega[1:nb.combi,1:nb.combi]

  param.ini <- param.beta

  print("=========================================================")
  print("=========================================================")
  print(" First pass: 'Ad hoc' Omega ")
  print("=========================================================")
  print("=========================================================")

  ini.value <- compute.distance.moments(param=param.ini,Z0,u,v,Omega,
                                        indic.3rd.4th.order.moments)
  print(paste("--- Initial value of the criteria: ",toString(ini.value),sep=""))

  last.best <- 10^10


  # (iv) Optimization  with initial weighting matrix:
  print("--- Starting Minimization of criteria ---")
  for(i in 1:nb.iterations.gmm){
    res.optim <- optim(param.ini,compute.distance.moments,
                       Z=Z0,
                       u=u,
                       v=v,
                       Omega=Omega,
                       indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                       gr = NULL,
                       #method="Nelder-Mead",
                       method="BFGS",
                       control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
    if(res.optim$value<last.best){
      param.ini <- res.optim$par
      last.best <- res.optim$value
    }
    res.optim <- optim(param.ini,compute.distance.moments,
                       Z=Z0,
                       u=u,
                       v=v,
                       Omega=Omega,
                       indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="BFGS",
                       control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
    if(res.optim$value<last.best){
      param.ini <- res.optim$par
      last.best <- res.optim$value
    }
    print(paste("Criteria at the end of iteration ",toString(i),
                " (Max eval:",toString(maxitNM),")",": ",toString(last.best),sep=""))
  }
  print("--- End of Minimization of criteria ---")

  C0.est <- matrix(res.optim$par[1:(n^2)],n,n)
  Theta.est <- matrix(res.optim$par[(n^2+1):(2*n^2)],n,n)
  C1.est <- - Theta.est %*% C0.est
  if(indic.3rd.4th.order.moments==1){
    kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
    kappa.4.est <- 1 + kappa.3.est^2 + abs(res.optim$par[(2*n^2+n+1):(2*n^2+2*n)]) - 3
  }else if(indic.3rd.4th.order.moments==3){
    kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
    kappa.4.est <- pmax(0,1 + kappa.3.est^2 - 3)
  }else if(indic.3rd.4th.order.moments==4){
    kappa.4.est <- 1 + abs(res.optim$par[(2*n^2+1):(2*n^2+n)]) - 3
    kappa.3.est <- 0*kappa.4.est
  }

  print(Theta.est)

  if(indic.Blaschke==1){
    all.MA <- compute.all.forms(Theta.est,C0.est)
    all.MA$all.Theta <- Re(all.MA$all.Theta)
    all.MA$all.C <- Re(all.MA$all.C)
    nb.cases <- dim(all.MA$all.Theta)[3]
    for(i in 1:(nb.cases-1)){#the last element of all.MA (last layers) correspond to the original estimation
      if(indic.print==1){
        print("======================================================")
        print(paste("End of estimation ",toString(i)," out of ",toString(nb.cases), " (using Blaschke transf.)",sep=""))
        print("======================================================")
      }
      if(indic.3rd.4th.order.moments==1){
        param.ini.i <- c(all.MA$all.C[,,i],all.MA$all.Theta[,,i],
                         kappa.3.ini,kappa.4.ini + 2 - kappa.3.ini^2)
      }else if(indic.3rd.4th.order.moments==3){
        param.ini.i <- c(all.MA$all.C[,,i],all.MA$all.Theta[,,i],
                         kappa.3.ini)
      }else if(indic.3rd.4th.order.moments==4){
        param.ini.i <- c(all.MA$all.C[,,i],all.MA$all.Theta[,,i],
                         kappa.4.ini + 2 - kappa.3.ini^2)
      }
      ini.value <- compute.distance.moments(param=param.ini.i,Z0,u,v,Omega,
                                            indic.3rd.4th.order.moments)
      print(paste("--- Initial value of the criteria: ",toString(ini.value),sep=""))

      last.best <- 10^10

      print("--- Starting Minimization of criteria ---")
      for(i in 1:nb.iterations.gmm){
        res.optim.i <- optim(param.ini.i,compute.distance.moments,
                             Z=Z0,
                             u=u,
                             v=v,
                             Omega=Omega,
                             indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                             gr = NULL,
                             #method="Nelder-Mead",
                             method="BFGS",
                             control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
        if(res.optim.i$value<last.best){
          param.ini.i <- res.optim.i$par
          last.best   <- res.optim.i$value
        }
        res.optim.i <- optim(param.ini.i,compute.distance.moments,
                             Z=Z0,
                             u=u,
                             v=v,
                             Omega=Omega,
                             indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                             gr = NULL,
                             method="Nelder-Mead",
                             #method="BFGS",
                             control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
        if(res.optim.i$value<last.best){
          param.ini.i <- res.optim.i$par
          last.best   <- res.optim.i$value
        }
        print(paste("Criteria at the end of iteration ",toString(i)," (Max eval:",
                    toString(maxitNM),")",": ",toString(last.best),sep=""))
      }
      print("--- End of Minimization of criteria ---")

      if(res.optim.i$value<res.optim$value){
        print("")
        print(" ++++ Improvement of criteria ++++")
        print("")
        res.optim <- res.optim.i
        C0.est <- matrix(res.optim$par[1:(n^2)],n,n)
        Theta.est <- matrix(res.optim$par[(n^2+1):(2*n^2)],n,n)
        C1.est <- - Theta.est %*% C0.est
        if(indic.3rd.4th.order.moments==1){
          kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
          kappa.4.est <- 1 + kappa.3.est^2 + abs(res.optim$par[(2*n^2+n+1):(2*n^2+2*n)]) - 3
        }else if(indic.3rd.4th.order.moments==3){
          kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
          kappa.4.est <- pmin(0,1 + kappa.3.est^2 - 3)
        }else if(indic.3rd.4th.order.moments==4){
          kappa.4.est <- 1 + abs(res.optim$par[(2*n^2+1):(2*n^2+n)]) - 3
          kappa.3.est <- 0*kappa.4.est
        }
        print(Theta.est)
      }
    }
  }

  print("=========================================================")
  print("=========================================================")
  print(" Second pass: 'Optimal' Omega ")
  print("=========================================================")
  print("=========================================================")

  # Recompute Omega:
  if(indic.3rd.4th.order.moments==1){
    param.beta <- c(C0.est,Theta.est,kappa.3.est,kappa.4.est + 2 - kappa.3.est^2)
  }else if(indic.3rd.4th.order.moments==3){
    param.beta <- c(C0.est,Theta.est,kappa.3.est)
  }else if(indic.3rd.4th.order.moments==4){
    param.beta <- c(C0.est,Theta.est,kappa.4.est + 2)
  }
  res.compute.asympt <- compute.asymptotic.distri.2SLS.GMM(Y,
                                                           matrix.g_t  = res.estim.IV$matrix.g_t,
                                                           Q.Pi        = res.estim.IV$Q.Pi,
                                                           param.alpha = res.estim.IV$param.alpha.est,
                                                           param.beta  = param.beta,
                                                           u,v,lag.NW,indic.estim.phi,
                                                           indic.3rd.4th.order.moments,p=p)

  Omega <- res.compute.asympt$Omega.star

  param.ini <- res.optim$par

  last.best <- 10^10

  for(i in 1:nb.iterations.gmm){
    res.optim <- optim(param.ini,compute.distance.moments,
                       Z=Z0,
                       u=u,
                       v=v,
                       Omega=Omega,
                       indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                       gr = NULL,
                       #method="Nelder-Mead",
                       method="BFGS",
                       control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
    if(res.optim$value<last.best){
      param.ini <- res.optim$par
      last.best <- res.optim$value
    }
    res.optim <- optim(param.ini,compute.distance.moments,
                       Z=Z0,
                       u=u,
                       v=v,
                       Omega=Omega,
                       indic.3rd.4th.order.moments=indic.3rd.4th.order.moments,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="BFGS",
                       control=list(trace=FALSE,maxit=maxitNM),hessian=FALSE)
    if(res.optim$value<last.best){
      param.ini <- res.optim$par
      last.best <- res.optim$value
    }
    print(paste("Criteria at the end of iteration ",toString(i)," (Max eval:",toString(maxitNM),")",": ",toString(last.best),sep=""))
  }
  print("--- End of Minimization of criteria ---")

  C0.est <- matrix(res.optim$par[1:(n^2)],n,n)
  Theta.est <- matrix(res.optim$par[(n^2+1):(2*n^2)],n,n)
  C1.est <- - Theta.est %*% C0.est
  if(indic.3rd.4th.order.moments==1){
    kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
    kappa.4.est <- 1 + kappa.3.est^2 + abs(res.optim$par[(2*n^2+n+1):(2*n^2+2*n)]) - 3
  }else if(indic.3rd.4th.order.moments==3){
    kappa.3.est <- res.optim$par[(2*n^2+1):(2*n^2+n)]
    kappa.4.est <- pmax(0,1 + kappa.3.est^2 - 3)
  }else if(indic.3rd.4th.order.moments==4){
    kappa.4.est <- 1 + abs(res.optim$par[(2*n^2+1):(2*n^2+n)]) - 3
    kappa.3.est <- 0*kappa.4.est
  }

  if(indic.3rd.4th.order.moments==1){
    param.beta <- c(C0.est,Theta.est,kappa.3.est,kappa.4.est + 2 - kappa.3.est^2)
  }else if(indic.3rd.4th.order.moments==3){
    param.beta <- c(C0.est,Theta.est,kappa.3.est)
  }else if(indic.3rd.4th.order.moments==4){
    param.beta <- c(C0.est,Theta.est,kappa.4.est + 2)
  }
  Asympt <- compute.asymptotic.distri.2SLS.GMM(Y,
                                               matrix.g_t  = res.estim.IV$matrix.g_t,
                                               Q.Pi        = res.estim.IV$Q.Pi,
                                               param.alpha = res.estim.IV$param.alpha.est,
                                               param.beta  = param.beta,
                                               u,v,lag.NW,indic.estim.phi,
                                               indic.3rd.4th.order.moments,p=p)
  return(
    list(
      mu.est = Mu.IV,
      Phi.est = Phi.IV,
      m.est = m.IV,
      C0.est = C0.est,
      C1.est = C1.est,
      Theta.est = Theta.est,
      kappa.3.est = kappa.3.est,
      kappa.4.est = kappa.4.est,
      Asympt = Asympt,
      param.alpha = ifelse(indic.estim.phi==1,res.estim.IV$param.alpha.est,NaN),
      param.beta = param.beta
    )
  )
}


compute.asymptotic.distri.2SLS.GMM <- function(Y,
                                               matrix.g_t,
                                               Q.Pi,
                                               param.alpha,
                                               param.beta,u,v,lag.NW=6,
                                               indic.estim.phi=1,
                                               indic.3rd.4th.order.moments=1,p=1){
  # This procedure computes various matrices entering the asymptotic distribution of GMM estimates
  # It can also be used to compute the "optimal" Omega matrix.
  # If indic.constant==1, then a constant is included in the VAR specification.

  n <- dim(Y)[2]
  T <- dim(Y)[1]

  if(indic.estim.phi==1){
    if(length(param.alpha)==n^2*p){
      indic.constant <- 0
      Model <- list(
        Mu = rep(0,n),
        Phi = array(param.alpha,c(n,n,p))
      )
    }else{
      indic.constant <- 1
      Model <- list(
        Mu = param.alpha[1:n],
        Phi = array(param.alpha[(n+1):(n+p*n^2)],c(n,n,p))
      )
    }
  }else{
    # In this case, estimation of a VMA model, so no Phi.
    Model <- list(
      Mu = rep(0,n),
      Phi = array(0,c(n,n,1))
    )
  }

  C0 <- matrix(param.beta[1:n^2],n,n)
  Theta <- matrix(param.beta[(n^2+1):(2*n^2)],n,n)
  C1 <- -Theta %*% C0

  if(indic.3rd.4th.order.moments==1){
    kappa.3 <- param.beta[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- param.beta[(2*n^2+n+1):(2*n^2+2*n)]
  }else if(indic.3rd.4th.order.moments==3){
    kappa.3 <- param.beta[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- 0*kappa.3
  }else if(indic.3rd.4th.order.moments==4){
    kappa.4 <- param.beta[(2*n^2+1):(2*n^2+n)]
    kappa.3 <- 0*kappa.4
  }

  h.0 <- compute.h.function(Y,Model,C0,C1,kappa.3,kappa.4,u,v,
                            indic.3rd.4th.order.moments)


  # ===========================
  # Computation of E.dh.dalpha.prime:
  eps <- .00000001
  #I_g <- NULL
  E.dh.dalpha.prime <- NULL

  if(indic.estim.phi==1){
    for(i in 1:length(param.alpha)){
      param.alpha.i <- param.alpha
      param.alpha.i[i] <- param.alpha.i[i] + eps

      if(indic.constant==0){
        Model.i <- list(
          Mu = rep(0,n),
          Phi = array(param.alpha.i,c(n,n,p))
        )
      }else{
        Model.i <- list(
          Mu = param.alpha.i[1:n],
          Phi = array(param.alpha.i[(n+1):(n+p*n^2)],c(n,n,p))
        )
      }

      h.i <- compute.h.function(Y,Model.i,C0,C1,kappa.3,kappa.4,u,v,
                                indic.3rd.4th.order.moments)
      d.h.i <- matrix((h.i$mean.h_t - h.0$mean.h_t)/eps,ncol=1)

      E.dh.dalpha.prime <- cbind(E.dh.dalpha.prime,d.h.i)
    }
  }else{
    # VMA case: no alpha estimation.
    E.dh.dalpha.prime <- NULL
  }


  # ===========================
  # Computation of Q and Q.0:
  dim.h <- length(h.0$mean.h_t)

  if(indic.estim.phi==1){
    matrix.g_t.h_t <- cbind(matrix.g_t,h.0$matrix.h_t.complete)

    # remove lines with NaNs:
    aux <- apply(matrix.g_t.h_t,1,sum)
    indic.not.NaN <- which(!is.na(aux))
    matrix.g_t.h_t <- matrix.g_t.h_t[indic.not.NaN,]

    Q <- NW.LongRunVariance(matrix.g_t.h_t,lag.NW)

    aux <- cbind(
      E.dh.dalpha.prime %*% Q.Pi, diag(dim.h)
    )
    Q.0 <- aux %*% Q %*% t(aux)

  }else{
    # no alpha estimation
    Q <- NW.LongRunVariance(h.0$matrix.h_t,lag.NW)
    Q.0 <- Q
  }


  # ===========================
  # Computation of E.dh.dbeta.prime:

  # Baseline value:
  if(indic.3rd.4th.order.moments==1){
    kappa.3 <- param.beta[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- 1 + kappa.3^2 + abs(param.beta[(2*n^2+n+1):(2*n^2+2*n)]) - 3
  }else if(indic.3rd.4th.order.moments==3){
    kappa.3 <- param.beta[(2*n^2+1):(2*n^2+n)]
    kappa.4 <- pmax(0,1 + kappa.3^2 - 3)
  }else if(indic.3rd.4th.order.moments==4){
    kappa.4 <- 1 + abs(param.beta[(2*n^2+1):(2*n^2+n)]) - 3
    kappa.3 <- 0*kappa.4
  }
  h.0 <- compute.h.function(Y,Model,C0,C1,kappa.3,kappa.4,u,v,
                            indic.3rd.4th.order.moments)

  E.dh.dbeta.prime <- NULL

  for(i in 1:length(param.beta)){
    param.beta.i <- param.beta
    param.beta.i[i] <- param.beta.i[i] + eps

    C0.i <- matrix(param.beta.i[1:n^2],n,n)
    Theta.i <- matrix(param.beta.i[(n^2+1):(2*n^2)],n,n)
    C1.i <- - Theta.i %*% C0.i

    if(indic.3rd.4th.order.moments==1){
      kappa.3.i <- param.beta.i[(2*n^2+1):(2*n^2+n)]
      kappa.4.i <- 1 + kappa.3.i^2 + abs(param.beta.i[(2*n^2+n+1):(2*n^2+2*n)]) - 3
    }else if(indic.3rd.4th.order.moments==3){
      kappa.3.i <- param.beta.i[(2*n^2+1):(2*n^2+n)]
      kappa.4.i <- pmax(0,1 + kappa.3.i^2 - 3)
    }else if(indic.3rd.4th.order.moments==4){
      kappa.4.i <- 1 + abs(param.beta.i[(2*n^2+1):(2*n^2+n)]) - 3
      kappa.3.i <- 0*kappa.4.i
    }

    h.i <- compute.h.function(Y,Model,C0.i,C1.i,kappa.3.i,kappa.4.i,u,v,
                              indic.3rd.4th.order.moments)
    d.h.i <- matrix((h.i$mean.h_t - h.0$mean.h_t)/eps,ncol=1)
    E.dh.dbeta.prime <- cbind(E.dh.dbeta.prime,d.h.i)
  }


  # ===========================
  # Computation of Var(alpha.hat):
  if(indic.estim.phi==1){
    Var.alpha.hat <- 1/T * Q.Pi %*% Q[1:dim(Q.Pi)[2],1:dim(Q.Pi)[2]] %*% t(Q.Pi)
  }else{
    Var.alpha.hat <- NULL
  }

  eig.Q.0 <- eigen(Q.0)
  Q.0_1 <- eig.Q.0$vectors %*% diag(1/eig.Q.0$values) %*% t(eig.Q.0$vectors)
  Omega.star <- Q.0_1

  AUX <- t(E.dh.dbeta.prime) %*% Omega.star %*% E.dh.dbeta.prime

  eig.AUX <- eigen(AUX)
  AUX_1 <- eig.AUX$vectors %*% diag(1/eig.AUX$values) %*% t(eig.AUX$vectors)
  Var.beta.hat <- 1/T * AUX_1

  return(list(
    Q = Q,
    Q.0 = Q.0,
    Var.alpha.hat = Var.alpha.hat,
    Var.beta.hat = Var.beta.hat,
    Omega.star = Omega.star,
    E.dh.dalpha.prime = E.dh.dalpha.prime
  ))
}


compute.kappa.3.4 <- function(p,mu,sigma){
  # compute kappa.3 and kappa.4 for a mixture of Gaussian with E(X)=0 and Var(X)=1

  mu.1 <- mu
  sigma.1 <- sigma

  mu.2 <- - p/(1-p)*mu.1
  sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )

  skewness <- p*(mu.1^3 + 3*mu.1*sigma.1^2) + (1-p)*(mu.2^3 + 3*mu.2*sigma.2^2)
  kurtosis <- p*(mu.1^4 + 6*mu.1^2*sigma.1^2 + 3*sigma.1^4) +
    (1-p)*(mu.2^4 + 6*mu.2^2*sigma.2^2 + 3*sigma.2^4)

  return(list(
    kappa.3 = skewness,
    kappa.4 = kurtosis - 3
  ))
}


gmm.estim.sigma.theta <- function(Cov,AutoCov,MAXIT.nlminb){
  # Cov is a sample cvariance matrix
  # AutoCov is the sample first-order covariance matrix
  n <- dim(Cov)[1]

  if(n>1){
    param.ini <- c(rep(.1,n^2),lower.tri(Cov),diag(Cov))

    maxlik <- optimx(param.ini,fff,
                     Cov=Cov,
                     AutoCov=AutoCov,
                     n=n,
                     method = "nlminb",
                     control=list(trace=0, maxit=MAXIT.nlminb, kkt = FALSE))
    param.est <- c(as.matrix(maxlik)[1:length(param.ini)])

    Theta <- matrix(param.est[1:n^2],n,n)
    Sigma <- matrix(0,n,n)
    Sigma[lower.tri(Sigma)] <- param.est[(n^2+1):(n^2+n*(n-1)/2)]
    diag(Sigma) <- abs(param.est[(n^2+n*(n-1)/2+1):(n^2+n*(n-1)/2+n)])
    Sigma <- Sigma %*% t(Sigma)
    Cov.est <- Sigma + Theta %*% Sigma %*% t(Theta)
    AutoCov.est <- - Theta %*% Sigma
  }else{
    Delta <- (Cov/AutoCov)^2 + 4
    solutions <- .5*c(-Cov/AutoCov - sqrt(Delta),-Cov/AutoCov + sqrt(Delta))
    Theta <- solutions[2]
    Sigma <- Cov/(1 + Theta^2)
    param.est <- NaN
    AutoCov.est <- NaN
    Cov.est <- NaN
  }

  return(list(param.est=param.est,
              Cov.est=Cov.est,
              AutoCov.est=AutoCov.est,
              Sigma.est=Sigma,
              Theta.est=Theta))
}


fff <- function(param,Cov,AutoCov,n){
  Theta <- matrix(param[1:n^2],n,n)
  Sigma <- matrix(0,n,n)
  Sigma[lower.tri(Sigma)] <- param[(n^2+1):(n^2+n*(n-1)/2)]
  Sigma <- Sigma + t(Sigma)
  diag(Sigma) <- abs(param[(n^2+n*(n-1)/2+1):(n^2+n*(n-1)/2+n)])
  Cov.est <- Sigma + Theta %*% Sigma %*% t(Theta)
  AutoCov.est <- - Theta %*% Sigma
  dev.Cov <- c(Cov - Cov.est)
  dev.AutoCov <- c(AutoCov - AutoCov.est)
  return(sum(dev.Cov^2+dev.AutoCov^2))
}


TSLS.estim <- function(Y,p,q,indic.constant=1,addit.IV=0,lag.NW=6){
  # IV estimation of the autoregressive part of a VARMA model.
  # if indic.constant = 1, then uses constant in the VAR model
  # If addit.IV=0 then the parameters are "just identified". Otherwise they are
  # over-identified, estinated by two-stage least squares.

  n <- dim(Y)[2]
  T <- dim(Y)[1]

  # Construcition of W and X:
  Y_i <- Y # k-lagged Y

  if(indic.constant==1){
    W <- matrix(1,T,1)
    X <- matrix(1,T,1)
  }else{
    W <- NULL
    X <- NULL
  }
  for(i in 1:(p+q+addit.IV)){
    Y_i <- matrix(rbind(matrix(NaN,1,n),matrix(Y_i[1:(T-1),],T-1,n)),T,n)
    if(i<=p){
      X <- cbind(X,Y_i)
    }
    if(i>q){
      W <- cbind(W,Y_i)
    }
  }
  aux <- apply(cbind(X,W),1,sum)
  indic.not.NaN <- which(!is.na(aux))
  W.no.NaN <- matrix(W[indic.not.NaN,],ncol=dim(W)[2])
  X.no.NaN <- matrix(X[indic.not.NaN,],ncol=dim(X)[2])
  Y.no.NaN <- matrix(Y[indic.not.NaN,],ncol=dim(Y)[2])

  X.hat.no.NaN <- W.no.NaN %*% solve(t(W.no.NaN)%*%W.no.NaN) %*% t(W.no.NaN) %*% X.no.NaN
  Pi.hat <- solve(t(X.hat.no.NaN)%*%X.hat.no.NaN) %*% t(X.hat.no.NaN) %*% Y.no.NaN

  Z.hat.no.NaN <- Y.no.NaN - X.hat.no.NaN %*% Pi.hat

  T.corrected <- dim(W.no.NaN)[1]
  Q.XW <- 1/T.corrected * t(X.no.NaN) %*% W.no.NaN
  Q.W  <- 1/T.corrected * t(W.no.NaN) %*% W.no.NaN

  Q.Pi <- (solve( Q.XW %*% solve(Q.W) %*% t(Q.XW) ) %*% Q.XW %*% solve(Q.W)) %x% diag(n)

  nb.col.Z <- dim(Z.hat.no.NaN)[2]
  nb.col.W <- dim(W.no.NaN)[2]

  AUX <- (matrix(1,1,nb.col.W) %x% Z.hat.no.NaN) * (W.no.NaN %x% matrix(1,1,nb.col.Z))
  matrix.g_t <- matrix(NaN,T,dim(AUX)[2])
  matrix.g_t[indic.not.NaN,] <- AUX

  param.alpha.est <- NULL
  if(indic.constant==1){
    Mu.est <- matrix(Pi.hat[1,],ncol=1)
    param.alpha.est <- c(Mu.est)
    Phi.est <- array(NaN,c(n,n,p))
    for(i in 1:p){
      Phi.est[,,i] <- t(Pi.hat[(1+n*(i-1)+1):(1+n*i),])
      param.alpha.est <- c(param.alpha.est,Phi.est[,,i])
    }
  }else{
    Mu.est <- matrix(0,n,1)
    Phi.est <- array(NaN,c(n,n,p))
    for(i in 1:p){
      Phi.est[,,i] <- t(Pi.hat[(n*(i-1)+1):(n*i),])
      param.alpha.est <- c(param.alpha.est,Phi.est[,,i])
    }
  }

  # Compute estimator of covariance matrix of alpha = vec(Pi.hat):
  LRV.g_t <- NW.LongRunVariance(AUX,q=lag.NW)
  Mat.covar <- 1/T* Q.Pi %*% LRV.g_t %*% t(Q.Pi)

  return(list(Pi.hat=Pi.hat,
              Mu.est  = Mu.est,
              Phi.est = Phi.est,
              Z.hat.no.NaN=Z.hat.no.NaN,
              W.no.NaN = W.no.NaN,
              Q.XW = Q.XW,
              Q.W = Q.W,
              Q.Pi = Q.Pi,
              matrix.g_t = matrix.g_t,
              indic.not.NaN = indic.not.NaN,
              Mat.covar = Mat.covar,
              param.alpha.est = param.alpha.est))
}









# ===============================================
# ===============================================
# Skewness, Kurtosis, and Mixture of Gaussian ditributions
# ===============================================
# ===============================================




density.mixt.gauss <- function(param,X){
  param.0 <- make.p.mu.sigma(param[1],param[2],param[3])
  p <- param.0$p
  mu.1 <- param.0$mu
  sigma.1 <- param.0$sigma

  mu.2 <- - p/(1-p)*mu.1
  sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )

  f <- p*1/sqrt(2*pi*sigma.1^2)*exp(-.5*(X-mu.1)^2/sigma.1^2) +
    (1-p)*1/sqrt(2*pi*sigma.2^2)*exp(-.5*(X-mu.2)^2/sigma.2^2)
  return(f)
}


loglik.mixt.Gauss <- function(param,X){
  f <- density.mixt.gauss(param,X)
  res <- -sum(log(f),na.rm=TRUE)  + 0*sum((abs(param)-10)^2*(abs(param) > 10))
  return(res)
}

estim.MLE.mixtures <- function(X,
                               maxit.NM=10000,
                               nb.loops=4){
  # This function estimates the three paramters of a mixture of Gaussian distri. by MLE.
  # X is a vector of dimension Tx1.

  # In a first step, X is demeaned and standardized.
  X.stand <- (X - mean(X,na.rm=TRUE))/sd(X, na.rm=TRUE)

  param.ini <- c(0,0,0)
  for(i in 1:nb.loops){
    res.optim <- optim(param.ini,loglik.mixt.Gauss,
                       X = X.stand,
                       gr = NULL,
                       method="Nelder-Mead",
                       #method="CG",
                       #method="BFGS",
                       control=list(trace=FALSE,maxit=maxit.NM),hessian=FALSE)
    param.ini <- res.optim$par
  }

  param.est <- param.ini

  param.0 <- make.p.mu.sigma(param.est[1],param.est[2],param.est[3])
  p <- param.0$p
  mu.1 <- param.0$mu
  sigma.1 <- param.0$sigma

  return(list(
    param.est = param.est,
    p.est = p,
    mu.est = mu.1,
    max.log.lik = -res.optim$value,
    sigma.est = sigma.1
  ))
}


get.Gaussian.mixture.from.kappa <- function(kappa.3,kappa.4,p,indic.plot.distri=TRUE){
  # Find parameterization that is consistent with kappa.3 and kappa.4:

  mu.1 <- seq(-30,30,by=.0001)

  kappa.4.vec <- kurtosis(mu.1,p,kappa.3) - 3
  kappa.4.vec[is.na(kappa.4.vec)] <- 10000

  indic.sol <- which((kappa.4.vec-kappa.4)^2==min((kappa.4.vec-kappa.4)^2,na.rm=TRUE))
  mu.1.est <- mu.1[indic.sol[1]]
  sigma.12 <- (1-p)/(3*p*mu.1.est)*(kappa.3 - mu.1.est^3*(p - p^3/(1-p)^2 + 3*p^2/(1-p)^2) +
                                      3*mu.1.est*p/(1-p))
  sigma.22 <- 1/(1-p)*(1 - p*sigma.12 - p/(1-p)*mu.1.est^2)
  mu.2 <- -p*mu.1.est/(1 - p)
  # Compute standard deviations:
  sigma.1 <- sqrt(sigma.12)
  sigma.2 <- sqrt(sigma.22)

  # Check results:
  step <- .001
  xx <- seq(-30,30,by=step)
  f.0 <- p * dnorm((xx-mu.1.est)/sigma.1)/sigma.1 + (1 - p) * dnorm((xx-mu.2)/sigma.2)/sigma.2
  if(indic.plot.distri){
    par(plt=c(.1,.9,.1,.9))
    plot(xx,f.0,type="l")
  }
  print(paste("Variance of distri.:",toString(sum(xx^2*f.0)*step),sep=""))
  print(paste("kappa_3 of distri.:",toString(sum(xx^3*f.0)*step),sep=""))
  print(paste("kappa_4 of distri.:",toString(sum(xx^4*f.0)*step-3),sep=""))
  #print(sum(xx^5*f.0)*step - 10*sum(xx^3*f.0)*step*sum(xx^2*f.0)*step)

  return(list(
    mu.1 = mu.1.est,
    mu.2 = mu.2,
    sigma.1 = sigma.1,
    sigma.2 = sigma.2
  ))
}



compute.mixture <- function(p,kappa.3,kappa.4,max.mu = 30,indic.print.distri = 1){
  # For a given p, this function computes mu.1 and sigma.1
  # that are such that the associated Gaussian mixture
  # p,N(mu.1,sigma.1^2),N(mu.2,sigma.2^2), of mean zero and unit variance
  # has as 3rd and 4th cumulants kappa.3 and kappa.4

  mu.1 <- seq(-max.mu,max.mu,by=.001)
  kappa.4.vec <- kurtosis(mu.1,p,kappa.3) - 3
  kappa.4.vec[is.na(kappa.4.vec)] <- 10000
  indic.sol <- which((kappa.4.vec-kappa.4)^2==min((kappa.4.vec-kappa.4)^2,na.rm=TRUE))
  mu.1 <- mu.1[indic.sol[1]]

  resulting.kappa.4 <- kurtosis(mu.1,p,kappa.3) - 3
  print(paste("Resulting kappa.4: ",toString(round(resulting.kappa.4,3)),sep=""))

  sigma.12 <- (1-p)/(3*p*mu.1)*(kappa.3 - mu.1^3*(p - p^3/(1-p)^2 + 3*p^2/(1-p)^2) +
                                  3*mu.1*p/(1-p))
  sigma.22 <- 1/(1-p)*(1 - p*sigma.12 - p/(1-p)*mu.1^2)
  mu.2 <- -p*mu.1/(1 - p)

  sigma.1 <- sqrt(sigma.12)
  sigma.2 <- sqrt(sigma.22)

  if(indic.print.distri==1){
    step <- .01
    xx <- seq(-6,6,by=step)
    f.0 <- p * dnorm((xx-mu.1)/sigma.1)/sigma.1 + (1 - p) * dnorm((xx-mu.2)/sigma.2)/sigma.2
    plot(xx,f.0,type="l")
  }

  print(paste("Check variance: ",toString(round(sum(xx^2*f.0)*step,3)),sep=""))
  print(paste("Check skewness: ",toString(round(sum(xx^3*f.0)*step,3)),sep=""))
  print(paste("Check kurtosis: ",toString(round(sum(xx^4*f.0)*step,3)),sep=""))

  return(
    list(
      mu.1 = mu.1,
      mu.2 = mu.2,
      sigma.1 = sigma.1,
      sigma.2 = sigma.2
    )
  )
}


Kurtosis <- function(Y){
  X <- Y[!is.na(Y)]
  MN <- mean(X)
  SD <- sd(X)
  if(class(X)=="numeric"){
    Y <- t(X-MN)/SD
  } else {
    Y <- (X-MN)/SD
  }
  k = 1/max(dim(Y))*sum(Y^4)
  return(k)
}

skewness <- function(Y){
  X <- Y[!is.na(Y)]
  MN <- mean(X)
  SD <- sd(X)
  if(class(X)=="numeric"){
    Y <- t(X-MN)/SD
  } else {
    Y <- (X-MN)/SD
  }
  k = 1/max(dim(Y))*sum(Y^3)
  return(k)
}


qmixt <- function(alpha,mu,sigma,p){
  # compute the quantiles (alpha) of mixtures of Gaussian
  # mu has two elements
  # sigma has two elements
  if(length(mu)<2){
    mu.1 <- mu
    sigma.1 <- sigma
    mu.2 <- - p/(1-p)*mu.1
    if((1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2))<0){
      sigma.2 <- 100000
    }else{
      sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
    }
  }else{
    mu.1 <- mu[1]
    mu.2 <- mu[2]
    sigma.1 <- sigma[1]
    sigma.2 <- sigma[2]
  }

  # First: apply a grid to find a good starting points.
  x.0 <- c(
    seq(mu.1 - 3*sigma.1,mu.1 + 3*sigma.1,by=(6*sigma.1)/20),
    seq(mu.2 - 3*sigma.2,mu.2 + 3*sigma.2,by=(6*sigma.2)/20)
  )
  f.0 <- p * pnorm((x.0-mu.1)/sigma.1) + (1 - p) * pnorm((x.0-mu.2)/sigma.2)

  x <- qnorm(alpha)
  matrix.alpha.f0 <- matrix(alpha,ncol=1) %*% matrix(1,1,length(f.0)) -
    matrix(1,length(alpha),ncol=1) %*% matrix(f.0,1,length(f.0))
  indic.min <- apply(matrix.alpha.f0^2,1,function(x){which(x==min(x))[1]})
  x <- x.0[indic.min]

  # compute min and max possible values (10 stdev)
  min.values <- pmin(mu.1 - 5*sigma.1,mu.2 - 5*sigma.2)
  max.values <- pmax(mu.1 + 5*sigma.1,mu.2 + 5*sigma.2)

  for(i in 1:5){
    phi.star <- p * pnorm((x-mu.1)/sigma.1) + (1 - p) * pnorm((x-mu.2)/sigma.2)
    f <- log(phi.star) - log(alpha)
    df <- (p/sigma.1 * dnorm((x-mu.1)/sigma.1) + (1 - p)/sigma.2 * dnorm((x-mu.2)/sigma.2))/phi.star
    x <- x - .5*f/df
    x <- pmin(pmax(x,min.values),max.values)
  }
  return(x)
}



make.p.mu.sigma <- function(theta1,theta2,theta3){
  p <- .001 + .998*.5*(1+atan(theta1)*2/pi)
  sigma <- sqrt(1/p*.5*(1+atan(theta2)*2/pi))
  mu <- atan(theta3)*2/pi*sqrt((1-p*sigma^2)*(1-p)/p)

  return(list(p=p,
              mu=mu,
              sigma=sigma))
}


make.theta1.2.3 <- function(p,mu,sigma){
  theta1 <- tan(((p-.001)/(.998*.5) - 1)*pi/2)
  theta2 <- tan((2*p*sigma^2 - 1)*pi/2)
  theta3 <- tan(
    mu*pi/2/sqrt((1-p*sigma^2)*(1-p)/p)
  )
  return(list(
    theta1 = theta1,
    theta2 = theta2,
    theta3 = theta3
  ))
}


kurtosis <- function(mu.1,p,S){
  # p and mu.1 are scalar
  # not appropriate for mu.1 = 0.
  s <- S
  s[s>2*(1-p)/sqrt(p*(2-p))] <- NaN
  s[s< -2*(1-p)/sqrt(p*(2-p))] <- NaN

  sigma.12 <- (1-p)/(3*p*mu.1)*(s - mu.1^3*(p - p^3/(1-p)^2 + 3*p^2/(1-p)^2) +
                                  3*mu.1*p/(1-p))
  indic.zero.skew <- which(S==0)
  sigma.12[indic.zero.skew] <- 1 - mu.1^2/3

  sigma.22 <- 1/(1-p)*(1 - p*sigma.12 - p/(1-p)*mu.1^2)

  mu.2 <- -p*mu.1/(1 - p)

  kurt <- p*(mu.1^4 + 6*mu.1^2*sigma.12 + 3*sigma.12^2) +
    (1-p)*(mu.2^4 + 6*mu.2^2*sigma.22 + 3*sigma.22^2)

  kurt[sigma.12 <0] <- NaN
  kurt[sigma.22 <0] <- NaN

  return(kurt)
}


find.mu.sigma.from.kappas <- function(kappa.3,kappa.4,p,max.mu.1=3){

  mu.1 <- seq(-max.mu.1,max.mu.1,by=.001)
  kappa.4.vec <- kurtosis(mu.1,p,kappa.3) - 3
  indic.sol <- which((kappa.4.vec-kappa.4)^2==min((kappa.4.vec-kappa.4)^2,na.rm=TRUE))
  mu.1.est <- mu.1[indic.sol[1]]

  sigma.12 <- (1-p)/(3*p*mu.1.est)*(kappa.3 - mu.1.est^3*(p - p^3/(1-p)^2 + 3*p^2/(1-p)^2) +
                                      3*mu.1.est*p/(1-p))
  sigma.22 <- 1/(1-p)*(1 - p*sigma.12 - p/(1-p)*mu.1.est^2)
  mu.2 <- -p*mu.1.est/(1 - p)

  sigma.1 <- sqrt(sigma.12)
  sigma.2 <- sqrt(sigma.22)

  return(list(mu.1 = mu.1.est,
              mu.2 = mu.2,
              sigma.1 = sigma.1,
              sigma.2 = sigma.2))
}





compute.skew.kurt.mixt.gaussian <- function(mu,sigma,p,step=.0025,x.min=-8,x.max=8){
  x <- seq(x.min,x.max,by=step)

  mu.1 <- mu
  sigma.1 <- sigma
  mu.2 <- - p/(1-p)*mu.1
  if((1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2))<0){
    sigma.2 <- 100000
  }else{
    sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
  }
  f <-p/sigma.1 * dnorm((x-mu.1)/sigma.1) + (1 - p)/sigma.2 * dnorm((x-mu.2)/sigma.2)

  Skew <- sum(x^3*f)*step
  Kurt <- sum(x^4*f)*step
  return(list(skewness=Skew,kurtosis=Kurt,f=f,x=x))
}


compute.dev.skew.kurt <- function(theta,Skew,Kurt){
  model.skew.kurt <- compute.skew.kurt.mixt.gaussian(mu=theta[1],sigma=abs(theta[2]),p=.5)
  return(
    (model.skew.kurt$skewness - Skew)^2 + (model.skew.kurt$kurtosis - Kurt)^2
  )
}


compute.dev.quantiles <- function(theta,vec.q,vec.p){
  model.q <- qmixt(vec.p,mu=theta[1],sigma=abs(theta[2]),p=.5)
  return(
    sum((model.q - vec.q)^2)
  )
}


fit.skew.kurt <- function(Skew,Kurt){
  theta.ini <- c(0,1)
  res.optim <- optim(theta.ini,compute.dev.skew.kurt,
                     Skew=Skew,Kurt=Kurt,
                     gr = NULL,
                     method="Nelder-Mead",
                     #method="CG",
                     #method="BFGS",
                     control=list(trace=FALSE,maxit=100),hessian=TRUE)
  res.optim$par[2] <- abs(res.optim$par[2])
  print(c(Skew,Kurt))
  res.aux <- compute.skew.kurt.mixt.gaussian(mu=res.optim$par[1],sigma=res.optim$par[2],p)
  print(c(res.aux$skewness,res.aux$kurtosis))
  return(res.optim$par)
}


fit.q <- function(vec.q,vec.p){
  theta.ini <- c(0,1)
  res.optim <- optim(theta.ini,compute.dev.quantiles,
                     vec.q=vec.q,vec.p=vec.p,
                     gr = NULL,
                     method="Nelder-Mead",
                     #method="CG",
                     #method="BFGS",
                     control=list(trace=FALSE,maxit=100),hessian=TRUE)
  res.optim$par[2] <- abs(res.optim$par[2])
  print(vec.q)
  res.aux <- qmixt(vec.p,mu=res.optim$par[1],sigma=abs(res.optim$par[2]),p=.5)
  print(c(res.aux))
  return(res.optim$par)
}







# ===============================================
# ===============================================
# Various SSVARMA procedures
# ===============================================
# ===============================================


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
  PHI <- make.PHI2(Model$Phi)
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

compute.V <- function(Model,Y){
  # Y is of dimension n x T
  n <- dim(Model$Phi)[1]
  p <- dim(Model$Phi)[3]
  q <- dim(Model$Theta)[3]
  T <- dim(Y)[2]

  MU <- c(Model$Mu,rep(0,n*(p-1)))
  PHI <- make.PHI2(Model$Phi)

  Y.augm <- NULL
  for(k in 1:p){
    Y.augm <- rbind(Y.augm,
                    cbind(matrix(NaN,n,k),matrix(Y[,1:(T-k)],n,T-k)))
  }

  Y.pred <- matrix(MU,n*p,T) + PHI %*% Y.augm
  V <- Y - Y.pred[1:n,]
  return(V)
}



make.Theta.C <- function(param,n){

  param.4.distri <- NULL # vector of parameters specifying the distribution of eta's. Remains NULL if indic.estim.mixtures=0

  if(length(param)==2*n^2){# Estimation of Theta and C only.
    Theta <- matrix(param[1:n^2],n,n)
    #Theta <- Model$Theta[,,1]
    size <- n^2
    C <- matrix(param[(size+1):(size+n^2)],n,n)
    aux <- matrix(0,n,n)
    diag(aux) <- sign(diag(C))
    if(n==1){
      C <- abs(C)
    }
  }else{ # Estimation od Theta, C + parameters governing the distri of eta
    # Specifically, we estimate Theta, C and gamma, where gamma are the parameters defining the distri. of eta's
    nb.param.per.distri <- 3 # 2 for Gaussian mixtures with E(eta)=0 and Var(eta)=1 and p=.5
    Theta <- matrix(param[1:n^2],n,n)
    size <- n^2
    C <- matrix(param[(size+1):(size+n^2)],n,n)
    aux <- matrix(0,n,n)
    diag(aux) <- sign(diag(C))
    #C <- C %*% aux
    if(n==1){
      C <- abs(C)
    }
    size <- size + n^2
    param.4.distri <- param[(size+1):(size+n*nb.param.per.distri)]
  }

  return(list(Theta=Theta,C=C,param.4.distri=param.4.distri))
}




estim.struct.shocks <- function(Y,Model){
  # estimate structural shocks, taking into account the fact that
  # Model may be non-fundamental.
  # Implemented for MA(1) only.
  # Y is of dimension T x n

  n <- dim(Y)[2]
  T <- dim(Y)[1]
  p <- dim(Model$Phi)[3]
  if(is.matrix(Model$Theta)){
    q <- 1
  }else{
    q <- dim(Model$Theta)[3]
  }

  if(q>1){
    print("Not implemented for q>1 yet")
    return(0)
  }

  Y.aug <- Y
  Y.aug_1 <- rbind(rep(NaN,n),matrix(Y[1:(T-1),],ncol=n))
  Y_k <- Y.aug_1
  if(p>1){
    for(i in 1:(p-1)){
      Y.aug <- cbind(Y.aug,Y_k)
      Y_k <- rbind(rep(NaN,n),Y_k[1:(T-1),])
      Y.aug_1 <- cbind(Y.aug_1,Y_k)
    }
  }

  vec.1 <- matrix(1,T,1)

  # First step: compute residuals:
  PHI <- make.PHI2(Model$Phi)
  V <- Y.aug - vec.1 %*% matrix(c(Model$Mu,rep(0,(p-1)*n)),nrow=1) - Y.aug_1 %*% t(PHI)

  # Second step: Schur decomposition:
  if(is.matrix(Model$Theta)){
    res.schur <- Schur(Model$Theta)
  }else{
    res.schur <- Schur(as.matrix(Model$Theta[,,1]))
  }
  U <- res.schur$T

  # Second step: diagonalize Theta:
  A <- res.schur$Q
  A_1 <- t(A)

  V.star <- V[,1:n] %*% t(A_1)
  V.star[is.na(V.star)] <- 0

  EPS.star.est <- NULL
  count <- n # will count number of variables treated
  while(count > 0){
    # detect if 1 (real eigenvalue) or 2 (complex eigenvalue) variables:
    if(count>1){
      if(U[count,count-1]==0){
        # In this case, real eigenvalue
        indices.columns <- count
        count <- count - 1
      }else{
        indices.columns <- c(count-1,count)
        count <- count - 2
      }
    }else{
      # In this case, real eigenvalue
      indices.columns <- count
      count <- count - 1
    }

    lambda <- as.matrix(U[indices.columns,indices.columns])

    #================================================
    #================================================
    #if(abs(lambda[1,1])<1){
    #================================================
    #================================================
    if(abs(eigen(lambda)$values)[1]<1){
      eps.t.star <- V.star[1,indices.columns] # first estimate of eps.t.star
      eps.star <- matrix(eps.t.star,1,length(indices.columns))
      for(t in 2:T){
        eps.t.star <- V.star[t,indices.columns] + lambda %*% c(eps.t.star)
        eps.star <- rbind(eps.star,c(eps.t.star)) # FORWARD DIRECTION
      }
    }else{
      eps.t.star <- rep(0,length(indices.columns)) # Last estimate of eps.t.star
      eps.star <- matrix(0,1,length(indices.columns))
      for(t in (T-1):1){
        U_k_1 <- solve(lambda)
        eps.t.star <- - U_k_1 %*% V.star[t+1,indices.columns] + U_k_1 %*% eps.t.star
        eps.star <- rbind(Re(c(eps.t.star)),eps.star) # BACKWARD DIRECTION
      }
    }

    EPS.star.est <- cbind(eps.star,EPS.star.est)

    if(count>0){
      V.star[,1:count] <- V.star[,1:count] +
        rbind(rep(0,length(indices.columns)),
              as.matrix(eps.star[1:(T-1),])) %*% t(matrix(U[1:count,indices.columns],count,length(indices.columns)))
    }
  }
  EPS.est <- EPS.star.est %*% t(A)

  if(indic.solve.C==0){
    ETA.est <- EPS.est %*% t(solve(Model$C))
  }else{
    C.eig <- eigen(Model$C)
    AAUUXX <- max(abs(C.eig$values))/min(abs(C.eig$values))
    AAUUXX2 <- abs(Im(C.eig$values)/Re(C.eig$values))
    if((AAUUXX>10^14)|(max(AAUUXX2)>10^5)){
      C_1 <- diag(dim(Model$C)[1])
    }else{
      C_1 <- solve(Model$C)
    }
    ETA.est <- EPS.est %*% t(C_1)
  }

  return(list(EPS.est=EPS.est,V.est=V,ETA.est=ETA.est))
}



build.othonormal.basis <- function(z){
  # z is a vector of dimension n
  # This function creates an orthonromal matrix K whose first column is the normalized value of z.
  n <- length(z)
  I <- diag(n)
  if(z[1]!=0){
    Z <- matrix(I[,2:n],n,n-1)
  }else{
    Z <- matrix(I[,1:(n-1)],n,n-1)
  }
  if(sum(z^2)==0){
    K <- NaN
    print("z should not be 0")
  }else{
    K <- matrix(z/sqrt(sum(z^2)),ncol=1)
    for(i in 2:n){
      y <- matrix(Z[,i-1],ncol=1)
      X <- K[,1:(i-1)]
      k <- (I - X %*% solve(t(X) %*% X) %*% t(X)) %*% y
      k <- k/sqrt(sum(k^2))
      K <- cbind(K,k)
    }
  }
  return(K)
}


compute.other.basic.form <- function(Theta,C,w){
  # w is a vector whose entries are 1 or -1 (see Lippi and Reichlin 1994)
  # w is of dimension n, where Theta is n x n.
  # If w_i=-1 then the output MA represention corresponds to an MA representation
  # with the same spetral density as the input one, but the i^th pole of det(I - Theta.L) would have been replaced
  # by its inverse conjugate.

  n <- dim(C)[1]

  poles.ini <- 1/eigen(Theta)$values

  Theta.k <- Theta
  C.k <- C

  for(nb.pole in 1:n){
    if(w[nb.pole]!=1){
      pole <- poles.ini[nb.pole]
      poles.k <- 1/eigen(Theta.k)$values
      indic.pole.k <- which(abs(poles.k-pole)==min(abs(poles.k-pole)))[1]

      k <- eigen(Theta.k)$vectors[,indic.pole.k]
      k.star <- solve(C.k) %*% k

      K <- build.othonormal.basis(k.star)

      CK.mod <- C.k %*% K
      a <- CK.mod[,1]
      CK.mod[,1] <- - a * (1/Conj(pole))

      Theta.mod <- Theta.k %*% C.k %*% K
      Theta.mod[,1] <- - a

      C.k <- CK.mod
      Theta.k <- Theta.mod %*% solve(CK.mod)
    }
  }
  return(list(Theta=Theta.k,C=C.k))
}


compute.all.forms <- function(Theta,C){

  n <- dim(C)[1]

  if(n==1){
    all.Theta <- array(c(1/Theta,Theta),c(1,1,2))
    all.C <- array(c(C*Theta,C),c(1,1,2))
  }else{
    combi <- matrix(0,2^n,n)
    basic <- matrix(c(0,1),ncol=1)
    for(i in 1:n){
      combi[,i] <- rep(basic,2^(i-1)) %x% rep(1,2^(n-i))
    }

    all.Theta <- array(NaN,c(n,n,2^n))
    all.C <- array(NaN,c(n,n,2^n))

    for(i in 1:(2^n)){
      aux <- compute.other.basic.form(Theta,C,combi[i,])
      all.Theta[,,i] <- aux$Theta
      all.C[,,i] <- aux$C
    }
  }

  return(list(all.Theta=all.Theta,all.C=all.C))
}







# ===============================================
# ===============================================
# Other various procedures
# ===============================================
# ===============================================


make.PHI2 <- function(Phi){
  p <- dim(Phi)[3]
  n <- dim(Phi)[1]
  PHI <- matrix(0,n*p,n*p)
  if(p>1){
    PHI[(n+1):(n*p),1:((p-1)*n)] <- diag((p-1)*n)
  }
  for(i in 1:p){
    PHI[1:n,((i-1)*n+1):(i*n)] <- Phi[,,i]
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


# log.g.gaussian <- function(x,mu,sigma){
#   # Gaussian distribution
#   log.g <- -1/2 * log(2*pi) - log(sigma) -(x - mu)^2/(2*sigma^2)
#   return(list(log.g=log.g))
# }
#
#
# log.g.student <- function(x,nu){
#   # Student distribution
#   log.g <- -(1 + nu)/2 * log(1 + (x * sqrt(nu/(nu-2)))^2/nu) +
#     log(gamma((nu+1)/2)) - .5*log(nu*pi) - log(gamma(nu/2)) + 1/2 * log(nu/(nu-2))
#   return(list(log.g=log.g))
# }
#
#
# log.g.mixt.gaussian <- function(x,mu,sigma,p){
#   # Mixture of Gaussian distributions
#   # mu is 2 x 1 vector; sigma is a 2 x1 vector, p is a scalar
#   mu_1 <- mu[1]
#   mu_2 <- mu[2]
#   sigma_1 <- sigma[1]
#   sigma_2 <- sigma[2]
#   aux <-     p / sqrt(2*pi) / sigma_1 * exp(- (x - mu_1)^2/(2 * sigma_1^2)) +
#     (1-p) / sqrt(2*pi) / sigma_2 * exp(- (x - mu_2)^2/(2 * sigma_2^2))
#   aux[aux==0] <- 10^(-100)
#   log.g <- log(aux)
#
#   return(list(log.g=log.g))
# }


g <- function(eps,Model){
  # evaluate the log pdf at the points given in eps, where eps is of dimension T x n
  n <- dim(eps)[2]
  distri <- Model$distri
  log.L <- NULL
  if(sum(abs(eigen(Model$C)$values)<10^(-10))>0){
    return(-100000)
  }
  if(indic.solve.C==0){
    eta <- eps %*% t(solve(Model$C))
  }else{
    C.eig <- eigen(Model$C)
    AAUUXX <- max(abs(C.eig$values))/min(abs(C.eig$values))
    AAUUXX2 <- abs(Im(C.eig$values)/Re(C.eig$values))
    if((AAUUXX>10^14)|(max(AAUUXX2)>10^5)){
      C_1 <- diag(dim(Model$C)[1])
    }else{
      C_1 <- solve(Model$C)
    }
    eta <- eps %*% t(C_1)
  }
  for(ii in 1:n){
    if(distri$type[ii] == "student"){
      g.aux <- log.g.student(eta[,ii],nu = distri$df[ii])
    }else if(distri$type[ii] == "gaussian"){
      g.aux <- log.g.gaussian(eta[,ii],mu = 0,sigma = 1)
    }else if(distri$type[ii] == "mixt.gaussian"){
      p = distri$p[ii]
      mu.1 <- distri$mu[ii]
      sigma.1 <- distri$sigma[ii]
      mu.2 <- - p/(1-p)*mu.1
      if(is.na(1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2))){
        print(p)
        print(mu.1)
        print(sigma.1)
      }
      if((1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2))<0){
        sigma.2 <- 100000
      }else{
        sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
      }
      g.aux <- log.g.mixt.gaussian(eta[,ii],mu = c(mu.1,mu.2),sigma = c(sigma.1,sigma.2),p)
    }
    log.L    <- cbind(log.L,g.aux$log.g)
  }
  logl <- log.L %*% matrix(1,n,1) - log(abs(det(Model$C)))
  return(logl)
}



make.nice.distri.name <- function(distri){
  vec.names <- NULL
  for(i in 1:length(distri$type)){
    if(distri$type[i]=="gaussian"){
      vec.names = c(vec.names,
                    "Gaussian")
    }else if(distri$type[i]=="student"){
      vec.names = c(vec.names,
                    paste("t(",toString(distri$df[i]),")",sep=""))
    }else if(distri$type[i]=="mixt.gaussian"){
      vec.names = c(vec.names,"mixt.-Gaussian")
    }
  }
  return(vec.names)
}


round.fixed.length <- function(X,n){
  # This procedure is used in the automatic creation of Latex tables
  # x is numeric. The output is a string with n numbers after ".", even if they are 0s.
  if(is.na(X)){
    string.x <- "NaN"
  }else{
    signX <- sign(X)
    if(signX>=0){
      signX <- NULL
    }else{
      signX <- "-"
    }
    x <- abs(X)
    string.integer <- toString(as.integer(x))
    string.decimal <- trunc(10^n*(x - as.integer(x,0)))/10^n
    # string.decimal <- toString(round(x - as.integer(x,0),n))
    while(nchar(string.decimal)<n+2){
      string.decimal <- paste(string.decimal,"0",sep="")
    }
    if(n==0){
      string.x <- paste(signX,string.integer,sep="")
    }else{
      string.x <- paste(signX,string.integer,".",str_replace(string.decimal,"0.",""),sep="")
    }
    if(as.numeric(string.x)==0){
      if(n==0){
        string.x <- paste(string.integer,sep="")
      }else{
        string.x <- paste(string.integer,".",str_replace(string.decimal,"0.",""),sep="")
      }
    }
  }
  return(string.x)
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


compute.cov.ratio <- function(All_res,vec.conf.levels,theta,sigma){

  # Keep all results:
  N <- dim(All_res)[1]
  quantiles <- qnorm(vec.conf.levels + (1 - vec.conf.levels)/2)
  theta.is.in.CI <- array(NaN,c(N,12,length(vec.conf.levels)))
  sigma.is.in.CI <- array(NaN,c(N,12,length(vec.conf.levels)))
  for(quant in 1:length(quantiles)){
    lower.bounds <- All_res[,c(1:4,21:24,41:44)] - All_res[,c(13:16,33:36,53:56)] * quantiles[quant]
    upper.bounds <- All_res[,c(1:4,21:24,41:44)] + All_res[,c(13:16,33:36,53:56)] * quantiles[quant]
    theta.is.in.CI[,,quant] <- (theta>lower.bounds)*(theta<upper.bounds)
    lower.bounds <- abs(All_res[,c(5:8,25:28,45:48)]) - All_res[,c(17:20,37:40,57:60)] * quantiles[quant]
    upper.bounds <- abs(All_res[,c(5:8,25:28,45:48)]) + All_res[,c(17:20,37:40,57:60)] * quantiles[quant]
    sigma.is.in.CI[,,quant] <- (sigma>lower.bounds)*(sigma<upper.bounds)
  }

  # results (lines = given distri,column = given confidence level)
  cov.ratio.res.theta <- apply(theta.is.in.CI,c(2,3),function(x){mean(x,na.rm=TRUE)})
  cov.ratio.res.sigma <- apply(sigma.is.in.CI,c(2,3),function(x){mean(x,na.rm=TRUE)})

  # Remove observations whose abs. value is on the wrong side:
  All_res_dupli <- All_res
  for(i in c(1:4,21:24,41:44)){
    if(abs(theta)<=1){
      All_res_dupli[abs(All_res_dupli[,i])>1,i] <- NaN
    }else{
      All_res_dupli[abs(All_res_dupli[,i])<=1,i] <- NaN
    }
  }
  quantiles <- qnorm(vec.conf.levels + (1 - vec.conf.levels)/2)
  theta.is.in.CI.right.side <- array(NaN,c(N,12,length(vec.conf.levels)))
  sigma.is.in.CI.right.side <- array(NaN,c(N,12,length(vec.conf.levels)))
  for(quant in 1:length(quantiles)){
    lower.bounds <- All_res_dupli[,c(1:4,21:24,41:44)] - All_res_dupli[,c(13:16,33:36,53:56)] * quantiles[quant]
    upper.bounds <- All_res_dupli[,c(1:4,21:24,41:44)] + All_res_dupli[,c(13:16,33:36,53:56)] * quantiles[quant]
    theta.is.in.CI.right.side[,,quant] <- (theta>lower.bounds)*(theta<upper.bounds)
    lower.bounds <- abs(All_res_dupli[,c(5:8,25:28,45:48)]) - All_res_dupli[,c(17:20,37:40,57:60)] * quantiles[quant]
    upper.bounds <- abs(All_res_dupli[,c(5:8,25:28,45:48)]) + All_res_dupli[,c(17:20,37:40,57:60)] * quantiles[quant]
    sigma.is.in.CI.right.side[,,quant] <- (sigma>lower.bounds)*(sigma<upper.bounds)
  }

  # results (lines = given distri,column = given confidence level)
  cov.ratio.res.theta.right.side <- apply(theta.is.in.CI.right.side,c(2,3),function(x){mean(x,na.rm=TRUE)})
  cov.ratio.res.sigma.right.side <- apply(sigma.is.in.CI.right.side,c(2,3),function(x){mean(x,na.rm=TRUE)})


  # compute share of theta's on the right side
  proportion.right.side <- NULL
  for(i in c(1:4,21:24,41:44)){
    aux <- All_res[!is.na(All_res[,i]),i]
    nb.total <- length(aux)
    aux <- All_res_dupli[!is.na(All_res_dupli[,i]),i]
    nb.right.side <- length(aux)
    proportion.right.side <- c(proportion.right.side,nb.right.side/nb.total)
  }

  return(
    list(
      cov.ratio.res.theta = cov.ratio.res.theta,
      cov.ratio.res.sigma = cov.ratio.res.sigma,
      cov.ratio.res.theta.right.side = cov.ratio.res.theta.right.side,
      cov.ratio.res.sigma.right.side = cov.ratio.res.sigma.right.side,
      proportion.right.side = proportion.right.side
    )
  )
}


compute.cov.ratio.n2 <- function(estimated.param,true.param,
                                 stdev.param,vec.conf.levels,
                                 indic.keep){
  # indic.keep is of the same length as estimated.param,
  # if different from 1, then not taken into account for "reduced" coverage ratios

  # Keep all results:
  N <- length(estimated.param)
  quantiles <- qnorm(vec.conf.levels + (1 - vec.conf.levels)/2)
  param.is.in.CI <- matrix(NaN,N,length(vec.conf.levels))
  for(quant in 1:length(quantiles)){
    lower.bounds <- estimated.param - stdev.param * quantiles[quant]
    upper.bounds <- estimated.param + stdev.param * quantiles[quant]
    param.is.in.CI[,quant] <- (true.param>lower.bounds)*(true.param<upper.bounds)
  }

  # results (lines = given distri,column = given confidence level)
  cov.ratio <- apply(param.is.in.CI,2,function(x){mean(x,na.rm=TRUE)})

  # Remove observations whose abs. value is on the wrong side:
  param.is.in.CI.right.reg <- param.is.in.CI
  param.is.in.CI.right.reg[!(indic.keep==1)] <- NaN

  # results (lines = given distri,column = given confidence level)
  cov.ratio.right.reg <- apply(param.is.in.CI.right.reg,2,function(x){mean(x,na.rm=TRUE)})


  # compute share of theta's on the right side
  nb.total <- sum(!is.na(estimated.param))
  indic.keep.aux <- indic.keep
  indic.keep.aux[is.na(indic.keep.aux)] <- 0
  nb.right.side <- sum(indic.keep.aux)
  proportion.right.reg <- nb.right.side/nb.total

  return(
    list(
      cov.ratio = cov.ratio,
      cov.ratio.right.reg = cov.ratio.right.reg,
      proportion.right.reg = proportion.right.reg
    )
  )
}



