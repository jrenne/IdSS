###########################################3
# Various procedures useful for
# VAR with sign restrictions
###########################################3


###################################################################3
# SVAR with sign and (short-run) zero restrictions

svar.signs <- function(Y,p,
                       nb.shocks,
                       nb.periods.IRF = 20,
                       bootstrap.replications = 1,
                       confidence.interval = 0.90, # expressed in pp.
                       indic.plot = 1, # Plots are displayed if = 1.
                       nb.draws = 1000, # number of draws.
                       sign.restrictions,
                       horizon,
                       recursive = 0,
                       SR.restrictions = list()
){
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # Uhlig (2005) and Arias et al. (2019) methodology
  # Inference according to Arias et al. (2009)
  
  names.of.variables <- colnames(Y)
  
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  
  # Matrices useful for simulated IRFs
  y0.star <- rep(0,n*p)

  # do Choleski (we need the FULL IRFs)
  cholesky.res <- svar.ordering.all(Y,p,
                                    nb.periods.IRF,
                                    n)
  # Store results
  IRFs <- cholesky.res$IRFs
  est.VAR <- cholesky.res$est.VAR
  Phi     <- Acoef(est.VAR)
  B.hat <- cholesky.res$B.hat
  cst     <- Bcoef(est.VAR)[,p*n+1]
  resids  <- residuals(est.VAR)
  Omega   <- var(resids)

# if no bootstrap then simply use point estimate
    simulated.IRFs  <- replicate(1,IRFs, simplify="array")
    simulated.B.hat <- replicate(1,B.hat, simplify="array")
    simulated.Phi   <- replicate(1,Phi, simplify="array")

# if bootstrap then generate and store simulated IRFs, B.hat and Phi
  if(bootstrap.replications>1){

    bootstrap.res <- nonparam.bootstrap(Y,p,
                                      nb.periods.IRF,
                                      n,
                                      bootstrap.replications,
                                      posit.of.shock = 0)
    simulated.IRFs  <- bootstrap.res$simulated.IRFs
    simulated.B.hat <- bootstrap.res$simulated.B.hat
    simulated.Phi   <- bootstrap.res$simulated.Phi
}

# IRFs satisfying sign restrictions (and possibly zero restrictions)
    if (length(SR.restrictions)==0 & recursive ==0){
      # sign restrictions only (non-recursive)
  signs.res <- svar.signs.aux(Y,p,
                                 nb.periods.IRF,
                                 nb.shocks,
                                 sign.restrictions,
                                 horizon,
                                 simulated.IRFs,
                                 simulated.B.hat,
                                 simulated.Phi,
                                 nb.draws,
                                 bootstrap.replications)
    }else if (length(SR.restrictions)==0 & recursive == 1){
      # sign restrictions only (recursive)
      signs.res <- svar.signs.aux.recursive(Y,p,
                                            nb.periods.IRF,
                                            nb.shocks,
                                            sign.restrictions,
                                            horizon,
                                            simulated.IRFs,
                                            simulated.B.hat,
                                            simulated.Phi,
                                            nb.draws,
                                            bootstrap.replications)      
    } else{
      # sign and zero restrictions
      signs.res <- svar.signs.aux.zeros(Y,p,
                                  nb.periods.IRF,
                                  nb.shocks,
                                  sign.restrictions,
                                  horizon,
                                  SR.restrictions,
                                  simulated.IRFs,
                                  simulated.B.hat,
                                  simulated.Phi,
                                  nb.draws,
                                  bootstrap.replications)
    }
  
  xx         <- signs.res$xx         # number of selected rotations
  IRFs.signs <- signs.res$IRFs.signs # Set of admissible IRFs

  for (j in 1:nb.shocks){
    # reshape the IRFs
    IRFs.signs.j <- IRFs.signs[,j,,]
    dim(IRFs.signs.j) <- c(n,nb.periods.IRF,xx)
   # compute the moments
    all.stdv.IRFs.j <- NULL
    all.CI.lower.bounds.j <- NULL
    all.CI.upper.bounds.j <- NULL
    all.CI.median.j <- NULL
    for(i in 1:n){
      all.IRFs.i <- matrix(IRFs.signs.j[i,,],nrow=nb.periods.IRF,ncol=xx)
      
      stdv.IRF <- apply(all.IRFs.i,1,sd)
      all.stdv.IRFs.j <- cbind(all.stdv.IRFs.j,stdv.IRF)
      
      all.CI.lower.bounds.j <- cbind(all.CI.lower.bounds.j,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,(1-confidence.interval)/2)}))
      all.CI.upper.bounds.j <- cbind(all.CI.upper.bounds.j,
                                   apply(all.IRFs.i,1,
                                         function(x){quantile(x,1-(1-confidence.interval)/2)}))
      all.CI.median.j <- cbind(all.CI.median.j,
                             apply(all.IRFs.i,1,
                                   function(x){quantile(x,0.5)}))
    }
    
  # Figure  
    if(indic.plot==1){
        par(mfrow=c(2,
                    ifelse(round(n/2)==n/2,n/2,(n+1)/2)))
        for(i in 1:n){
          plot(all.CI.median.j[,i],type="l",lwd=2,xlab="",ylab="",
               ylim=c(min(all.CI.lower.bounds.j[,i]),
                      max(all.CI.upper.bounds.j[,i])),
               main=paste("Effect of shock ",j," on ",names.of.variables[i],sep=""))
          abline(h=0,col="grey")
          lines(all.CI.lower.bounds.j[,i],col="red",lty=2,lwd=2)
          lines(all.CI.upper.bounds.j[,i],col="red",lty=2,lwd=2)
      }
    }
  }
  all.CI.lower.bounds <- list()
  all.CI.upper.bounds <- list()
  all.CI.median       <- list()
  all.stdv.IRFs       <- list()
  for (j in 1:nb.shocks){
    all.CI.lower.bounds[[j]] <- all.CI.lower.bounds.j
    all.CI.upper.bounds[[j]] <- all.CI.upper.bounds.j
    all.CI.median[[j]] <- all.CI.median.j
    all.stdv.IRFs[[j]] <- all.stdv.IRFs.j
  }  
  
  return(
    list(
      xx = xx,
      IRFs.signs = IRFs.signs,
      simulated.IRFs = simulated.IRFs,
      all.CI.median = all.CI.median,
      all.stdv.IRFs = all.stdv.IRFs,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds
    ))
}


###################################################################3
# Sign restrictions (non-recursive)

svar.signs.aux <- function(Y,p,
                              nb.periods.IRF,
                              nb.shocks,
                              sign.restrictions,
                              horizon,
                              simulated.IRFs,
                              simulated.B.hat,
                              simulated.Phi,
                              nb.draws,
                              bootstrap.replications){
  
  considered.variables <- colnames(Y)
  n <- dim(Y)[2]
  
  # Draw orthonormal matrices and parameters
  xx <- 0
  admissible.rotation <- list()
  IRFs.draw <- list()
  B.hat.draw <- list()
  Phi.draw <- list()
  for (l in bootstrap.replications){
    for(i in 1:nb.draws){
      A <- matrix(rnorm(n*n), nrow = n, ncol = n) # draw a random matrix
      AA <- qr(A)
      Q <- qr.Q(AA, complete = FALSE) # use QR decomposition of A to compute Q
      ff <- 0
        for (shock in 1:nb.shocks){
          k <- NULL
          for (h in horizon[[shock]]){
            K <- sign.restrictions[[shock]]%*%simulated.IRFs[,,h,l]%*%Q[,shock]
            kk <- matrix(K[K >0], nrow=1)
            k <- cbind(k,kk)
          }
          # check if sign restrictions are satisfied
        if (length(k) == rankMatrix(sign.restrictions[[shock]])[1]*length(horizon[[shock]]))
          {ff <- ff+1}
      }
      if (ff==nb.shocks)
        {xx <- xx+1
        admissible.rotation[[xx]] <- Q
        IRFs.draw[[xx]] <- simulated.IRFs[,,,l]
        B.hat.draw[[xx]] <- simulated.B.hat[,,l]
        Phi.draw[[xx]] <- simulated.Phi[,l]
      }
    }
  }
  

  IRFs.signs <- array(NaN,c(n,n,nb.periods.IRF,xx))
   for (k in 1:xx){
     for (t in 1:nb.periods.IRF){
       IRFs.signs[,,t,k] <- IRFs.draw[[k]][,,t]%*%admissible.rotation[[k]]
    }
  }

  
  return(list(
    xx = xx,
    IRFs.signs = IRFs.signs, # IRFs.draw*Q
    B.hat.draw = B.hat.draw,
    Phi.draw = Phi.draw,
    admissible.rotation = admissible.rotation # Q
  ))
}

###################################################################3
# Sign restrictions (recursive)

svar.signs.aux.recursive <- function(Y,p,
                                     nb.periods.IRF,
                                     nb.shocks,
                                     sign.restrictions,
                                     horizon,
                                     simulated.IRFs,
                                     simulated.B.hat,
                                     simulated.Phi,
                                     nb.draws,
                                     bootstrap.replications){
  
  considered.variables <- colnames(Y)
  n <- dim(Y)[2]
  
  # Draw orthonormal matrices and parameters
  xx <- 0
  admissible.rotation <- list()
  IRFs.draw <- list()
  B.hat.draw <- list()
  Phi.draw <- list()
  for (l in bootstrap.replications){
    for(i in 1:nb.draws){
      Q <- NULL # initialize Q
      ff <- 0
      for (shock in 1:nb.shocks){
        k <- NULL
        a <- matrix(rnorm(n-shock+1,1), nrow = n-shock+1, ncol = 1)
        a <- a%*%solve((t(a)%*%a)^0.5)
        if (shock ==1){q <-a}else{
        q <- Null(Q)%*%a} # q is of length 1 and orthogonal to the vectors of Q
        for (h in horizon[[shock]]){
          K <- sign.restrictions[[shock]]%*%simulated.IRFs[,,h,l]%*%q
          kk <- matrix(K[K >0], nrow=1)
          k <- cbind(k,kk)
        }
        # check if sign restrictions are satisfied
        if (length(k) == rankMatrix(sign.restrictions[[shock]])[1]*length(horizon[[shock]])){
          ff <- ff+1
          Q <- cbind(Q,q) # update Q with the new vector q
        }
        if (is.null(Q)){break}
      }
      if (ff==nb.shocks)
        {xx <- xx+1
        Q <- cbind(Q,Null(Q))
        admissible.rotation[[xx]] <- Q
        IRFs.draw[[xx]] <- simulated.IRFs[,,,l]
        B.hat.draw[[xx]] <- simulated.B.hat[,,l]
        Phi.draw[[xx]] <- simulated.Phi[,l]
        }
      }
  }
  
  
  IRFs.signs <- array(NaN,c(n,n,nb.periods.IRF,xx))
  for (k in 1:xx){
    for (t in 1:nb.periods.IRF){
      IRFs.signs[,,t,k] <- IRFs.draw[[k]][,,t]%*%admissible.rotation[[k]]
    }
  }
  
  
  return(list(
    xx = xx,
    IRFs.signs = IRFs.signs, # IRFs.draw*Q
    B.hat.draw = B.hat.draw,
    Phi.draw = Phi.draw,
    admissible.rotation = admissible.rotation # Q
  ))
}

###################################################################3
# Sign and zero restrictions

svar.signs.aux.zeros <- function(Y,p,
                                     nb.periods.IRF,
                                     nb.shocks,
                                     sign.restrictions,
                                     horizon,
                                     SR.restrictions,
                                     simulated.IRFs,
                                     simulated.B.hat,
                                     simulated.Phi,
                                     nb.draws,
                                     bootstrap.replications){
  
  considered.variables <- colnames(Y)
  n <- dim(Y)[2]
  
  # Draw orthonormal matrices and parameters
  xx <- 0
  admissible.rotation <- list()
  IRFs.draw <- list()
  B.hat.draw <- list()
  Phi.draw <- list()
  for (l in bootstrap.replications){
    for(i in 1:nb.draws){
      Q <- NULL # initialize Q
      ff <- 0
      for (shock in 1:nb.shocks){
        k <- NULL
          J <- cbind(Q,t(SR.restrictions[[shock]]%*%simulated.B.hat[,,l]))
          a <- matrix(rnorm(n-dim(J)[2],1), nrow = n-dim(J)[2], ncol = 1)
          a <- a%*%solve((t(a)%*%a)^0.5)
          q <- Null(J)%*%a # q is orthogonal to Q and satisfies the zero restrictions
        for (h in horizon[[shock]]){
          K <- sign.restrictions[[shock]]%*%simulated.IRFs[,,h,l]%*%q
          kk <- matrix(K[K >0], nrow=1)
          k <- cbind(k,kk)
        }
          # check if sign restrictions are satisfied
        if (length(k) == rankMatrix(sign.restrictions[[shock]])[1]*length(horizon[[shock]])){
          ff <- ff+1
          Q <- cbind(Q,q) # update Q
        }
          if (is.null(Q)){break}
      }
      if (ff==nb.shocks)
        {xx <- xx+1
        Q <- cbind(Q,Null(Q))
        admissible.rotation[[xx]] <- Q
        IRFs.draw[[xx]] <- simulated.IRFs[,,,l]
        B.hat.draw[[xx]] <- simulated.B.hat[,,l]
        Phi.draw[[xx]] <- simulated.Phi[,l]
        }
      }
  }
  
  
  IRFs.signs <- array(NaN,c(n,n,nb.periods.IRF,xx))
  for (k in 1:xx){
    for (t in 1:nb.periods.IRF){
      IRFs.signs[,,t,k] <- IRFs.draw[[k]][,,t]%*%admissible.rotation[[k]]
    }
  }
  
  
  return(list(
    xx = xx,
    IRFs.signs = IRFs.signs, # IRFs.draw*Q
    B.hat.draw = B.hat.draw,
    Phi.draw = Phi.draw,
    admissible.rotation = admissible.rotation # Q
  ))
}



