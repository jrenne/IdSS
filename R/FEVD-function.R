
#' Identify structural VAR shocks via FEVD maximization (Uhlig-style)
#'
#' Identifies one or more structural shocks in a VAR by
#' maximizing the (windowed) forecast-error variance (FEV) of user-selected
#' variables over user-selected horizons. For each requested shock, it finds the
#' orthonormal rotation that maximizes the FEV contribution of the target
#' variable on a given horizon window, applies a sign normalization, enforces
#' orthogonality with previously identified shocks, and returns impulse
#' responses (IRFs). Optional parametric bootstrap yields IRF uncertainty bands.
#'
#' @param Y A numeric matrix of size \eqn{T \times n} with the VAR variables.
#' @param p VAR lag order (integer).
#' @param nb.shocks Number of structural shocks to identify
#'   sequentially.
#' @param names.of.shocks Character vector of length `nb.shocks`. Labels used
#'   in plots.
#' @param H1 Integer vectors (length `nb.shocks`). Start horizons
#' for the FEV-maximization window of each shock.
#' @param H2 Integer vectors (length `nb.shocks`). End horizons (inclusive)
#' for the FEV-maximization window of each shock.
#' @param variable Integer vector. 1-based indices of the
#'   target variables whose FEV is maximized for each shock.
#' @param norm Integer vector of length `nb.shocks`. For each shock, the horizon
#'   IRF array index at which its response is normalized to be positive. Index
#'   1 is impact (horizon zero). H1 and H2 instead count forecast steps: a
#'   forecast of H steps uses response indices 1 through H.
#' @param nb.periods.IRF Number of IRF horizons to compute.
#' @param bootstrap.replications Number of parametric bootstrap
#'   replications. If 1, returns point estimates only (no uncertainty).
#' @param confidence.interval Numeric in (0,1). Central probability mass for IRF
#'   bands (default 0.90).
#' @param indic.plot If 1 (default), plots median IRFs with
#'   upper/lower quantile bands for each identified shock.
#'
#' @return A list with components:
#' \describe{
#'   \item{simulated.IRFs}{Array of identified IRFs across bootstrap draws.}
#'   \item{stdv.IRFs}{Array of pointwise bootstrap standard deviations.}
#'   \item{CI.lower.bounds}{Array of lower pointwise quantiles.}
#'   \item{CI.upper.bounds}{Array of upper pointwise quantiles.}
#' }
#' @examples
#' set.seed(123)
#' Y <- cbind(output = rnorm(160), inflation = rnorm(160))
#' # Identify the shock explaining most of output's FEV over horizons 1 to 4.
#' # One replication requests point estimates, without bootstrap inference.
#' fit <- svar.fevmax(Y, p = 1, nb.shocks = 1,
#'                   names.of.shocks = "Output shock", H1 = 1, H2 = 4,
#'                   variable = 1, norm = 1, nb.periods.IRF = 8,
#'                   bootstrap.replications = 1, indic.plot = 0)
#' dim(fit$simulated.IRFs)
#' fit$simulated.IRFs[, 1, 1, 1]  # impact responses to the identified shock
#' @export
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
  simulated.Phi   <- matrix(rep(Phi,1L),nrow=p,ncol=1L)
  # if bootstrap then generate and store simulated IRFs, B.hat and Phi
  if(bootstrap.replications>1){
    bootstrap.res <- param.bootstrap(Y,p,nb.periods.IRF,n,bootstrap.replications,
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
      for (h in seq_len(H2[i])){
        V99 <- simulated.IRFs[,,h,l]%*%Q[,i:n,l]
        # Notice that we use only columns 2 to n of Q:
        #  the first column selects the TFP surprise shock,
        #  which is the first shock in the Cholesky
        #  decomposition where TFP is ordered first.
        JJ <- (H2[i]+1-max(H1[i],h))*t(V99)%*%E[[i]]%*%V99
        WWW <- WWW+JJ}
      r <- eigen(WWW)
      # Take the eigenvector with the highest eigenvalue
      eigvec <- matrix(r$vectors[,1],n-i+1,1)
      # We might need to adjust the sign
      if (simulated.IRFs[variable[i],,norm[i],l]%*%Q[,i:n,l]%*%eigvec>0){
        Q[,i,l] <- Q[,i:n,l]%*%eigvec
      }else{
        Q[,i,l] <- - Q[,i:n,l]%*%eigvec}
      if (i < n) {
        Q[,(i+1):n,l] <- Null(Q[,1:i,l])
      }
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
        plot(seq_len(nb.periods.IRF)-1L,CI.median[i,j,],type="l",lwd=2,xlab="",ylab="",
             ylim=c(min(CI.lower.bounds[i,j,]),
                    max(CI.upper.bounds[i,j,])),
             main=paste0("(",letters[i],") ",names.of.variables[i]))
        abline(h=0,col="grey")
        lines(seq_len(nb.periods.IRF)-1L,CI.lower.bounds[i,j,],col="grey45",lty=2,lwd=2)
        lines(seq_len(nb.periods.IRF)-1L,CI.upper.bounds[i,j,],col="grey45",lty=2,lwd=2)}
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
