
#' Simulate a (S)VAR Process
#'
#' This function has two modes, controlled by \code{indic.IRF}:
#' \itemize{
#'  \item 1) Simulation mode (default, \code{indic.IRF == 0}):
#' starting from the initial state \code{y0.star}, simulates
#' \code{nb.sim} periods of a VAR(p).
#' If \code{eta} is supplied (an \code{nb.sim x n} matrix), its rows are used as
#' the structural shocks. Otherwise shocks are generated as i.i.d. standard
#' normal and then mapped by \code{B}.
#' The result is an \code{nb.sim x n} matrix of simulated observations.
#' \item 2) IRF mode (\code{indic.IRF != 0}):
#' compute a finite-horizon impulse response to a one-time structural shock
#' (length n). At t = 1, the state is set to the effect
#' of \code{u.shock} through \code{B}; for t >= 2 the system evolves
#' deterministically using the companion dynamics (i.e., no further shocks).
#' The output is an \code{nb.sim x n} matrix containing the IRF
#' path. In IRF mode, \code{eta} is ignored. }
#' @param c Numeric vector of length \eqn{n}; the intercept term.
#' @param Phi A list, each element of which is a \eqn{Phi_i} matrix. Hence it has p elements if we consider a VAR(p).
#' @param B Numeric \eqn{n \times n} matrix; contemporaneous impact matrix
#'   mapping structural shocks \eqn{u_t}to the reduced-form residuals.
#' @param nb.sim Integer; number of simulated periods to produce (horizon).
#' @param y0.star Numeric vector of length \eqn{n p}; initial companion-state
#'  . Required when
#'   \code{indic.IRF == 0}. Default is \code{NaN}.
#' @param indic.IRF Integer flag. If \code{0} (default), run a forward
#'   simulation from \code{y0.star}. If non-zero, compute an IRF.
#' @param u.shock Numeric vector of length n; one-time structural shock used only in IRF mode.
#' @param eta  Optional matrix of dimension \code{nb.sim x n} containing the
#'   structural shocks \eqn{u_t} to use in simulation mode. If supplied, these are
#'   used; otherwise \eqn{u_t \sim \mathcal{N}(0, I_n)} is generated.
#'
#' @return A numeric matrix of size \code{nb.sim x n}:
#'   simulated series in simulation mode, or the IRF path in IRF mode.
#' @examples
#' ##  Example with (n = 2, p = 2)
#' Phi <- list(
#'   matrix(c(0.5, 0.1, 0.0, 0.3), 2, 2, byrow = TRUE),
#'   matrix(c(0.2, 0.0, 0.0, 0.1), 2, 2, byrow = TRUE)
#' )
#' B  <- diag(2)
#' c0 <- c(0, 0)
#' y0 <- rep(0, length(Phi) * nrow(Phi[[1]]))
#'
#' ## 1) Simulation
#' set.seed(1)
#' Y  <- simul.VAR(c0, Phi, B, nb.sim = 10, y0.star = y0)
#'
#' ## 2) IRF (unit shock in variable 1)
#' irf <- simul.VAR(c0, Phi, B, nb.sim = 12, y0.star = y0,
#'                  indic.IRF = 1, u.shock = c(1, 0))
#'
#' @export
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




#' Simulate a (S)VARMA Process
#'
#'  This function has two modes, controlled by \code{indic.IRF}:
#'  \itemize{
#'    \item 1) Simulation mode (default, \code{indic.IRF == 0}):
#'      generates simulated time series from a Vector Autoregressive Moving-Average (VARMA) model.
#'      The model is defined by autoregressive and moving-average coefficient matrices, a constant,
#'      and a mixing (impulse) matrix for the shocks. The simulation
#'      starts from \code{Y0} when specified; otherwise, from the
#'      process’s unconditional mean.
#'      The result is an \code{nb.sim x n} matrix of simulated observations.
#'    \item 2) IRF mode (\code{indic.IRF != 0}):
#' computes a finite-horizon impulse response to a one-time structural shock \code{eta0}
#' (length n). At t = 1, the system is hit by the shock \code{eta0} (applied through B or C,
#' depending on the model specification).; For t ≥ 2, the system evolves deterministically under the companion dynamics (i.e., with no further shocks).
#' The output is an \code{nb.sim × n} matrix containing the impulse response path.}
#'
#' @param Model A list specifying the VARMA model. Must contain:
#' \itemize{
#'     \item{\code{c}}{ Vector of constants.}
#'     \item{\code{Phi}}{ Array of autoregressive coefficient matrices.}
#'     \item{\code{Theta}}{ Array of moving-average coefficient matrices.}
#'     \item{\code{C or B}}{ Mixing/impulse matrix of dimension \eqn{n \times n}. Captures the contemporaneous impact of
#'     structural shock on \eqn{y_t}.}
#'     \item{\code{distri}}{ list which characterizes the distribution of the shocks. Default is standard normal distributions.}
#'   }
#'
#' @param nb.sim Integer. Number of simulation periods.
#' @param Y0 Vector containing the initial values of the endogenous variables \eqn{Y}.
#'   Must be of dimension \eqn{(p \times n) \times 1}, where n is the dimension of \eqn{y_t} (concatenation of lags \eqn{Y_1, \dots, Y_p}).
#' @param eta0 Vector containing the initial values of the shocks \eqn{\eta}.
#'   Must be of dimension \eqn{(q \times n) \times 1} (concatenation of lags \eqn{\eta_1, \dots, \eta_q}).
#' @param indic.IRF Logical or numeric (0/1). If 1, the function produces IRFs by setting shocks to zero
#'   after the first period. Default is 0.
#'@details
#'In Model, if B is specified, and not C, then the SVARMA specification is:
#' \eqn{y_t = c + \Phi_1·y_{t-1} + ... + \Phi_p·y_{t-p} +
#'       B·\eta_t + \Theta_1·B·\eta_{t-1} + ... + \Theta_q·B·\eta_{t-q}}.
#'        Otherwise, if C is specified, the VARMA specification is:
#'       \eqn{ y_t = c + \Phi_1·y_{t-1} + ... + \Phi_p·y_{t-p} +
#'              C·\eta_t - \Theta_1·C·\eta_{t-1} - ... - \Theta_q·C·\eta_{t-q}}
#'              (Notice the change in the signs in the MA part).
#' @return A list with components:
#' \describe{
#'   \item{\code{Y}}{ Simulated endogenous variables.}
#'   \item{\code{EPS}}{ Simulated reduced-form shocks.}
#'   \item{\code{ETA}}{ Simulated structural shocks.}
#'   \item{\code{V}}{ Simulated moving-average terms.}
#' }
#' @examples
#' Simulate a simple 2-variable VARMA(1,1) process
#' n <- 2
#' p <- 1
#' q <- 1
#'
#' Model <- list(
#'   c = rep(0, n),
#'   Phi = array(0.2, dim = c(n, n, p)),
#'   Theta = array(0.1, dim = c(n, n, q)),
#'   C = diag(n),
#'   distri = list(type = rep("gaussian", n))
#' )
#'
#' Y0 <- rep(0, n * p)
#' eta0 <- rep(0, n * q)
#'
#' sim <- simul.VARMA(Model, nb.sim = 50, Y0 = Y0, eta0 = eta0)
#' irf <- simul.VARMA(Model, nb.sim = 50, Y0 = Y0, eta0 = eta0, indic.IRF = 1)
#' @export
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

#'
#'Build the VAR companion matrix
#'
#' Constructs the block companion matrix \eqn{\Phi} for a VAR(p) from the
#' coefficient matrices \eqn{\Phi_1,\ldots,\Phi_p}. The input can be either
#' a list of \eqn{n \times n} matrices or a 3-D array of dimension
#' \eqn{n \times n \times p}. The returned matrix has size \eqn{(np) \times (np)}
#' with the first block row \eqn{[\Phi_1 \ \Phi_2 \ \cdots \ \Phi_p]}, the
#' sub-diagonal equal to an identity matrix \eqn{I_{n(p-1)}}, and zeros elsewhere.
#' When \eqn{p = 1}, the result is simply \eqn{\Phi_1}.
#'
#' @param Phi: a list object or a 3-D array with coefficient matrices for the lagged endogenous variables.
#'
#' @returns  A numeric matrix representing the
#'   VAR(p) companion matrix.
#' @examples
#' ## Example with a list of coefficient matrices (n = 2, p = 2)
#'Phi_list <- list(
#'  matrix(c(0.5, 0.8,
#'           0.3, 0.2), 2, 2, byrow = TRUE),
#'  matrix(c(0.1, 0.0,
#'           0.0, 0.1), 2, 2, byrow = TRUE))
#'
#'## Example to check for stationarity
#'  library(vars)
#'  data <- US3var[,c("y.gdp.gap","infl")]
#'  estimated.var <- VAR(data,p=3)
#'  Phi <- Acoef(estimated.var)
#'  PHI <- make.PHI(Phi) # autoregressive matrix of companion form.
#'  print(abs(eigen(PHI)$values)) # check stationarity
#'
#' @export
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

#' Simulate independent shocks from user-specified distributions
#'
#' @description
#' Generates \code{nb.sim} i.i.d. draws for each shock described in \code{distri}
#' and returns them as a matrix whose columns correspond to shocks.
#'
#' The distribution of each column is selected by \code{distri$type[i]}.
#' Supported options are:
#' \itemize{
#'   \item \code{"gaussian"}: Standard normal \eqn{\mathcal{N}(0,1)}.
#'   \item \code{"mixt.gaussian"}: Two-component Gaussian mixture with mixing
#'   probability \code{p}. The first component has mean \code{mu} and standard
#'   deviation \code{sigma}. The second component parameters are chosen so that
#'   the mixture has mean 0 and variance 1 (under the constraints implied by the
#'   provided \code{p}, \code{mu}, and \code{sigma}).
#'   \item \code{"student"}: Student-\eqn{t} with \code{df} degrees of freedom,
#'   rescaled to have variance 1 when \code{df > 2}.
#'   \item \code{"laplace"}: Laplace(0, b) draws generated by inverse CDF, with
#'   \code{b = 1/sqrt(2)} so that the variance equals 1.
#'   \item \code{"hyper.sec"}: Hyperbolic secant distribution generated by inverse
#'   transform sampling.
#' }
#'
#' @details
#' The function loops over \code{i = 1, ..., length(distri$type)} and simulates
#' one vector \code{eps.i} of length \code{nb.sim} per distribution specification,
#' then column-binds these vectors.
#'
#' @param distri A list or data frame describing each shock distribution. It must
#' contain at least:
#' \itemize{
#'   \item \code{type}: character vector of length \eqn{k} indicating each distribution.
#' }
#' Additional fields are required depending on \code{type}:
#' \itemize{
#'   \item If \code{type == "mixt.gaussian"}: \code{p}, \code{mu}, \code{sigma}.
#'   \item If \code{type == "student"}: \code{df}.
#' }
#' All fields are expected to be indexable as \code{distri$<field>[i]}.
#'
#' @param nb.sim Positive integer. Number of simulations (rows) to generate for
#' each shock.
#'
#' @return A numeric matrix of dimension \code{nb.sim x k}, where \code{k =
#' length(distri$type)}. Column \code{i} contains the simulated draws for the
#' \code{i}-th shock specification.
#'
#' @examples
#' # Two shocks: standard normal and Student-t with df = 5 (variance-normalized)
#' distri <- list(
#'   type = c("gaussian", "student"),
#'   df   = c(NA, 5)
#' )
#' eps <- simul.distri(distri, nb.sim = 1000)
#' dim(eps)  # 1000 x 2
#'
#' # Gaussian mixture shock (requires p, mu, sigma for that component)
#' distri2 <- list(
#'   type  = c("mixt.gaussian"),
#'   p     = c(0.3),
#'   mu    = c(1.0),
#'   sigma = c(0.5)
#' )
#' eps2 <- simul.distri(distri2, nb.sim = 500)
#'
#' @seealso \code{\link[stats]{rnorm}}, \code{\link[stats]{rt}}, \code{\link[stats]{runif}}
#'
#' @export
simul.distri <- function(distri,nb.sim){
  # Simulation of independent shocks
  eps <- NULL
  for(i in 1:length(distri$type)){
    if(distri$type[i]=="gaussian"){
      eps.i <- rnorm(nb.sim)
    }else if(distri$type[i]=="mixt.gaussian"){
      p <- distri$p[i]
      mu.1 <- distri$mu[i]
      sigma.1 <- distri$sigma[i]
      mu.2 <- - p/(1-p)*mu.1
      sigma.2 <- sqrt( 1/(1-p) * (1 - p*sigma.1^2 - p/(1-p)*mu.1^2) )
      B <- (runif(nb.sim)<p)*1
      eps.i <- B*rnorm(nb.sim,mean = mu.1,sd = sigma.1) +
        (1-B)*rnorm(nb.sim,mean = mu.2,sd = sigma.2)
    }else if(distri$type[i]=="student"){
      nu <- distri$df[i]
      eps.i <- rt(nb.sim,df = nu)/sqrt(nu/(nu-2))
    }else if(distri$type[i]=="laplace"){
      U <- runif(nb.sim) - .5
      b <- 1/sqrt(2)
      eps.i <- - b * sign(U) * log(1 - 2 * abs(U))
    }else if(distri$type[i]=="hyper.sec"){
      U <- runif(nb.sim)
      eps.i <- 2/pi * log(tan(pi/2*U))
    }
    eps <- cbind(eps,eps.i)
  }
  return(eps)
}


#' Variance Decomposition for SVAR Models
#'
#'This function produces a plot showing the variance decomposition
#'associated with a VAR model defined by autoregressive matrices Phi
#' and an impact matrix B.
#'
#' @details
#' The function first constructs impulse responses by simulating the VAR
#' (one unit shock at a time). Then, it computes the variance decomposition.
#' Each panel corresponds to one endogenous variable; within a panel, the stacked
#' colored bands sum to 1 at each horizon, showing the fraction of variance
#' attributable to each shock.
#'
#' @param Phi Array of autoregressive coefficient matrices.
#' @param B Numeric \eqn{n \times n} matrix; contemporaneous impact matrix
#'   mapping structural shocks \eqn{u_t}to the reduced-form residuals.
#' @param maxHorizon Integer; maximum horizon for the variance decomposition.
#' @param mfrow Optional length-2 integer vector passed to \code{par(mfrow = ...)}
#' to arrange multiple panels. Defaults internally to \code{c(1, n)} when not specified (i.e., when \code{mfrow = NA}).
#' @param names.var Optional character vector of length \eqn{n} with labels for variables.
#'   Defaults to \code{"Variable 1"}, \code{"Variable 2"}, … .
#' @param names.shock Optional character vector of length \eqn{n} with labels for shocks.
#'   Defaults to \code{"Shock 1"}, \code{"Shock 2"}, … .
#'
#' @returns the function output is a set of plots.
#' @examples
#' Phi <- array(NaN,c(2,2,2)) # (2,2,2) for (n,n,p)
#' Phi[,,1] <- matrix(c(.6,0,.2,.5),2,2)
#' Phi[,,2] <- matrix(c(-.1,.2,.1,.3),2,2)
#' B <- matrix(c(.5,-1,1.5,.8),2,2)
#' make.variance.decompo(Phi,B,maxHorizon=20)
#' @export
make.variance.decompo <- function(Phi, B, maxHorizon,
                                  mfrow = NA,
                                  names.var = NA,
                                  names.shock = NA) {
  # This function produces a plot showing the variance decomposition
  # associated with a VAR model defined by autoregressive matrices Phi
  # and an impact matrix B.

  n <- dim(B)[1]
  if (is.na(names.var[1])) names.var <- paste("Variable", 1:n)
  if (is.na(names.shock[1])) names.shock <- paste("Shock", 1:n)

  IRFs <- array(NA, c(n, n, maxHorizon))
  for (i in 1:n) {
    u.shock <- rep(0, n)
    u.shock[i] <- 1
    Y <- simul.VAR(c = matrix(0, n, 1), Phi = Phi, B = B,
                   nb.sim = maxHorizon, indic.IRF = 1, u.shock = u.shock)
    IRFs[, i, ] <- t(Y)
  }
  res <- variance.decomp(IRFs)$vardecomp

  # Handle layout flexibly
  if (is.na(mfrow[1])) {
    mfrow <- c(1, n)
  }
  par(mfrow = mfrow)

  # Adjust margins dynamically depending on layout
  nrows <- mfrow[1]
  ncols <- mfrow[2]

  par(plt=c(.15,.95,.2,.8))

  for (variable in 1:n) {
    res_variable <- res[variable, variable,,]
    res_variable <- apply(res_variable, 1, cumsum)
    res_variable <- rbind(0, res_variable)

    plot(1:maxHorizon, rep(0, maxHorizon),
         ylim = c(0, 1), col = "white", las = 1,
         xlab = "", ylab = "Variance share",
         main = names.var[variable])

    for (k in 1:n) {
      polygon(c(1:maxHorizon, maxHorizon:1),
              c(res_variable[k,], rev(res_variable[k + 1,])),
              col = k + 1, border = NA)
    }

    # Add custom "Horizon" label inside the plot
    usr <- par("usr")
    text(x = usr[2] - 0.3 * (usr[2] - usr[1]),
         y = usr[3] + 0.07 * (usr[4] - usr[3]),
         labels = "Horizon",
         adj = c(1, 0), xpd = TRUE)

    # Legend adjustment: fit better inside multi-plot layouts
    legend("topleft", horiz = FALSE, bty = "n",
           pch = 22,pt.bg = 2:(n + 1),
           col = "black", legend = names.shock,
           xpd = TRUE, pt.cex = 1.2, cex = 1)
  }

  invisible(NULL)
}

