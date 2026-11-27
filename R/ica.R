# Independent-component identification and the Cayley parametrization.

#' Estimate an SVAR model identified by Independent Component Analysis (ICA)
#'
#' This function estimates a Structural Vector Autoregression (SVAR) using
#' Independent Component Analysis (ICA) for shock identification.
#' The reduced-form VAR residuals are first pre-whitened, and the structural
#' impact matrix is then identified by maximizing a pseudo-maximum likelihood
#' criterion under non-Gaussianity assumptions on the structural shocks.
#'
#' The identification follows a parametrization based on a skew-symmetric
#' matrix, ensuring orthogonality of the resulting transformation.
#'
#' @param endog.var A \eqn{T \times n} numeric matrix of endogenous variables,
#' ordered consistently with the VAR specification.
#'
#' @param distri A list specifying the marginal distributions of the
#' structural shocks, as accepted by [simul.distri()] (for example,
#' `list(type = c("student", "student"), df = c(5, 7))`).
#'
#' @param p An integer specifying the number of lags in the VAR model.
#'
#' @details
#' The procedure consists of the following steps:
#' \enumerate{
#'   \item Estimation of a reduced-form VAR(\eqn{p}) equation-by-equation
#'   using OLS.
#'   \item Construction of the reduced-form residual matrix.
#'   \item Pre-whitening of residuals using the Cholesky decomposition of
#'   their covariance matrix.
#'   \item ICA-based identification via optimization of a pseudo-likelihood
#'   objective function under the assumed non-Gaussian distribution.
#'   \item Recovery of the contemporaneous impact matrix and asymptotic
#'   covariance matrix of the ICA parameters.
#' }
#'
#' The ICA transformation matrix is parameterized using a skew-symmetric matrix
#' whose free parameters are stacked in a vector of dimension \eqn{n(n-1)/2}.
#'
#' @return A list containing:
#' \describe{
#'   \item{C.PML}{An \eqn{n \times n} ICA rotation matrix estimated by
#'   pseudo-maximum likelihood.}
#'   \item{B.PML}{The estimated contemporaneous impact matrix of structural
#'   shocks.}
#'   \item{eta.PML}{A \eqn{T-p \times n} matrix of estimated structural shocks.}
#'   \item{V}{The asymptotic covariance matrix of the ICA parameters.}
#' }
#'
#' @examples
#' First.date <- "1959-04-01"
#' Last.date <- "2015-01-01"
#' data <- US3var
#' data <- data[(data$Date >= First.date) & (data$Date <= Last.date), ]
#' Y <- as.matrix(data[c("infl", "y.gdp.gap", "r")])
#' distri <- list(
#'   type = c("mixt.gaussian", "mixt.gaussian", "mixt.gaussian"),
#'   df = c(NaN, NaN, NaN),
#'   p = c(0.5, .5, .5), mu = c(.1, .1, .1), sigma = c(.5, .7, 1.3)
#' )
#' ICA.res.no.commo <- estim.SVAR.ICA(Y, distri = distri, p = 6)
#' @export
estim.SVAR.ICA <- function(endog.var, distri, p) {
  n <- dim(endog.var)[2]
  T <- dim(endog.var)[1]

  YY <- NULL
  for (i in 1:p) {
    YY <- cbind(YY, endog.var[(1 + p - i):(T - i), ])
  }

  eps <- NULL
  for (j in 1:n) {
    eq <- lm(endog.var[(1 + p):T, j] ~ YY)
    eps <- cbind(eps, eq$residuals)
  }

  # Pre-whiten Y:
  P <- t(chol(var(eps)))
  Y <- eps %*% solve(t(P))

  AA.0 <- rep(0, n * (n - 1) / 2)
  res.optim <- optim(AA.0, func.2.minimize,
    Y = Y,
    distri = distri,
    gr = d.func.2.minimize,
    method = "Nelder-Mead",
    control = list(trace = FALSE, maxit = 500)
  )
  AA.0 <- res.optim$par
  res.optim <- optim(AA.0, func.2.minimize,
    Y = Y,
    distri = distri,
    gr = d.func.2.minimize,
    method = "BFGS",
    control = list(trace = FALSE)
  )
  AA.0 <- res.optim$par
  res.optim <- optim(AA.0, func.2.minimize,
    Y = Y,
    distri = distri,
    gr = d.func.2.minimize,
    method = "Nelder-Mead",
    control = list(trace = FALSE, maxit = 500)
  )
  AA.0 <- res.optim$par

  AA.est <- AA.0
  M <- make.M(n)
  A.est <- matrix(M %*% AA.est, n, n)
  C.PML <- (diag(n) + A.est) %*% solve(diag(n) - A.est)

  B.PML <- P %*% C.PML

  eta.PML <- Y %*% C.PML
  V <- make.Asympt.Cov.delta(eta.PML, distri, C.PML)

  return(list(
    C.PML = C.PML,
    B.PML = B.PML,
    eta.PML = eta.PML,
    V = V
  ))
}

#' Construct the ICA Skew-Symmetric Basis
#'
#' Constructs the linear map used to transform the free lower-triangular
#' parameters of a skew-symmetric matrix into its vectorized form.
#'
#' @param n Positive integer giving the number of structural shocks.
#' @return A numeric matrix with `n^2` rows and `n * (n - 1) / 2` columns.
#' @examples
#' M <- make.M(3)
#' dim(M)
#' @export
make.M <- function(n) {
  # function that builds M_n, the n2 x ((n-1)n/2) matrix that is such that Vec(A) = M_n a where a is the vectorized entries of the
  # lower triangular matrix of A and A is such that t(A) = -A.
  aux1 <- 1:(n^2)
  aux2 <- matrix(0, n, n)
  aux2[lower.tri(aux2)] <- 1:((n - 1) * n / 2)
  aux2 <- matrix(aux2, ncol = 1)
  aux3 <- matrix(0, n, n)
  aux3[lower.tri(aux3)] <- 1:((n - 1) * n / 2)
  aux3 <- t(aux3)
  aux3 <- matrix(aux3, ncol = 1)

  M <- matrix(0, n^2, n * (n - 1) / 2)
  indic.row <- which(aux2 > 0)
  M[indic.row, ] <- diag(n * (n - 1) / 2)
  indic.row <- which(aux3 > 0)
  M[indic.row, ] <- -diag(n * (n - 1) / 2)
  return(M)
}

#' ICA Score Contributions
#'
#' Evaluates marginal log pseudo-densities and, optionally, their first two
#' derivatives with respect to the residuals.
#'
#' @param eps Numeric matrix whose rows are observations and columns are shocks.
#' @param distri Distribution specification list accepted by
#'   [simul.distri()].
#' @param indic.deriv Integer flag. Zero returns score contributions; a non-zero
#'   value requests the corresponding derivative contribution.
#' @return A list containing `log.L`, a matrix with the dimensions of `eps`,
#' and `d.log.L` and `d2.log.L`, the corresponding derivative matrices.
#' The derivative components are `NULL` when `indic.deriv = 0`.
#' @examples
#' eps <- matrix(rnorm(20), 10, 2)
#' distri <- list(
#'   type = c("gaussian", "gaussian"),
#'   mu = c(0, 0), sigma = c(1, 1),
#'   df = c(NA, NA), p = c(NA, NA)
#' )
#' gICA(eps, distri)
#' @export
gICA <- function(eps, distri, indic.deriv = 0) {
  # Y is of dimension T x n
  # distri is a list:
  #       - distri$type contains the types of distributions ("student" or "gaussian")
  #       - distri$df contains the degrees of freedom for student distri (NaN for non student distri)
  n <- dim(eps)[2]
  log.L <- NULL
  d.log.L <- NULL
  d2.log.L <- NULL
  for (ii in 1:n) {
    if (distri$type[ii] == "student") {
      g.aux <- ica_log_student(eps[, ii], nu = distri$df[ii], indic.deriv)
    } else if (distri$type[ii] == "gaussian") {
      g.aux <- ica_log_gaussian(eps[, ii], mu = 0, sigma = 1, indic.deriv)
    } else if (distri$type[ii] == "mixt.gaussian") {
      p <- distri$p[ii]
      mu.1 <- distri$mu[ii]
      sigma.1 <- distri$sigma[ii]
      mu.2 <- -p / (1 - p) * mu.1
      sigma.2 <- sqrt(1 / (1 - p) * (1 - p * sigma.1^2 - p / (1 - p) * mu.1^2))
      g.aux <- ica_log_gaussian_mixture(eps[, ii],
        mu = c(mu.1, mu.2), sigma = c(sigma.1, sigma.2),
        p, indic.deriv
      )
    } else if (distri$type[ii] == "laplace") {
      g.aux <- ica_log_laplace(eps[, ii], indic.deriv)
    } else if (distri$type[ii] == "hyper.sec") {
      g.aux <- ica_log_hyperbolic_secant(eps[, ii], indic.deriv)
    }
    log.L <- cbind(log.L, g.aux$log.g)
    d.log.L <- cbind(d.log.L, g.aux$d.log.g)
    d2.log.L <- cbind(d2.log.L, g.aux$d2.log.g)
  }
  return(list(log.L = log.L, d.log.L = d.log.L, d2.log.L = d2.log.L))
}

#' ICA pseudo-likelihood objective
#'
#' Evaluates the negative pseudo-log-likelihood used to estimate the ICA
#' rotation in a structural VAR. It is exposed because the Identification
#' course uses it directly with [stats::optim()]. Most users will normally
#' call [estim.SVAR.ICA()] instead.
#'
#' @param theta Numeric vector containing the free parameters of the
#'   skew-symmetric ICA parametrization.
#' @param Y Numeric matrix whose rows are pre-whitened reduced-form residuals.
#' @param distri List describing the marginal pseudo-distributions; see
#'   [gICA()] for the supported components.
#'
#' @return A numeric scalar equal to the negative pseudo-log-likelihood.
#' @seealso [d.func.2.minimize()], [estim.SVAR.ICA()]
#' @examples
#' set.seed(123)
#' errors <- cbind(rt(160, df = 5), rt(160, df = 7))
#' # Center and whiten the two series before evaluating an ICA rotation.
#' Y <- scale(errors, center = TRUE, scale = FALSE) %*% solve(chol(cov(errors)))
#' distri <- list(type = c("student", "student"), df = c(5, 7))
#' objective <- func.2.minimize(theta = 0.1, Y = Y, distri = distri)
#' objective
#' fit <- optim(0, func.2.minimize,
#'   gr = d.func.2.minimize,
#'   Y = Y, distri = distri, method = "BFGS"
#' )
#' fit$convergence
#' @export
func.2.minimize <- function(theta, Y, distri) {
  aux <- pseudo.log.L(Y, theta, distri)
  return(-aux$log.L)
}

#' Gradient of the ICA pseudo-likelihood objective
#'
#' Computes the gradient of [func.2.minimize()] with respect to the free ICA
#' parameters. It can be supplied as the `gr` argument of [stats::optim()].
#'
#' @inheritParams func.2.minimize
#'
#' @return A numeric vector with the same length as `theta`.
#' @seealso [func.2.minimize()], [estim.SVAR.ICA()]
#' @examples
#' set.seed(123)
#' errors <- cbind(rt(160, df = 5), rt(160, df = 7))
#' # Center and whiten the two series before evaluating an ICA rotation.
#' Y <- scale(errors, center = TRUE, scale = FALSE) %*% solve(chol(cov(errors)))
#' distri <- list(type = c("student", "student"), df = c(5, 7))
#' theta <- 0.1
#' gradient <- d.func.2.minimize(theta, Y, distri)
#' # Compare the analytic gradient with a central finite difference.
#' h <- 1e-6
#' finite_difference <- (func.2.minimize(theta + h, Y, distri) -
#'   func.2.minimize(theta - h, Y, distri)) / (2 * h)
#' c(analytic = unname(gradient), finite_difference = finite_difference)
#' stopifnot(isTRUE(all.equal(unname(gradient), finite_difference,
#'   tolerance = 1e-5
#' )))
#' @export
d.func.2.minimize <- function(theta, Y, distri) {
  aux <- pseudo.log.L(Y, theta, distri, indic.Jacobian = 1)
  return(-aux$Jacobian)
}

#' @keywords internal
#' @noRd
make.C <- function(AA, n) {
  A <- matrix(0, n, n)
  A[lower.tri(A)] <- AA
  A[t(lower.tri(A))] <- -AA
  if (prod(A == (0 * A)) == 1) {
    C <- diag(n)
  } else {
    C <- (diag(n) + A) %*% solve(diag(n) - A)
  }
  return(C)
}

#' @keywords internal
#' @noRd
pseudo.log.L <- function(Y, AA, distri, indic.Jacobian = 0) {
  # AA is \mathcal{A}
  # Y is the matrix of observations
  # g is a likelihood function
  n <- ncol(Y)
  M <- make.M(n)
  I <- diag(n)
  A <- matrix(M %*% AA, n, n)
  C <- make.C(AA, n)
  eps <- Y %*% C
  g.eps <- gICA(eps, distri, indic.deriv = indic.Jacobian)
  vector.log.L <- g.eps$log.L
  log.L <- sum(vector.log.L)
  if (indic.Jacobian == 0) {
    Jacobian <- NULL
  } else {
    # C = (I + A) (I - A)^(-1) = 2 (I - A)^(-1) - I.
    # Hence d vec(C) = 2 ((I-A)^(-T) %x% (I-A)^(-1)) M d theta.
    inverse <- solve(I - A)
    aux1 <- 2 * (t(inverse) %x% inverse) %*% M
    d.log.L <- g.eps$d.log.L # This is a T x n matrix
    aux2.1 <- matrix(1, 1, n) %x% Y
    aux2.2 <- d.log.L %x% matrix(1, 1, n)
    aux2 <- aux2.1 * aux2.2
    Jacobian <- apply(aux2 %*% aux1, 2, sum)
  }
  return(list(log.L = log.L, Jacobian = Jacobian))
}
