# Newey-West covariance estimation and two-stage least squares.

#' Newey–West Long-Run Variance Estimator
#'
#' Computes the Newey–West estimator of the long-run variance
#' (or long-run covariance matrix) of a univariate or
#' multivariate time series.
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
NW.LongRunVariance <- function(X, q) {
  gamma0 <- autocov(X, 0)
  LRV <- gamma0
  if (q > 0) {
    weights <- 1 - (1:q) / (q + 1)
    for (i in 1:q) {
      LRV <- LRV + weights[i] * (autocov(X, i) + t(autocov(X, i)))
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
#' Z <- cbind(1, rnorm(T))
#' X <- cbind(1, Z[, 2] + 0.5 * rnorm(T))
#' Y <- X %*% c(1, 2) + rnorm(T)
#'
#' res <- tsls(Y = as.matrix(Y), X = X, Z = Z, q = 4)
#' res$b.iv
#' res$covmat.b.iv
#'
#' @export
tsls <- function(Y, X, Z, q) {
  # Regress y on x, using z as an instrument.
  # Y, X and Z have to be matrices
  # The covariance matrix of parameter estimates is corrected for HAC (Newey-West)

  if (!is.matrix(Y)) {
    print("Error: Y should be a matrix.")
    return(0)
  }
  if (!is.matrix(X)) {
    print("Error: X should be a matrix.")
    return(0)
  }
  if (!is.matrix(Z)) {
    print("Error: Z should be a matrix.")
    return(0)
  }

  # Determine sample size:
  T <- length(Y)

  # Compute the IV estimator:
  PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  Y.hat <- PZ %*% Y
  X.hat <- PZ %*% X
  PX.hat <- solve(t(X.hat) %*% X.hat) %*% t(X.hat)
  b.iv <- PX.hat %*% Y.hat

  # Residual estimates:
  eps <- Y - X %*% b.iv

  # Approximation of the covariance matrix of b.iv:
  Q <- T * solve(t(X) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*%
    t(solve(t(Z) %*% Z) %*% t(Z) %*% X)

  Z.times.eps <- Z * matrix(eps, dim(Z)[1], dim(Z)[2])
  S <- NW.LongRunVariance(Z.times.eps, q)

  covmat.b.iv <- 1 / T * Q %*% S %*% t(Q)

  return(list(
    b.iv = b.iv,
    covmat.b.iv = covmat.b.iv
  ))
}

#' @keywords internal
#' @noRd
autocov <- function(X, n) {
  if (!is.matrix(X)) {
    T <- length(X)
    X.1 <- X[1:(T - n)] - mean(X)
    X.2 <- X[(n + 1):T] - mean(X)
    return(1 / T * sum(X.1 * X.2))
  } else {
    T <- dim(X)[1]
    k <- dim(X)[2]
    vec.1 <- matrix(1, 1, k)
    mean.X <- apply(X, 2, mean)
    X.1 <- X[1:(T - n), ] - t(matrix(mean.X, k, T - n))
    X.2 <- X[(n + 1):T, ] - t(matrix(mean.X, k, T - n))
    return(
      matrix(1 / T * apply((X.1 %x% vec.1) * (vec.1 %x% X.2), 2, sum), k, k)
    )
  }
}
