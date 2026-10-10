#' Simulate an ARMA or ARMAX Process
#'
#' Simulates one or more paths from an ARMA model, optionally with exogenous
#' regressors, or computes the response to a one-standard-deviation innovation.
#' The function retains the interface historically used in the companion
#' courses.
#'
#' @param c Numeric scalar. Intercept in the observation equation.
#' @param phi Numeric vector of autoregressive coefficients. For an AR(0)
#'   model, use `0` and set `y.0 = 0`.
#' @param theta Numeric vector whose first entry multiplies the contemporaneous
#'   innovation. For an AR model without an MA component, use `theta = 1`.
#' @param sigma Non-negative numeric scalar giving the innovation standard
#'   deviation. In impulse-response mode this is the size of the initial shock.
#' @param T Positive integer giving the number of simulated periods or response
#'   horizons.
#' @param y.0 Numeric vector of initial observations, with the same length as
#'   `phi`, ordered from the most recent lag to the oldest lag.
#' @param nb.sim Positive integer giving the number of simulated paths.
#' @param make.IRF Numeric flag. If zero (the default), simulate random paths.
#'   Otherwise, set all innovations to zero except for an initial shock of size
#'   `sigma`.
#' @param X Optional numeric vector or matrix with `T` rows containing exogenous
#'   regressors. The historical sentinel `NaN` denotes no regressors.
#' @param beta Optional coefficient vector conformable with `X`. The historical
#'   sentinel `NaN` is retained for backwards compatibility.
#'
#' @return A numeric matrix with `T` rows and `nb.sim` columns. Each column is a
#'   simulated path or, when `make.IRF != 0`, an impulse response.
#'
#' @details The simulated model is
#'   \deqn{y_t = c + \sum_{i=1}^p \phi_i y_{t-i}
#'   + \sum_{j=0}^q \theta_j \varepsilon_{t-j} + X_t\beta,}
#'   where the innovations have standard deviation `sigma`.
#'
#' @examples
#' set.seed(123)
#' paths <- sim.arma(c = 0, phi = 0.7, theta = 1, sigma = 1,
#'                   T = 20, y.0 = 0, nb.sim = 2)
#' irf <- sim.arma(c = 0, phi = 0.7, theta = 1, sigma = 1,
#'                 T = 12, y.0 = 0, make.IRF = 1)
#' @export
sim.arma <- function(c, phi, theta, sigma, T, y.0, nb.sim = 1,
                     make.IRF = 0, X = NaN, beta = NaN) {
  p <- length(phi)
  q <- length(theta)
  if (length(y.0) != p) {
    stop("`y.0` must have the same length as `phi`.", call. = FALSE)
  }
  if (length(T) != 1L || !is.finite(T) || T < 1 || T != as.integer(T)) {
    stop("`T` must be a positive integer.", call. = FALSE)
  }
  if (length(nb.sim) != 1L || !is.finite(nb.sim) || nb.sim < 1 ||
      nb.sim != as.integer(nb.sim)) {
    stop("`nb.sim` must be a positive integer.", call. = FALSE)
  }

  has_x <- !(length(X) == 1L && is.na(X[1]))
  if (has_x) {
    X <- as.matrix(X)
    if (nrow(X) != T) {
      stop("`X` must have `T` rows.", call. = FALSE)
    }
    if (length(beta) != ncol(X)) {
      stop("`beta` must have one entry per column of `X`.", call. = FALSE)
    }
    X.beta <- as.numeric(X %*% beta)
  } else {
    X.beta <- rep(0, T)
  }

  if (make.IRF == 0) {
    eps <- sigma * matrix(rnorm(nb.sim * T), T, nb.sim)
  } else {
    eps <- matrix(0, T, nb.sim)
    eps[1, ] <- sigma
  }
  eps.t <- matrix(0, q, nb.sim)
  F <- make_ar_companion(phi)
  Y <- matrix(NA_real_, T, nb.sim)
  y.t <- matrix(y.0, p, nb.sim)

  for (t in seq_len(T)) {
    if (q > 1) {
      eps.t <- rbind(eps[t, ], matrix(eps.t[seq_len(q - 1), ],
                                     q - 1, nb.sim))
    } else {
      eps.t <- matrix(eps[t, ], nrow = 1)
    }
    innovation <- matrix(theta, nrow = 1) %*% eps.t
    innovation <- rbind(innovation, matrix(0, p - 1, nb.sim))
    y.t <- c(c, rep(0, p - 1)) + F %*% y.t + innovation + X.beta[t]
    Y[t, ] <- y.t[1, ]
  }
  Y
}

#' Estimate an ARMAX Model by Conditional Gaussian Maximum Likelihood
#'
#' Fits an ARMA(`p`, `q`) model, optionally augmented with contemporaneous
#' exogenous regressors. Starting values are obtained from an autoregression and
#' the conditional Gaussian likelihood is optimized numerically.
#'
#' @param Y Numeric vector containing the dependent time series.
#' @param p Non-negative integer giving the autoregressive order.
#' @param q Non-negative integer giving the moving-average order.
#' @param X Optional numeric vector or matrix of exogenous regressors with one
#'   row per observation in `Y`. The historical sentinel `NaN` denotes an ARMA
#'   model without exogenous regressors.
#'
#' @return A list with components `c`, `phi`, `theta`, `sigma`, and `beta` for
#'   the estimated coefficients; `max.logl` for the maximized conditional
#'   log-likelihood; `I` for the inverse observed Hessian; `st.dev` for the
#'   associated standard errors; and `THETA` for the complete parameter vector.
#'
#' @details The first element of `theta` is fixed at one. The likelihood is
#'   conditional on zero pre-sample innovations. This routine is primarily
#'   intended to reproduce the instructional examples in the companion books;
#'   for general-purpose ARIMA estimation, see [stats::arima()].
#'
#' @examples
#' set.seed(42)
#' y <- as.numeric(sim.arma(0.2, 0.5, 1, 0.4, 80, 0))
#' fit <- estim.armax(y, p = 1, q = 0)
#' fit$phi
#' @export
estim.armax <- function(Y, p, q, X = NaN) {
  if (!is.numeric(Y) || anyNA(Y) || length(Y) <= p) {
    stop("`Y` must be a complete numeric vector longer than `p`.", call. = FALSE)
  }
  if (p < 0 || q < 0 || p != as.integer(p) || q != as.integer(q)) {
    stop("`p` and `q` must be non-negative integers.", call. = FALSE)
  }
  T <- length(Y)
  has_x <- !(length(X) == 1L && is.na(X[1]))
  beta.0 <- NULL

  if (has_x) {
    X <- as.matrix(X)
    if (nrow(X) != T) stop("`X` must have one row per observation in `Y`.")
  }

  if (p >= 1) {
    y <- Y[(p + 1):T]
    Y_1 <- embed(Y, p + 1)[, -1, drop = FALSE]
    if (has_x) {
      eq <- lm(y ~ Y_1 + X[(p + 1):T, , drop = FALSE])
      beta.0 <- tail(coef(eq), ncol(X))
    } else {
      eq <- lm(y ~ Y_1)
    }
    c.0 <- coef(eq)[1]
    phi.0 <- coef(eq)[seq_len(p) + 1]
    sigma.0 <- sd(residuals(eq))
  } else if (has_x) {
    eq <- lm(Y ~ X)
    c.0 <- coef(eq)[1]
    phi.0 <- numeric(0)
    beta.0 <- coef(eq)[seq_len(ncol(X)) + 1]
    sigma.0 <- sd(residuals(eq))
  } else {
    c.0 <- mean(Y)
    phi.0 <- numeric(0)
    sigma.0 <- sd(Y)
  }

  THETA.0 <- c(c.0, phi.0, rep(0, q), sigma.0, beta.0)
  maxit_nm <- 100 * length(THETA.0)
  result <- NULL
  for (i in seq_len(4)) {
    result <- optim(par = THETA.0, fn = armax_objective,
                    Y = Y, p = p, q = q, X = X,
                    method = "Nelder-Mead", control = list(maxit = maxit_nm))
    THETA.0 <- result$par
    result <- optim(par = THETA.0, fn = armax_objective,
                    Y = Y, p = p, q = q, X = X,
                    method = "BFGS", control = list(maxit = 10),
                    hessian = i == 4)
    THETA.0 <- result$par
  }

  THETA <- result$par
  information_inverse <- tryCatch(solve(result$hessian),
                                  error = function(e) matrix(NA_real_,
                                    length(THETA), length(THETA)))
  st.dev <- sqrt(diag(information_inverse))
  intercept <- THETA[1]
  phi <- if (p > 0) THETA[seq_len(p) + 1] else numeric(0)
  theta <- if (q > 0) c(1, THETA[seq_len(q) + p + 1]) else 1
  sigma <- abs(THETA[p + q + 2])
  beta <- if (has_x) tail(THETA, ncol(X)) else NULL

  list(c = intercept, phi = phi, theta = theta, sigma = sigma, beta = beta,
       max.logl = -result$value, I = information_inverse,
       st.dev = st.dev, THETA = THETA)
}

#' Evaluate an ARCH/GARCH Conditional Variance Recursion
#'
#' Computes the conditional variance sequence and the negative Gaussian
#' log-likelihood for an ARCH(`m0`) or GARCH(`m0`, `r0`) specification.
#'
#' @param theta Numeric parameter vector ordered as `zeta`, followed by `m0`
#'   ARCH coefficients and `r0` GARCH coefficients. Absolute values are used for
#'   these parameters, preserving the historical course implementation.
#' @param x Numeric vector of mean-zero observations or residuals.
#' @param m0 Non-negative integer giving the number of ARCH terms.
#' @param r0 Non-negative integer giving the number of lagged conditional
#'   variances.
#'
#' @return A list with `logf`, the negative Gaussian log-likelihood (or a large
#'   penalty for a non-stationary parameterization); `h`, the conditional
#'   variance sequence padded with initial `NaN` values; and `names.param`, the
#'   parameter labels.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(100)
#' arch <- compute.garch(c(1, 0.1), x, m0 = 1, r0 = 0)
#' length(arch$h)
#' @export
compute.garch <- function(theta, x, m0, r0) {
  if (!is.numeric(x) || anyNA(x)) stop("`x` must be a complete numeric vector.")
  if (m0 < 0 || r0 < 0 || m0 != as.integer(m0) || r0 != as.integer(r0)) {
    stop("`m0` and `r0` must be non-negative integers.", call. = FALSE)
  }
  expected <- 1 + m0 + r0
  if (length(theta) != expected) {
    stop("`theta` must contain 1 + m0 + r0 parameters.", call. = FALSE)
  }
  zeta <- abs(theta[1])
  alpha <- if (m0 > 0) abs(theta[seq_len(m0) + 1]) else 0
  delta <- if (r0 > 0) abs(theta[seq_len(r0) + m0 + 1]) else 0
  m <- max(1L, m0)
  r <- max(1L, r0)
  p <- max(m, r)
  if (length(x) <= p) stop("`x` is too short for the requested lag orders.")

  kappa <- zeta * (1 - sum(delta))
  persistence <- sum(alpha) + sum(delta)
  if (!is.finite(kappa) || persistence >= 1) {
    return(list(logf = 1e20 * max(persistence - 1, .Machine$double.eps),
                h = rep(NaN, length(x)),
                names.param = garch_parameter_names(m0, r0)))
  }
  Eu2 <- kappa / (1 - persistence)
  vec.u2 <- if (m0 > 0) x[p:(p - m0 + 1)]^2 else 0
  vec.h <- rep(Eu2, r)
  h <- numeric(length(x) - p)

  for (i in (p + 1):length(x)) {
    last.h <- kappa + sum(alpha * vec.u2) + sum(delta * vec.h)
    h[i - p] <- last.h
    vec.u2 <- if (m0 > 1) c(x[i]^2, vec.u2[seq_len(m0 - 1)]) else x[i]^2
    vec.h <- if (r0 > 1) c(last.h, vec.h[seq_len(r0 - 1)]) else last.h
  }
  if (any(!is.finite(h)) || any(h <= 0)) {
    logf <- 1e20
  } else {
    logf <- -sum(-0.5 * log(2 * pi * h) - x[(p + 1):length(x)]^2 / (2 * h))
  }
  list(logf = logf, h = c(rep(NaN, p), h),
       names.param = garch_parameter_names(m0, r0))
}

make_ar_companion <- function(phi) {
  p <- length(phi)
  F <- matrix(0, p, p)
  F[1, ] <- phi
  if (p > 1) F[2:p, seq_len(p - 1)] <- diag(p - 1)
  F
}

armax_objective <- function(THETA, Y, p, q, X = NaN) {
  intercept <- THETA[1]
  phi <- if (p > 0) THETA[seq_len(p) + 1] else numeric(0)
  theta <- if (q > 0) c(1, THETA[seq_len(q) + p + 1]) else 1
  sigma <- abs(THETA[p + q + 2])
  has_x <- !(length(X) == 1L && is.na(X[1]))
  beta <- if (has_x) tail(THETA, ncol(as.matrix(X))) else NULL
  if (!is.finite(sigma) || sigma <= .Machine$double.eps) return(.Machine$double.xmax)
  residuals <- armax_residuals(Y, intercept, phi, theta, X, beta)
  -sum(-0.5 * log(2 * pi * sigma^2) - residuals^2 / (2 * sigma^2))
}

armax_residuals <- function(Y, intercept, phi, theta, X = NaN, beta = NULL) {
  T <- length(Y)
  p <- length(phi)
  q <- length(theta) - 1
  has_x <- !(length(X) == 1L && is.na(X[1]))
  X.beta <- if (has_x) as.numeric(as.matrix(X) %*% beta) else rep(0, T)
  eps <- numeric(T)
  for (t in seq_len(T)) {
    ar_lags <- if (p > 0 && t > 1) {
      available <- seq_len(min(p, t - 1))
      sum(phi[available] * Y[t - available])
    } else 0
    ma_lags <- if (q > 0 && t > 1) {
      available <- seq_len(min(q, t - 1))
      sum(theta[available + 1] * eps[t - available])
    } else 0
    eps[t] <- Y[t] - intercept - ar_lags - ma_lags - X.beta[t]
  }
  eps
}

garch_parameter_names <- function(m0, r0) {
  c("zeta", if (m0 > 0) paste0("alpha", seq_len(m0)),
    if (r0 > 0) paste0("delta", seq_len(r0)))
}
