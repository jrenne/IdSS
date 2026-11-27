# ICA log densities, scores and second derivatives.

#' @keywords internal
#' @noRd
ica_log_gaussian <- function(x, mu, sigma, indic.deriv = 0) {
  # Gaussian distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- -(x - mu)^2 / (2 * sigma^2)
  if (indic.deriv == 1) {
    d.log.g <- -(x - mu) / sigma^2
    d2.log.g <- -1 / sigma^2 * rep(1, length(x))
  } else {
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g = log.g, d.log.g = d.log.g, d2.log.g = d2.log.g))
}

#' @keywords internal
#' @noRd
ica_log_student <- function(x, nu, indic.deriv = 0) {
  # Student distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- -(1 + nu) / 2 * log(1 + (x * sqrt(nu / (nu - 2)))^2 / nu)
  if (indic.deriv == 1) {
    d.log.g <- -x * (1 + nu) / (nu - 2 + x^2)
    d2.log.g <- -(1 + nu) * (nu - 2 - x^2) / (nu - 2 + x^2)^2
  } else {
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g = log.g, d.log.g = d.log.g, d2.log.g = d2.log.g))
}

#' @keywords internal
#' @noRd
ica_log_laplace <- function(x, indic.deriv = 0) {
  # Laplace distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  b <- 1 / sqrt(2)
  log.g <- -abs(x) / b
  if (indic.deriv == 1) {
    d.log.g <- -sign(x) / b
    d2.log.g <- 0
  } else {
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g = log.g, d.log.g = d.log.g, d2.log.g = d2.log.g))
}

#' @keywords internal
#' @noRd
ica_log_hyperbolic_secant <- function(x, indic.deriv = 0) {
  # Hyperbolic secant distribution
  # if indic.deriv = 1 -> computes the derivatives (1st and 2nd) of log.g w.r.t. x
  log.g <- -log(exp(pi / 2 * x) + exp(-pi / 2 * x))
  if (indic.deriv == 1) {
    d.log.g <- -(1 / (exp(pi / 2 * x) + exp(-pi / 2 * x)) *
      (exp(pi / 2 * x) * pi / 2 + exp(-pi / 2 * x) * (-pi / 2)))
    d2.log.g <- ((-pi^2 * exp(-1 / 2 * pi * x)^2) / (exp(-1 / 2 * pi * x)^4 + 2 * exp(-1 / 2 * pi * x)^2 + 1))
  } else {
    d.log.g <- NULL
    d2.log.g <- NULL
  }
  return(list(log.g = log.g, d.log.g = d.log.g, d2.log.g = d2.log.g))
}

#' @keywords internal
#' @noRd
ica_log_gaussian_mixture <- function(x, mu, sigma, p, indic.deriv = 0) {
  density_1 <- p * dnorm(x, mean = mu[1], sd = sigma[1])
  density_2 <- (1 - p) * dnorm(x, mean = mu[2], sd = sigma[2])
  density <- density_1 + density_2
  log.g <- log(density)
  d.log.g <- d2.log.g <- NULL

  if (indic.deriv == 1) {
    weight <- density_1 / density
    score_1 <- -(x - mu[1]) / sigma[1]^2
    score_2 <- -(x - mu[2]) / sigma[2]^2
    d.log.g <- weight * score_1 + (1 - weight) * score_2
    # E[component curvature | x] + Var[component score | x].
    d2.log.g <- -weight / sigma[1]^2 - (1 - weight) / sigma[2]^2 +
      weight * (1 - weight) * (score_1 - score_2)^2
  }
  list(log.g = log.g, d.log.g = d.log.g, d2.log.g = d2.log.g)
}
