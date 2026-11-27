test_that("Gaussian-mixture scores and curvature match numerical derivatives", {
  x <- seq(-2, 2, length.out = 19)
  mu <- c(-0.3, 0.4)
  sigma <- c(0.4, 1.2)
  h <- 1e-4
  for (p in c(0.1, 0.5, 0.9)) {
    log_density <- function(z) {
      log(p * dnorm(z, mu[1], sigma[1]) + (1 - p) * dnorm(z, mu[2], sigma[2]))
    }
    actual <- IdSS:::ica_log_gaussian_mixture(x, mu, sigma, p, indic.deriv = 1)
    score <- (log_density(x + h) - log_density(x - h)) / (2 * h)
    curvature <- (log_density(x + h) - 2 * log_density(x) + log_density(x - h)) / h^2
    expect_equal(actual$log.g, log_density(x))
    expect_equal(actual$d.log.g, score, tolerance = 1e-6)
    expect_equal(actual$d2.log.g, curvature, tolerance = 1e-5)
  }
})

test_that("a single Gaussian component has its usual score and curvature", {
  x <- seq(-1, 1, length.out = 9)
  for (p in c(0, 1)) {
    mu <- c(-0.3, 0.4)
    sigma <- c(0.4, 1.2)
    i <- if (p == 1) 1L else 2L
    actual <- IdSS:::ica_log_gaussian_mixture(x, mu, sigma, p, indic.deriv = 1)
    expect_equal(actual$log.g, dnorm(x, mu[i], sigma[i], log = TRUE))
    expect_equal(actual$d.log.g, -(x - mu[i]) / sigma[i]^2)
    expect_equal(actual$d2.log.g, rep(-1 / sigma[i]^2, length(x)))
  }
})
