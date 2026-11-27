test_that("Independent regime-switching identities hold", {
  # Independent checks for the chapter's filter, smoother, EM and forecasts.
  near <- function(a, b, tol = 1e-9) expect_equal(as.numeric(a), as.numeric(b), tolerance = tol)
  P <- matrix(c(.92, .08, .13, .87), 2, byrow = TRUE)
  c0 <- c(-.6, 1.1)
  phi <- c(.3, -.2)
  sigma <- c(.7, 1.2)
  q <- ms_stationary(P)
  near(as.vector(crossprod(P, q)), q)
  near(q, c(.13, .08) / .21)
  # Enumerate every state path to validate filtering, smoothing, likelihood and
  # expected transition counts independently of the recursion (2^6 paths).
  y <- c(-.8, .2, 1.4, .7, -.5, -.1)
  y0 <- .15
  prior <- c(.4, .6)
  paths <- as.matrix(expand.grid(rep(list(1:2), length(y))))
  path_weights <- function(values) {
    n <- length(values)
    mat <- as.matrix(expand.grid(rep(list(1:2), n)))
    weights <- apply(mat, 1, function(s) {
      density <- prod(dnorm(values, c0[s] + phi[s] * c(y0, head(values, -1)), sigma[s]))
      prob <- prior[s[1]]
      if (n > 1) for (t in 2:n) prob <- prob * P[s[t - 1], s[t]]
      density * prob
    })
    list(paths = mat, weight = weights)
  }
  f <- ms_filter(y, c0, phi, sigma, P, y0, prior)
  for (t in seq_along(y)) {
    e <- path_weights(y[1:t])
    w <- e$weight / sum(e$weight)
    near(f$filtered[t, ], sapply(1:2, function(j) sum(w[e$paths[, t] == j])))
    near(sum(f$increments[1:t]), log(sum(e$weight)))
  }
  e <- path_weights(y)
  w <- e$weight / sum(e$weight)
  for (t in seq_along(y)) near(f$smoothed[t, ], sapply(1:2, function(j) sum(w[paths[, t] == j])))
  counts <- matrix(0, 2, 2)
  for (i in 1:2) {
    for (j in 1:2) {
      for (t in 2:length(y)) {
        counts[i, j] <- counts[i, j] + sum(w[paths[, t - 1] == i & paths[, t] == j])
      }
    }
  }
  near(f$transition_counts, counts)
  near(rowSums(f$smoothed), rep(1, length(y)))
  near(f$smoothed[nrow(f$smoothed), ], tail(f$filtered, 1)[1, ])
  # Relabelling must leave likelihood unchanged and permute posterior columns.
  p <- 2:1
  g <- ms_filter(y, c0[p], phi[p], sigma[p], P[p, p], y0, prior[p])
  near(g$loglik, f$loglik)
  near(g$smoothed, f$smoothed[, p])
  # Forward probabilities at t cannot depend on observations after t.
  changed <- ms_filter(c(y[1:3], 10, -10, 8), c0, phi, sigma, P, y0, prior)
  near(changed$filtered[1:3, ], f$filtered[1:3, ])
  stopifnot(max(abs(changed$smoothed[1:3, ] - f$smoothed[1:3, ])) > .05)
  # An i.i.d. state sequence reduces to a static Gaussian mixture.
  Q <- matrix(rep(q, each = 2), 2)
  z <- ms_filter(y, c0, c(0, 0), sigma, Q, prior = q)
  near(z$loglik, sum(log(q[1] * dnorm(y, c0[1], sigma[1]) + q[2] * dnorm(y, c0[2], sigma[2]))))
  # Log-sum-exp must not underflow for unlikely but finite observations.
  rare <- ms_filter(c(50, -60, 45), c0, phi, sigma, P, y0, prior)
  stopifnot(is.finite(rare$loglik))
  near(rowSums(rare$filtered), rep(1, 3))
  # Exact forecasts checked against the explicit state-path mixture at each horizon.
  last_y <- .4
  a0 <- c(.25, .75)
  pred <- ms_forecast_moments(last_y, a0, c0, phi, sigma, P, horizon = 5)
  for (h in 1:5) {
    s <- as.matrix(expand.grid(rep(list(1:2), h + 1L)))
    weights <- means <- variances <- numeric(nrow(s))
    for (k in 1:nrow(s)) {
      path <- s[k, ]
      weights[k] <- a0[path[1]]
      m <- last_y
      v <- 0
      for (i in 2:length(path)) {
        weights[k] <- weights[k] * P[path[i - 1], path[i]]
        m <- c0[path[i]] + phi[path[i]] * m
        v <- sigma[path[i]]^2 + phi[path[i]]^2 * v
      }
      means[k] <- m
      variances[k] <- v
    }
    near(pred$mean[h], sum(weights * means))
    near(pred$variance[h], sum(weights * (variances + means^2)) - sum(weights * means)^2)
  }
  # Long-horizon forecasts must approach independently solved stationary moments.
  m <- solve(diag(2) - diag(phi) %*% t(P), q * c0)
  r <- solve(
    diag(2) - diag(phi^2) %*% t(P),
    q * (c0^2 + sigma^2) + 2 * c0 * phi * as.vector(crossprod(P, m))
  )
  long <- ms_forecast_moments(last_y, a0, c0, phi, sigma, P, horizon = 400)
  near(tail(long$mean, 1), sum(m))
  near(tail(long$variance, 1), sum(r) - sum(m)^2)

  # EM monotonicity and admissibility on a fixed, identifiable simulated sample.
  set.seed(8128)
  sim <- ms_simulate(350, c(-1, 1), c(.25, .25), c(.5, .5), P)
  fit <- fit_msar(sim$y, starts = 3, maxit = 600)
  stopifnot(
    fit$converged, !fit$boundary, min(diff(fit$history)) > -1e-6,
    max(abs(fit$phi)) < 1, all(fit$sigma > 0), fit$c[1] < fit$c[2]
  )
  near(rowSums(fit$P), c(1, 1))
})
