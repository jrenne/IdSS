test_that("Independent kalman identities hold", {
  # Independent Gaussian conditioning checks; no additional package.
  near <- function(x, y, tol = 1e-8) expect_equal(as.numeric(x), as.numeric(y), tolerance = tol)
  block_diag <- function(blocks) {
    sizes <- vapply(blocks, nrow, integer(1))
    ends <- cumsum(sizes)
    out <- matrix(0, sum(sizes), sum(sizes))
    for (i in seq_along(blocks)) {
      k <- (ends[i] - sizes[i] + 1L):ends[i]
      out[k, k] <- blocks[[i]]
    }
    out
  }
  n <- 7L
  m <- p <- 2L
  F <- Z <- Q <- R <- array(0, c(2, 2, n))
  for (t in 1:n) {
    F[, , t] <- matrix(c(.7, .08, -.05, .6), 2) + diag(c(.01, -.01) * t)
    Z[, , t] <- matrix(c(1, .4, -.2, 1), 2) + diag(c(.01, .02) * t)
    Q[, , t] <- matrix(c(.3, .06, .06, .2), 2)
    R[, , t] <- matrix(c(.5, .12, .12, .4), 2)
  }
  a0 <- c(.2, -.1)
  P0 <- matrix(c(.8, .1, .1, .6), 2)
  c0 <- cbind(seq(.01, .07, length.out = n), rep(-.03, n))
  d0 <- cbind(rep(.4, n), seq(0, .12, length.out = n))
  y <- cbind(c(.1, NA, NA, .3, -.4, .8, .2), c(.7, .2, NA, NA, .1, .6, -.3))
  f <- kf_filter(y, Z, F, Q, R, a0, P0, d0, c0)
  # Construct the complete joint Gaussian law using explicit shock-loading matrices.
  # This calculation does not use innovation or smoothing recursions.
  W <- block_diag(c(list(P0), lapply(1:n, function(t) Q[, , t])))
  L <- matrix(0, m, m * (n + 1L))
  L[, 1:m] <- diag(m)
  A <- matrix(0, n * m, ncol(L))
  mean_state <- numeric(n * m)
  previous <- a0
  for (t in 1:n) {
    L <- F[, , t] %*% L
    L[, (m * t + 1L):(m * (t + 1L))] <- diag(m)
    k <- ((t - 1L) * m + 1L):(t * m)
    A[k, ] <- L
    previous <- as.vector(c0[t, ] + F[, , t] %*% previous)
    mean_state[k] <- previous
  }
  Vstate <- A %*% W %*% t(A)
  D <- block_diag(lapply(1:n, function(t) Z[, , t]))
  mean_y <- as.vector(t(d0)) + as.vector(D %*% mean_state)
  Vy <- D %*% Vstate %*% t(D) + block_diag(lapply(1:n, function(t) R[, , t]))
  C <- Vstate %*% t(D)
  flat_y <- as.vector(t(y))
  all_obs <- which(!is.na(flat_y))
  conditional <- function(k, obs) {
    B <- C[k, obs, drop = FALSE] %*% solve(Vy[obs, obs, drop = FALSE])
    list(
      mean = mean_state[k] + as.vector(B %*% (flat_y[obs] - mean_y[obs])),
      cov = Vstate[k, k, drop = FALSE] - B %*% t(C[k, obs, drop = FALSE])
    )
  }
  for (t in 1:n) {
    k <- ((t - 1L) * m + 1L):(t * m)
    obs <- all_obs[all_obs <= t * p]
    dense <- conditional(k, obs)
    near(f$filtered[t, ], dense$mean)
    near(f$filt_cov[, , t], dense$cov)
    dense <- conditional(k, all_obs)
    near(f$smoothed[t, ], dense$mean)
    near(f$smooth_cov[, , t], dense$cov)
    stopifnot(min(eigen(f$filt_cov[, , t] - f$smooth_cov[, , t], symmetric = TRUE)$values) > -1e-9)
  }
  U <- chol(Vy[all_obs, all_obs])
  v <- flat_y[all_obs] - mean_y[all_obs]
  ll <- -.5 * (length(v) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(forwardsolve(t(U), v)^2))
  near(ll, f$loglik)
  near(f$filtered[3, ], f$predicted[3, ])
  near(f$filt_cov[, , 3], f$pred_cov[, , 3])
  near(f$increments[3], 0)
  near(tail(f$smoothed, 1), tail(f$filtered, 1))
  near(f$smooth_cov[, , n], f$filt_cov[, , n])
  # Future observations revise smoothers, never earlier filtered values.
  y2 <- y
  y2[6:7, ] <- 10
  changed <- kf_filter(y2, Z, F, Q, R, a0, P0, d0, c0)
  near(f$filtered[1:5, ], changed$filtered[1:5, ])
  stopifnot(max(abs(f$smoothed[1:5, ] - changed$smoothed[1:5, ])) > .1)
  # Scalar update and one-step smoothing worked example.
  z <- kf_filter(c(2, 0), matrix(1), matrix(1), matrix(1), matrix(1), 0, matrix(1))
  near(z$filtered[, 1], c(4 / 3, .5))
  near(z$filt_cov[1, 1, ], c(2 / 3, 5 / 8))
  near(z$smoothed[1, 1], 1)
  near(z$smooth_cov[1, 1, 1], .5)
  # ARMA(1,1) singular state noise, exact measurement and stationary initialization.
  set.seed(9731)
  phi <- .6
  theta <- .4
  sd <- .8
  mu <- 1.5
  x <- as.numeric(arima.sim(list(ar = phi, ma = theta), n = 100, sd = sd)) + mu
  Fa <- matrix(c(phi, 0, theta, 0), 2)
  Qa <- sd^2 * matrix(1, 2, 2)
  Pa <- kf_stationary_cov(Fa, Qa)
  near(Pa, Fa %*% Pa %*% t(Fa) + Qa)
  a <- kf_filter(x - mu, matrix(c(1, 0), 1), Fa, Qa, matrix(0), c(0, 0), Pa, smooth = FALSE)
  a_ref <- arima(x, order = c(1, 0, 1), fixed = c(phi, theta, mu), transform.pars = FALSE, method = "ML", SSinit = "Rossignol2011")
  # stats::arima profiles out the innovation variance: compare at its fitted scale.
  Qref <- a_ref$sigma2 * matrix(1, 2, 2)
  a_scaled <- kf_filter(x - mu, matrix(c(1, 0), 1), Fa, Qref, matrix(0), c(0, 0),
    kf_stationary_cov(Fa, Qref),
    smooth = FALSE
  )
  near(a_scaled$loglik, as.numeric(logLik(a_ref)), 1e-7)
  # Observation forecast variance must include both propagated state risk and R.
  F0 <- matrix(c(.8, .1, 0, .7), 2)
  Q0 <- diag(c(.2, .1))
  Z0 <- matrix(c(1, .5), 1)
  R0 <- matrix(.4)
  a_last <- c(.2, -.3)
  P_last <- diag(c(.6, .2))
  pred <- kf_forecast(a_last, P_last, Z0, F0, Q0, R0, horizon = 5)
  Fpower <- diag(2)
  Vsum <- matrix(0, 2, 2)
  for (h in 1:5) {
    Fpower <- F0 %*% Fpower
    G <- diag(2)
    Vsum <- matrix(0, 2, 2)
    for (j in 0:(h - 1L)) {
      Vsum <- Vsum + G %*% Q0 %*% t(G)
      G <- F0 %*% G
    }
    near(pred$mean[h, 1], as.vector(Z0 %*% Fpower %*% a_last))
    near(pred$covariance[1, 1, h], as.vector(Z0 %*% (Fpower %*% P_last %*% t(Fpower) + Vsum) %*% t(Z0) + R0))
  }
  # A noiseless observation is allowed, but redundant singular observations must fail.
  exact <- kf_filter(c(.1, .4), matrix(1), matrix(.8), matrix(.2), matrix(0), 0, matrix(1), smooth = FALSE)
  near(exact$filtered[, 1], c(.1, .4))
  near(exact$filt_cov, array(0, c(1, 1, 2)))
  bad <- try(kf_filter(matrix(1, 2, 2), matrix(c(1, 1), 2), matrix(1), matrix(1), matrix(0, 2, 2), 0, matrix(1)), silent = TRUE)
  stopifnot(inherits(bad, "try-error"))
})
