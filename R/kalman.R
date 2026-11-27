# Linear Gaussian state-space procedures; no additional dependency.
# Timing: alpha_t = c_t + F_t alpha_(t-1) + eta_t; y_t = d_t + Z_t alpha_t + eps_t.
# a0 and P0 describe alpha_0, BEFORE the first transition. Missing values are NA.

kf_sym <- function(x) (x + t(x)) / 2
kf_matrix_at <- function(x, t) if (length(dim(x)) == 3L) x[, , t, drop = FALSE][, , 1] else x
kf_offset <- function(x, n, k) {
  if (is.null(x)) {
    return(matrix(0, n, k))
  }
  if (.idss_real(x) && is.null(dim(x)) && length(x) == k && all(is.finite(x))) {
    return(matrix(rep(x, each = n), n, k))
  }
  if (is.matrix(x) && .idss_real(x) && identical(dim(x), c(n, k)) && all(is.finite(x))) {
    return(x)
  }
  stop("An offset must be a vector of the appropriate dimension or an n-by-k matrix")
}
kf_psd <- function(x) {
  all(is.finite(x)) && max(abs(x - t(x))) < 1e-9 * (1 + max(abs(x))) &&
    min(eigen(kf_sym(x), symmetric = TRUE, only.values = TRUE)$values) >= -1e-9 * (1 + max(abs(x)))
}

# Gaussian filtering, observed-component likelihood and RTS smoothing.
# Z, F, Q and R are constant matrices, or arrays with time as third dimension.
# The smoother requires positive-definite one-step state prediction covariances.
#' Kalman filtering, observed-data likelihood and RTS smoothing
#'
#' Filter and optionally smooth a linear Gaussian state-space model, allowing time-varying
#'   matrices and partially or wholly missing measurement vectors.
#'
#' @param Z Measurement matrix, with observation dimension by state dimension, or an array with
#'   these first two dimensions and time as the third.
#' @param F State transition matrix or time-indexed array. F at time t maps the state at t-1
#'   into the state at t.
#' @param Q Positive-semidefinite state disturbance covariance matrix or time-indexed array.
#' @param R Positive-semidefinite measurement disturbance covariance matrix or time-indexed
#'   array. State and measurement disturbances are independent.
#' @param d Optional measurement intercept vector, or matrix with time in rows. NULL means zero.
#' @param c Optional state-transition intercept vector, or matrix with time in rows. NULL means
#'   zero.
#' @param y Numeric vector for a scalar measurement, or numeric matrix with time in rows and
#'   observations in columns. NA entries are missing measurements.
#' @param a0 Finite numeric vector giving the mean of the state at time zero, before the first
#'   transition.
#' @param P0 Positive-semidefinite covariance matrix of the same initial state.
#' @param smooth Logical scalar. If TRUE, also run the Rauch-Tung-Striebel backward smoother.
#' @return A list with n-by-m predicted and filtered state means, m-by-m-by-n pred_cov and
#'   filt_cov arrays, n-by-p innovation matrix (NA for absent coordinates), innovation_cov (a
#'   list of observed-subspace covariance matrices; NULL at wholly missing dates), loglik,
#'   per-date increments, observed component counts, and the model inputs. With smooth=TRUE it
#'   also contains smoothed and smooth_cov, with the same dimensions as filtered and filt_cov.
#' @details The timing is alpha[t] = c[t] + F[t] alpha[t-1] + eta[t] and y[t] = d[t] + Z[t]
#'   alpha[t] + epsilon[t]. Disturbances are mutually independent, independent over time and
#'   independent of the Gaussian initial state. Matrices are deterministic or conditioned-on
#'   exogenous paths. A time-varying array must have exactly n time slices. The filter uses the
#'   observed rows of Z and the corresponding principal submatrix of R. A wholly absent
#'   observation skips updating and contributes zero to the observed-data likelihood. This
#'   conditions on an ignorable availability pattern.
#' The implementation uses Cholesky solves and Joseph-form covariance updates, without adding
#'   arbitrary jitter. Observed-subspace innovation covariances must be positive definite. RTS
#'   smoothing additionally requires positive-definite predicted state covariances; use a
#'   nonredundant state representation for degenerate cases. Covariances may otherwise be
#'   positive semidefinite, including exact measurements. All state covariances condition on
#'   fixed parameters and exclude parameter-estimation uncertainty. A finite large prior
#'   covariance is not exact diffuse initialization.
#' @examples
#' y <- c(1, NA, .5, -.2)
#' f <- kf_filter(y,
#'   Z = matrix(1), F = matrix(.9), Q = matrix(.2), R = matrix(.5),
#'   a0 = 0, P0 = matrix(1)
#' )
#' cbind(predicted = f$predicted[, 1], filtered = f$filtered[, 1], smoothed = f$smoothed[, 1])
#' stopifnot(f$increments[2] == 0)
#' @family Kalman state-space models
#' @export
kf_filter <- function(y, Z, F, Q, R, a0, P0, d = NULL, c = NULL, smooth = TRUE) {
  if (is.null(dim(y))) y <- matrix(y, ncol = 1L)
  if (!is.matrix(y) || !.idss_real(y) || !nrow(y) || !ncol(y) || any(is.infinite(y))) {
    stop("y must be a nonempty numeric matrix or vector, with NA for missing data")
  }
  if (!.idss_flag(smooth)) stop("smooth must be a logical scalar")
  n <- nrow(y)
  p <- ncol(y)
  m <- length(a0)
  if (!m || !.idss_real(a0) || !is.null(dim(a0)) || any(!is.finite(a0)) ||
    !is.matrix(P0) || !.idss_real(P0) ||
    !identical(dim(P0), c(m, m)) || !kf_psd(P0)) {
    stop("Invalid prior")
  }
  dims <- list(Z = c(p, m), F = c(m, m), Q = c(m, m), R = c(p, p))
  matrices <- list(Z = Z, F = F, Q = Q, R = R)
  for (name in names(matrices)) {
    x <- matrices[[name]]
    dx <- dim(x)
    if (!.idss_real(x) || !length(dx) || !length(dx) %in% c(2L, 3L) ||
      !identical(dx[1:2], dims[[name]]) ||
      (length(dx) == 3L && dx[3] != n) || any(!is.finite(x))) {
      stop("Invalid dimensions or values in ", name)
    }
    if (name %in% c("Q", "R")) {
      for (t in seq_len(if (length(dx) == 3L) n else 1L)) {
        if (!kf_psd(matrix(kf_matrix_at(x, t), dx[1], dx[2]))) stop(name, " must be positive semidefinite")
      }
    }
  }
  d <- kf_offset(d, n, p)
  c <- kf_offset(c, n, m)
  predicted <- filtered <- matrix(NA_real_, n, m)
  pred_cov <- filt_cov <- array(NA_real_, c(m, m, n))
  innovation <- matrix(NA_real_, n, p)
  innovation_cov <- vector("list", n)
  increments <- numeric(n)
  observed <- rowSums(!is.na(y))
  previous_a <- as.numeric(a0)
  previous_P <- P0
  for (t in seq_len(n)) {
    Ft <- matrix(kf_matrix_at(F, t), m, m)
    Qt <- matrix(kf_matrix_at(Q, t), m, m)
    Zt <- matrix(kf_matrix_at(Z, t), p, m)
    Rt <- matrix(kf_matrix_at(R, t), p, p)
    a <- as.vector(c[t, ] + Ft %*% previous_a)
    V <- kf_sym(Ft %*% previous_P %*% t(Ft) + Qt)
    predicted[t, ] <- a
    pred_cov[, , t] <- V
    keep <- which(!is.na(y[t, ]))
    if (length(keep)) {
      Zo <- Zt[keep, , drop = FALSE]
      Ro <- Rt[keep, keep, drop = FALSE]
      v <- y[t, keep] - d[t, keep] - as.vector(Zo %*% a)
      S <- kf_sym(Zo %*% V %*% t(Zo) + Ro)
      U <- tryCatch(chol(S), error = function(e) stop("Singular observation innovation covariance at t=", t))
      solve_S <- function(b) backsolve(U, forwardsolve(t(U), b))
      K <- t(solve_S(Zo %*% V))
      a <- as.vector(a + K %*% v)
      A <- diag(m) - K %*% Zo
      V <- kf_sym(A %*% V %*% t(A) + K %*% Ro %*% t(K)) # Joseph form
      innovation[t, keep] <- v
      innovation_cov[[t]] <- S
      increments[t] <- -.5 * (length(keep) * log(2 * pi) + 2 * sum(log(diag(U))) +
        sum(v * solve_S(v)))
    }
    filtered[t, ] <- a
    filt_cov[, , t] <- V
    previous_a <- a
    previous_P <- V
  }
  result <- list(
    predicted = predicted, filtered = filtered, pred_cov = pred_cov,
    filt_cov = filt_cov, innovation = innovation, innovation_cov = innovation_cov,
    loglik = sum(increments), increments = increments, observed = observed,
    y = y, Z = Z, F = F, Q = Q, R = R, a0 = a0, P0 = P0, d = d, c = c
  )
  if (smooth) {
    smoothed <- filtered
    smooth_cov <- filt_cov
    if (n > 1L) {
      for (t in (n - 1L):1L) {
        Ft <- matrix(kf_matrix_at(F, t + 1L), m, m)
        V <- matrix(pred_cov[, , t + 1L], m, m)
        U <- tryCatch(
          chol(V),
          error = function(e) stop("RTS requires nonsingular state prediction covariance")
        )
        B <- matrix(filt_cov[, , t], m, m) %*% t(Ft)
        J <- t(backsolve(U, forwardsolve(t(U), t(B))))
        smoothed[t, ] <- filtered[t, ] + as.vector(J %*% (smoothed[t + 1L, ] - predicted[t + 1L, ]))
        smooth_cov[, , t] <- kf_sym(matrix(filt_cov[, , t], m, m) +
          J %*% (matrix(smooth_cov[, , t + 1L], m, m) - V) %*% t(J))
      }
    }
    result$smoothed <- smoothed
    result$smooth_cov <- smooth_cov
  }
  result
}

# Stationary state covariance for spectral radius(F)<1, including singular Q.
#' Stationary covariance of a stable linear Gaussian state
#'
#' Solve the discrete Lyapunov equation for a stable state transition, including lower-rank
#'   state disturbance covariances.
#'
#' @param F Finite square transition matrix with spectral radius strictly below one.
#' @param Q Compatible positive-semidefinite state disturbance covariance matrix.
#' @return Symmetric state covariance matrix satisfying P = F P t(F) + Q.
#' @details The calculation solves the vectorized Lyapunov equation. It does not apply to random
#'   walks or local linear trends with unit roots. It supplies a covariance, not the stationary
#'   mean. For a transition intercept c, the stationary mean is solve(diag(nrow(F))-F,c).
#' @examples
#' F <- diag(c(.8, .6))
#' Q <- matrix(c(.2, .05, .05, .1), 2)
#' P <- kf_stationary_cov(F, Q)
#' stopifnot(max(abs(P - F %*% P %*% t(F) - Q)) < 1e-10)
#' @family Kalman state-space models
#' @export
kf_stationary_cov <- function(F, Q) {
  m <- nrow(F)
  if (!is.matrix(F) || !.idss_real(F) || !nrow(F) || !identical(dim(F), c(m, m)) ||
    !is.matrix(Q) || !.idss_real(Q) || !identical(dim(Q), c(m, m)) || !kf_psd(Q) ||
    any(!is.finite(F)) || max(Mod(eigen(F, only.values = TRUE)$values)) >= 1) {
    stop("A stable square transition and compatible positive semidefinite Q are required")
  }
  kf_sym(matrix(solve(diag(m * m) - kronecker(F, F), as.vector(Q)), m, m))
}

# Forecasts for constant matrices or explicitly supplied future arrays.
# These are observation forecasts, so their covariance includes R.
#' Forecast observations and states from a Kalman origin
#'
#' Propagate a filtered state distribution without future measurement updates, and compute
#'   Gaussian observation forecasts.
#'
#' @param Z Measurement matrix, with observation dimension by state dimension, or an array with
#'   these first two dimensions and time as the third.
#' @param F State transition matrix or time-indexed array. F at time t maps the state at t-1
#'   into the state at t.
#' @param Q Positive-semidefinite state disturbance covariance matrix or time-indexed array.
#' @param R Positive-semidefinite measurement disturbance covariance matrix or time-indexed
#'   array. State and measurement disturbances are independent.
#' @param d Optional measurement intercept vector, or matrix with time in rows. NULL means zero.
#' @param c Optional state-transition intercept vector, or matrix with time in rows. NULL means
#'   zero.
#' @param a Finite numeric filtered state mean vector at the forecast origin.
#' @param P Positive-semidefinite filtered state covariance at that origin.
#' @param horizon Positive integer forecast horizon. Time-indexed inputs must have this many
#'   future slices; offset matrices must have this many rows.
#' @return A list with horizon-by-p observation mean, p-by-p-by-horizon observation covariance,
#'   horizon-by-m state_mean, and m-by-m-by-horizon state_cov.
#' @details Future state disturbances and observation noise are included. In particular,
#'   observation covariance contains R; noise-free signal covariance omits R. Parameters and any
#'   future regressors/paths are fixed, so their uncertainty is not included. For a future
#'   array, its first slice is the first transition or measurement after the origin, not the
#'   last in-sample slice.
#' @examples
#' f <- kf_forecast(
#'   a = 0, P = matrix(.2), Z = matrix(1), F = matrix(.9),
#'   Q = matrix(.1), R = matrix(.5), horizon = 5
#' )
#' cbind(mean = f$mean[, 1], sd = sqrt(f$covariance[1, 1, ]))
#' stopifnot(all(f$covariance[1, 1, ] > .5))
#' @family Kalman state-space models
#' @export
kf_forecast <- function(a, P, Z, F, Q, R, horizon = 12L, d = NULL, c = NULL) {
  if (!.idss_integer(horizon)) {
    stop("horizon must be a positive integer")
  }
  if (!.idss_real(Z) || !length(dim(Z)) || !length(dim(Z)) %in% c(2L, 3L) ||
    dim(Z)[1] < 1L) {
    stop("Z must be a numeric matrix or time-indexed array")
  }
  p <- dim(Z)[1]
  m <- length(a)
  f <- kf_filter(matrix(NA_real_, horizon, p), Z, F, Q, R, a, P, d, c, smooth = FALSE)
  mean <- matrix(0, horizon, p)
  covariance <- array(0, c(p, p, horizon))
  for (t in seq_len(horizon)) {
    Zt <- matrix(kf_matrix_at(Z, t), p, m)
    Rt <- matrix(kf_matrix_at(R, t), p, p)
    mean[t, ] <- f$d[t, ] + as.vector(Zt %*% f$predicted[t, ])
    covariance[, , t] <- kf_sym(Zt %*% matrix(f$pred_cov[, , t], m, m) %*% t(Zt) + Rt)
  }
  list(mean = mean, covariance = covariance, state_mean = f$predicted, state_cov = f$pred_cov)
}

# Local-level ML conditional on fixed a0/P0. Positive variance bounds are explicit.
# Both parameters are variances; returns sd as well for convenient reporting.
#' Maximum-likelihood estimation of local-level variances
#'
#' Estimate state and measurement variances of a Gaussian local-level model using a fixed proper
#'   initial-state prior and three optimization starts.
#'
#' @param y Finite numeric observation vector with at least 20 entries; missing values are not
#'   supported by this estimator.
#' @param a0 Finite scalar initial-state mean, before the first transition.
#' @param P0 Finite non-negative scalar initial-state variance.
#' @param variance_floor Positive scalar lower bound on both variances.
#' @return A kf_filter result with smoothing, plus named variance and sd vectors (State and
#'   Measurement), converged, boundary, variance_floor and variance_upper. The best converged
#'   optimization is retained.
#' @details The model is alpha[t] = alpha[t-1] + eta[t], y[t] = alpha[t] + epsilon[t], with
#'   independent Gaussian disturbances. Log variances are optimized between the stated lower
#'   bound and max(100*var(diff(y)),100*variance_floor). The initial prior is held fixed in all
#'   evaluations. A boundary flag signals an active lower or upper constraint. No standard
#'   errors are reported. A local level has no stationary state distribution; a proper finite
#'   prior is distinct from an exact diffuse treatment.
#' @examples
#' set.seed(123)
#' n <- 100
#' state <- cumsum(c(rnorm(1), rnorm(n, sd = .3)))
#' y <- state[-1] + rnorm(n, sd = .6)
#' f <- kf_fit_local_level(y, a0 = 0, P0 = 1)
#' f$variance
#' f$boundary
#' stopifnot(f$converged)
#' @family Kalman state-space models
#' @export
#' @importFrom stats setNames
kf_fit_local_level <- function(y, a0 = 0, P0 = 1, variance_floor = 1e-6) {
  if (!.idss_real(y) || !is.null(dim(y))) stop("y must be a numeric vector")
  y <- as.numeric(y)
  if (length(y) < 20L || any(!is.finite(y)) || !.idss_scalar(a0) || !.idss_scalar(P0) || P0 < 0 ||
    !.idss_scalar(variance_floor) || variance_floor <= 0) {
    stop("Invalid local-level inputs")
  }
  scale <- var(diff(y))
  if (!is.finite(scale)) stop("The variance of first differences must be finite")
  upper <- max(100 * scale, 100 * variance_floor)
  objective <- function(par) {
    variance <- exp(par)
    -kf_filter(y, matrix(1), matrix(1), matrix(variance[1]), matrix(variance[2]),
      a0, matrix(P0),
      smooth = FALSE
    )$loglik
  }
  starts <- list(c(.1, .5), c(.5, .1), c(.25, .25))
  fits <- lapply(starts, function(x) {
    optim(log(pmax(variance_floor, scale * x)), objective,
      method = "L-BFGS-B", lower = rep(log(variance_floor), 2), upper = rep(log(upper), 2)
    )
  })
  valid <- which(vapply(fits, function(x) x$convergence == 0, logical(1)))
  if (!length(valid)) stop("Local-level optimization did not converge")
  best <- fits[[valid[which.min(vapply(fits[valid], function(x) x$value, numeric(1)))]]]
  v <- exp(best$par)
  f <- kf_filter(y, matrix(1), matrix(1), matrix(v[1]), matrix(v[2]), a0, matrix(P0))
  f$variance <- setNames(v, c("State", "Measurement"))
  f$sd <- sqrt(f$variance)
  f$converged <- TRUE
  f$boundary <- any(v <= variance_floor * (1 + 1e-4) | v >= upper * (1 - 1e-4))
  f$variance_floor <- variance_floor
  f$variance_upper <- upper
  f
}
