# Local projections and instrumental-variable local projections.

#' Local-Projection IRFs à la Jordà (2005)
#'
#' Computes impulse response functions (IRFs) using Jordà’s local-projection
#' method. For each horizon \eqn{h = 0, \dots, H}, the routine runs a separate
#' regression of \eqn{y_{i,t+h}} on the contemporaneous shock series (chosen by
#' \code{posit.of.shock}) and optional lags of the endogenous variables, then
#' stacks the coefficients across horizons. Pointwise standard errors are
#' computed via a Newey–West estimator. Responses are multiplied by the
#' standard deviation of the impact-date policy innovation, estimated once
#' using the full impact regression and held fixed across horizons.
#'
#' @param Y Numeric matrix of size \eqn{T \times n} containing the endogenous
#'   variables (one column per variable). Column names are used for plot titles.
#' @param posit.of.shock An integer giving the position (column index in \code{Y})
#' of the structural shock of interest. Default is 1.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param nb.lags.endog.var.4.control Integer; number of lags of all endogenous
#'   variables to include as controls in each local projection (default 0).
#' @param confidence.interval A numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#'
#' @return A list with:
#' \item{IRFs}{Numeric matrix \eqn{(H+1) \times n} of responses to a one-standard-deviation
#'   impact innovation. Row \code{h+1} is horizon \code{h}.}
#' \item{all.stdv.IRFs}{Numeric matrix \eqn{(H+1) \times n} of pointwise
#'   (Newey–West) standard errors corresponding to \code{IRFs}.}
#' @export
#'
#' @references
#' Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
#' Projections. \emph{American Economic Review}, 95(1), 161–182.
#' https://doi.org/10.1257/0002828053828474
#'
#'
#' @examples
#' data("USmonthlyExample")
#' considered.variables <- c("LIP", "UNEMP", "LCPI", "LPCOM", "FFR", "NBR", "TTR", "M1")
#' y <- as.matrix(USmonthlyExample[considered.variables])
#' res.jorda <- make.jorda.irf(y,
#'   posit.of.shock = 5,
#'   nb.periods.IRF = 12, nb.lags.endog.var.4.control = 3,
#'   indic.plot = 1, # Plots are displayed if = 1.
#'   confidence.interval = 0.90
#' )
#' @param hac.lag Optional non-negative integer for the Bartlett/Newey--West
#'   bandwidth. The default uses the larger of horizon plus one and
#'   floor(4*(N/100)^(2/9)), capped at N-1 usable observations.
#' @details The reported standard errors are pointwise normal approximations
#'   under valid identifying moments, adequate controls, moment/CLT conditions
#'   and a suitable HAC bandwidth. They are not weak-instrument robust or
#'   simultaneous path bands. For `make.jorda.irf`, the estimated impact-shock
#'   standard deviation is treated as fixed in the reported standard errors.
#' @export
make.jorda.irf <- function(Y,
                           posit.of.shock = 1,
                           nb.periods.IRF = 20,
                           nb.lags.endog.var.4.control = 0,
                           indic.plot = 1,
                           confidence.interval = 0.90,
                           hac.lag = NULL) {
  Y <- .lp_input(Y, nb.periods.IRF, nb.lags.endog.var.4.control)
  n <- ncol(Y)
  T <- nrow(Y)
  H <- nb.periods.IRF
  p <- nb.lags.endog.var.4.control
  k <- posit.of.shock
  if (length(k) != 1L || !is.finite(k) || k != as.integer(k) || k < 1 || k > n) {
    stop("posit.of.shock must select one column of Y")
  }
  .lp_hac_lag(hac.lag, 0, T - p)
  design <- function(h) {
    dates <- seq.int(p + 1L, T - h)
    X <- cbind(1, Y[dates, seq_len(k), drop = FALSE])
    if (p > 0) for (j in seq_len(p)) X <- cbind(X, Y[dates - j, , drop = FALSE])
    list(X = X, dates = dates)
  }
  first <- design(0)
  controls <- first$X[, -(k + 1L), drop = FALSE]
  innovation <- stats::lm.fit(controls, Y[first$dates, k])$residuals
  shock.sd <- stats::sd(innovation)
  if (!is.finite(shock.sd) || shock.sd <= 0) stop("The policy innovation has zero variance")
  IRFs <- SE <- matrix(NA_real_, H + 1L, n)
  for (h in 0:H) {
    message("Jorda's approach, Currently working on horizon h=", h)
    d <- design(h)
    X <- d$X
    N <- nrow(X)
    if (N <= ncol(X) || qr(X)$rank < ncol(X)) stop("Insufficient observations or collinear LP controls")
    inverse <- solve(crossprod(X))
    for (i in seq_len(n)) {
      outcome <- Y[d$dates + h, i]
      b <- inverse %*% crossprod(X, outcome)
      e <- as.numeric(outcome - X %*% b)
      scores <- X * e
      S <- NW.LongRunVariance(scores, .lp_hac_lag(hac.lag, h, N))
      V <- N * inverse %*% S %*% inverse
      IRFs[h + 1L, i] <- b[k + 1L] * shock.sd
      SE[h + 1L, i] <- sqrt(max(0, V[k + 1L, k + 1L])) * shock.sd
    }
  }
  if (indic.plot == 1) .lp_plot(IRFs, SE, colnames(Y), confidence.interval)
  list(IRFs = IRFs, all.stdv.IRFs = SE)
}

#' Local Projection IV impulse response functions (LP-IV)
#'
#' @description
#' Computes impulse response functions (IRFs) using the Local Projection
#' Instrumental Variables (LP-IV) approach. For each horizon \eqn{h}, the function
#' estimates a separate IV regression of future outcomes on a shock proxy,
#' instrumented by \code{Z}, optionally controlling for lags of the endogenous
#' variables and the instrument.
#'
#' By convention, the shock is normalized to have a unit effect on the first
#' column of \code{Y}.
#' @details
#' Let \eqn{Y_t} be an \eqn{n}-dimensional vector of endogenous variables and
#' \eqn{Z_t} an instrument. For each horizon \eqn{h = 0, \dots, H} and each variable
#' \eqn{i = 1, \dots, n}, the function estimates:
#' \deqn{
#' Y_{i,t+h} = \alpha_{i,h} + \beta_{i,h} Y_{1,t} + \Gamma_{i,h}' W_t + u_{i,t+h},
#' }
#' using two-stage least squares, where:
#' \itemize{
#'   \item \eqn{Y_{1,t}} is the shock variable (first column of \code{Y}),
#'   \item \eqn{Z_t} is the instrument for \eqn{Y_{1,t}},
#'   \item \eqn{W_t} includes optional lags of \code{Y} and \code{Z},
#'   \item a constant is always included.
#' }
#'
#' The coefficient \eqn{\beta_{i,h}} is interpreted as the impulse response of
#' variable \eqn{i} at horizon \eqn{h}. Standard errors are obtained from the
#' HAC-corrected covariance matrix returned by \code{\link{tsls}}, using a
#' Newey–West lag length of \code{h + 1}.
#'
#' @param Y Numeric matrix of dimension \eqn{T \times n} containing the endogenous
#' variables. The first column is interpreted as the shock variable whose impulse
#' responses are traced out.
#'
#' @param Z Numeric matrix of dimension \eqn{T \times m} containing the instrument(s)
#' for the shock variable.
#'
#' @param nb.periods.IRF Non-negative integer specifying the maximum horizon
#' \eqn{H} for which impulse responses are computed.
#'
#' @param nb.lags.Y.4.control Non-negative integer specifying the number of lags
#' of \code{Y} included as control variables in each local projection.
#'
#' @param nb.lags.Z.4.control Non-negative integer specifying the number of lags
#' of \code{Z} included as control variables in each local projection.
#'
#' @param indic.plot Logical or integer (1/0). If equal to 1, impulse responses
#' with confidence intervals are plotted.
#'
#' @param confidence.interval Numeric scalar in \eqn{(0,1)} specifying the
#' confidence level used to construct pointwise confidence intervals.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{IRFs}: A matrix of dimension \eqn{(H+1) \times n} containing the
#'   estimated impulse responses for each variable and horizon.
#'   \item \code{all.stdv.IRFs}: A matrix of dimension \eqn{(H+1) \times n}
#'   containing the corresponding standard errors.
#' }
#'
#' @examples
#' data("USmonthly_1990_to_2012")
#' Z <- as.matrix(USmonthly_1990_to_2012[, "BRW", drop = FALSE])
#' considered.variables <- c("GS1", "LIP", "LCPI", "EBP")
#' Y <- as.matrix(USmonthly_1990_to_2012[, considered.variables])
#' keep <- complete.cases(Y, Z)
#' Y <- Y[keep, , drop = FALSE]
#' Z <- Z[keep, , drop = FALSE]
#' res.LP.IV <- make.LPIV.irf(Y, Z,
#'   nb.periods.IRF = 20,
#'   nb.lags.Y.4.control = 4,
#'   nb.lags.Z.4.control = 4,
#'   indic.plot = 1, # Plots are displayed if = 1.
#'   confidence.interval = 0.90
#' )
#'
#' @param hac.lag Optional non-negative integer for the Bartlett/Newey--West
#'   bandwidth. The default uses the larger of horizon plus one and
#'   floor(4*(N/100)^(2/9)), capped at N-1 usable observations.
#' @details The reported standard errors are pointwise normal approximations
#'   under valid identifying moments, adequate controls, moment/CLT conditions
#'   and a suitable HAC bandwidth. They are not weak-instrument robust or
#'   simultaneous path bands. For `make.jorda.irf`, the estimated impact-shock
#'   standard deviation is treated as fixed in the reported standard errors.
#' @export
make.LPIV.irf <- function(Y, Z,
                          # posit.of.shock = 1,
                          nb.periods.IRF = 20,
                          nb.lags.Y.4.control = 0,
                          nb.lags.Z.4.control = 0,
                          indic.plot = 1, # Plots are displayed if = 1.
                          confidence.interval = 0.90,
                          hac.lag = NULL) {
  # This function implements the LP-IV approach.
  # By convention, the shock has a unit effect on the first component of Y
  # Z is the instrumental variable.

  Y <- .lp_input(Y, nb.periods.IRF, max(nb.lags.Y.4.control, nb.lags.Z.4.control))
  if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != nrow(Y) || any(!is.finite(Z))) {
    stop("Z must be a finite numeric matrix with the same rows as Y")
  }
  .lp_hac_lag(hac.lag, 0, nrow(Y))
  considered.variables <- colnames(Y)
  T <- dim(Y)[1]
  n <- dim(Y)[2]

  max.nb.lags.control <- max(nb.lags.Y.4.control, nb.lags.Z.4.control)

  IRFs <- matrix(NaN, nb.periods.IRF + 1, n)
  all.stdv.IRFs <- matrix(NaN, nb.periods.IRF + 1, n)

  for (h in 0:nb.periods.IRF) {
    message(paste("LP-IV approach, Currently working on horizon h=", toString(h), sep = ""))

    for (i in 1:n) {
      ZZ <- matrix(Z[(1 + max.nb.lags.control):(T - h), ], ncol = dim(Z)[2])
      X <- matrix(Y[(1 + max.nb.lags.control):(T - h), 1], ncol = 1)

      y.i.t.p.h <- matrix(Y[(h + 1 + max.nb.lags.control):T, i], ncol = 1)
      if (nb.lags.Y.4.control > 0) {
        for (k in 1:nb.lags.Y.4.control) {
          X <- cbind(X, Y[(1 + max.nb.lags.control - k):(T - h - k), ])
          ZZ <- cbind(ZZ, Y[(1 + max.nb.lags.control - k):(T - h - k), ])
        }
      }
      if (nb.lags.Z.4.control > 0) {
        for (k in 1:nb.lags.Z.4.control) {
          X <- cbind(X, Z[(1 + max.nb.lags.control - k):(T - h - k), ])
          ZZ <- cbind(ZZ, Z[(1 + max.nb.lags.control - k):(T - h - k), ])
        }
      }
      X <- cbind(1, X) # add a constant in the regression
      ZZ <- cbind(ZZ, 1) # add a constant in the regression
      res.iv.est <- tsls(y.i.t.p.h, X, ZZ, .lp_hac_lag(hac.lag, h, nrow(X)))
      IRFs[h + 1, i] <- res.iv.est$b.iv[2] # Note: res.iv.est$b.iv[1] corresponds to the constant
      all.stdv.IRFs[h + 1, i] <- sqrt(max(0, res.iv.est$covmat.b.iv[2, 2]))
    }
  }

  if (indic.plot == 1) .lp_plot(IRFs, all.stdv.IRFs, considered.variables, confidence.interval)

  return(
    list(
      IRFs = IRFs,
      all.stdv.IRFs = all.stdv.IRFs
    )
  )
}

#' @keywords internal
#' @noRd
.lp_input <- function(Y, H, p) {
  if (!is.matrix(Y) || !is.numeric(Y) || any(!is.finite(Y))) stop("Y must be a finite numeric matrix")
  if (length(H) != 1L || !is.finite(H) || H < 0 || H != as.integer(H) ||
    length(p) != 1L || !is.finite(p) || p < 0 || p != as.integer(p)) {
    stop("Horizons and lag counts must be non-negative integers")
  }
  if (nrow(Y) - p - H < 3L || ncol(Y) < 1L) {
    stop("Insufficient observations for the requested horizons and lags")
  }
  Y
}

#' @keywords internal
#' @noRd
.lp_hac_lag <- function(q, h, N) {
  if (is.null(q)) {
    return(min(N - 1L, max(h + 1L, floor(4 * (N / 100)^(2 / 9)))))
  }
  if (length(q) != 1L || !is.finite(q) || q < 0 || q != as.integer(q)) {
    stop("hac.lag must be NULL or a non-negative integer")
  }
  min(q, N - 1L)
}

#' @keywords internal
#' @noRd
.lp_plot <- function(IRFs, SE, labels, confidence.interval) {
  n <- ncol(IRFs)
  H <- nrow(IRFs) - 1L
  if (is.null(labels)) labels <- paste("Variable", seq_len(n))
  par(mfrow = c(2, ceiling(n / 2)))
  critical <- stats::qnorm((1 + confidence.interval) / 2)
  for (i in seq_len(n)) {
    lo <- IRFs[, i] - critical * SE[, i]
    hi <- IRFs[, i] + critical * SE[, i]
    plot(0:H, IRFs[, i],
      type = "l", lwd = 2, xlab = "Horizon", ylab = "",
      ylim = range(c(lo, hi)), main = paste0("(", letters[i], ") ", labels[i])
    )
    abline(h = 0, col = "grey")
    lines(0:H, lo, col = "grey45", lty = 2, lwd = 2)
    lines(0:H, hi, col = "grey45", lty = 2, lwd = 2)
  }
}
