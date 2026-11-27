# Recursive (Cholesky) identification of structural shocks.

#' Identification of SVAR - Estimation of IRF
#'
#' This functions identifies a SVAR model and computes IRF. The estimation of structural shocks
#'   is done by Cholesky decomposition. Confidence intervals are obtained by boostrapping the
#'   estimated VAR model.
#' Optionally, the function performs inference on the IRFs using parametric bootstrap,
#'   non-parametric
#' bootstrap, Monte Carlo, or bootstrap-after-bootstrap procedures, and plots
#' confidence intervals if requested.
#' @param Y     Numeric matrix of size \code{T×n} containing the endogenous variables (one
#'   column per variable).
#' @param p Integer for the lag order.
#' @param posit.of.shock An integer giving the position (column index in \code{Y})
#' of the structural shock of interest. Default is 1.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param inference Integer. Method used for inference:
#'   \itemize{
#'     \item 0 = no inference (default)
#'     \item 1 = parametric bootstrap
#'     \item 2 = non-parametric bootstrap
#'     \item 3 = Monte Carlo
#'     \item 4 = bootstrap-after-bootstrap
#'   }
#' @param nb.draws Positive integer giving the number of Monte Carlo or
#'   bootstrap draws used when `inference` is non-zero.
#' @param nb.bootstrap.replications Deprecated compatibility alias for
#'   `nb.draws`. When supplied while `inference = 0`, parametric bootstrap
#'   inference is selected, matching the historical course interface.
#' @param confidence.interval A numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{IRFs}}{Matrix of estimated impulse response functions.}
#'   \item{\code{simulated.IRFs}}{Array of simulated IRFs across bootstrap or
#'     Monte Carlo replications (if \code{inference > 0}).}
#'   \item{\code{all.stdv.IRFs}}{Matrix of standard deviations of the simulated IRFs
#'     (if \code{inference > 0}).}
#'   \item{\code{all.CI.lower.bounds}}{Matrix of lower bounds of the confidence intervals.}
#'   \item{\code{all.CI.upper.bounds}}{Matrix of upper bounds of the confidence intervals.}
#'   \item{\code{all.CI.median}}{Matrix of medians of the simulated IRFs.}
#' }
#' @examples
#' library(vars)
#' library(IdSS)
#' data("USmonthlyExample")
#' considered.variables <- c("LIP", "UNEMP", "LCPI", "LPCOM", "FFR", "NBR", "TTR", "M1")
#' y <- as.matrix(USmonthlyExample[considered.variables])
#' ## Monte-Carlo
#' res.svar.ordering <-
#'   svar.ordering(y,
#'     p = 3, posit.of.shock = 5,
#'     nb.periods.IRF = 20,
#'     inference = 3,
#'     nb.draws = 200,
#'     confidence.interval = 0.90,
#'     indic.plot = 1
#'   )
#' @export
svar.ordering <- function(Y, p,
                          posit.of.shock = 1,
                          nb.periods.IRF = 20,
                          inference = 0,
                          nb.draws = 100, # This is used if inference >0
                          nb.bootstrap.replications = NULL,
                          confidence.interval = 0.90, # expressed in pp.
                          indic.plot = 1 # Plots are displayed if = 1.
) {
  if (!is.null(nb.bootstrap.replications)) {
    nb.draws <- nb.bootstrap.replications
    if (inference == 0) inference <- 1
  }

  names.of.variables <- colnames(Y)

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  baseline.res <- svar.ordering.aux(
    Y, p,
    posit.of.shock,
    nb.periods.IRF
  )
  IRFs <- baseline.res$IRFs

  if (inference > 0) {
    # All four samplers return the same simulated.IRFs array.
    sampler <- switch(as.character(inference),
      `1` = param.bootstrap,
      `2` = nonparam.bootstrap,
      `3` = montecarlo,
      `4` = bootstrap.after.bootstrap,
      stop("inference must be 0, 1, 2, 3 or 4")
    )
    simulated.res <- sampler(Y, p, nb.periods.IRF, n, nb.draws, posit.of.shock)
    all.simulated.IRFs.res <- simulated.res$simulated.IRFs

    # compute some key moments of the simulated IRFs
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    all.CI.median <- NULL
    for (i in 1:n) {
      all.IRFs.i <- matrix(all.simulated.IRFs.res[, i, ], ncol = nb.draws)

      stdv.IRF <- apply(all.IRFs.i, 1, sd)
      all.stdv.IRFs <- cbind(all.stdv.IRFs, stdv.IRF)

      all.CI.lower.bounds <- cbind(
        all.CI.lower.bounds,
        apply(
          all.IRFs.i, 1,
          function(x) {
            quantile(x, (1 - confidence.interval) / 2)
          }
        )
      )
      all.CI.upper.bounds <- cbind(
        all.CI.upper.bounds,
        apply(
          all.IRFs.i, 1,
          function(x) {
            quantile(x, 1 - (1 - confidence.interval) / 2)
          }
        )
      )
      all.CI.median <- cbind(
        all.CI.median,
        apply(
          all.IRFs.i, 1,
          function(x) {
            quantile(x, 0.5)
          }
        )
      )
    }
  } else {
    # if no inference, simply use the estimated IRFs

    all.simulated.IRFs.res <- NULL
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- IRFs
    all.CI.upper.bounds <- IRFs
    all.CI.median <- IRFs
  }

  # here we plot the IRFs (if asked)
  if (indic.plot == 1) {
    par(mfrow = c(2, ifelse(round(n / 2) == n / 2, n / 2, (n + 1) / 2)))
    for (i in 1:n) {
      plot(seq_len(nb.periods.IRF) - 1L, IRFs[, i],
        type = "l", lwd = 2, xlab = "", ylab = "",
        ylim = c(
          min(all.CI.lower.bounds[, i]),
          max(all.CI.upper.bounds[, i])
        ),
        main = paste0("(", letters[i], ") ", names.of.variables[i])
      )
      abline(h = 0, col = "grey")
      if (inference > 0) {
        lines(seq_len(nb.periods.IRF) - 1L, all.CI.lower.bounds[, i], col = "grey45", lty = 2, lwd = 2)
        lines(seq_len(nb.periods.IRF) - 1L, all.CI.upper.bounds[, i], col = "grey45", lty = 2, lwd = 2)
      }
    }
  }

  return(
    list(
      IRFs = IRFs,
      simulated.IRFs = all.simulated.IRFs.res,
      all.stdv.IRFs = all.stdv.IRFs,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds,
      all.CI.median = all.CI.median
    )
  )
}

#' @keywords internal
#' @noRd
svar.ordering.aux <- function(Y, p,
                              posit.of.shock,
                              nb.periods.IRF) {
  # This functions computes IRF (potentially with Confidence Intervals) using the
  # Christiano, Eichembaum and Evans methodology

  considered.variables <- colnames(Y)
  n <- dim(Y)[2]

  # Select number of lags in VAR models:

  est.VAR <- VAR(Y, p = p)

  # get estimated Phi matrices:
  Phi <- Acoef(est.VAR)

  # Compute the covariance matrix of VAR residuals:
  resids <- residuals(est.VAR)
  Omega <- var(resids)

  # Short-Run Restrictions ~ Cholesky:
  B.hat <- t(chol(Omega))

  y0.star <- rep(0, n * p)

  u.shock <- rep(0, n)
  u.shock[posit.of.shock] <- 1

  par(mfrow = c(2, 4))

  IRFs <- simul.VAR(
    c = rep(0, dim(Y)[2]),
    Phi,
    B.hat,
    nb.periods.IRF,
    y0.star,
    indic.IRF = 1,
    u.shock = u.shock
  )

  return(list(
    IRFs = IRFs,
    est.VAR = est.VAR,
    B.hat = B.hat
  ))
}

#' @keywords internal
#' @noRd
svar.ordering.all <- function(Y, p,
                              nb.periods.IRF,
                              n) {
  IRFs <- array(NaN, c(n, n, nb.periods.IRF))
  # Call svar.ordering.aux for all shocks
  for (j in 1:n) {
    baseline.res <- svar.ordering.aux(
      Y, p,
      j,
      nb.periods.IRF
    )
    IRFs.j <- baseline.res$IRFs
    dim(IRFs.j) <- c(nb.periods.IRF, n, 1)
    IRFs[, j, ] <- aperm(IRFs.j, c(2, 3, 1))
  }

  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[, p * n + 1]
  resids <- residuals(est.VAR)
  Omega <- var(resids)

  return(list(
    IRFs = IRFs,
    est.VAR = est.VAR,
    Phi = Phi, # IRFs.draw*Q
    B.hat = B.hat,
    cst = cst,
    resids = resids,
    Omega = Omega # Q
  ))
}

#' @keywords internal
#' @noRd
simul.VAR.all <- function(c = rep(0, n), n, Phi,
                          B.hat,
                          nb.periods.IRF,
                          y0.star,
                          indic.IRF = 1,
                          U.shocks = diag(n)) {
  IRFs <- array(NaN, c(n, n, nb.periods.IRF))
  # Call simul.VAR for all shocks
  for (j in 1:n) {
    IRFs.j <- simul.VAR(
      c = rep(0, n),
      Phi,
      B.hat,
      nb.periods.IRF,
      y0.star,
      indic.IRF = 1,
      U.shocks[, j]
    )
    IRFs.j <- t(IRFs.j)
    dim(IRFs.j) <- c(n, 1, nb.periods.IRF)
    IRFs[, j, ] <- IRFs.j
  }
  return(IRFs = IRFs)
}
