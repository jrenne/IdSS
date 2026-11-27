# Identification with external instruments.

#' Implementation of the IV-SVAR approach
#'
#' This function implements the IV-SVAR approach.
#'  We are interested in the IRF of Y following a structural shock (\eqn{\eta_1}),
#'  Z being an instrument for \eqn{\eta_1}
#'  Y is modelled as a VAR(p) model.
#'  The IRF are scaled in such a way that the contemporaneous impact of \eqn{\eta_1} on
#'   \eqn{Y_1} is 1.
#' @param Y A \code{T x n} numeric matrix
#'   of endogenous variables, with rows as time and columns as series.
#' @param Z A T × k (k ≥ 1) numeric matrix with external instrument(s) for the
#'   targeted structural shock (default: the first shock). If k > 1, columns are used
#'   jointly as instruments. Must be time-aligned with Y.
#' @param p Integer, the VAR lag order for \code{Y}.
#' @param names.of.variables Character vector of length \eqn{n} with the name of the endogenous
#'   variables.
#' @param nb.periods.IRF An integer specifying the number of periods for which
#'   impulse response functions are computed. Default is 20.
#' @param z.AR.order Integer \eqn{\ge 0}. AR order for \code{Z}. Used in the parametric
#'   bootstrap only.
#'  Default is \code{3}.
#' @param nb.bootstrap.replications An integer giving the number of bootstrap
#'   replications used to compute confidence intervals. Applies only when a
#'   parametric bootstrap is used. Default is 0 (no bootstrap).
#' @param confidence.interval  numeric value between 0 and 1 indicating the
#'   confidence level for bootstrap intervals (e.g. 0.90 for 90% intervals).
#'   Default is 0.90.
#' @param indic.plot Plots are displayed if = 1.
#'
#' @return
#' \item{IRFs}{A matrix of IRFs (horizons in rows).}
#' \item{all.simulated.IRFs.res}{Array  of
#'   bootstrap IRFs; \code{NULL} if no bootstrap.}
#' \item{all.stdv.IRFs}{Matrix of bootstrap standard
#'   deviations; \code{NULL} if no bootstrap.}
#' \item{Sigma.IV}{Empirical covariance matrix of the stacked innovation vector
#'   \eqn{(\hat{\varepsilon}_Y^\top,\; \hat{u}_Z)^\top} used in the bootstrap;
#'   \code{NULL} if no bootstrap.}
#' \item{cor.Sigma.IV}{Correlation matrix corresponding to \code{Sigma.IV};
#'   \code{NULL} if no bootstrap.}
#' \item{all.CI.lower.bounds, all.CI.upper.bounds}{
#'   Matrices of bootstrap percentile bands; \code{NULL} if no bootstrap.}
#' \item{all.B.tilde.1}{Matrix collecting the normalized impact vectors from
#'   each bootstrap draw (one column per replication); \code{NULL} if no
#'   bootstrap.}
#'
#'
#' @examples
#' library(vars)
#' data("USmonthly_1990_to_2012")
#' Z <- as.matrix(USmonthly_1990_to_2012[, "BRW", drop = FALSE])
#' considered.variables <- c("GS1", "LIP", "LCPI", "EBP")
#' Y <- as.matrix(USmonthly_1990_to_2012[, considered.variables])
#' keep <- complete.cases(Y, Z)
#' Y <- Y[keep, , drop = FALSE]
#' Z <- Z[keep, , drop = FALSE]
#' svar.iv(Y, Z,
#'   p = 4, names.of.variables = considered.variables,
#'   nb.periods.IRF = 20,
#'   z.AR.order = 1,
#'   nb.bootstrap.replications = 100,
#'   confidence.interval = 0.90,
#'   indic.plot = 1
#' )
#' @export
svar.iv <- function(Y, Z, p,
                    names.of.variables,
                    nb.periods.IRF = 20,
                    z.AR.order = 3, # This is used in the parametric bootstrap only
                    nb.bootstrap.replications = 0, # This is used in the parametric bootstrap only
                    confidence.interval = 0.90, # expressed in pp.
                    indic.plot = 1 # Plots are displayed if = 1.
) {
  # This function implements the IV-SVAR approach.
  # We are interested in the IRF of Y following a shock on a structural shock (eta_1), Z being an instrument for eta_1
  # Y is modelled as a VAR(p) model.
  # The IRF are scaled in such a way that the contemporaneous impact of eta_1 on Y_1 is 1.

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  baseline.res <- svar.iv.aux(
    Y, Z, p,
    names.of.variables,
    nb.periods.IRF
  )
  IRFs <- baseline.res$IRFs

  if (nb.bootstrap.replications > 0) {
    all.simulated.IRFs.res <- array(NaN, c(nb.periods.IRF, n, nb.bootstrap.replications))

    # Computation of std dev of IRF by parametric Gaussian bootstrap
    # (see Stock and Watson, appendix A.2)

    # step 0: fit an AR process for z:
    est.AR.z <- arima(c(Z), order = c(z.AR.order, 0, 0))

    # step 1: parameterize a VAR for [y',z]':
    est.VAR <- baseline.res$est.VAR
    Phi <- Acoef(est.VAR)
    resids <- residuals(est.VAR)

    c.IV <- Bcoef(est.VAR)[, p * n + 1]
    c.IV <- c(c.IV, 0) # add one term for the instrument
    Phi.IV <- list()
    for (k in 1:max(p, z.AR.order)) {
      aux <- matrix(0, n + 1, n + 1)
      if (k <= p) {
        aux[1:n, 1:n] <- Phi[[k]]
      }
      if (k <= z.AR.order) {
        aux[n + 1, n + 1] <- est.AR.z$coef[k]
      }
      Phi.IV[[k]] <- aux
    }
    RESIDS.eps.plus.z <- cbind(resids, est.AR.z$residuals[(1 + p):T])
    Sigma.IV <- var(RESIDS.eps.plus.z)
    cor.Sigma.IV <- cor(RESIDS.eps.plus.z)
    B.IV <- t(chol(Sigma.IV))
    y0.star <- NULL
    for (k in max(p, z.AR.order):1) {
      y0.star <- c(y0.star, Y[k, ], 0)
    }

    all.B.tilde.1 <- NULL
    for (i in 1:nb.bootstrap.replications) {
      simulated.y.plus.z <- simul.VAR(c.IV, Phi.IV, B.IV, nb.sim = dim(Y)[1], y0.star)
      simulated.Y <- simulated.y.plus.z[, 1:n]
      simulated.Z <- matrix(simulated.y.plus.z[, n + 1], ncol = 1)

      simulated.res <- svar.iv.aux(
        simulated.Y, simulated.Z, p,
        names.of.variables,
        nb.periods.IRF
      )

      all.simulated.IRFs.res[, , i] <- simulated.res$IRFs
      all.B.tilde.1 <- cbind(all.B.tilde.1, simulated.res$B.tilde.1)
    }

    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    for (i in 1:n) {
      all.IRFs.i <- matrix(all.simulated.IRFs.res[, i, ], ncol = nb.bootstrap.replications)

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
    }
  } else {
    all.simulated.IRFs.res <- NULL
    all.stdv.IRFs <- NULL
    all.CI.lower.bounds <- NULL
    all.CI.upper.bounds <- NULL
    all.B.tilde.1 <- NULL
    Sigma.IV <- NULL
    cor.Sigma.IV <- NULL
  }

  if (indic.plot == 1) {
    par(mfrow = c(
      2,
      ifelse(round(n / 2) == n / 2, n / 2, (n + 1) / 2)
    ))
    for (i in 1:n) {
      ylim <- if (is.null(all.CI.lower.bounds)) {
        range(IRFs[, i])
      } else {
        range(all.CI.lower.bounds[, i], all.CI.upper.bounds[, i])
      }
      plot(seq_len(nb.periods.IRF) - 1L, IRFs[, i],
        type = "l", lwd = 2, xlab = "", ylab = "", ylim = ylim,
        main = paste0("(", letters[i], ") ", names.of.variables[i])
      )
      abline(h = 0, col = "grey")
      if (!is.null(all.CI.lower.bounds)) {
        lines(seq_len(nb.periods.IRF) - 1L, all.CI.lower.bounds[, i], col = "grey45", lty = 2, lwd = 2)
        lines(seq_len(nb.periods.IRF) - 1L, all.CI.upper.bounds[, i], col = "grey45", lty = 2, lwd = 2)
      }
    }
  }

  return(
    list(
      IRFs = IRFs,
      all.simulated.IRFs.res = all.simulated.IRFs.res,
      all.stdv.IRFs = all.stdv.IRFs,
      Sigma.IV = Sigma.IV,
      cor.Sigma.IV = cor.Sigma.IV,
      all.CI.lower.bounds = all.CI.lower.bounds,
      all.CI.upper.bounds = all.CI.upper.bounds,
      all.B.tilde.1 = all.B.tilde.1
    )
  )
}

#' @keywords internal
#' @noRd
svar.iv.aux <- function(Y, Z, p,
                        names.of.variables,
                        nb.periods.IRF) {
  # This function implements the IV-SVAR approach.
  # We are interested in the IRF of Y following a shock on a structural shock (eta_1), Z being an instrument for eta_1
  # Y is modelled as a VAR(p) model.
  # The IRF are scaled in such a way that the contemporaneous impact of eta_1 on Y_1 is 1.

  if (!is.matrix(Y)) {
    print("Error: Y should be a matrix.")
    return(0)
  }
  if (!is.matrix(Z)) {
    print("Error: Z should be a matrix.")
    return(0)
  }

  T <- dim(Y)[1]
  n <- dim(Y)[2]

  # Estimate a VAR(p) model for Y:
  colnames(Y) <- names.of.variables
  est.VAR <- VAR(Y, p = p)

  # get estimated Phi matrices:
  Phi <- Acoef(est.VAR)

  # Compute the covariance matrix of VAR residuals:
  resids <- residuals(est.VAR)

  # Run TSLS estimations (eps[i] on eps[1], using Z as an instrument):
  B.tilde.1 <- NULL
  for (i in 1:n) {
    YY <- matrix(resids[, i], ncol = 1)
    XX <- matrix(resids[, 1], ncol = 1)
    ZZ <- matrix(Z[(1 + p):dim(Z)[1], ], ncol = dim(Z)[2])
    eq.iv <- tsls(YY, XX, ZZ, q = 3)
    B.tilde.1[i] <- eq.iv$b.iv
  }

  # Compute IRFs
  B.tilde <- matrix(0, n, n)
  B.tilde[, 1] <- B.tilde.1

  y0.star <- rep(0, dim(Y)[2] * p)

  nb.sim <- nb.periods.IRF
  IRFs <- simul.VAR(
    c = rep(0, dim(Y)[2]),
    Phi,
    B.tilde,
    nb.sim,
    y0.star,
    indic.IRF = 1,
    u.shock = c(1, rep(0, n - 1))
  )

  return(list(
    B.tilde.1 = B.tilde.1,
    IRFs = IRFs,
    est.VAR = est.VAR
  ))
}
