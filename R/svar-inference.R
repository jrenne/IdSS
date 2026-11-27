# Sampling uncertainty for recursive SVAR impulse responses.

#' @keywords internal
#' @noRd
VAR.regressors <- function(Y, p) {
  T <- dim(Y)[1]
  n <- dim(Y)[2]
  X <- array(1, c(T - p, 1))
  for (i in 1:p) {
    X <- cbind(X, Y[(1 + p - i):(T - i), ])
  }
  return(X = X)
}

#' @keywords internal
#' @noRd
VEC.coef <- function(Phi, cst, n, p) {
  phi <- cst
  for (i in 1:p) {
    phi <- cbind(phi, Phi[[i]])
  }
  phi <- t(phi)
  phi <- matrix(phi, ncol = 1)
  return(phi = phi)
}

#' @keywords internal
#' @noRd
distrib.phi <- function(Y, p, Phi, cst, Omega, n) {
  phi <- VEC.coef(Phi, cst, n, p)
  X <- VAR.regressors(Y, p)
  XX <- t(X) %*% X
  var.phi <- Omega %x% solve(XX)
  return(list(
    var.phi = var.phi,
    phi = phi
  ))
}

#' @keywords internal
#' @noRd
distrib.omega <- function(Omega, n, T) {
  omega <- NULL
  for (i in 1:n) {
    omega <- rbind(omega, matrix(Omega[i:n, i], ncol = 1))
  }
  var.omega <- array(NaN, c(dim(omega)[1], dim(omega)[1]))
  ind.r <- 1 # keep track of row index
  ind.c <- 1 # keep track of column index
  for (i in 1:n) {
    for (j in i:n) {
      for (k in 1:n) {
        for (l in k:n) {
          # correlation between sigma_ij and sigma_kl
          var.omega[ind.r, ind.c] <- (Omega[i, k] * Omega[j, l] + Omega[i, l] * Omega[j, k]) / T
          ind.c <- ind.c + 1
        }
      }
      ind.r <- ind.r + 1
      ind.c <- 1
    }
  }
  return(list(
    omega = omega,
    var.omega = var.omega
  ))
}

#' @keywords internal
#' @noRd
montecarlo <- function(Y, p,
                       nb.periods.IRF,
                       n,
                       nb.draws,
                       posit.of.shock = 0 # =0 to simulate for all shocks
) {
  # Estimate IRFs for all shocks or just one
  T <- dim(Y)[1]
  if (posit.of.shock != 0) {
    baseline.res <- svar.ordering.aux(
      Y, p,
      posit.of.shock,
      nb.periods.IRF
    )
  } else {
    baseline.res <- svar.ordering.all(
      Y, p,
      nb.periods.IRF,
      n
    )
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[, p * n + 1]
  resids <- residuals(est.VAR)
  Omega <- var(resids)

  simulated.IRFs <- replicate(nb.draws, IRFs, simplify = "array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify = "array")
  simulated.Phi <- matrix(rep(Phi, nb.draws), nrow = p, ncol = nb.draws)

  # Get the asymptotic distribution of estimated coefficients
  distrib.phi.res <- distrib.phi(Y, p, Phi, cst, Omega, n)
  phi <- distrib.phi.res$phi
  var.phi <- distrib.phi.res$var.phi
  Phi.simul <- list()
  distrib.omega.res <- distrib.omega(Omega, n, nrow(resids))
  omega <- distrib.omega.res$omega
  var.omega <- distrib.omega.res$var.omega
  Omega.simul <- array(0, c(n, n))

  i <- 0
  xx <- 0
  # Perform the simulations
  while (xx < nb.draws) {
    i <- i + 1
    if (i > max(1000L, 100L * nb.draws)) stop("Unable to obtain enough positive-definite covariance draws")
    # here draw B.hat
    omega.simul <- mvrnorm(1, omega, var.omega)
    j <- 1
    for (k in 1:n) {
      Omega.simul[k:n, k] <- omega.simul[j:(j + n - k)]
      Omega.simul[k, k:n] <- t(omega.simul[j:(j + n - k)])
      j <- j + n - k + 1
    }
    r <- eigen(Omega.simul) # use this if some Omega.simul are not positive definite
    eig <- r$values
    if (length(eig[eig > 0]) != n) {
      next
    } else {
      xx <- xx + 1
      # here draw Phi
      phi.simul <- mvrnorm(1, phi, var.phi)
      B.hat.simul <- t(chol(Omega.simul))
      phi.simul.temp <- t(matrix(phi.simul, ncol = n, nrow = 1 + n * p))
      for (j in 1:p) {
        Phi.simul[[j]] <- phi.simul.temp[, (2 + (j - 1) * n):(1 + j * n)]
      }
      # Simulate IRFs for all shocks or just one
      if (posit.of.shock != 0) {
        IRFs.simul <- simul.VAR(
          c = rep(0, n), Phi.simul,
          B.hat.simul,
          nb.periods.IRF,
          y0.star = rep(0, n * p),
          indic.IRF = 1,
          diag(n)[, posit.of.shock]
        )

        simulated.IRFs[, , xx] <- IRFs.simul
      } else {
        IRFs.simul <- simul.VAR.all(
          c = rep(0, n), n, Phi.simul,
          B.hat.simul,
          nb.periods.IRF,
          y0.star = rep(0, n * p),
          indic.IRF = 1,
          diag(n)
        )

        simulated.IRFs[, , , xx] <- IRFs.simul
      }
      simulated.B.hat[, , xx] <- B.hat.simul
      simulated.Phi[, xx] <- Phi.simul
    }
  }
  return(list(
    simulated.IRFs = simulated.IRFs,
    simulated.B.hat = simulated.B.hat,
    simulated.Phi = simulated.Phi
  ))
}

#' @keywords internal
#' @noRd
param.bootstrap <- function(Y, p,
                            nb.periods.IRF,
                            n,
                            nb.draws,
                            posit.of.shock = 0 # =0 to simulate for all shocks
) {
  # Estimate IRFs for all shocks or just one
  if (posit.of.shock != 0) {
    baseline.res <- svar.ordering.aux(
      Y, p,
      posit.of.shock,
      nb.periods.IRF
    )
  } else {
    baseline.res <- svar.ordering.all(
      Y, p,
      nb.periods.IRF,
      n
    )
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[, p * n + 1]

  simulated.IRFs <- replicate(nb.draws, IRFs, simplify = "array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify = "array")
  simulated.Phi <- matrix(rep(Phi, nb.draws), nrow = p, ncol = nb.draws)

  y0.star <- NULL
  for (k in p:1) {
    y0.star <- c(y0.star, Y[k, ])
  }

  # Perform the simulations
  for (i in 1:nb.draws) {
    # Simulate T observations
    simulated.Y <- simul.VAR(cst, Phi, B.hat, nb.sim = dim(Y)[1], y0.star)
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock != 0) {
      simulated.res <- svar.ordering.aux(
        simulated.Y, p,
        posit.of.shock,
        nb.periods.IRF
      )
      simulated.IRFs[, , i] <- simulated.res$IRFs
    } else {
      simulated.res <- svar.ordering.all(
        simulated.Y, p,
        nb.periods.IRF,
        n
      )
      simulated.IRFs[, , , i] <- simulated.res$IRFs
    }

    simulated.B.hat[, , i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi[, i] <- Acoef(simulated.est.VAR)
  }

  return(list(
    simulated.IRFs = simulated.IRFs,
    simulated.B.hat = simulated.B.hat,
    simulated.Phi = simulated.Phi
  ))
}

#' @keywords internal
#' @noRd
nonparam.bootstrap <- function(Y, p,
                               nb.periods.IRF,
                               n,
                               nb.draws,
                               posit.of.shock) {
  # Estimate IRFs for all shocks or just one
  if (posit.of.shock != 0) {
    baseline.res <- svar.ordering.aux(
      Y, p,
      posit.of.shock,
      nb.periods.IRF
    )
  } else {
    baseline.res <- svar.ordering.all(
      Y, p,
      nb.periods.IRF,
      n
    )
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[, p * n + 1]
  resids <- residuals(est.VAR)

  simulated.IRFs <- replicate(nb.draws, IRFs, simplify = "array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify = "array")
  simulated.Phi <- matrix(rep(Phi, nb.draws), nrow = p, ncol = nb.draws)

  y0.star <- NULL
  for (k in p:1) {
    y0.star <- c(y0.star, Y[k, ])
  }

  # Perform the simulations
  for (i in 1:nb.draws) {
    # draw set residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1], dim(Y)[1], replace = TRUE)
    simul.residuals <- resids[simul.index, ]
    # Simulate T observations
    simulated.Y <- simul.VAR(cst, Phi,
      B = diag(n), nb.sim = dim(Y)[1],
      y0.star, indic.IRF = 0, eta = simul.residuals
    )
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock != 0) {
      simulated.res <- svar.ordering.aux(
        simulated.Y, p,
        posit.of.shock,
        nb.periods.IRF
      )
      simulated.IRFs[, , i] <- simulated.res$IRFs
    } else {
      simulated.res <- svar.ordering.all(
        simulated.Y, p,
        nb.periods.IRF,
        n
      )
      simulated.IRFs[, , , i] <- simulated.res$IRFs
    }

    simulated.B.hat[, , i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi[, i] <- Acoef(simulated.est.VAR)
  }

  return(list(
    simulated.IRFs = simulated.IRFs,
    simulated.B.hat = simulated.B.hat,
    simulated.Phi = simulated.Phi
  ))
}

#' @keywords internal
#' @noRd
bootstrap.after.bootstrap <- function(Y, p,
                                      nb.periods.IRF,
                                      n,
                                      nb.draws,
                                      posit.of.shock) {
  # Estimate IRFs for all shocks or just one
  if (posit.of.shock != 0) {
    baseline.res <- svar.ordering.aux(
      Y, p,
      posit.of.shock,
      nb.periods.IRF
    )
  } else {
    baseline.res <- svar.ordering.all(
      Y, p,
      nb.periods.IRF,
      n
    )
  }
  # Store results
  IRFs <- baseline.res$IRFs
  est.VAR <- baseline.res$est.VAR
  Phi <- Acoef(est.VAR)
  B.hat <- baseline.res$B.hat
  cst <- Bcoef(est.VAR)[, p * n + 1]
  resids <- residuals(est.VAR)
  phi <- VEC.coef(Phi, cst, n, p)

  simulated.IRFs <- replicate(nb.draws, IRFs, simplify = "array")
  simulated.B.hat <- replicate(nb.draws, B.hat, simplify = "array")
  simulated.Phi <- matrix(rep(Phi, nb.draws), nrow = p, ncol = nb.draws)
  simulated.phi <- array(replicate(nb.draws, phi, simplify = "array"), c(dim(phi)[1], nb.draws))

  y0.star <- NULL
  for (k in p:1) {
    y0.star <- c(y0.star, Y[k, ])
  }

  # FIRST BOOTSTRAP
  # Perform the simulations
  for (i in 1:nb.draws) {
    # draw set of residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1], dim(Y)[1], replace = TRUE)
    simul.residuals <- resids[simul.index, ]
    # Simulate T observations
    simulated.Y <- simul.VAR(cst, Phi,
      B = diag(n), nb.sim = dim(Y)[1],
      y0.star, indic.IRF = 0, eta = simul.residuals
    )
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock != 0) {
      simulated.res <- svar.ordering.aux(
        simulated.Y, p,
        posit.of.shock,
        nb.periods.IRF
      )
      simulated.IRFs[, , i] <- simulated.res$IRFs
    } else {
      simulated.res <- svar.ordering.all(
        simulated.Y, p,
        nb.periods.IRF,
        n
      )
      simulated.IRFs[, , , i] <- simulated.res$IRFs
    }

    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi.temp <- Acoef(simulated.est.VAR)
    simulated.phi[, i] <- VEC.coef(simulated.Phi.temp, Bcoef(simulated.est.VAR)[, p * n + 1], n, p)
  }

  # Estimated mean bias; intercept entries are not part of the stability check.
  bias <- rowMeans(simulated.phi) - as.numeric(phi)
  phi.list <- function(v) {
    coefficients <- t(matrix(v, nrow = 1 + n * p, ncol = n))
    lapply(seq_len(p), function(j) coefficients[, (2 + (j - 1) * n):(1 + j * n), drop = FALSE])
  }
  stable <- function(v) max(Mod(eigen(make.PHI(phi.list(v)), only.values = TRUE)$values)) < 1
  shrink <- function(v, b) {
    # No stationarity-based correction is applied to an unstable OLS draw.
    if (!stable(v)) {
      return(list(coef = v, bias = rep(0, length(b))))
    }
    for (weight in seq(1, 0, by = -.01)) {
      candidate <- v - weight * b
      if (stable(candidate)) {
        return(list(coef = candidate, bias = weight * b))
      }
    }
    stop("Could not shrink the bias correction to a stable estimate")
  }
  if (!stable(phi)) stop("Bootstrap bias correction requires a stable fitted VAR")
  corrected <- shrink(as.numeric(phi), bias)
  bias <- corrected$bias
  unbiased.Phi <- phi.list(corrected$coef)
  fitted.mean <- solve(diag(n) - Reduce(`+`, Phi), cst)
  corrected.cst <- (diag(n) - Reduce(`+`, unbiased.Phi)) %*% fitted.mean
  # SECOND BOOTSTRAP
  # Perform the simulations
  for (i in 1:nb.draws) {
    # draw set of residuals from estimated ones
    simul.index <- sample(1:dim(resids)[1], dim(Y)[1], replace = TRUE)
    simul.residuals <- resids[simul.index, ]
    # Simulate T observations
    simulated.Y <- simul.VAR(as.numeric(corrected.cst), unbiased.Phi,
      B = diag(n), nb.sim = dim(Y)[1],
      y0.star, indic.IRF = 0, eta = simul.residuals
    )
    colnames(simulated.Y) <- colnames(Y)

    # Estimate IRFs on simulated data for all shocks or just one
    if (posit.of.shock != 0) {
      simulated.res <- svar.ordering.aux(
        simulated.Y, p,
        posit.of.shock,
        nb.periods.IRF
      )
    } else {
      simulated.res <- svar.ordering.all(
        simulated.Y, p,
        nb.periods.IRF,
        n
      )
    }

    simulated.B.hat[, , i] <- simulated.res$B.hat
    simulated.est.VAR <- simulated.res$est.VAR
    simulated.Phi.temp <- Acoef(simulated.est.VAR)
    simulated.phi <- VEC.coef(simulated.Phi.temp, Bcoef(simulated.est.VAR)[, p * n + 1], n, p)
    # Compute the unbiased Phi
    unbiased.simulated.phi <- shrink(as.numeric(simulated.phi), bias)$coef
    unbiased.simulated.Phi.temp <- t(matrix(unbiased.simulated.phi, ncol = n, nrow = 1 + n * p))
    unbiased.simulated.Phi <- list()
    for (j in 1:p) {
      unbiased.simulated.Phi[[j]] <- unbiased.simulated.Phi.temp[, (2 + (j - 1) * n):(1 + j * n)]
    }
    simulated.Phi[, i] <- unbiased.simulated.Phi

    if (posit.of.shock != 0) {
      IRFs.simul <- simul.VAR(
        c = rep(0, n), simulated.Phi[, i],
        simulated.B.hat[, , i],
        nb.periods.IRF,
        y0.star = rep(0, n * p),
        indic.IRF = 1,
        diag(n)[, posit.of.shock]
      )

      simulated.IRFs[, , i] <- IRFs.simul
    } else {
      IRFs.simul <- simul.VAR.all(
        c = rep(0, n), n, simulated.Phi[, i],
        simulated.B.hat[, , i],
        nb.periods.IRF,
        y0.star = rep(0, n * p),
        indic.IRF = 1,
        diag(n)
      )

      simulated.IRFs[, , , i] <- IRFs.simul
    }
  }

  return(list(
    simulated.IRFs = simulated.IRFs,
    simulated.B.hat = simulated.B.hat,
    simulated.Phi = simulated.Phi
  ))
}
