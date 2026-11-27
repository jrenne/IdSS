# Forecast-error variance decomposition and plots.

#' Variance Decomposition from IRFs
#'
#' This functions computes the variance decomposition associated to a set
#' of IRFs organised as variables x shocks x horizons x IRF draw.
#'
#' @details
#' We must have that the number of variables is equal to the number of shocks. Hence,
#' IRFs has 4 dimensions. If there is only one draw (e.g. estimated model)
#' , then IRFs may have only 3 dimensions. For a single draw, the returned
#' arrays omit the singleton draw dimension, even when the input is 4-dimensional.
#' @param IRFs An array of IRFs with dimensions variables x shocks x horizons x IRF draws.
#'
#' @returns A list with components:
#' \item{variance}{Array where variance(i,j,h,k,l) is the contribution
#' of shock l to the forecast error covariance between variable i and j
#' for IRF draw k at horizon h (if i=j then it's the contribution to the variance of i).}
#' \item{Variance}{Array where Variance(i,j,h,k) is the forecast error covariance
#' between variable i and j for IRF draw k at horizon h (if i=j then it's the variance of i).}
#' \item{vardecomp}{Array where vardecomp(i,j,h,k,l) = variance(i,j,h,k,l)/Variance(i,j,h,k).}
#' @examples
#' # A two-variable model with two independent structural shocks.
#' B <- matrix(c(1, 0.3, 0.2, 1), nrow = 2)
#' IRFs <- array(0, dim = c(2, 2, 6))
#' for (h in 1:6) IRFs[, , h] <- diag(c(0.6, 0.4)^(h - 1)) %*% B
#' decomposition <- variance.decomp(IRFs)
#' # Shares of variable 1's forecast-error variance at the final horizon.
#' shares <- decomposition$vardecomp[1, 1, 6, ]
#' shares
#' sum(shares) # the two shock contributions sum to one
#' @export
variance.decomp <- function(IRFs) {
  # This functions computes the variance decomposition associated to a set
  # of IRFs organised as variables x shocks x horizons x IRF draws (we must
  # have that the number of variables is equal to the number of shocks)
  # Hence, IRFs has 4 dimensions. If there is only one draw (e.g. estimated model),
  #   then IRFs has only 3 dimensions.

  n <- dim(IRFs)[1]
  H <- dim(IRFs)[3]

  if (length(dim(IRFs)) < 4) {
    # In that case, add the fourth dimension (remains 1)
    N <- 1
    IRFs <- array(IRFs, c(n, n, H, N))
  }

  N <- dim(IRFs)[4]

  Variance <- array(0, c(n, n, H, N))
  variance <- array(0, c(n, n, H, N, n))
  Select <- array(0, c(n, n, n))
  for (k in 1:N) {
    Variance[, , 1, k] <- IRFs[, , 1, k] %*% t(IRFs[, , 1, k])
    for (i in 1:n) {
      Select[i, i, i] <- 1
      A <- Select[, , i]
      variance[, , 1, k, i] <- IRFs[, , 1, k] %*% A %*% t(IRFs[, , 1, k])
    }
    for (h in 2:H) {
      Variance[, , h, k] <- IRFs[, , h, k] %*% t(IRFs[, , h, k]) + Variance[, , h - 1, k]
      for (i in 1:n) {
        A <- Select[, , i]
        variance[, , h, k, i] <- IRFs[, , h, k] %*% A %*% t(IRFs[, , h, k]) + variance[, , h - 1, k, i]
      }
    }
  }
  Variance_rep <- replicate(n, Variance, simplify = "array")
  vardecomp <- variance / Variance_rep

  if (N == 1) {
    variance <- array(variance, c(n, n, H, n))
    Variance <- array(Variance, c(n, n, H))
    vardecomp <- array(vardecomp, c(n, n, H, n))
  }

  return(
    list(
      variance = variance,
      # variance(i,j,h,k,l): contribution of shock l to the forecast error
      # covariance between variable i and j for IRF draw k at horizon h
      # (if i=j then it's the contribution to the variance of i)
      Variance = Variance,
      # Variance(i,j,h,k): forecast error covariance between variable i
      # and j for IRF draw k at horizon h (if i=j then it's the
      # variance of i)
      vardecomp = vardecomp
      # vardecomp(i,j,h,k,l) = variance(i,j,h,k,l)/Variance(i,j,h,k)
    )
  )
}

#' Variance Decomposition for SVAR Models
#'
#' This function produces a plot showing the variance decomposition
#' associated with a VAR model defined by autoregressive matrices Phi
#' and an impact matrix B.
#'
#' @details
#' The function first constructs impulse responses by simulating the VAR
#' (one unit shock at a time). Then, it computes the variance decomposition.
#' Each panel corresponds to one endogenous variable; within a panel, the stacked
#' grayscale bands sum to 1 at each horizon, showing the fraction of variance
#' attributable to each shock.
#'
#' @param Phi Array of autoregressive coefficient matrices.
#' @param B Numeric \eqn{n \times n} matrix; contemporaneous impact matrix
#'   mapping structural shocks \eqn{u_t}to the reduced-form residuals.
#' @param maxHorizon Integer; maximum horizon for the variance decomposition.
#' @param mfrow Optional length-2 integer vector passed to \code{par(mfrow = ...)}
#' to arrange multiple panels. Defaults internally to \code{c(1, n)} when not specified (i.e.,
#'   when \code{mfrow = NA}).
#' @param names.var Optional character vector of length \eqn{n} with labels for variables.
#'   Defaults to \code{"Variable 1"}, \code{"Variable 2"}, … .
#' @param names.shock Optional character vector of length \eqn{n} with labels for shocks.
#'   Defaults to \code{"Shock 1"}, \code{"Shock 2"}, … .
#'
#' @returns the function output is a set of plots.
#' @examples
#' Phi <- array(NaN, c(2, 2, 2)) # (2,2,2) for (n,n,p)
#' Phi[, , 1] <- matrix(c(.6, 0, .2, .5), 2, 2)
#' Phi[, , 2] <- matrix(c(-.1, .2, .1, .3), 2, 2)
#' B <- matrix(c(.5, -1, 1.5, .8), 2, 2)
#' make.variance.decompo(Phi, B, maxHorizon = 20)
#' @export
make.variance.decompo <- function(Phi, B, maxHorizon,
                                  mfrow = NA,
                                  names.var = NA,
                                  names.shock = NA) {
  # This function produces a plot showing the variance decomposition
  # associated with a VAR model defined by autoregressive matrices Phi
  # and an impact matrix B.

  n <- dim(B)[1]
  if (is.na(names.var[1])) names.var <- paste("Variable", 1:n)
  if (is.na(names.shock[1])) names.shock <- paste("Shock", 1:n)

  IRFs <- array(NA, c(n, n, maxHorizon))
  for (i in 1:n) {
    u.shock <- rep(0, n)
    u.shock[i] <- 1
    Y <- simul.VAR(
      c = matrix(0, n, 1), Phi = Phi, B = B,
      nb.sim = maxHorizon, indic.IRF = 1, u.shock = u.shock
    )
    IRFs[, i, ] <- t(Y)
  }
  res <- variance.decomp(IRFs)$vardecomp

  # Handle layout flexibly
  if (is.na(mfrow[1])) {
    mfrow <- c(1, n)
  }
  par(mfrow = mfrow)

  # Adjust margins dynamically depending on layout
  nrows <- mfrow[1]
  ncols <- mfrow[2]

  par(plt = c(.15, .95, .2, .8))

  for (variable in 1:n) {
    res_variable <- res[variable, variable, , ]
    res_variable <- apply(res_variable, 1, cumsum)
    res_variable <- rbind(0, res_variable)

    plot(1:maxHorizon, rep(0, maxHorizon),
      ylim = c(0, 1), col = "white", las = 1,
      xlab = "", ylab = "Variance share",
      main = paste0("(", letters[variable], ") ", names.var[variable])
    )

    fills <- grDevices::gray.colors(n, start = 0.25, end = 0.85)
    for (k in 1:n) {
      polygon(c(1:maxHorizon, maxHorizon:1),
        c(res_variable[k, ], rev(res_variable[k + 1, ])),
        col = fills[k], border = NA
      )
    }

    # Add custom "Horizon" label inside the plot
    usr <- par("usr")
    text(
      x = usr[2] - 0.3 * (usr[2] - usr[1]),
      y = usr[3] + 0.07 * (usr[4] - usr[3]),
      labels = "Horizon",
      adj = c(1, 0), xpd = TRUE
    )

    # Legend adjustment: fit better inside multi-plot layouts
    legend("topleft",
      horiz = FALSE, bty = "n",
      pch = 22, pt.bg = fills,
      col = "black", legend = names.shock,
      xpd = TRUE, pt.cex = 1.2, cex = 1
    )
  }

  invisible(NULL)
}
