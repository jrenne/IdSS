# Asymptotic covariance of the ICA rotation.

#' Asymptotic Covariance for ICA Rotation Parameters
#'
#' Computes the sandwich covariance matrix used for the free rotation
#' parameters in the non-Gaussian SVAR identification procedure.
#'
#' @param eps Numeric matrix of estimated structural shocks.
#' @param distri Distribution specification list used by [gICA()].
#' @param C Numeric square ICA rotation matrix.
#' @return A symmetric covariance matrix of dimension
#'   `n * (n - 1) / 2`, where `n = ncol(eps)`.
#' @examples
#' set.seed(1)
#' eps <- matrix(rnorm(300), 100, 3)
#' distri <- list(
#'   type = rep("gaussian", 3), mu = rep(0, 3),
#'   sigma = rep(1, 3), df = rep(NA, 3), p = rep(NA, 3)
#' )
#' V <- make.Asympt.Cov.delta(eps, distri, diag(3))
#' dim(V)
#' @export
make.Asympt.Cov.delta <- function(eps, distri, C) {
  n <- dim(eps)[2]
  T <- dim(eps)[1]
  A <- make.A.matrix(eps, distri, C)
  Omega <- make.Omega(eps, distri)
  Omega.augment <- matrix(0, n^2, n^2)
  Omega.augment[1:(n * (n - 1) / 2), 1:(n * (n - 1) / 2)] <- Omega
  A_1 <- solve(A)
  mat.cov <- 1 / T * A_1 %*% Omega.augment %*% t(A_1)
  return(mat.cov)
}

#' @keywords internal
#' @noRd
make.g.stars <- function(eps, distri) {
  res.g <- gICA(eps, distri, indic.deriv = 1)
  g.1.star <- apply(res.g$d.log.L, 2, mean)
  g.2.star <- apply(res.g$d.log.L^2, 2, mean)
  g.3.star <- apply(eps * res.g$d.log.L, 2, mean)
  g.4.star <- apply(res.g$d2.log.L, 2, mean)
  return(list(g.1.star = g.1.star, g.2.star = g.2.star, g.3.star = g.3.star, g.4.star = g.4.star))
}

#' @keywords internal
#' @noRd
make.Omega <- function(eps, distri) {
  n <- dim(eps)[2]

  all.g.stars <- make.g.stars(eps, distri)
  g.1 <- matrix(all.g.stars$g.1.star, ncol = 1)
  g.2 <- matrix(all.g.stars$g.2.star, ncol = 1)
  g.3 <- matrix(all.g.stars$g.3.star, ncol = 1)
  vec.1 <- matrix(1, n, 1)

  g.1.g.1 <- g.1 %*% t(g.1)

  Omega <- g.1.g.1 %x% diag(n)

  for (i in 1:(n - 1)) {
    aux <- diag(n) %x% g.1.g.1[, i]
    Omega[(i * n + 1):n^2, ((i - 1) * n + 1):(i * n)] <-
      Omega[(i * n + 1):n^2, ((i - 1) * n + 1):(i * n)] - aux[(i * n + 1):n^2, ]
  }

  # Add lower-tri parts of block diagonal:
  Omega <- Omega + diag(n) %x% g.1.g.1

  # Remove all entries on and above diag:
  diag(Omega) <- 0
  Omega[upper.tri(Omega)] <- 0

  # Make it symmetric:
  Omega <- Omega + t(Omega)

  # Fill diagonal:
  aux <- g.2 %*% t(vec.1) + vec.1 %*% t(g.2) - 2 * (g.3 %*% t(vec.1)) * (vec.1 %*% t(g.3))
  diag(Omega) <- c(aux)

  # Remove rows and columns to get a [n(n-1)/2] x [n(n-1)/2] matrix:
  aux <- upper.tri(diag(n))
  diag(aux) <- TRUE
  indic.rows.to.keep <- which(!aux)

  Omega <- Omega[indic.rows.to.keep, ]
  if (is.matrix(Omega)) {
    Omega <- Omega[, indic.rows.to.keep]
  } else {
    Omega <- Omega[indic.rows.to.keep]
  }

  return(Omega)
}

#' @keywords internal
#' @noRd
make.A.matrix <- function(eps, distri, C) {
  n <- dim(eps)[2]

  all.g.stars <- make.g.stars(eps, distri)
  g.1 <- matrix(all.g.stars$g.1.star, ncol = 1)
  g.2 <- matrix(all.g.stars$g.2.star, ncol = 1)
  g.3 <- matrix(all.g.stars$g.3.star, ncol = 1)
  g.4 <- matrix(all.g.stars$g.4.star, ncol = 1)
  vec.1 <- matrix(1, n, 1)

  A.1 <- matrix(0, n^2, n^2)
  A.2 <- matrix(0, n^2, n^2)
  for (i in 1:n) {
    # compute the a_{i,.}:
    A_i <- (-g.4[i] + g.3 %*% t(vec.1)) * t(C)
    # Put these in A:
    A.1[((i - 1) * n + 1):(i * n), ((i - 1) * n + 1):(i * n)] <- A_i
    A.2[((i - 1) * n + 1):(i * n), ((i - 1) * n + 1):(i * n)] <- t(C)
  }
  for (j in 1:n) {
    # compute the a_{.,j}:
    A_j <- (-g.4 %*% t(vec.1) + g.3[j]) * (vec.1 %*% matrix(C[, j], nrow = 1))
    # Put these in A (diagonal blocks):
    for (i in 1:n) {
      A.1[((j - 1) * n + i), ((i - 1) * n + 1):(i * n)] <- -A_j[i, ]
      A.2[((j - 1) * n + i), ((i - 1) * n + 1):(i * n)] <- c(C[, j])
    }
  }
  # remove rows from A.1
  aux <- upper.tri(diag(n))
  diag(aux) <- TRUE
  indic.rows.to.keep <- which(!aux)
  A.1 <- A.1[indic.rows.to.keep, ]
  A.2 <- A.2[indic.rows.to.keep, ]

  A.3 <- matrix(0, n, n^2)
  for (i in 1:n) {
    A.3[i, ((i - 1) * n + 1):(i * n)] <- c(C[, i])
  }

  A <- rbind(A.1, A.2, A.3)

  return(A)
}
