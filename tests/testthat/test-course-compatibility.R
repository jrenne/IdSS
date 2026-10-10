test_that("the public interface required by both books is exported", {
  required <- c(
    "compute.garch", "d.func.2.minimize", "estim.SVAR.ICA", "estim.armax",
    "func.2.minimize", "gICA",
    "make.Asympt.Cov.delta", "make.LPIV.irf", "make.M", "make.PHI",
    "make.jorda.irf", "make.variance.decompo", "NW.LongRunVariance",
    "sim.arma", "simul.distri", "simul.VAR", "simul.VARMA",
    "svar.fevmax", "svar.iv", "svar.ordering", "svar.signs", "tsls",
    "variance.decomp"
  )
  expect_setequal(getNamespaceExports("IdSS"), required)
})

test_that("all course datasets are installed", {
  required <- c(
    "BQ", "international", "JST", "levpan", "QuarterlyUS", "Ramey",
    "smi", "US_CH_quarterly", "US3var", "USmonthly",
    "USmonthly_1990_to_2012", "USmonthlyExample", "USyields"
  )
  available <- data(package = "IdSS")$results[, "Item"]
  expect_true(all(required %in% available))
})

test_that("ARMA simulation reproduces elementary impulse responses", {
  irf <- sim.arma(c = 0, phi = 0.5, theta = 1, sigma = 1,
                  T = 5, y.0 = 0, make.IRF = 1)
  expect_equal(as.numeric(irf), 0.5^(0:4), tolerance = 1e-12)

  set.seed(10)
  paths <- sim.arma(c = 0, phi = 0.5, theta = 1, sigma = 1,
                    T = 8, y.0 = 0, nb.sim = 3)
  expect_identical(dim(paths), c(8L, 3L))
})

test_that("ARMAX fitting returns the documented components", {
  set.seed(11)
  y <- as.numeric(sim.arma(0.1, 0.4, 1, 0.5, 80, 0))
  fit <- estim.armax(y, 1, 0)
  expect_named(fit, c("c", "phi", "theta", "sigma", "beta", "max.logl",
                      "I", "st.dev", "THETA"))
  expect_length(fit$phi, 1)
  expect_gt(abs(fit$sigma), 0)
})

test_that("GARCH recursion returns a positive conditional variance", {
  set.seed(12)
  x <- rnorm(100)
  result <- compute.garch(c(1, 0.1, 0.8), x, m0 = 1, r0 = 1)
  expect_length(result$h, length(x))
  expect_true(all(result$h[is.finite(result$h)] > 0))
  expect_true(is.finite(result$logf))
})

test_that("VAR companion construction has the expected block form", {
  phi <- list(diag(c(0.5, 0.4)), diag(c(0.1, 0.2)))
  companion <- make.PHI(phi)
  expect_identical(dim(companion), c(4L, 4L))
  expect_equal(companion[1:2, ], cbind(phi[[1]], phi[[2]]))
  expect_equal(companion[3:4, 1:2], diag(2))
})

test_that("VARMA simulation accepts the historical Mu intercept alias", {
  model <- list(
    Mu = c(0, 0),
    Phi = array(diag(c(0.5, 0.4)), dim = c(2, 2, 1)),
    Theta = array(diag(c(0.1, 0.2)), dim = c(2, 2, 1)),
    C = diag(2),
    distri = list(type = c("gaussian", "gaussian"))
  )
  result <- simul.VARMA(model, nb.sim = 5, Y0 = c(0, 0),
                        eta0 = c(1, 0), indic.IRF = 1)
  expect_identical(dim(result$Y), c(2L, 5L))
  expect_equal(result$Y[, 1], c(1, 0))
})

test_that("SVAR ordering accepts the historical bootstrap argument", {
  set.seed(13)
  y <- matrix(rnorm(80), ncol = 2)
  result <- svar.ordering(
    y, p = 1, nb.periods.IRF = 3,
    nb.bootstrap.replications = 2, indic.plot = 0
  )
  expect_equal(dim(result$simulated.IRFs), c(3L, 2L, 2L))
})

test_that("every exported object has an installed help alias", {
  documented <- vapply(getNamespaceExports("IdSS"), function(topic) {
    length(utils::help(topic, package = "IdSS")) > 0
  }, logical(1))
  expect_true(all(documented), info = paste(names(documented)[!documented],
                                            collapse = ", "))
})
