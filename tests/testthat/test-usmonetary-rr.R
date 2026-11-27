test_that("narrative-shock data retain their calendar, schema and source evidence", {
  d <- USmonetaryRR
  expect_identical(names(d), c(
    "date", "policy_rate", "unemployment",
    "inflation", "rr_shock"
  ))
  expect_identical(dim(d), c(466L, 5L))
  expect_identical(d$date, seq(as.Date("1969-03-01"),
    as.Date("2007-12-01"),
    by = "month"
  ))
  expect_false(anyNA(d))
  expect_true(all(d$unemployment >= 0 & d$unemployment <= 100))
  expect_true(all(vapply(d[-1], is.numeric, logical(1))))
  evidence <- system.file("data-sources", "USmonetaryRR", package = "IdSS")
  expect_true(nzchar(evidence))
  for (file in c("CC0-1.0.txt", "upstream-README.txt", "provenance.json")) {
    expect_true(file.exists(file.path(evidence, file)), info = file)
  }
  expect_match(paste(readLines(file.path(evidence, "upstream-README.txt")), collapse = "\n"),
    "RRCGShock", fixed = TRUE)
  licence <- paste(readLines(file.path(evidence, "CC0-1.0.txt")), collapse = "\n")
  expect_match(licence, "CC0 1.0 Universal", fixed = TRUE)
  expect_match(paste(readLines(file.path(evidence, "provenance.json")), collapse = ""),
    "655696c1c576b7537c5a939d2c261f0a111ae663",
    fixed = TRUE
  )
})

test_that("the 24-month LP-IV teaching benchmark is reproduced", {
  d <- USmonetaryRR
  Y <- as.matrix(d[c("policy_rate", "unemployment", "inflation")])
  p <- 12L
  h <- 24L
  dates <- seq.int(p + 1L, nrow(Y) - h)
  W <- cbind(1, do.call(cbind, lapply(seq_len(p), function(j) Y[dates - j, ])))
  xr <- stats::lm.fit(W, Y[dates, 1])$residuals
  zr <- stats::lm.fit(W, d$rr_shock[dates])$residuals
  yr <- stats::lm.fit(W, Y[dates + h, 2])$residuals
  beta <- sum(zr * yr) / sum(zr * xr)

  # Independent scalar Bartlett meat, not the package's HAC implementation.
  scores <- zr * (yr - beta * xr)
  scores <- scores - mean(scores)
  N <- length(scores)
  q <- max(h + 1L, floor(4 * (N / 100)^(2 / 9)))
  meat <- sum(scores^2)
  for (lag in seq_len(q)) {
    meat <- meat + 2 * (1 - lag / (q + 1)) *
      sum(scores[(lag + 1L):N] * scores[seq_len(N - lag)])
  }
  se <- sqrt(meat) / abs(sum(zr * xr))
  X <- cbind(W[, 1], Y[dates, 1], W[, -1, drop = FALSE])
  Z <- cbind(d$rr_shock[dates], W[, -1, drop = FALSE], 1)
  fit <- tsls(as.matrix(Y[dates + h, 2]), X, Z, q)
  expect_equal(unname(fit$b.iv[2]), beta, tolerance = 1e-9)
  expect_equal(sqrt(fit$covmat.b.iv[2, 2]), se, tolerance = 1e-9)
  expect_equal(beta, 0.49542502, tolerance = 1e-7)
  expect_lt(abs(beta - stats::qnorm(0.975) * se - 0.1551499), 1e-7)
  expect_lt(abs(beta + stats::qnorm(0.975) * se - 0.8357001), 1e-7)
})
