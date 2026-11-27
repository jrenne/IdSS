# Offline reconstruction from the reviewed compact extraction; run at package root.
input <- "data-raw/USmonetaryRR/rr-monetary.csv"
expected_md5 <- "eed207f3cf12a5a8afce4e0742c3b5da"
stopifnot(identical(unname(tools::md5sum(input)), expected_md5))
USmonetaryRR <- read.csv(input, stringsAsFactors = FALSE)
USmonetaryRR$date <- as.Date(USmonetaryRR$date)
stopifnot(
  identical(names(USmonetaryRR), c(
    "date", "policy_rate", "unemployment",
    "inflation", "rr_shock"
  )),
  nrow(USmonetaryRR) == 466L,
  identical(USmonetaryRR$date, seq(as.Date("1969-03-01"),
    as.Date("2007-12-01"),
    by = "month"
  )),
  !anyNA(USmonetaryRR), all(USmonetaryRR$unemployment >= 0),
  all(USmonetaryRR$unemployment <= 100)
)
save(USmonetaryRR, file = "data/USmonetaryRR.rda", version = 2, compress = "xz")
