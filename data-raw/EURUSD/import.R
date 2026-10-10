# Rebuild the public-domain EUR/USD dataset from the archived FRED download.
raw <- read.csv("data-raw/EURUSD/DEXUSEU.csv", na.strings = ".")
eurusd <- data.frame(
  Date = as.Date(raw$observation_date),
  USD_per_EUR = raw$DEXUSEU
)
eurusd <- eurusd[is.finite(eurusd$USD_per_EUR), ]
rownames(eurusd) <- NULL
save(eurusd, file = "data/eurusd.rda", compress = "xz")
