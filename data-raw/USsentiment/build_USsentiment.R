# Rebuild the redistributable teaching dataset used in the sentiment-shock
# example. Run from the IdSS package root. Downloaded files are kept in a
# temporary directory and are not committed to the package repository.

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install the 'readxl' package before rebuilding USsentiment.")
}

cache <- tempfile("IdSS-USsentiment-")
dir.create(cache)
on.exit(unlink(cache, recursive = TRUE), add = TRUE)

download <- function(url, name, mode = "wb") {
  destination <- file.path(cache, name)
  utils::download.file(url, destination, mode = mode, quiet = TRUE)
  destination
}

qdate <- function(year, quarter) {
  as.Date(sprintf("%04d-%02d-15", year, 3L * quarter))
}

quarter_number <- function(date) {
  year <- as.integer(format(date, "%Y"))
  quarter <- (as.integer(format(date, "%m")) - 1L) %/% 3L + 1L
  year + (quarter - 1L) / 4
}

fred_quarterly <- download(
  paste0(
    "https://fred.stlouisfed.org/graph/fredgraph.csv?id=",
    "GDPC1,TOTLQ,DNDGRG3Q086SBEA,DSERRG3Q086SBEA,PCND,PCESV"
  ),
  "fred_quarterly.csv", "w"
)
fred_population <- download(
  "https://fred.stlouisfed.org/graph/fredgraph.csv?id=CNP16OV",
  "fred_population.csv", "w"
)
fernald_file <- download(
  "https://www.frbsf.org/wp-content/uploads/quarterly_tfp.xlsx",
  "quarterly_tfp.xlsx"
)
tealbook_zip <- download(
  paste0(
    "https://www.philadelphiafed.org/-/media/FRBP/Assets/Surveys-And-Data/",
    "greenbook-data/GBweb/Gbweb_gRGDP_Column_Format.zip"
  ),
  "tealbook_grgdp.zip"
)
utils::unzip(tealbook_zip, exdir = file.path(cache, "tealbook"))

# Current-vintage utilization-adjusted TFP. The growth rates are annualized;
# division by 400 converts them into quarterly log changes. The log index is
# normalized to zero in 1947 Q1.
tfp <- readxl::read_excel(fernald_file, sheet = "quarterly", skip = 1)
tfp <- tfp[grepl("^[0-9]{4}:Q[1-4]$", tfp$date), ]
tfp_parts <- do.call(rbind, strsplit(tfp$date, ":Q", fixed = TRUE))
tfp$date <- qdate(as.integer(tfp_parts[, 1]), as.integer(tfp_parts[, 2]))
tfp$tfp <- cumsum(replace(tfp$dtfp_util, is.na(tfp$dtfp_util), 0)) / 400

# Current-vintage BEA/BLS series distributed through FRED.
quarterly <- utils::read.csv(fred_quarterly, na.strings = ".")
quarterly$date0 <- as.Date(quarterly$observation_date)
quarterly$date <- qdate(
  as.integer(format(quarterly$date0, "%Y")),
  (as.integer(format(quarterly$date0, "%m")) - 1L) %/% 3L + 1L
)

monthly <- utils::read.csv(fred_population, na.strings = ".")
monthly$date0 <- as.Date(monthly$observation_date)
monthly$year <- as.integer(format(monthly$date0, "%Y"))
monthly$quarter <- (as.integer(format(monthly$date0, "%m")) - 1L) %/% 3L + 1L
population <- stats::aggregate(CNP16OV ~ year + quarter, monthly, mean,
                               na.rm = TRUE)
population$date <- qdate(population$year, population$quarter)

macro <- merge(quarterly, population[, c("date", "CNP16OV")], by = "date")
macro <- macro[order(macro$date), ]
macro$ln_gdp_pc <- log(macro$GDPC1 * 1e6 / macro$CNP16OV)
macro$ln_hours_pc <- log(macro$TOTLQ * 1e6 / macro$CNP16OV)

# Törnqvist aggregation of nondurable goods and services. Nominal expenditure
# shares weight log changes in real quantities (nominal spending divided by
# its component price index); population growth is then removed. The resulting
# per-capita log index is normalized to zero in the first available quarter.
share_nd <- macro$PCND / (macro$PCND + macro$PCESV)
quantity_nd <- macro$PCND / macro$DNDGRG3Q086SBEA
quantity_sv <- macro$PCESV / macro$DSERRG3Q086SBEA
average_share_nd <- (share_nd[-1] + share_nd[-length(share_nd)]) / 2
dlog_consumption_pc <- c(
  0,
  average_share_nd * diff(log(quantity_nd)) +
    (1 - average_share_nd) * diff(log(quantity_sv)) -
    diff(log(macro$CNP16OV))
)
macro$ln_consumption_pc <- cumsum(dlog_consumption_pc)

# One-quarter-ahead Greenbook/Tealbook real-GDP-growth forecast. To obtain one
# quarterly observation, use the publication closest to the 15th day of the
# middle month of the quarter, which approximates the SPF survey timing.
tealbook_files <- list.files(file.path(cache, "tealbook"), pattern = "[.]xlsx$",
                             full.names = TRUE)
read_tealbook <- function(path) {
  x <- as.data.frame(readxl::read_excel(path), check.names = FALSE)
  names(x)[1] <- "target"
  x
}
tealbook_parts <- lapply(tealbook_files, read_tealbook)
tealbook <- Reduce(function(x, y) merge(x, y, by = "target", all = TRUE),
                   tealbook_parts)
release_names <- names(tealbook)[-1]
release_dates <- as.Date(sub("^gRGDP_", "", release_names), format = "%Y%m%d")
release_quarters <- quarter_number(release_dates)
target_number <- as.numeric(tealbook$target)
target_quarters <- floor(target_number) +
  (round(10 * (target_number - floor(target_number))) - 1) / 4

sample_quarters <- seq(1968.75, 2010.50, by = 0.25)
expectations <- lapply(sample_quarters, function(q) {
  release_index <- which(abs(release_quarters - q) < 1e-8)
  year <- floor(q)
  quarter <- round(4 * (q - year)) + 1L
  anchor <- as.Date(sprintf("%04d-%02d-15", year, 3L * quarter - 1L))
  release_index <- release_index[
    which.min(abs(as.numeric(release_dates[release_index] - anchor)))
  ]
  target_index <- which(abs(target_quarters - (q + 0.25)) < 1e-8)
  data.frame(
    date = qdate(year, quarter),
    forecast = as.numeric(tealbook[target_index, release_index + 1L])
  )
})
expectations <- do.call(rbind, expectations)
expectations$expectation <- as.numeric(scale(expectations$forecast))

USsentiment <- Reduce(function(x, y) merge(x, y, by = "date"), list(
  tfp[, c("date", "tfp")],
  macro[, c("date", "ln_consumption_pc", "ln_gdp_pc", "ln_hours_pc")],
  expectations[, c("date", "expectation")]
))
USsentiment <- USsentiment[
  USsentiment$date >= as.Date("1968-12-15") &
    USsentiment$date <= as.Date("2010-09-15"),
  c("date", "tfp", "ln_consumption_pc", "ln_gdp_pc", "ln_hours_pc",
    "expectation")
]
row.names(USsentiment) <- NULL

stopifnot(nrow(USsentiment) == 168L, !anyNA(USsentiment))
save(USsentiment, file = "data/USsentiment.rda", compress = "xz")
