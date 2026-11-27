# Optional extraction from the immutable author archive; no download is performed.
# Usage at package root: Rscript --vanilla data-raw/USmonetaryRR/extract-from-upstream.R /path/to/LP_JEL_Replication.zip
# haven is needed for this maintenance operation only, not to install/use IdSS.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !requireNamespace("haven", quietly = TRUE)) {
  stop("Supply the original ZIP path and install haven for extraction")
}
archive <- normalizePath(args[1], mustWork = TRUE)
stopifnot(identical(unname(tools::md5sum(archive)), "2ec51f1712573bc5870afb817c7ab1a3"))
folder <- tempfile("idss-rr-source-")
dir.create(folder)
member <- "LP_JEL_Replication/Example6_JointInference/data_fred.dta"
utils::unzip(archive, files = member, exdir = folder)
source_data <- as.data.frame(haven::read_dta(file.path(folder, member)))
keep <- source_data$daten >= as.Date("1969-03-01") &
  source_data$daten <= as.Date("2007-12-01")
x <- source_data[keep, ]
out <- data.frame(
  date = as.character(x$daten), policy_rate = x$ffr,
  unemployment = x$urate, inflation = x$infl, rr_shock = x$RRCGShock
)
stopifnot(nrow(out) == 466L, !anyNA(out))
write.csv(out, "data-raw/USmonetaryRR/rr-monetary.csv", row.names = FALSE)
unlink(folder, recursive = TRUE)
stopifnot(identical(
  unname(tools::md5sum("data-raw/USmonetaryRR/rr-monetary.csv")),
  "eed207f3cf12a5a8afce4e0742c3b5da"
))
