# Build the redistributable teaching datasets used by IdSS.
#
# Raw inputs in this directory are snapshots of public U.S. government data:
# - usmonthly_fred.csv: INDPRO, UNRATE, CPIAUCSL, PPIACO, FEDFUNDS,
#   BOGNONBR, RESBALNS, and M1SL, retrieved through FRED;
# - gs1_fred.csv: GS1, retrieved through FRED;
# - bq_gnp_quarterly.csv and bq_unrate_monthly.csv: GNPC96 and UNRATE,
#   retrieved through FRED;
# - brw_shock_series.csv: Bu, Rogers, and Wu monetary-policy shocks,
#   downloaded from the Board of Governors of the Federal Reserve System.
# - ebp_csv.csv: the Favara--Gilchrist--Lewis--Zakrajsek excess bond premium,
#   downloaded from the Board of Governors of the Federal Reserve System.

monthly <- read.csv("data-raw/OpenTeachingData/usmonthly_fred.csv",
                    na.strings = c("", "."), check.names = FALSE)
monthly$DATES <- as.Date(monthly$observation_date)

gs1 <- read.csv("data-raw/OpenTeachingData/gs1_fred.csv",
                na.strings = c("", "."), check.names = FALSE)
gs1$DATES <- as.Date(gs1$observation_date)
monthly <- merge(monthly, gs1[c("DATES", "GS1")], by = "DATES", all.x = TRUE)

# The official CSV contains the revised series on the left and the original
# release on the right, separated by a blank column. Retain the revised monthly
# measure and convert labels such as 1994m1 to first-of-month dates.
brw.raw <- read.csv("data-raw/OpenTeachingData/brw_shock_series.csv",
                    na.strings = c("", "."), check.names = FALSE,
                    skip = 0, fill = TRUE)
brw.label <- trimws(as.character(brw.raw[[1L]]))
keep <- grepl("^[0-9]{4}m[0-9]{1,2}$", brw.label)
brw.label <- brw.label[keep]
brw.value <- suppressWarnings(as.numeric(brw.raw[[2L]][keep]))
parts <- strsplit(brw.label, "m", fixed = TRUE)
brw <- data.frame(
  DATES = as.Date(sprintf("%04d-%02d-01",
                          as.integer(vapply(parts, `[`, "", 1L)),
                          as.integer(vapply(parts, `[`, "", 2L)))),
  BRW = brw.value
)

monthly <- merge(monthly, brw, by = "DATES", all.x = TRUE)
ebp <- read.csv("data-raw/OpenTeachingData/ebp_csv.csv",
                na.strings = c("", "."), check.names = FALSE)
ebp$DATES <- as.Date(ebp$date)
monthly <- merge(monthly, ebp[c("DATES", "ebp")], by = "DATES", all.x = TRUE)
monthly <- monthly[monthly$DATES >= as.Date("1959-04-01") &
                     monthly$DATES <= as.Date("2015-01-01"), ]

positive.log <- function(x) {
  x[x <= 0] <- NA_real_
  log(x)
}

USmonthly <- data.frame(
  DATES = monthly$DATES,
  LIP = positive.log(monthly$INDPRO),
  UNEMP = monthly$UNRATE,
  LCPI = positive.log(monthly$CPIAUCSL),
  LPCOM = positive.log(monthly$PPIACO),
  FFR = monthly$FEDFUNDS,
  NBR = positive.log(monthly$BOGNONBR),
  TTR = positive.log(monthly$RESBALNS),
  M1 = positive.log(monthly$M1SL),
  GS1 = monthly$GS1,
  EBP = monthly$ebp,
  BRW = monthly$BRW
)

first_difference <- function(x) c(NA_real_, diff(x))
USmonthly$d.LIP <- first_difference(USmonthly$LIP)
USmonthly$d.LCPI <- first_difference(USmonthly$LCPI)
USmonthly$d.LPCOM <- first_difference(USmonthly$LPCOM)
USmonthly$d.GS1 <- first_difference(USmonthly$GS1)

USmonthlyExample <- USmonthly[
  USmonthly$DATES >= as.Date("1965-01-01") &
    USmonthly$DATES <= as.Date("1995-06-01"), ]
USmonthly_1990_to_2012 <- USmonthly[
  USmonthly$DATES >= as.Date("1990-05-01") &
    USmonthly$DATES <= as.Date("2012-06-01"), ]

save(USmonthly, file = "data/USmonthly.rda", compress = "xz")
save(USmonthlyExample, file = "data/USmonthlyExample.rda", compress = "xz")
save(USmonthly_1990_to_2012,
     file = "data/USmonthly_1990_to_2012.rda", compress = "xz")

# Open reconstruction of the Blanchard--Quah teaching example.
gnp <- read.csv("data-raw/OpenTeachingData/bq_gnp_quarterly.csv",
                na.strings = c("", "."), check.names = FALSE)
gnp$Date <- as.Date(gnp$observation_date)

unemp.m <- read.csv("data-raw/OpenTeachingData/bq_unrate_monthly.csv",
                    na.strings = c("", "."), check.names = FALSE)
unemp.m$Date <- as.Date(unemp.m$observation_date)
unemp.m$quarter <- paste(format(unemp.m$Date, "%Y"),
                         (as.integer(format(unemp.m$Date, "%m")) - 1L) %/% 3L + 1L,
                         sep = "Q")
unemp.q <- aggregate(UNRATE ~ quarter, data = unemp.m, FUN = mean, na.rm = TRUE)
unemp.q$Date <- as.Date(paste0(substr(unemp.q$quarter, 1L, 4L), "-",
                              sprintf("%02d", 1L +
                                        3L * (as.integer(substr(unemp.q$quarter, 6L, 6L)) - 1L)),
                              "-01"))

bq <- merge(gnp[c("Date", "GNPC96")], unemp.q[c("Date", "UNRATE")],
            by = "Date", all = FALSE)
bq <- bq[bq$Date >= as.Date("1948-01-01") &
           bq$Date <= as.Date("1987-10-01"), ]
bq$Dgdp <- 100 * c(NA_real_, diff(log(bq$GNPC96)))
break_at <- 104L
bq$Dgdp[seq_len(break_at)] <- bq$Dgdp[seq_len(break_at)] -
  mean(bq$Dgdp[seq_len(break_at)], na.rm = TRUE)
if (nrow(bq) > break_at) {
  after <- (break_at + 1L):nrow(bq)
  bq$Dgdp[after] <- bq$Dgdp[after] - mean(bq$Dgdp[after], na.rm = TRUE)
}
bq$unemp <- residuals(lm(UNRATE ~ seq_len(nrow(bq)), data = bq))
BQ <- bq[!is.na(bq$Dgdp), c("Date", "Dgdp", "unemp")]

save(BQ, file = "data/BQ.rda", compress = "xz")
