

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/USyields")

library(fredr)
fred_api_key <- Sys.getenv("FRED_API_KEY")
if (!nzchar(fred_api_key)) {
  stop("Set FRED_API_KEY in the environment before refreshing USyields.")
}
fredr_set_key(fred_api_key)
start_date <- as.Date("1990-01-01"); end_date <- as.Date("2022-01-01")
f <- function(ticker){
  fredr(series_id = ticker,
        observation_start = start_date,observation_end = end_date,
        frequency = "m",aggregation_method = "avg")
}

yd6m <- f("DGS6MO")
yd1 <- f("DGS1")
yd2 <- f("DGS2")
yd3 <- f("DGS3")
yd5 <- f("DGS5")
yd7 <- f("DGS7")
yd10 <- f("DGS10")
yd20 <- f("DGS20")
yd30 <- f("DGS30")

USyields <- data.frame(date=yd1$date,
                       Y1 = yd1$value,
                       Y2 = yd2$value,
                       Y3 = yd3$value,
                       Y5 = yd5$value,
                       Y7 = yd7$value,
                       Y10 = yd10$value,
                       Y20 = yd20$value,
                       Y30 = yd30$value)

save(USyields,file="../../data/USyields.rda")



