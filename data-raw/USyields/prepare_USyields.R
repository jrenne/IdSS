

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/USyields")

library(fredr)
fredr_set_key("df65e14c054697a52b4511e77fcfa1f3")
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




