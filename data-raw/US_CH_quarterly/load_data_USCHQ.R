
US_CH_quarterly <- read.csv2("data-raw/US_CH_quarterly/US_CH_quarterly.csv")

US_CH_quarterly$observation_date <- as.Date(US_CH_quarterly$observation_date)
names(US_CH_quarterly)[1] <- "date"

save(US_CH_quarterly,file="data/US_CH_quarterly.rda")

