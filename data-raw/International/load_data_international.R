
international <- read.csv2("data-raw/International/CPI_GDP.csv")

international$date <- as.Date(international$date)

save(international,file="data/international.rda")

