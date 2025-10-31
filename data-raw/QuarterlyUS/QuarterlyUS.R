
QuarterlyUS <- read.csv2("data-raw/QuarterlyUS/QVAR.csv")

QuarterlyUS$date <- as.Date(QuarterlyUS$date)

save(QuarterlyUS,file="data/QuarterlyUS.rda")
