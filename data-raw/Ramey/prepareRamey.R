

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/Ramey")

Ramey <- read.csv("Ramey2.csv")

Ramey$DATES <- as.Date(Ramey$DATES,"%m/%d/%Y")

save(Ramey,file="../../data/Ramey.rda")
