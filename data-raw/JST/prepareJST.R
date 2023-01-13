

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/JST")

JST <- read.csv("JSTdatasetR6.csv")

save(JST,file="../../data/JST.rda")
