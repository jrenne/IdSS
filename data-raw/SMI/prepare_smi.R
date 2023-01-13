


setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/SMI")

smi <- read.csv("SMI.csv")

#smi[,2:7] <- as.numeric(as.matrix(smi[,2:7]))

smi$Date <- as.Date(smi$Date,"%m/%d/%y")

save(smi,file="../../data/smi.rda")


