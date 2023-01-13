
# Load BQ data.
# Source of the data: https://econpapers.repec.org/software/bocbocode/rtz00017.htm

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/BQ")

bqdata <- read.csv("bqdata.csv", sep=";")
bqdata$DATE <- as.Date(bqdata$DATE,"%Y-%m-%d")

plot(bqdata$DATE,bqdata$LHMUR,type="l")

real.GNP <- bqdata$GNP/bqdata$GD87
T <- length(real.GNP)

gnp.growth <- real.GNP/(c(NaN,real.GNP[1:(T-1)])) - 1

y <- gnp.growth
y[1:104] <- y[1:104] - mean(y[1:104],na.rm=TRUE)
y[105:T] <- y[105:T] - mean(y[105:T],na.rm=TRUE)

u <- bqdata$LHMUR
trend <- 1:T

eq <- lm(u ~ trend)

u <- eq$residuals

y <- cbind(100*y,u)
y <- y[2:T,]

BQ <- data.frame(Date=bqdata$DATE[2:T],Dgdp=y[,1],unemp=y[,2])
save(BQ,file="../../data/BQ.rda")


