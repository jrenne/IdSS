# Load data from txt files (extracted fromn the FFED database)

library(mFilter)

indic.detrend <- 2 # 1 for linear, 2 for quadratic and 3 for cubic trend (GDP)

First.date <- "1959-04-01"
Last.date  <- "2015-01-01"

setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/US3var")

QVAR_Monthly <- read.delim("QVAR_new/QVAR_Monthly.txt")
QVAR_Quarterly <- read.delim("QVAR_new/QVAR_Quarterly.txt")

QVAR_Monthly$DATE <- as.Date(QVAR_Monthly$DATE,"%Y-%m-%d")
QVAR_Quarterly$DATE <- as.Date(QVAR_Quarterly$DATE,"%Y-%m-%d")

QVAR_Monthly$month <- as.integer(format(QVAR_Monthly$DATE,"%m"))
QVAR_Monthly <- subset(QVAR_Monthly,month %in% c(1,4,7,10))

DATA <- merge(QVAR_Monthly,QVAR_Quarterly,by="DATE",all=TRUE)

# Select first and last dates of the estimation sample:


first.date <- which(DATA$DATE==as.Date(First.date,"%Y-%m-%d"))
last.date <- which(DATA$DATE==as.Date(Last.date,"%Y-%m-%d"))

# Prepare data:

# Output:
gdp <- log(DATA$GDPC1[first.date:last.date])
trend <- 1:length(gdp)
trend2 <- trend^2
trend3 <- trend^3
if(indic.detrend == 1){
  eq <- lm(gdp ~ trend)
}else if(indic.detrend == 2){
  eq <- lm(gdp ~ trend + trend2)
}else if(indic.detrend == 3){
  eq <- lm(gdp ~ trend + trend2 + trend3)
}
y.gdp <- 100*eq$residuals

gdp.pot <- log(DATA$GDPPOT[first.date:last.date])
y.gdp.gap <- 100*(gdp - gdp.pot)

y.gdp.gap.hp <- 100*hpfilter(gdp,freq=1600)$cycle
y.gdp.gap.bp <- 10*cffilter(gdp,pl=6,pu=40)$cycle
plot(y.gdp.gap,type="l")
lines(y.gdp.gap.hp,col="red")
lines(y.gdp.gap.bp,col="blue")

# Unemployment
u.gap <- DATA$UNRATE[first.date:last.date] - DATA$NROU[first.date:last.date]
y.u.gap <- u.gap

# Inflation:
p <- log(DATA$GDPDEF)
infl <- 400 * (p - c(NaN,p[1:(length(p)-1)]))
#infl <- 100 * (p - c(rep(NaN,3),p[1:(length(p)-3)]))
infl <- infl[first.date:last.date]

# Commodities:
DATA$OILPRICE[is.na(DATA$OILPRICE)] <- DATA$MCOILWTICO[is.na(DATA$OILPRICE)]
aux <- log(DATA$OILPRICE)
commo <- 400 * (aux - c(rep(NaN,1),aux[1:(length(aux)-1)]))
commo <- 100 * (aux - c(rep(NaN,4),aux[1:(length(aux)-4)]))
#commo <- aux
commo_1 <- commo[(first.date-1):(last.date-1)]
commo_2 <- commo[(first.date-2):(last.date-2)]
commo <- commo[first.date:last.date]


# Fed fund rate:
r <- DATA$FEDFUNDS[first.date:last.date]

vec.dates <- DATA$DATE[first.date:last.date]

par(mfrow=c(2,2))
par(plt=c(.1,.9,.15,.7))
plot(vec.dates,y.gdp.gap,type="l",main="Real activity")
lines(vec.dates,y.u.gap,col="red")
abline(h=0,col="grey")
plot(vec.dates,infl,type="l",main="Inflation")
abline(h=0,col="grey")
plot(vec.dates,r,type="l",main="Federal fund rate")
plot(vec.dates,commo,type="l",main="Commodity prices")
abline(h=0,col="grey")

US3var <- data.frame(Date = vec.dates,
                     y.gdp.gap = y.gdp.gap,
                     y.u.gap = y.u.gap,
                     infl = infl,
                     r=r,commo=commo)

save(US3var,file="../../data/US3var.rda")



