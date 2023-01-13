# =========================================
# Indentification of Structural Shocks
# -----------------------------------------
# Course
# Kenza Benhima and Jean-Paul Renne
# 2019
# =========================================


setwd("~/Dropbox/Teaching/MyRpackages/AEC/data-raw/USmonthly")

DATA <- read.delim("CEEmonthly_Monthly.txt")
DATA$DATE <- as.Date(DATA$DATE,"%Y-%m-%d")

# Declare data from FRED:
LIP <- log(DATA$INDPRO)
UNE <- DATA$UNRATE
LPI <- log(DATA$CPIAUCSL)
LCO <- log(DATA$PPIACO)
FFR <- DATA$FEDFUNDS
DATA$BOGNONBR[DATA$BOGNONBR<=0] <- 1
NBR <- log(DATA$BOGNONBR)
TTR <- log(DATA$RESBALNS)
M1  <- log(DATA$M1SL)

names.of.variables <- c("LIP","UNE","LPI","LCO","FFR","NBR","TTR","M1")
y <- cbind(LIP,UNE,LPI,LCO,FFR,NBR,TTR,M1)

names.of.variables <- c("NBR","TTR","M1")
y <- cbind(NBR,TTR,M1)

colnames(y)  <- names.of.variables

DATAFRAME <- as.data.frame(y)

DATAFRAME$DATES <- DATA$DATE

# Add Ramey's data:
RameyDATA <- read.csv2("Ramey_Monetarydat.csv")
RameyDATA$DATES <- as.Date(RameyDATA$DATES,"%d.%m.%Y")


DATAFRAME <- merge(DATAFRAME,RameyDATA,by="DATES",all=TRUE)
T.all <- dim(DATAFRAME)[1]


DATAFRAME$d.LIP <- NaN
DATAFRAME$d.LIP[2:T.all] <- DATAFRAME$LIP[2:T.all] - DATAFRAME$LIP[1:(T.all-1)]
DATAFRAME$d.LCPI[2:T.all] <- DATAFRAME$LCPI[2:T.all] - DATAFRAME$LCPI[1:(T.all-1)]
DATAFRAME$d.LPCOM[2:T.all] <- DATAFRAME$LPCOM[2:T.all] - DATAFRAME$LPCOM[1:(T.all-1)]
DATAFRAME$d.GS1[2:T.all] <- DATAFRAME$GS1[2:T.all] - DATAFRAME$GS1[1:(T.all-1)]

indic.first <- which(DATAFRAME$DATES==First.date)
indic.last  <- which(DATAFRAME$DATES==Last.date)

DATAFRAME <- DATAFRAME[indic.first:indic.last,]
vector.of.dates <- DATAFRAME$DATES

USmonthly <- DATAFRAME

save(USmonthly,file="../../data/USmonthly.rda")


