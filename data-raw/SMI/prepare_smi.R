# Rebuild the packaged SMI data from the archived CSV snapshot.
smi <- utils::read.csv("data-raw/SMI/SMI.csv")
smi$Date <- as.Date(smi$Date, "%m/%d/%y")
save(smi, file = "data/smi.rda", version = 3)
