
library(readxl)
levpan <- read_excel("data-raw/LevPan/levpan_replicationdata.xls")
plot(levpan$date,levpan$tfp_lev)

save(levpan,file="data/levpan.rda")

