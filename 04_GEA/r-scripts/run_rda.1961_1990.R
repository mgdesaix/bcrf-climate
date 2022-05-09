#!/projects/mgdesaix@colostate.edu/mambaforge/envs/lea/bin/Rscript

library(vegan)
library(data.table)

# read in large genetics file
resp.file <- "./data/bcrf.nokinerror.snps.9replaced.space.transposed.txt"
resp <- fread(resp.file)

# read in environmental and mem predictors
pred.file <- "./data/pred.env.1961_1990.mems.csv"
pred_and_meta <- read.csv(pred.file)
pred <- pred_and_meta[,-c(1,2)]


############ run conditional RDA
#

rda.con2 <- rda(resp ~ Condition(MEM1) + Elev + MWMT + PAS + SHM,data = pred, scale=T)
# rda.con2 <- readRDS("./output/1961_1990/rda.conditionMEM1.MWMT.PAS.Elev.SHM.rds")
saveRDS(rda.con2, "./output/1961_1990/rda.conditionMEM1.MWMT.PAS.Elev.SHM.rds")

rda.con2.sum <- summary(rda.con2)$concont
write.csv(rda.con2.sum, "./output/1961_1990/rda.conditionMEM1.MWMT.PAS.Elev.SHM.sum.csv")

load.rda1 <- summary(rda.con2)$species[,1:2]
cand1 <- outliers(load.rda1[,1],3)
cand2 <- outliers(load.rda1[,2],3)
rda.cand <- c(names(cand1), names(cand2))
rda.cand <- rda.cand[!duplicated(rda.cand)]

saveRDS(rda.cand, "./output/1961_1990/rda.conditionMEM1.MWMT.PAS.Elev.SHM.candSNPs.rds")