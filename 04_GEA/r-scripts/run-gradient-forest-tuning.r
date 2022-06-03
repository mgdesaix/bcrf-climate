library(gradientForest)
library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly = T)

print(args)

mtry <- as.numeric(args[1])
ntree <- as.numeric(args[2])
Env.pred.file <- "./data/pred.env.1991.csv"
year <- 1991

## allele frequency data, remove EBMO site and pop names

Gen.file <- paste0("../allele-freq/rda-lfmm-intersect/rda.lfmm.cand.",year,".af")
# ex. rda.lfmm.cand.1961.af
Gen.resp <- fread(Gen.file) %>%
  as.data.frame() %>%
  filter(V1 != "EBMO") %>%
  arrange(V1) %>%
  dplyr::select(-c(V1)) %>%
  as.matrix()
colnames(Gen.resp) <- sub("^", "V", (1:ncol(Gen.resp)))

## environmental data

Env.pred <- read.csv(Env.pred.file)

nSites <- dim(Gen.resp)[1]
nSpecs <- dim(Gen.resp)[2]
# set depth of conditional permutation
lev <- floor(log2(nSites*0.368/2))

print(mtry)
print(ntree)

bcrfForest <- gradientForest(cbind(Env.pred, Gen.resp),
                                 predictor.vars = colnames(Env.pred[,-1]),
                                 response.vars = colnames(Gen.resp),
                                 ntree = ntree,
                                 mtry = mtry,
                                 transform = NULL,
                                 compact = T,
                                 maxLevel = lev,
                                 corr.threshold = 0.5,
                                 trace = F)

out.name <- paste0("./tuning-params/bcrf.rda_lfmm_cand500.env.",year,".ntree.",ntree,".mtry.",mtry,".rds")
saveRDS(bcrfForest, file = out.name)
