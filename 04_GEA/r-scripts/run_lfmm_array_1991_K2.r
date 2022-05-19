#!/projects/mgdesaix@colostate.edu/mambaforge/envs/lea/bin/Rscript

# library(vegan)
library(LEA)

myargs <- commandArgs(trailingOnly=TRUE)
pred.file <- myargs[1]

lfmm.file <- "./data/1991_2020/bcrf.nokinerror.QC.01maf.ld02.sexFst001.imp.lfmm"

# lfmm.k1 <- lfmm(lfmm.file, pred.file, K = 1, CPU = 4, repetitions=5, project = "new")
lfmm.k2 <- lfmm(lfmm.file, pred.file, K = 2, CPU = 4, repetitions=5, project = "new")
