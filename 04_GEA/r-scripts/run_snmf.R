#!/projects/mgdesaix@colostate.edu/mambaforge/envs/lea/bin/Rscript

library(LEA)

vcf <- "../snps/recalibrated/gatherVcfs/bcrf/norelate/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.vcf"
output <- "./output/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.geno"

vcf2geno(input.file = vcf, output.file = output)

snmf2 <- LEA::snmf(output, K=1:6, ploidy=2, entropy=T, alpha=100, project="new")
saveRDS(snmf2, "./output/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.geno.snmf_k1_6.rds")