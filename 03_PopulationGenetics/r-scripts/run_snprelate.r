#!/projects/mgdesaix@colostate.edu/mambaforge/envs/snprelate/bin/Rscript

library(gdsfmt)
library(SNPRelate)

vcf.fn <- "/scratch/summit/mgdesaix@colostate.edu/ROFI/snps/recalibrated/gatherVcfs/bcrf/norelate/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.vcf.gz"

gds.file <- "./data/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.gatk.gds"
snpgdsVCF2GDS(vcf.fn, gds.file, method="biallelic.only")

genofile <- snpgdsOpen(gds.file)

pop_code <- scan("bcrf.vcf.populations.txt", what = character())

snpset_ld02 <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = FALSE)
snpset_ld04 <- snpgdsLDpruning(genofile, ld.threshold=0.4, autosome.only = FALSE)

snpset_ld02.id <- unlist(unname(snpset_ld02))
snpset_ld04.id <- unlist(unname(snpset_ld04))

pca_ld02 <- snpgdsPCA(genofile, snp.id=snpset_ld02.id, autosome.only = FALSE)
pca_ld04 <- snpgdsPCA(genofile, snp.id=snpset_ld04.id, autosome.only = FALSE)

saveRDS(pca_ld02, "./output/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.gatk.pca_ld02.rds")
saveRDS(pca_ld04, "./output/bcrf.nokinerror.snps.10missing.05maf.5_4depth12.het.gatk.pca_ld04.rds")
