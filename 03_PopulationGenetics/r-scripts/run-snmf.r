library(LEA)

vcf <- "../variants/bcrf.nokinerror.QC.01maf.ld02.sexFst001.vcf"
output <- "../variants/bcrf.nokinerror.QC.01maf.ld02.sexFst001.geno"

vcf2geno(input.file = vcf, output.file = output)

snmf2 <- LEA::snmf(output, K=1:6, ploidy=2, entropy=T, alpha=100, project="new")
saveRDS(snmf2, "./output/bcrf.nokinerror.QC.01maf.ld02.sexFst001.geno.snmf_k1_6.rds")
