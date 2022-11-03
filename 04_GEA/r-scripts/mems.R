## Spatial predictors with Moran Eigenvector's Maps (MEMs)

# libraries
library(codep)
library(adespatial)
library(adegraphics)

#############################################
# Brown-capped rosy-finch sampling data
# using subset of samples that aren't related

# meta data with lat/long for individual's sampling site
bcrf.individuals <- read_csv("../../data/ROFI-meta-tidy.csv") %>%
  dplyr::filter(Alpha_Code == "BCRF") %>%
  dplyr::select(ShortName, Lat, Long, Group_ID)

# this file is read in that has individuals in the same order as they are in the vcf
# easy to get sample order with: bcftools query -l input.vcf
bcrf.nokinerror <- read_csv("../data/snprelate/bcrf.vcf.nokinerror.individuals.txt") %>%
  left_join(bcrf.individuals, by = c("id" = "ShortName"))

#############################################
# MEMs

# individual coordinates in order of lat and lon
bcrf.coordinates <- as.matrix(bcrf.individuals[,c(2,3)])

# compuate spatial distances
bcrf.spatial <- gcd.hf(bcrf.coordinates)
bcrf.dbmem <- dbmem(bcrf.spatial)

# first mem explains ~48% of the geographic variation

## check importance of the different components

bcrf.moran.randtest <- moran.randtest(bcrf.dbmem)


