# Identifying putative adaptive variation through genetic-environment association (GEA) analyses

To identify putative adaptive variation, I identified candidate SNPs with redundancy analysis (RDA) and latent factor mixed models (LFMM). The intersection of these candidate SNP sets were used for subsequent analysis with gradient forest.

## Environmental variable selection

Current and projected climate data for North America were obtained from the [AdaptWest Project](https://adaptwest.databasin.org/pages/adaptwest-climatena/)(AdaptWest Project 2021). The data set included 33 bioclimate variables (described [here](http://s3-us-west-2.amazonaws.com/www.cacpd.org/CMIP6/ReadMe_ClimateNA_CMIP6.txt)) and elevation.

Correlated environmental variables were identified by $|r| < 0.75$ based on the environment values extracted for an individual's sampling location. Correlated variables were reduced by selecting environmental variables based on the biological relevance to the Rosy-Finch (see [environmental variable script](./r-scripts/environmental_variable_selection.R)). This resulted in a final data set of 4 variables:

- Mean temperature of the warmest month (MWMT)
- Precipitation as snow (PAS)
- Elevation (Elev)
- Summer heat moisture index (SHM)

The filtered genotypes of the 104 unrelated BCRF individuals were used as the genetic data for these analyses. This data set had 2,367,581 SNPs.

## Redundancy analysis (RDA)

First, I created spatial variables to use to condition the RDA model on, as a method to account for population structure. Then, I performed RDA using individual genotypes as the response variable and the environmental variables identified above as the predictor variables.

### Moran Eigenvector Maps (MEMs)
To account for population structure in the RDA, I created Moran Eigenvector Maps (MEMs) as spatial variables in the R package `adespatial`. In short, I used the lat/lon coordinates of the BCRF sampling locations to create spatial distances, which, in turn, were used to create the MEMs. Additional explanation of MEMs can be found [here by Stéphane Dray](https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html). My code for creating MEMs for Brown-capped Rosy-Finch is [here](./r-scripts/mems.R).

### Running RDA

RDA was performed on the cluster in R (see [RDA R script](./r-scripts/run_rda.1961_1990.R)). The general model used was:

```
Genotypes ~ Condition(MEM1) + Elev + MWMT + PAS + SHM
```

## Latent factor mixed models (LFMM)

To correct for population structure in LFMM, I calculated ancestry values with the `snmf` function from the `LEA` R package ([see script](./r-scripts/run_snmf.R)). The results of this showed that K=1 was the best option, therefore latent factors were set as K=1 in the subsequent runs of LFMM.








## References

AdaptWest Project. 2021. Gridded current and projected climate data for North America at 1km resolution, generated using the ClimateNA v7.01 software (T. Wang et al., 2021). Available at adaptwest.databasin.org.