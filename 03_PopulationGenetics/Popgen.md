# Population Genetics

## Linkage disequilibrium

SNPs with LD > 0.2 were pruned from the data set using [Plink2](https://www.cog-genomics.org/plink/2.0/ld)

## Principal components analysis

PCA was initially used to visualize individual genetic variation (see [PCA script](./r-scripts/run_snprelate.r)) and was implemented with [SNPRelate](https://github.com/zhengxwen/SNPRelate). The results showed minimal genetic population structure, and highlighted the potential for **related** individuals by principal components that had only a pair of individuals being separated from the rest of the individuals.

## Identify related individuals

I used the [KING software](https://www.kingrelatedness.com/manual.shtml) to identify related individuals up to the 2nd degree (see [relatedness script](./slurm-scripts/get_related.sh)). To do so, I used all SNPs that passed quality filters, but were not filtered by minor allele frequency, based on the recommendations of the developers. The results of this analysis resulted in the removal of 12 individuals due to relatedness up to the 2nd degree.


## Sex-linked variation

Interestingly, PCA showed evidence of separation of individuals by the sexes. This was made more apparent when ancestry values were computed for K=2 with `snmf` from the `LEA` R package (see [snmf script](./r-scripts/run_snmf.r). No geographic structure was evident with the K=2 bar plot, but comparing the ancestry Q values with sex revealed clear separation of male vs female.

![](images/bcrf.k2.barplot.png?raw=true)

![](images/bcrf.sex.k2.png?raw=true)

To remove sex-linked variation (autosomal or sex chromosomes), I used the ancestry values from `snmf` to classify individuals as male or female, and then calculated pairwise Fst values for all SNPs in `vcftools` (see [script](./slurm-scripts/get_pairwise_fsts.sh). SNPs were filtered out that had > 0.01 FST between males and females (see [filtering script](./slurm-scripts/get_sex_linked/removed.sh))

## Population structure

After removing sex-linked variation, I checked population structure with principal components analysis (again with `SNPRelate` in R), ancestry values (again with `snmf` from `LEA` package in R), and also checked isolation by distance.

For isolation by distance, I used pairwise FSTs to derive the linearized FST (i.e. $\frac{\text{FST}}{1-\text{FST}}$) among sampling sites, and compared to geographic distance. I formally tested the association between genetic and geographic distances with a Mantel test in the R package `ade4`. The results did not suggest a pattern of isolation by distance, but the plotting of sites revealed that the Pike's Peak site (site ID = PIPE) had elevated genetic differentiation with other sites.

![](images/ibd_PIPE.png?raw=true)

PCA also revealed weak clustering of PIPE individuals from other sampling sites. Ancestry values had the lowest cross-validation error for K=1, but K=2 had only slightly higher error and showed PIPE individuals as separate. These results provided evidence that some weak population structure should be accounted for in the subsequent [genetic environment association analyses](https://github.com/mgdesaix/bcrf-climate/blob/main/04_GEA/GEA.md)


