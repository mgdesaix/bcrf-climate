# bcrf-climate

Bioinformatics and analysis scripts for Brown-capped Rosy-Finch (*Leucosticte australis*) climate change research that is published and open access [DeSaix et al. 2022. Diversity and Distributions](https://onlinelibrary.wiley.com/doi/full/10.1111/ddi.13628):

1.  [Data preprocessing](https://github.com/mgdesaix/bcrf-climate/blob/master/01_Preprocessing/Preprocessing.md): Includes trimming FASTQs, mapping to reference, merging, marking PCR duplicates

2.  [Variant calling](https://github.com/mgdesaix/bcrf-climate/blob/master/02_VariantCalling/Variants.md): Base quality recalibration, filtering high-quality SNPs

3.  [Population genetics things](https://github.com/mgdesaix/bcrf-climate/blob/main/03_PopulationGenetics/Popgen.md): Principal components analysis, relatedness, linkage disequilibrium, ancestry values, and more!

4.  [Identifying putative adaptive variants](https://github.com/mgdesaix/bcrf-climate/blob/main/04_GEA/GEA.md): Outlier analyses using **redundancy analysis** (RDA) and **latent factor mixed models** (LFMM)

5. [Landscape turnover of putative adaptive variation](https://github.com/mgdesaix/bcrf-climate/blob/main/05_GradientForest/gradient-forest.md): Using **gradient forest** to predict genetic composition across the breeding range

6. [Ecological niche models](https://github.com/mgdesaix/bcrf-climate/blob/main/06_EcologicalNicheModeling/enm.md): Modeling current and future habitat suitability

7. [Genomic offset](https://github.com/mgdesaix/bcrf-climate/blob/main/07_GenomicOffset/genomic-offset.md): Calculating genomic offset from gradient forest predictions of genetic composition

8. [Niche margin index](https://github.com/mgdesaix/bcrf-climate/blob/main/08_NicheMarginIndex/nmi.md) Measuring environmental distance from trained values in gradient forest to assess extent of model extrapolation


### Supercomputing acknowledgements:

“This work utilized the Summit supercomputer, which is supported by the National Science Foundation (awards ACI-1532235 and ACI-1532236), the University of Colorado Boulder, and Colorado State University. The Summit supercomputer is a joint effort of the University of Colorado Boulder and Colorado State University.”