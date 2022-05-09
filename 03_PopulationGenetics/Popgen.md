# Population genetics



## Principal components analysis

PCA was initially used to visualize individual genetic variation (see [PCA script](./r-scripts/run_snprelate.r)) and was implemented with [SNPRelate](https://github.com/zhengxwen/SNPRelate). The results showed minimal genetic population structure, and highlighted the potential for **related** individuals by principal components that had only a pair of individuals being separated from the rest of the individuals.

## Identify related individuals

I used the [KING software](https://www.kingrelatedness.com/manual.shtml) to identify related individuals up to the 2nd degree (see [relatedness script](./slurm-scripts/get_related.sh)). To do so, I used all SNPs that passed quality filters, but were not filtered by minor allele frequency, based on the recommendations of the developers. The results of this analysis resulted in the removal of 12 individuals due to relatedness up to the 2nd degree.