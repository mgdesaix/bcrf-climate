# Variant calling in GATK

Software I use:

- GATK4 (ver 4.1.4.0)

- bcftools (ver 1.9)

- samtools (ver 1.9)

- htslib (ver 1.9)

Below I outline my workflow for obtaining a final VCF file.

**1.** Preparing genomic intervals

**2.** Base quality recalibration

**3.** HaplotypeCaller

**4.** GatherVCFs


## 1) Preparing genomic intervals

The [HaplotypeCaller function from GATK](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) calls SNPs. This process can be very slow, so the best way to do this is to run parallel jobs for different regions of the genome, instead of the entire genome at once. Some reference genomes can have extremely large scaffolds (>100 MB) and these will cause Haplotypecaller to take much longer than 24 hours. I like to parallelize jobs so they take less than 24 hours and to do this I provide HaplotypeCaller with genomic region files of 2-4 MB (This can be changed easily to your own needs). I will demonstrate how to provide HaplotypeCaller with genomic regions of ~2MB.

I first create a file `scaffold-lengths` that has two columns: 1) Scaffold names, 2) Scaffold lengths. This is generally straightforward to do by using the FASTA reference genome file. For example, my ROFI reference genome headers look like this:

```sh
% head -n 1 leucosticte_australis_final_assembly.fasta
>Scaffold_1__1_contigs__length_151072562
```

The `faidx` function in [samtools](http://www.htslib.org/doc/samtools-faidx.html) can be used to get the scaffold length of the reference genome (Note: Even though the scaffold names have a length listed, this doesn't necessarily align with the actual scaffold length produced from the index `faidx` function. This will produce downstream errors if you try to use the length from the scaffold name).

After extracting the scaffold name and length, I then can get a file like this:

```sh
% head -n 3 scaffold-lengths.txt
>Scaffold_1__1_contigs__length_151072562	151072562
>Scaffold_2__1_contigs__length_114774641	114774641
>Scaffold_3__1_contigs__length_112408489	112408489
```

Then in `awk`, I can create a new file that breaks up large scaffolds into smaller regions. In this example, I break them up into 1 MB chunks (*side note:* I will ultimately be analyzing regions in 2 MB chunks as mentioned above, but this avoids potentially combining chunks slightly less than 2 MB which could result in a chunk of almost 4 MB).

```sh
cat scaffold-lengths.txt | awk '{tot=0; while(tot<$2) {start=tot+1; tot+=1e6; if(tot>$2) tot=$2; printf("%d\t%s:%d-%d\n",++n, $1,start,tot);}}' > ./intervals/intervals1mb.txt
```

Now I have a file with the large scaffolds broken up into regions:

```sh
--% head -n 3 intervals1mb.txt
Scaffold_1__1_contigs__length_151072562:1-1000000
Scaffold_1__1_contigs__length_151072562:1000001-2000000
Scaffold_1__1_contigs__length_151072562:2000001-3000000
```

Regions of a scaffold are specified by `SCAFFOLD_NAME:START-END`, where START and END are the genomic positions within the given scaffold.

But the later scaffolds are really small and there are lots of these genomic regions, and I don't want to run 11,488 jobs when I could combine the smaller scaffolds to total closer to my 2 MB length objective:

```sh
--% tail -n 3 intervals1mb.txt
Scaffold_10476__1_contigs__length_666:1-666
Scaffold_10477__1_contigs__length_638:1-638
Scaffold_10478__1_contigs__length_603:1-603
--% wc -l intervals1mb.txt
11488 intervals1mb.txt
```

So, again using the magic of awk, the goal is to read down the list of broken up scaffolds and their lengths, and save regions that total to 2 MB in a file. Initial files will just have 2 of the 1 MB genomic regions from a scaffold, but later files may have lots of the small scaffolds to total to 2 MB.

By calculating length from the start and end of each region, I can produce a scaffold of the genomic regions and the length:

```sh
--% head -n 1 rofi-scaffold-interval-lengths-1mb.txt
Scaffold_1__1_contigs__length_151072562:1-1000000	1000000
Scaffold_1__1_contigs__length_151072562:1000001-2000000	1000000
Scaffold_1__1_contigs__length_151072562:2000001-3000000	1000000
```

Cool. So now all I need is the awk code to read through the lengths of the second column, and spit out the scaffold names/regions into files until that file totals to 2 MB, and then it will spit out scaffold names/regions to a new file, on and on.

```sh
 awk 'BEGIN {id=1; tot=0} {if(tot<2e6) {tot+=$2; k=sprintf("%03d",id); print $1 >> "interval2mb_" k ".list"} else {id+=1; tot=$2; k=sprintf("%03d", id); print $1 >> "interval2mb_" k ".list";}}' rofi-scaffold-interval-lengths-1mb.txt
```

**Note:** The *.list* suffix is necessary for GATK - otherwise headaches insue.

```sh
--% ls interval2mb_* | tail -n 3
interval2mb_517.list
interval2mb_518.list
interval2mb_519.list
```

The first files are short, with just the two 1 MB regions, and the last files are long with lots of short scaffolds:

```sh
--% wc -l interval2mb_001.list
2 interval2mb_001.txt
--% wc -l interval2mb_519.list
531 interval2mb_519.txt
```

## 2) Base quality score recalibration

Base quality score recalibration (BQSR) is used to detect systematic errors during sequencing, and is highly recommended by the folks at GATK. For more details, see their [technical documentation for BQSR](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-).

For the Brown-capped Rosy-Finch data, I first used two algorithms to call variants (GATK and samtools) and took the intersection of these variant sets. I then used this set of variants to perform base quality score recalibration of the original BAM files. Finally, using the recalibrated BAM files, I then proceeded to create the final variant data set.

To determine a set of high confidence variants, I used the `mpileup` function from `bcftools` to implement the Samtools algorithm for variant calling (see [mpileup script](./slurm-scripts/get_mpileup_regions.sh)). I also used GATK's algorithm implemented in the `HaplotypeCaller` function (see details below). Both of these variant data sets were filted with `bcftools` by only including genotypes with <25% missing data and filtered by depth, specifically twice the average depth and 3/4s of the average depth (between 4.5 and 12; see example [filtering script](./get_filtered.sh)). I then used the `bcftools isec` function to find the intersection of the 2 variant data sets, which produced a VCF file of 396,021 variants (see [intersection script](./slurm-scripts/get_intersection.sh)).

I then followed the [GATK workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-) for base recalibration, which uses the two functions `BaseRecalibrator` and `ApplyBQSR` (see [BQSR script](./slurm-scripts/get_baserecalibrated.sh)).

## 3) HaplotypeCaller

To determine my final set of genomic variants, I used GATK's [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller). I ran `HaplotypeCaller` across the genomic regions described earlier (see [HaplotypeCaller script](./slurm-scripts/get_haplotypecaller.sh). This can take a while.

## 4) GatherVCFs

I then used GATK's [GatherVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037422071-GatherVcfs-Picard-) to combine all the individual VCF files of the different genomic intervals (see [GatherVCFs script](./slurm-scripts/get_gatherVcfs.sh)). I conducted final filtering for SNPs in this last data set using the following parameters in `bcftools view`:

- 2 alleles (-m 2 -M 2)
- SNPs (--type snps)
- Missing genotypes < 10% (-i 'F_MISSING < 0.1')
- Minor allele frequency of < 0.1 (--min-af 0.1 --max-af 0.1)
- Depth of less than 2x the average depth and greater than average depth based on the number of missing individuals allowed. With mean depth of 6 and 10% missing allowed, the minimum was 6*0.9=5.4 (-i 'AVG(FORMAT/DP)>5.4 & AVG(FORMAT/DP)<12')


