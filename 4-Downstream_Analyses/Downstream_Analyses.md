---
title: "Downstream_Analyses"
author: "Enrico"
date: "2022-10-17"
output: html_document
---

### Generate raw Files

To perform downstream analysis in R, I first of all need to extract the data from candidate SNPs in [PLINK](https://linkinghub.elsevier.com/retrieve/pii/S0002929707613524)'s RAW format.

Based on different combinations of predictors and conditioning used during the genome scans I run a custom script for creating [neutral](code/get_neutral_vcf.sh) and [candidate](code/get_candidate_vcf.sh) VCF files, and another [script](code/create_candidate_and_neutral_raw.sh) for creating RAW files

```{bash}
# For vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC1+PC2
4-Downstream_Analyses/code/get_neutral_vcf.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC1+PC2

4-Downstream_Analyses/code/get_candidate_vcf.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC1+PC2

4-Downstream_Analyses/code/create_candidate_and_neutral_raw.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC1+PC2

###

# For vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC2
4-Downstream_Analyses/code/get_neutral_vcf.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC2

4-Downstream_Analyses/code/get_candidate_vcf.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC2

4-Downstream_Analyses/code/create_candidate_and_neutral_raw.sh \
  vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC2

###

# For vars_bio9+bio2+bio16+jan_mean_depth.cond_PC2

4-Downstream_Analyses/code/get_neutral_vcf.sh \
  vars_bio9+bio2+bio16+jan_mean_depth.cond_PC2

4-Downstream_Analyses/code/get_candidate_vcf.sh \
  vars_bio9+bio2+bio16+jan_mean_depth.cond_PC2

4-Downstream_Analyses/code/create_candidate_and_neutral_raw.sh \
  vars_bio9+bio2+bio16+jan_mean_depth.cond_PC2
```

### Gene Enrichment

I use a custom [bash script](code/get_gene_list_from_bed.sh) to extract the list of genes found within or around (± 20 kbp) the candidate windows using [betools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).

```{bash}
4-Downstream_Analyses/code/get_gene_list_from_bed.sh \
  3-Identifying_Candidate_Loci/tables/vars_bio9+bio2+bio15+bio16+mean_snow_days.cond_PC1+PC2.overlap_sd_win_nodup.bed \
  > 4-Downstream_Analyses/tables/candidate_genes_list.txt
```

### PCA and Distance calculations

Principal component analysis and calculation of pairwise genetic distances between samples based on neutral and adaptive SNPs are performed in a custom [R script](code/run_pca_distances.R) (some remnants of data explorations with different methods are commented out).

The scripts outputs PC1-PC2 and PC3-PC4 biplots of neutral and adaptive PCAs, as well as a plot of the per-population distributions of pairwise genetic distances.

### Generalized Dissimilarity Modeling

I ran Generalized Dissimilarity Modeling (GDM) using a custom [R script](code/run_gdm.R) that uses the package [gdm](https://CRAN.R-project.org/package=gdm).

I followed guidelines found in the [original article](https://onlinelibrary.wiley.com/doi/10.1111/ele.12376) and a [follow up guide](https://onlinelibrary.wiley.com/doi/10.1111/ele.12376) on its applications.

It outputs maps of neutral and adaptive genetic turnover based on the first three principal components of the GDM transformed raster layers, the associated biplots, a zoom of the Latvia and NE-Poland populations, and the residuals of Procrustes superimposition of neutral vs adaptive turnovers.
