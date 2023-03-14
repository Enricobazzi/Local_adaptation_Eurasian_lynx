---
title: "Identifying Candidate Loci"
author: "Enrico"
date: "2022-10-17"
output: html_document
---

GENERAL STRATEGY : Identifying candidate loci using the combination of BayPass univariate analysis of correlation between population allele frequencies and individual variables and RDA multivariate analysis of individual's genotypes and all of the selected variables together.

### step 1 - Variable Selection

Following the [tutorial](https://github.com/Capblancq/RDA-landscape-genomics) associated with the work by [Capblancq & Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13722) we use an approach based on forward model selection.

This [script](code/decide_variables_to_analyze.R) was used to run the model selection and store some of its results. In general two sets of variables make most sense to me:

- bio2 + bio9 + bio16 + jan_mean_depth : approach focused on minimizing collinearity between variables -> discards bio15 and mean_snow_days in favor of just jan_mean_depth to reduce correlation between variables, with r < 0.8 and VIF < 4
- bio2 + bio9 + bio15 + bio16 + mean_snow_days : approach focused on maximizing variance explained by the model, which increases with the inclusion of bio15 and with mean_snow_days compared to jan_mean_depth, while adding a bit of collinearity r < 0.85 and VIF < 6 (?? double check)

### step 2 - Variance Partitioning

In order to analyze the way variance in the RDA model is partitioned between genetic structure, geographic structure and environmental structure, we followed the same [tutorial](https://github.com/Capblancq/RDA-landscape-genomics) to create summary tables of the variance explained by:
- full model => selected variables + genetic structure + geography
- pure climate model => selected variables - genetic structure - geography
- pure genetic structure model => genetic structure - selected variables - geography
- pure geographic structure model => geography - selected variables - genetic structure
using this [script](code/variance_partitioning_in_rda_model.R)

### step 3 - Univariate Analysis with BayPass

Univariate GEA was performed with the software [Baypass](https://academic.oup.com/genetics/article/201/4/1555/5930067). This software uses Allele Frequency data and Bayesian Hierarchical Models to generate a distribution of differentiation coefficients of SNPs across the genome. First a scaled population covariance matrix based on population allele frequencies is generated, to evaluate the relationships between populations. Based on this covariance matrix, and the supposed ancestral allele frequency (inferred from weighted mean of the reference allele frequency), a CORE model is generated.

SNPs that present values of differentiation (XtX, a SNP-specific Fst explicitly corrected for the scaled covariance of population allele frequencies) that exceed the amount expected under the core model can be identified as candidate selection loci.

Furthermore, this software allows the evaluation of the association of particular SNPs to some environmental Covariate. With one measurement for population for each covariate, the software will evaluate the data under to additional models:

(1) The Standard Covariate model (STD): which adds an association covariable (given by the correlation coefficient between the covariate measurement and XtX) to the CORE model. (2) The Auxiliary Variable Covariate model (AUX): which further builds on the STD model by attaching a binary variable (0 or 1 -> association or no association) to each locus' regression coefficient. The posterior mean of this variable will indicate the posterior probability of the association of that variable with a particular SNP.

In order to run the software we will need:

– Allele Count data for all the considered populations: in the form of a space delimited file, with one row for each SNP and two columns for each population, with the allele counts for the reference and alternative alleles (obtained using this [script](code/Allele_Count_generation.sh)).

– Environmental variable files: in the form of different files with one column per population and their respective measurement for that particular variable (see 2-Preparing_Environmental_Data for how it was prepared).

To avoid possible bias due to consecutive SNPs being non-independant observations because of Linkage Disequilibrium, I will divide my dataset into 50 different dataset, made of one SNP every 50 the following way:

```{bash}
# For the awk script to work I need to iterate from 0 to 49
for n in {0..49}
 do
  echo ${n}
  awk -v number="$n" 'NR % 50 == 0+number' all.allelecounts \
  > all.allelecounts.${n}
done

# But the 0 iteration actually starts from the 50th line of the file
# so I will change the name of the file:
mv all.allelecounts.0 all.allelecounts.50
```

Finally I can run baypass for each of the environmental variables using a custom made [script](code/baypass_aux_v2.sh) in the following way:

```{bash}
cd /home/ebazzicalupo/Local_adaptation_Eurasian_lynx

for var in bio9 bio2 bio15 bio16 mean_snow_days jan_mean_depth
 do
  echo "analyzing association with ${var}"
  screen -dmS aux_${var}  sh -c "3-Identifying_Candidate_Loci/code/baypass_aux_v2.sh ${var}; exec /bin/bash"
done
```

### step 4 - Genomic Windows with GenWin

GenWin is a program that trys to detect inflection points in summary statistics calculated in a locus per locus manner (e.g. Fst). Using these inflection points it will divide the genome in windows which are not of arbitrary size, but based on properties of the summary statistic values. The program will output the number of SNPs inside each window and give a summary value based on the spline (Wstat). Wstat is a value that depends on the mean of the summary statistic, weighted against the mean and standard deviation of the dataset and the number of SNPs inside each window.

Wstat can be then used to calculate outlier windows with a quantile criteria (e.g. above 99.9% as done in the GenWin paper).

Using a custom [R script](code/run_genwin.R) for each selected variable, I will use it to divide the Lynx lynx genome in windows based on the Bayes Factor summary statistic from BayPass. Outlier windows are saved as TSV and BED for downstream analyses.

```{bash}
for var in bio9 bio2 bio15 bio16 mean_snow_days jan_mean_depth
 do
  echo "generating genwin windows from baypass results of ${var}"
  Rscript 3-Identifying_Candidate_Loci/code/run_genwin.R ${var}
done
```

### step 5 - Multivariate Analysis ReDundancy Analysis (RDA)

To run multivariate RDA using the pruned SNPs dataset as the response variables, the selected environmental predictors as explanatory variables and conditioning based on neutral genetic structure, I used the package vegan in a custom [R script](code/run_rda_exploration.R). 

This script allows the use of multiple combinations of environmental predictors and genetic structure conditioning which were saved in a formula file.

The combinations I wanted to try found in the formula file were:
 - line 1 : bio9+bio2+bio15+bio16+mean_snow_days	PC1+PC2
 - line 2 : bio9+bio2+bio15+bio16+mean_snow_days	PC2
 - line 8 : bio9+bio2+bio16+jan_mean_depth	PC2

Outlier SNPs are saved in both TSV and BED formats. RDA biplots and Manhattan plots are also produced.

```{bash}
for line in 1 2 8
 do
  echo "running RDA with line $line of formula file"
  Rscript 3-Identifying_Candidate_Loci/code/run_rda_exploration.R $line
done
```

### step 6 - Overlap of Multivariate and Univariate results

In order to find overlap between the Baypass/Genwin windows and RDA SNPs I wrote an ambitious, poorly written and computationally inefficient [R script](code/overlap_baypass_rda.R).

As the Previous RDA script it works on different combinations of variables and conditioning based on the same formula file.

It generates:
 - a summary table of the overlap between the univariate and multivariate outliers
 - tables in TSV and BED format listing overlapping candidate windows (windows from univariate method with at least 1 SNP from multivariate)
 - Manhattan plots of the univariate windows of each environmental predictor, showing which are outliers, and which of the outliers have at least one SNP from univariate analysis

```{bash}
for line in 1 2 8
 do
  echo "overlap between RDA and BayPass with line $line of formula file"
  Rscript 3-Identifying_Candidate_Loci/code/overlap_baypass_rda.R $line
done
```

