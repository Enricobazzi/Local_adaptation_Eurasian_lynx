#!/bin/bash

# to get the vcf of candidate SNPs I can intersect the vcf with the bed of 
# candidate snps
vcf=1-Preparing_Genetic_Data/tables/ll_wholegenome_LyCa_ref.sorted.filter7.vcf

sd_cand_bed=3-Identifying_Candidate_Loci/tables/${1}.overlap_sd_snps_nodup.bed
sd_cand_vcf=4-Downstream_Analyses/tables/ll_wholegenome_LyCa_ref.sorted.filter7.${1}.sd_cand.vcf

qvals_cand_bed=3-Identifying_Candidate_Loci/tables/${1}.overlap_qvals_snps_nodup.bed
qvals_cand_vcf=4-Downstream_Analyses/tables/ll_wholegenome_LyCa_ref.sorted.filter7.${1}.qvals_cand.vcf

# sd candidates:
echo "extracting sd candidate vcf for ${1}"
bedtools intersect -header \
    -a $vcf \
    -b $sd_cand_bed \
    > $sd_cand_vcf

# qvals candidates:
echo "extracting qvals candidate vcf for ${1}"
bedtools intersect -header \
    -a $vcf \
    -b $qvals_cand_bed \
    > $qvals_cand_vcf
