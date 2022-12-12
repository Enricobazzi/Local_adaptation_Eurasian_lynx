#!/bin/bash

# vcf with all the samples
vcf=1-Preparing_Genetic_Data/tables/ll_wholegenome_LyCa_ref.sorted.filter7.vcf
# vcf with all the samples minus the ones to be filtered out
# (c_ll_ba_0216 c_ll_ba_0233 c_ll_cr_0211 h_ll_ba_0214 h_ll_ba_0215)
vcf_final=1-Preparing_Genetic_Data/tables/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.vcf
# list of samples to be kept and their total number
sample_list=($(grep -m1 "#CHR" ${vcf} | tr '\t' '\n' | grep "_" | grep -vE "c_ll_ba_0216|c_ll_ba_0233|c_ll_cr_0211|h_ll_ba_0214|h_ll_ba_0215"))
nsamples=($(grep -m1 "#CHR" ${vcf} | tr '\t' '\n' | grep "_" | grep -vE "c_ll_ba_0216|c_ll_ba_0233|c_ll_cr_0211|h_ll_ba_0214|h_ll_ba_0215" | wc -l))

# remove samples that will be filtered out in downstream analyses
/opt/gatk-4.1.0.0/gatk SelectVariants \
 -R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
 -V ${vcf} \
 $(for j in ${sample_list[@]}; do echo "-sn ${j}";done) \
 -O ${vcf_final}

# get the chromosome, position and alternative allele columns (for SNP IDs)
# grep -v "#" ${vcf_final} | cut -f 1,2,5 | sed 's/\t/:/' | sed 's/\t/_/' > chr_pos.tmp
# get the chromosome and position (for SNP IDs)
grep -v "#" ${vcf_final} | cut -f 1,2 | sed 's/\t/:/' > chr_pos.tmp

# get the number of missing genotypes for each SNP
grep -v "#" ${vcf_final} | cut -f8 | cut -d';' -f3 | 
   cut -d'=' -f2 | awk -v nsam="${nsamples}" '{print ((nsam*2)-$1)/2}' \
   > miss.tmp

# paste SNP IDs and Missing data
paste chr_pos.tmp miss.tmp > 1-Preparing_Genetic_Data/tables/finalset_missing_gts.tsv

# remove tmp files
rm *.tmp
