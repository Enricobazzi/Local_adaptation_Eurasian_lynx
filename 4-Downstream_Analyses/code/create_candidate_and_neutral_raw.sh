#!/bin/bash

# on genomics-a
cd /home/ebazzicalupo/Local_adaptation_Eurasian_lynx

# define vcf from which to extract candidate and neutral SNPs
sd_cand_vcf=4-Downstream_Analyses/tables/ll_wholegenome_LyCa_ref.sorted.filter7.${1}.sd_cand.vcf
qvals_cand_vcf=4-Downstream_Analyses/tables/ll_wholegenome_LyCa_ref.sorted.filter7.${1}.qvals_cand.vcf
vcf_neutral=4-Downstream_Analyses/tables/ll_wholegenome_LyCa_ref.sorted.filter7.${1}.neutral.vcf

# I have a samplestoremove.txt that includes samples I want to remove from the analysis:
# cat samplestoremove.txt
# # c_ll_ba_0216 c_ll_ba_0216 0 0 0 -9
# # c_ll_ba_0233 c_ll_ba_0233 0 0 0 -9
# # c_ll_cr_0211 c_ll_cr_0211 0 0 0 -9
# # h_ll_ba_0214 h_ll_ba_0214 0 0 0 -9
# # h_ll_ba_0215 h_ll_ba_0215 0 0 0 -9

# sd candidate snps:
echo "extracting sd candidate raw for ${1}"
plink_1.9 --vcf $sd_cand_vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove 4-Downstream_Analyses/tables/samplestoremove.txt --geno 0 \
 --recode A --out 4-Downstream_Analyses/tables/${1}.sd_candidate_snps

# qvals candidate snps:
echo "extracting qvals candidate raw for ${1}"
plink_1.9 --vcf $qvals_cand_vcf \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove 4-Downstream_Analyses/tables/samplestoremove.txt --geno 0 \
 --recode A --out 4-Downstream_Analyses/tables/${1}.qvals_candidate_snps
 
# neutral snps - prune first:
echo "extracting neutral raw for ${1}"
plink_1.9 --vcf $vcf_neutral \
--double-id --allow-extra-chr --set-missing-var-ids @:# --geno 0 --maf 0.05 \
--remove 4-Downstream_Analyses/tables/samplestoremove.txt --indep 100 10 2 \
--out 4-Downstream_Analyses/tables/${1}.neutral_pruned
# neutral snps:
plink_1.9 --vcf $vcf_neutral \
 --double-id --allow-extra-chr --set-missing-var-ids @:# \
 --remove 4-Downstream_Analyses/tables/samplestoremove.txt \
 --extract <(shuf --random-source=<(yes 42) 4-Downstream_Analyses/tables/neutral_pruned.prune.in | head -n 10000) \
 --recode A --out 4-Downstream_Analyses/tables/${1}.neutral_snps

