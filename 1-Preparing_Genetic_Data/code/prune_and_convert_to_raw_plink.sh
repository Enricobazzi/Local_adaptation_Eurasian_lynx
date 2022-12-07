# script to prune snps from the vcf and convert it into raw format for analysis in R
cd /home/ebazzicalupo/Local_adaptation_Eurasian_lynx/1-Preparing_Genetic_Data

## For the whole genome :

# vcf to be converted
vcf=tables/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.vcf

# prune snps 
  # handle weird chromosome names and no snp names
  # no missing data is allowed by RDA
  # prune 5kb windows, step of 5 variants, max r2 of 0.8

plink_1.9 --vcf $vcf \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --geno 0 \
  --indep-pairwise 5kb 1 0.7 \
  --out tables/finalset.maf5pc.nomiss_indep_5kb_1_07

# snps to be kept are stored in tables/finalset.maf5pc.nomiss_indep_5kb_1_07.prune.in

# generate raw file of pruned dataset
  # handle weird chromosome names and no snp names
  # use only pruned in
  # generate raw

plink_1.9 --vcf $vcf \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract tables/finalset.maf5pc.nomiss_indep_5kb_1_07.prune.in \
  --recode A --out tables/finalset.maf5pc.pruned

## For the neutral genome (removing the genes) :

# Genes and 5kbp flanking regions were defined as follows by Dani:
# cd /GRUPOS/grupolince/reference_genomes/lynx_canadensis
# 
# awk -F"\t" '$3 == "gene" {printf ("%s\t%s\t%s\n", $1, $4-5001, $5+5000)}' lc4.NCBI.nr_main.gff3 \
#  > lc4.NCBI.nr_main.genes.plus5000.temp_bed
# 
# join -1 1 -2 1 <(LANG=en_EN sort -k1,1 -k2,2n -k3,3n lc_ref_all_the_genome.bed) \
#  <(LANG=en_EN sort -k1,1 -k2,2n -k3,3n lc4.NCBI.nr_main.genes.plus5000.temp_bed) | 
#  awk -v OFS='\t' '{if ($4<0) {$4="1"} else {$4=$4}; print}' | 
#  awk -v OFS='\t' '{if ($5>$3) {$5=$3} else {$5=$5}; print $1,$4,$5}' | 
#  bedtools merge -i stdin -d 1 > lc4.NCBI.nr_main.genes.plus5000.bed
# 
# rm lc4.NCBI.nr_main.genes.plus5000.temp_bed

# vcf without genes to be converted
vcf_neutral=tables/ll_wholegenome_LyCa_ref.sorted.filter7.finalset.maf5pc.nogenes.vcf

# remove the genes from the VCF using bedtools
bedtools subtract -header -a $vcf \
 -b /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.genes.plus5000.bed \
 > $vcf_neutral

# generate raw file of pruned neutral dataset
  # handle weird chromosome names and no snp names
  # use only pruned in
  # generate raw

plink_1.9 --vcf $vcf_neutral \
  --double-id --allow-extra-chr --set-missing-var-ids @:# \
  --extract tables/finalset.maf5pc.nomiss_indep_5kb_1_07.prune.in \
  --recode A --out tables/finalset.maf5pc.pruned.nogenes
