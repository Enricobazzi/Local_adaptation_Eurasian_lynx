#!/bin/bash

# Create a Genome region file for Bedtools:
# A file with the list of chromosomes as col1 and their length as col2,
# tab separated
# Basically the first two columns of a FAI file:
cut -f1,2 $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.fa.fai > \
$LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.genome

# Bedtools random to generate file of 200 random segments of 100000 bp
# Output a BED file:
bedtools random -l 100000 -n 200 -g $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.genome | \
sort > $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.200x100kbp.genome.bed

# I also want to remove the repetetive regions from this subset,
# as depth at those positions doesn't matter and might bias 
# our final calculations
bedtools subtract -a $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.200x100kbp.genome.bed \
-b $LUSTRE/LL_selection/Ref_Genome_LyCa/repetitive_regions/lc_rep_ALL_scaffold_coord.bed \
> $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.200x100kbp.masked.genome.bed

# Array of populations except mongolia
nomopopARRAY=($(ls $LUSTRE/LL_selection/LyCaRef_bams/*.namelist | rev | cut -d'/' -f1 | rev | grep -v "all-samples" | cut -d '.' -f1 | sort -u | grep -v "mo"))
# Generate population bamlist except mongolia
for pop in ${nomopopARRAY[@]}
  do
  echo "generating bamlist of ${pop}"
  ls $LUSTRE/LL_selection/LyCaRef_bams/*${pop}*_indelrealigner.bam |
  grep -vE "c_ll_cr_0212|c_ll_ki_0090|c_ll_vl_0112|c_ll_ya_0146|c_ll_ba_0224|c_ll_ca_0240|c_ll_ca_0249|c_ll_ca_0253" \
  > $LUSTRE/LL_selection/LyCaRef_bams/${pop}.depth.bamlist
done

# Generate mongolia bamlist
mopopARRAY=($(cat $LUSTRE/LL_selection/LyCaRef_bams/mo.namelist | cut -d'_' -f3 | sort -u))
for pop in ${mopopARRAY[@]}
  do
  echo "generating bamlist of mongolia - adding ${pop}"
  ls $LUSTRE/LL_selection/LyCaRef_bams/*${pop}*_indelrealigner.bam |
  grep -vE "c_ll_cr_0212|c_ll_ki_0090|c_ll_vl_0112|c_ll_ya_0146|c_ll_ba_0224|c_ll_ca_0240" \
  >> $LUSTRE/LL_selection/LyCaRef_bams/mo.depth.bamlist
done

# Array of bamlist files
BAMlistArray=($(ls $LUSTRE/LL_selection/LyCaRef_bams/*.depth.bamlist | rev | cut -d'/' -f 1 | rev))

# Loop of Samtools depth calculations for each bamlist
for bamlist in ${BAMlistArray[@]}
  do
  echo "Calculating depth for ${bamlist}"
  samtools depth -a -b $LUSTRE/LL_selection/Ref_Genome_LyCa/lc4.200x100kbp.masked.genome.bed \
  -f $LUSTRE/LL_selection/LyCaRef_bams/${bamlist} \
  > $LUSTRE/LL_selection/SamTools_Depth/${bamlist}.200x100kbp.masked.depth &
done
