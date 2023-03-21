#!/bin/bash

# script to get the supplementary table with:
# adaptive windows coordinates,
# the variable from which they were calculated as outliers in BayPass+Genwin,
# their support (in terms of Wstat and number of RDA snps)
# the genes found within

# overlap is calculated 20kbp up- and down- stream of windows

genes_gff3=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3

for var in bio9 bio2 bio15 bio16 mean_snow_days
 do
  echo "variable $var"
  
  echo "genes" > tmp_${var}_genes
  
  bed=4-Downstream_Analyses/tables/${var}_top_adaptive_wins.bed
  tsv=4-Downstream_Analyses/tables/${var}_top_adaptive_wins.tsv
  
  for i in $(seq 1 $(wc -l < $bed))
   do
   
   if [ $i -eq $(wc -l < $bed) ]
    then
     echo "window : ${i}"
   else
     echo -ne "window : ${i} \r"
   fi
   
   bedtools intersect -wa -wb \
    -a <(sed "${i}q;d" $bed | awk '{if ($2-20000 < 0) print $1, 0, $3+20000; else print $1, $2-20000, $3+20000}' | tr ' ' '\t') \
    -b $genes_gff3 > tmp_hits
   
   if [ $(wc -l < tmp_hits) -gt 0 ]
    then
     grep -oE 'gene=[[:alnum:]]+;' tmp_hits | sort -u | sed 's/gene=//' |
      sed 's/;//' | paste -sd'|' - >> tmp_${var}_genes
   else
     echo "NA" >> tmp_${var}_genes
   fi
   
  done
  
  paste $tsv tmp_${var}_genes > tmp && mv tmp $tsv
  
  rm tmp_*
  
done

cat 4-Downstream_Analyses/tables/*_top_adaptive_wins.tsv |
 cut -f1,2,3,5,8,10,11 | head -1 \
 > 4-Downstream_Analyses/tables/adaptive_windows_with_genes.tsv

cat 4-Downstream_Analyses/tables/*_top_adaptive_wins.tsv |
 cut -f1,2,3,5,8,10,11 | grep -v "overlap_sd" \
 >> 4-Downstream_Analyses/tables/adaptive_windows_with_genes.tsv
