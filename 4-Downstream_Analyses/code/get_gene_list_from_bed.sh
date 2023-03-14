#!/bin/bash

# script to get the list of genes from the Canada lynx reference genes which have
# any overlap with a list of genomic coordinates from a BED file

# overlap is calculated 20kbp up- and down- stream of windows

# first argument is the bed file
bed=${1}

# genes gff3
genes_gff3=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.NCBI.nr_main.gff3

# print genes
bedtools intersect \
 -a $genes_gff3 \
 -b <(awk '{if ($2-20000 < 0) print $1, 0, $3+20000; else print $1, $2-20000, $3+20000}' $bed | tr ' ' '\t') |
 | grep -oE 'gene=[[:alnum:]]+;' | uniq |
 tr ';' ',' | sed 's/gene=//' | sed '$s/,$//'
