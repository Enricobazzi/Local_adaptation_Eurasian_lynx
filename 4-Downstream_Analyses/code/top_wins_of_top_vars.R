## yea why not

library(tidyverse)
library(RColorBrewer)
library(GenomicRanges)

### ###  ### ### ### ### 
### parse arguments ####

args = commandArgs(trailingOnly=TRUE)
formula_file_line = args[1]

# read formula table
formula_file = "3-Identifying_Candidate_Loci/tables/rda_exploration_formulas.txt"
form <- read.table(formula_file, header = T, sep = "\t")

# calculate K based on PC condition
if (form[formula_file_line,2] == "PC1+PC2"){
  K = 2
} else if (form[formula_file_line,2] == "PC2"){
  K = 3
}

# read two components of the formula - [1] variables and [2] condition
v_form <- c(form[formula_file_line,1], form[formula_file_line,2])

# sd candidates table name
sd_candidates_table = paste0("vars_", v_form[1],
                             ".cond_", v_form[2],
                             ".cand_sd.K", K,
                             ".candidates.tsv")

## folder where results are stored
results_folder = "3-Identifying_Candidate_Loci/tables/"

## missing genotypes table
miss_gts_tsv = "1-Preparing_Genetic_Data/tables/finalset_missing_gts.tsv"

## genes
# genes_gff3 = "4-Downstream_Analyses/tables/lc4.genes.gff3"
# 
# genes = read.table(genes_gff3, header = FALSE, sep = "\t",
#                    nrows = 19340)[,c(1,4,5,9)]
# 
# colnames(genes) = c("chromosome", "start", "end", "gene")
# 
# genes$gene = sapply(str_split(sapply(str_split(genes$gene, ";"), "[[", 5), "="), "[[", 2)

## outlier SNPs with standard deviation (loose)
outlier_snps_sd = read.table(paste0(results_folder,
                                    sd_candidates_table),
                             header = T)
# changing snp-name values because alt allele don't always match
snps = c()
poss = c()
chrs = c()
for (snp in outlier_snps_sd$snp){
  # snp name = chr:pos
  snp = strsplit(snp, "_")
  snp = rev(rev(snp[[1]])[-1])
  snp = paste(snp, collapse = "_")
  snps = c(snps, snp)
  # position
  snp_pos = unlist(strsplit(snp, ":"))[2]
  poss = c(poss, snp_pos)
  # chromosome
  snp_chr = unlist(strsplit(snp, ":"))[1]
  chrs = c(chrs, snp_chr)
}
outlier_snps_sd$snp = snps
outlier_snps_sd$chromosome = chrs
outlier_snps_sd$position = as.numeric(poss)

## read missing data
miss_gts = read.table(miss_gts_tsv, header = F, sep = "\t",
                      comment.char = "",
                      col.names = c("snp", "missing"))

# add missing data
outlier_snps_sd = merge(outlier_snps_sd, miss_gts, by = "snp")

# => variables who's BayPass result windows should be merged
vars = strsplit(v_form[1], split = "+", fixed = T)[[1]]

# iterate through variables 
for (var in vars){
  print(paste("reading windows of", var))
  # outliers for variable
  var_wins = read.table(paste0(results_folder, var,
                               "_genwin_windows.tsv"), header = T) %>%
    filter(outlier == "yes")
  
  # number of overlapping snps
  n_snps_sd = c()
  
  # iterate through windows
  for (w in 1:nrow(var_wins)){
    message('\r', round((w / nrow(var_wins) * 100), 2), appendLF = FALSE)
    # extract window
    window = var_wins[w,]
    # overlapping sd snps
    ol_sd = outlier_snps_sd %>% 
      filter(chromosome == window$scaffold & position >=  window$WindowStart & position <= window$WindowStop) %>%
      nrow()
    n_snps_sd = c(n_snps_sd, ol_sd)
  }
  var_wins$overlap_sd <- n_snps_sd
  
  top_adaptive_wins <- filter(var_wins, overlap_sd > 0)
  
  write.table(top_adaptive_wins,
              paste0("4-Downstream_Analyses/tables/",
                     var, "_top_adaptive_wins.tsv"),
              row.names = F, sep = "\t", col.names = T, quote = F)
  
  # genes_vec <- c()
  # for (w in 1:nrow(top_adaptive_wins)){
  #   message('\r', round((w / nrow(top_adaptive_wins) * 100), 2), appendLF = FALSE)
  #   window = top_adaptive_wins[w,]
  #   # overlapping genes (± 20kbp)
  #   ol_genes = genes %>%
  #     filter(chromosome == window$scaffold & start >=  window$WindowStart - 20000 & end <= window$WindowStop + 20000)
  #   # add gene
  #   genes_vec <- c(genes_vec, ifelse(nrow(ol_genes) > 0, ol_genes$gene, NA))
  # }
  # 
  # top_adaptive_wins$gene <- genes_vec
  
  top_adaptive_wins_bed <- data.frame(
    chr = top_adaptive_wins$scaffold,
    start = top_adaptive_wins$WindowStart,
    end = top_adaptive_wins$WindowStop
  )
  
  write.table(top_adaptive_wins_bed,
              paste0("4-Downstream_Analyses/tables/",
                     var, "_top_adaptive_wins.bed"),
              row.names = F, sep = "\t", col.names = F, quote = F)
}

# mean_snow_days
# COP1:
# https://onlinelibrary.wiley.com/doi/full/10.1111/jpi.12340
# PDZRN4:
# https://www.sciencedirect.com/science/article/pii/S1751731121001841

# bio2
# CADM2:
# https://www.sciencedirect.com/science/article/pii/S221287781730772X
# PDZRN4:
# https://www.sciencedirect.com/science/article/pii/S1751731121001841
# REV3L:
# https://www.sciencedirect.com/science/article/pii/S0960982200007259
# PTPRM:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8977644/
