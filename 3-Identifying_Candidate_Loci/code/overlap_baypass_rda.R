### load libraries ###
library(tidyverse)
library(GenomicRanges)

# => variables who's baypass result windows should be merged
vars = c("bio9", "bio2", "bio16", "jan_mean_depth")

### load outlier windows ###
results_folder = "3-Identifying_Candidate_Loci/tables/"
outlier_windows = data.frame()
all_windows = data.frame()
for (var in vars){
  # outliers
  var_wins = read.table(paste0(results_folder, var,
                               "_genwin_windows_outliers.bed"),
                        col.names = c("chromosome", "start", "end"))
  var_wins$var = var
  outlier_windows = rbind(outlier_windows, var_wins)
  # all windows
  var_wins = read.table(paste0(results_folder, var,
                               "_genwin_windows.tsv"),
                        header = T)
  all_windows = rbind(all_windows, var_wins)
}

# sorted bed
sorted_outlier_windows = outlier_windows[order(outlier_windows[,1],
                                               outlier_windows[,2],
                                               outlier_windows[,3]), ]
# merge overlapping windows
gr = reduce(GRanges(
  seqnames = Rle(sorted_outlier_windows$chromosome),
  ranges = IRanges(sorted_outlier_windows$start, end = sorted_outlier_windows$end),
  score = c(sorted_outlier_windows$chromosome)
  ))

a = seqnames(gr)
chrs = c()
for (i in 1:length(a@lengths)){
  chrs = c(chrs, rep(as.character(a@values)[i], a@lengths[i]))
}

merged_outlier_windows = data.frame(chromosome = chrs, 
                                    start = b@start,
                                    end = b@start + b@width - 1) 

sorted_all_windows = all_windows[order(all_windows[,1],
                                       all_windows[,2],
                                       all_windows[,3]), ]

### outlier SNPs SD ###
outlier_snps_sd = read.table(paste0(results_folder,
      "vars_bio9-bio2-bio16-janmeandepth.cond_PC2.cand_sd.K3.candidates.tsv"),
      header = T)

outlier_snps_sd_bed = data.frame()
for (s in 1:nrow(outlier_snps_sd)){
  # snp info
  snp = outlier_snps_sd[s,"snp"]
  snp_chr = unlist(strsplit(snp,":"))[1]
  snp_pos = as.numeric(unlist(strsplit(unlist(strsplit(snp,":"))[2], "_"))[1])
  snp_load = outlier_snps_sd[s, "loading"]
  # overlap
  ol = merged_outlier_windows %>%
    filter(chromosome == snp_chr & start - 5000 <= snp_pos & end + 5000 >= snp_pos)
  # check number of overlapping windows
  if (nrow(ol) == 0){
    snp_ol = "none"
  } else {
    snp_ol = paste(snp_chr, min(ol$start), max(ol$end), sep = "_")
  }
  snp_bed = data.frame(chr = snp_chr, position = snp_pos,
                       load = snp_load, overlap = snp_ol)
  outlier_snps_sd_bed = rbind(outlier_snps_sd_bed, snp_bed)
}

overlap_sd_snps = outlier_snps_sd_bed %>%
  filter(overlap != "none")
overlap_sd_snps_nodup = overlap_sd_snps[order(overlap_sd_snps[,'overlap'],
                                              -overlap_sd_snps[,'load']),]
overlap_sd_snps_nodup = overlap_sd_snps_nodup[!duplicated(overlap_sd_snps_nodup$overlap),]

# save table
write.table(overlap_sd_snps_nodup,
            "3-Identifying_Candidate_Loci/tables/overlap_sd_snps_nodup.tsv",
            row.names = F, sep = "\t")

# save bed
overlap_sd_snps_nodup_bed = data.frame(chromosome = overlap_sd_snps_nodup$chr,
                                       start = overlap_sd_snps_nodup$position - 1,
                                       end = overlap_sd_snps_nodup$position)

write.table(overlap_sd_snps_nodup_bed,
            "3-Identifying_Candidate_Loci/tables/overlap_sd_snps_nodup.bed",
            row.names = F, sep = "\t", col.names = F, quote = F)


### outlier SNPs qvals ###
outlier_snps_qvals = read.table(paste0(results_folder,
   "vars_bio9-bio2-bio16-janmeandepth.cond_PC2.cand_qvals.K3.candidates.tsv"),
   header = T)

outlier_snps_qvals_bed = data.frame()
for (s in 1:nrow(outlier_snps_qvals)){
  # snp info
  snp = outlier_snps_qvals[s,'snp']
  snp_chr = unlist(strsplit(snp,':'))[1]
  snp_pos = as.numeric(unlist(strsplit(unlist(strsplit(snp,":"))[2], "_"))[1])
  snp_load = outlier_snps_qvals[s, 'q.values']
  # overlap
  ol = merged_outlier_windows %>%
    filter(chromosome == snp_chr & start - 5000 <= snp_pos & end + 5000 >= snp_pos)
  # check number of overlapping windows
  if (nrow(ol) == 0){
    snp_ol = "none"
  } else {
    snp_ol = paste(snp_chr, min(ol$start), max(ol$end), sep = "_")
  }
  snp_bed = data.frame(chr = snp_chr, position = snp_pos,
                       load = snp_load, overlap = snp_ol)
  outlier_snps_qvals_bed = rbind(outlier_snps_qvals_bed, snp_bed)
}

overlap_qvals_snps = outlier_snps_qvals_bed %>%
  filter(overlap != "none")
overlap_qvals_snps_nodup = overlap_qvals_snps[order(overlap_qvals_snps[,'overlap'],
                                              -overlap_qvals_snps[,'load']),]
overlap_qvals_snps_nodup = overlap_qvals_snps_nodup[!duplicated(overlap_qvals_snps_nodup$overlap),]

# save table
write.table(overlap_qvals_snps_nodup,
            "3-Identifying_Candidate_Loci/tables/overlap_qvals_snps_nodup.tsv",
            row.names = F, sep = "\t")

# save bed
overlap_qvals_snps_nodup_bed = data.frame(chromosome = overlap_qvals_snps_nodup$chr,
                                       start = overlap_qvals_snps_nodup$position - 1,
                                       end = overlap_qvals_snps_nodup$position)

write.table(overlap_qvals_snps_nodup_bed,
            "3-Identifying_Candidate_Loci/tables/overlap_qvals_snps_nodup.bed",
            row.names = F, sep = "\t", col.names = F, quote = F)


## PLOTS ##

sorted_all_windows
merged_outlier_windows
overlap_sd_snps_nodup
overlap_qvals_snps_nodup