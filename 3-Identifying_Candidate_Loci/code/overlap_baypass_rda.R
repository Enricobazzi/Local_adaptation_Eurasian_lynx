### ### ### ### ### ###
### load libraries ####

library(tidyverse)
library(GenomicRanges)
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html

### ### ### ### ### ### ### #
### prepare the datasets ####

print("preparing the datasets")

## folder where results are stored
results_folder = "3-Identifying_Candidate_Loci/tables/"

## missing genotypes table
miss_gts_tsv = "1-Preparing_Genetic_Data/tables/finalset_missing_gts.tsv"

## outlier SNPs with standard deviation (loose)
outlier_snps_sd = read.table(paste0(results_folder,
                                    "vars_bio9-bio2-bio16-janmeandepth.cond_PC2.cand_sd.K3.candidates.tsv"),
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
outlier_snps_sd$position = poss

## outlier SNPs with q-values (strict)
outlier_snps_qvals = read.table(paste0(results_folder,
                                       "vars_bio9-bio2-bio16-janmeandepth.cond_PC2.cand_qvals.K3.candidates.tsv"),
                                header = T)
# changing snp-name values because alt allele don't always match
snps = c()
poss = c()
chrs = c()
for (snp in outlier_snps_qvals$snp){
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
outlier_snps_qvals$snp = snps
outlier_snps_qvals$chromosome = chrs
outlier_snps_qvals$position = poss

## add missing data
miss_gts = read.table(miss_gts_tsv)

outlier_snps_sd$missing = miss_gts[which(
  miss_gts$V1 %in% outlier_snps_sd$snp),2]

outlier_snps_qvals$missing = miss_gts[which(
  miss_gts$V1 %in% outlier_snps_qvals$snp),2]


### ### ### ### ### ### ### ### ### ### ### ### ###
### create unified genwin windows data frames ####

print("create unified genwin windows data frames")

# => variables who's baypass result windows should be merged
vars = c("bio9", "bio2", "bio16", "jan_mean_depth")

# data frames for outlier windows and all windows
outlier_windows = data.frame()
all_windows = data.frame()

# iterate through variables 
for (var in vars){
  print(paste("reading windows of", var))
  # outliers for variable
  var_wins = read.table(paste0(results_folder, var,
                               "_genwin_windows_outliers.bed"),
                        col.names = c("chromosome", "start", "end"))
  var_wins$var = var
  outlier_windows = rbind(outlier_windows, var_wins)
  
  # all windows for variable
  var_wins = read.table(paste0(results_folder, var,
                               "_genwin_windows.tsv"),
                        header = T)
  all_windows = rbind(all_windows, var_wins)
}

# sorted beds
sorted_outlier_windows = outlier_windows[order(outlier_windows[,1],
                                               outlier_windows[,2],
                                               outlier_windows[,3]), ]

sorted_all_windows = all_windows[order(all_windows[,1],
                                       all_windows[,2],
                                       all_windows[,3]), ]

### ### ### ### ### ### ### ### 
### generate summary table ####

print("generating summary table")

# data frame for summary of overlap
overlap_summary = data.frame()

# iterate through variables 
for (vari in vars){
  print(paste("summarizing", vari))
  
  # extract variable's windows
  var_wins = outlier_windows %>% filter(var == vari)
  
  # window length and number of overlapping snps
  n_snps_sd = c()
  n_snps_qvals = c()
  
  # iterate through windows
  for (w in 1:nrow(var_wins)){
    # extract window
    window = var_wins[w,]
    
    # overlapping sd snps
    ol_sd = outlier_snps_sd %>% 
      filter(chromosome == window$chromosome & position >=  window$start & position <= window$end) %>%
      nrow()
    n_snps_sd = c(n_snps_sd, ol_sd)
    
    # overlapping qvals snps
    ol_qvals = outlier_snps_qvals %>% 
      filter(chromosome == window$chromosome & position >=  window$start & position <= window$end) %>%
      nrow()
    n_snps_qvals = c(n_snps_qvals, ol_qvals)
  }
  # summary of variable:
  var_summary = data.frame(
    variable = vari,
    mean_w_len = mean(var_wins$end - var_wins$start),
    min_w_len = min(var_wins$end - var_wins$start),
    max_w_len = max(var_wins$end - var_wins$start),
    sd_w_len = sd(var_wins$end - var_wins$start),
    mean_n_snps_sd = mean(n_snps_sd),
    windows_with_snps_sd = length(n_snps_sd[n_snps_sd != 0]),
    max_n_snps_sd = max(n_snps_sd),
    sd_n_snps_sd = sd(n_snps_sd),
    mean_n_snps_qvals = mean(n_snps_qvals),
    windows_with_snps_qvals = length(n_snps_qvals[n_snps_qvals != 0]),
    max_n_snps_qvals = max(n_snps_qvals),
    sd_n_snps_qvals = sd(n_snps_qvals)
  )
  # add to the general summary table
  overlap_summary = rbind(overlap_summary, var_summary)
}
# write table
write.table(overlap_summary,
            "3-Identifying_Candidate_Loci/tables/overlap_summary_table.tsv",
            row.names = F, sep = "\t", col.names = T, quote = F)


### ### ### ### ### ### ### ### ### ### ### ### ### 
### venn diagram - overlap between variables ####

# no idea how :(


### ### ### ### ### ### ### ### ### ###
### merge outlier windows together ####

print("merging outlier windows together")

## bed-like df of windows merged outlier windows
# merge overlapping windows with reduce
gr = reduce(GRanges(
  seqnames = Rle(sorted_outlier_windows$chromosome),
  ranges = IRanges(sorted_outlier_windows$start, end = sorted_outlier_windows$end),
  score = c(sorted_outlier_windows$var)
))

# call chromosome names
a = seqnames(gr)
chrs = c()
for (i in 1:length(a@lengths)){
  chrs = c(chrs, rep(as.character(a@values)[i], a@lengths[i]))
}
# call ranges
b = ranges(gr)

# build bed-like df of merged outlier windows
merged_outlier_windows = data.frame(chromosome = chrs, 
                                    start = b@start,
                                    end = b@start + b@width - 1)

# add a variable column
merged_outlier_windows$variable = NA

# add what variable it was an outlier for
for (vari in vars){
  print(paste("marking windows from", vari))
  var_wins = outlier_windows %>% filter(var == vari)
  gr_vari = reduce(GRanges(
    seqnames = Rle(var_wins$chromosome),
    ranges = IRanges(var_wins$start, end = var_wins$end),
    score = c(var_wins$var)
  ))
  overlaps = findOverlaps(gr, gr_vari)
  merged_outlier_windows[overlaps@from, "variable"] = paste0(
    merged_outlier_windows[overlaps@from, "variable"], ".", vari)
}
# remove NA at start of string
merged_outlier_windows$variable = str_remove(merged_outlier_windows$variable,
                                             "NA.")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### sd outliers - generate overlap window and snps beds #### 

candidate_snps_sd_bed = data.frame()
candidate_win_sd_bed = data.frame()

for (s in 1:nrow(outlier_snps_sd)){
  # snp info
  snp = outlier_snps_sd[s,"snp"]
  snp_chr = unlist(strsplit(snp,":"))[1]
  snp_pos = as.numeric(unlist(strsplit(unlist(strsplit(snp,":"))[2], "_"))[1])
  snp_load = outlier_snps_sd[s, "loading"]
  snp_miss = outlier_snps_sd[s, "missing"]
  # overlap
  ol = merged_outlier_windows %>%
    filter(chromosome == snp_chr & start - 5000 <= snp_pos & end + 5000 >= snp_pos)
  # check number of overlapping windows
  if (nrow(ol) > 0){
    # candidate snp
    snp_ol = paste(snp_chr, min(ol$start), max(ol$end), sep = "_")
    snp_bed = data.frame(chr = snp_chr, position = snp_pos,
                         load = abs(snp_load), overlap = snp_ol,
                         missing = snp_miss)
    candidate_snps_sd_bed = rbind(candidate_snps_sd_bed, snp_bed)
    # candidate window
    win_bed = data.frame(chr = snp_chr, start = min(ol$start),
                         end = max(ol$end), 
                         variable = paste(ol$variable, collapse = "."))
    candidate_win_sd_bed = rbind(candidate_win_sd_bed, win_bed)
  }
}

# remove duplicates windows 
# (have more than one snp so they appear multiple times in table)
overlap_sd_win_nodup = candidate_win_sd_bed[which(!duplicated(candidate_win_sd_bed)),]

# save bed
write.table(overlap_sd_win_nodup,
            "3-Identifying_Candidate_Loci/tables/overlap_sd_win_nodup.bed",
            row.names = F, sep = "\t", col.names = F, quote = F)


# remove "duplicate" snps
# (are in same window and we keep the highest scoring without missing data)
overlap_sd_snps_nodup = candidate_snps_sd_bed[order(candidate_snps_sd_bed[,'overlap'],
                                                    candidate_snps_sd_bed[,'missing'],
                                                   -candidate_snps_sd_bed[,'load']),]
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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### qvals outliers - generate overlap window and snps beds #### 

candidate_snps_qvals_bed = data.frame()
candidate_win_qvals_bed = data.frame()

for (s in 1:nrow(outlier_snps_qvals)){
  # snp info
  snp = outlier_snps_qvals[s,"snp"]
  snp_chr = unlist(strsplit(snp,":"))[1]
  snp_pos = as.numeric(unlist(strsplit(unlist(strsplit(snp,":"))[2], "_"))[1])
  snp_load = outlier_snps_qvals[s, "q.values"]
  snp_miss = outlier_snps_qvals[s, "missing"]
  # overlap
  ol = merged_outlier_windows %>%
    filter(chromosome == snp_chr & start - 5000 <= snp_pos & end + 5000 >= snp_pos)
  # check number of overlapping windows
  if (nrow(ol) > 0){
    # candidate snp
    snp_ol = paste(snp_chr, min(ol$start), max(ol$end), sep = "_")
    snp_bed = data.frame(chr = snp_chr, position = snp_pos,
                         load = abs(snp_load), overlap = snp_ol,
                         missing = snp_miss)
    candidate_snps_qvals_bed = rbind(candidate_snps_qvals_bed, snp_bed)
    # candidate window
    win_bed = data.frame(chr = snp_chr, start = min(ol$start),
                         end = max(ol$end), 
                         variable = paste(ol$variable, collapse = "."))
    candidate_win_qvals_bed = rbind(candidate_win_qvals_bed, win_bed)
  }
}

# remove duplicates windows 
# (have more than one snp so they appear multiple times in table)
overlap_qvals_win_nodup = candidate_win_qvals_bed[which(!duplicated(candidate_win_qvals_bed)),]

# save bed
write.table(overlap_qvals_win_nodup,
            "3-Identifying_Candidate_Loci/tables/overlap_qvals_win_nodup.bed",
            row.names = F, sep = "\t", col.names = F, quote = F)


# remove "duplicate" snps
# (are in same window and we keep the highest scoring without missing data)
overlap_qvals_snps_nodup = candidate_snps_qvals_bed[order(candidate_snps_qvals_bed[,'overlap'],
                                                          candidate_snps_qvals_bed[,'missing'],
                                                         -candidate_snps_qvals_bed[,'load']),]
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


 ### ### ### ### ### ### 
### manhattan plots ####

