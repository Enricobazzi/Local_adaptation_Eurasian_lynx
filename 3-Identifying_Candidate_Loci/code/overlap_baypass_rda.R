### ### ### ### ### ###
### load libraries ####

library(tidyverse)
library(RColorBrewer)
library(GenomicRanges)
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html

# load a dictionary for the variables
var_dict <- c("bio1" = "T_mean_year",
              "bio2" = "T_range_day",
              "bio3" = "Iso_T",
              "bio4" = "T_seasonalit",
              "bio5" = "T_max_warm",
              "bio6" = "T_min_cold",
              "bio7" = "T_range_year",
              "bio8" = "T_wet_quart",
              "bio9" = "T_dry_quart",
              "bio10" = "T_warm_quart",
              "bio11" = "T_cold_quart",
              "bio12" = "P_annual",
              "bio13" = "P_wet_month",
              "bio14" = "P_dry_month",
              "bio15" = "P_seasonality",
              "bio16" = "P_wet_quart",
              "bio17" = "P_dry_quart",
              "bio18" = "P_warm_quart",
              "bio19" = "P_cold_quart",
              "jan_mean_depth" = "Jan_mean_depth",
              "mean_snow_days" = "Mean_snow_days",
              "Geographic" = "Geographic",
              "xCoord" = "x-coord",
              "yCoord" = "y-coord")

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

# qvals candidates table name
qvals_candidates_table = paste0("vars_", v_form[1],
                             ".cond_", v_form[2],
                             ".cand_qvals.K", K,
                             ".candidates.tsv")

### ### ### ### ### ### ### ###
###  prepare the datasets  ####

print("preparing the datasets")

## folder where results are stored
results_folder = "3-Identifying_Candidate_Loci/tables/"

## missing genotypes table
miss_gts_tsv = "1-Preparing_Genetic_Data/tables/finalset_missing_gts.tsv"

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

## outlier SNPs with q-values (strict)
outlier_snps_qvals = read.table(paste0(results_folder,
                                       qvals_candidates_table),
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
outlier_snps_qvals$position = as.numeric(poss)

## read missing data
miss_gts = read.table(miss_gts_tsv, header = F, sep = "\t",
                      comment.char = "",
                      col.names = c("snp", "missing"))

# add missing data
outlier_snps_sd = merge(outlier_snps_sd, miss_gts, by = "snp")
outlier_snps_qvals = merge(outlier_snps_qvals, miss_gts, by = "snp")


### ### ### ### ### ### ### ### ### ### ### ### ###
### create unified genwin windows data frames ####

print("create unified genwin windows data frames")

# => variables who's BayPass result windows should be merged
vars = strsplit(v_form[1], split = "+", fixed = T)[[1]]

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
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_summary_table.tsv"),
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
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_sd_win_nodup.bed"),
            row.names = F, sep = "\t", col.names = F, quote = F)

# remove "duplicate" snps
# (are in same window and we keep the highest scoring without missing data)
overlap_sd_snps_nodup = candidate_snps_sd_bed[order(candidate_snps_sd_bed[,'overlap'],
                                                    candidate_snps_sd_bed[,'missing'],
                                                   -candidate_snps_sd_bed[,'load']),]
overlap_sd_snps_nodup = overlap_sd_snps_nodup[!duplicated(overlap_sd_snps_nodup$overlap),]

# save table
write.table(overlap_sd_snps_nodup,
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_sd_snps_nodup.tsv"),
            row.names = F, sep = "\t")

# save bed
overlap_sd_snps_nodup_bed = data.frame(chromosome = overlap_sd_snps_nodup$chr,
                                       start = overlap_sd_snps_nodup$position - 1,
                                       end = overlap_sd_snps_nodup$position)

write.table(overlap_sd_snps_nodup_bed,
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_sd_snps_nodup.bed"),
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
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_qvals_win_nodup.bed"),
            row.names = F, sep = "\t", col.names = F, quote = F)


# remove "duplicate" snps
# (are in same window and we keep the highest scoring without missing data)
overlap_qvals_snps_nodup = candidate_snps_qvals_bed[order(candidate_snps_qvals_bed[,'overlap'],
                                                          candidate_snps_qvals_bed[,'missing'],
                                                         -candidate_snps_qvals_bed[,'load']),]
overlap_qvals_snps_nodup = overlap_qvals_snps_nodup[!duplicated(overlap_qvals_snps_nodup$overlap),]

# save table
write.table(overlap_qvals_snps_nodup,
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_qvals_snps_nodup.tsv"),
            row.names = F, sep = "\t")

# save bed
overlap_qvals_snps_nodup_bed = data.frame(chromosome = overlap_qvals_snps_nodup$chr,
                                          start = overlap_qvals_snps_nodup$position - 1,
                                          end = overlap_qvals_snps_nodup$position)

write.table(overlap_qvals_snps_nodup_bed,
            paste0("3-Identifying_Candidate_Loci/tables/",
                   "vars_", v_form[1],
                   ".cond_", v_form[2],
                   ".overlap_qvals_snps_nodup.bed"),
            row.names = F, sep = "\t", col.names = F, quote = F)


 ### ### ### ### ### ### 
### manhattan plots ####

for (var in vars) {
  print(var)
  
  # subset variable and prepare vectors for storing overlap
  var_wins <- sorted_all_windows %>% filter(variable == var)
  
  var_wins_overlaps <- data.frame()
  
  # loop through windows registering overlaps between
  # outlier windows of baypass and outlier snps of rda
  for (w in 1:nrow(var_wins)){
    message('\r', round((w / nrow(var_wins) * 100), 2), appendLF = FALSE)
    
    # extract window
    window = var_wins[w,]
    
    # overlapping sd snps
    ol_sd = outlier_snps_sd[outlier_snps_sd$chromosome == window$scaffold,]
    ol_sd = ol_sd[ol_sd$position >= window$WindowStart,]
    ol_sd = ol_sd[ol_sd$position <= window$WindowStop,]
    ol_sd = nrow(ol_sd)
    
    # add yes or no based on if there is overlap
    if (ol_sd > 0){
      cand_sd = "yes"
    } else {
      cand_sd = "no"
    }
    
    # add number of overlapping snps
    # n_snps_sd[w] = ol_sd
    
    # overlapping qvals snps
    ol_qvals = outlier_snps_qvals[outlier_snps_qvals$chromosome == window$scaffold,]
    ol_qvals = ol_qvals[ol_qvals$position >= window$WindowStart,]
    ol_qvals = ol_qvals[ol_qvals$position <= window$WindowStop,]
    ol_qvals = nrow(ol_qvals)
    
    # add yes or no based on if there is overlap
    if (ol_qvals > 0){
      cand_qvals = "yes"
    } else {
      cand_qvals = "no"
    }
    
    win_ol <- data.frame(sd_cand = cand_sd, n_snps_sd = ol_sd,
                         qvals_cand = cand_qvals, n_snps_qvals = ol_qvals)
    
    var_wins_overlaps <- rbind(var_wins_overlaps, win_ol)
    
  }
  
  var_wins <- cbind(var_wins, var_wins_overlaps)
  
  # add overlaps windownumbers and colors to dataframe for plotting
  var_wins$WindowNumber <- 1:nrow(var_wins)
  var_wins$color <- ifelse(var_wins$n_snps_sd == 0, "grey32",
                           ifelse(var_wins$n_snps_sd == 1, brewer.pal(11,"RdYlBu")[10],
                                  ifelse(var_wins$n_snps_sd == 2, brewer.pal(11,"RdYlBu")[8],
                                         ifelse(var_wins$n_snps_sd == 3, brewer.pal(11,"RdYlBu")[6],
                                                ifelse(var_wins$n_snps_sd == 4, brewer.pal(11,"RdYlBu")[4],
                                                       ifelse(var_wins$n_snps_sd == 5, brewer.pal(11,"RdYlBu")[2],
                                                              ifelse(var_wins$n_snps_sd > 5, brewer.pal(8,"Dark2")[4], 
                                                                     brewer.pal(8,"Dark2")[8])))))))
  display.brewer.pal(8,"Dark2")
  # windows that are not outliers of any method
  crapwins <- var_wins %>%
    filter(outlier == "no" & sd_cand == "no" & qvals_cand == "no")
  
  # windows that are outliers of baypass but not of rda
  bpout_nosnp <- var_wins %>%
    filter(outlier == "yes" & sd_cand == "no" & qvals_cand == "no")
  
  # windows that are outliers of baypass and of rda SD but not rda qvals
  bpout_sdout_noqval <- var_wins %>%
    filter(outlier == "yes" & sd_cand == "yes" & qvals_cand == "no")
  
  # windows that are outliers of baypass and of rda SD and rda qvals
  bpout_sdout_qvalout <- var_wins %>%
    filter(outlier == "yes" & sd_cand == "yes" & qvals_cand == "yes")
  
  # plot
  p <- ggplot() +
    geom_point(data=crapwins, aes(x=WindowNumber, y=Wstat),
               color="grey32", shape=4, size=1.5) +
    geom_point(data=bpout_nosnp, aes(x=WindowNumber, y=Wstat),
               color="gold2",shape=4, size=1.5) +
    geom_point(data=bpout_sdout_noqval, aes(x=WindowNumber, y=Wstat),
               color="grey32", fill="darkorange", shape=21, size=3) +
    geom_point(data=bpout_sdout_qvalout, aes(x=WindowNumber, y=Wstat),
               color="grey32", fill="darkorange", shape=21, size=3) +
    xlab("Window Number") +
    ggtitle(as.vector(var_dict[var])) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  p
  # facet_wrap(~ scaffold, scales = "free_x")
  
  ggsave(filename = paste0("3-Identifying_Candidate_Loci/plots/", var,
                           "_windows_overlap_",
                           "vars_", v_form[1],
                           ".cond_", v_form[2], 
                           "_manhattan.pdf"),
         plot = p, width = 6, height = 3, units = "in")
}
