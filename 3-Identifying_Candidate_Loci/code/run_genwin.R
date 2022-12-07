source("3-Identifying_Candidate_Loci/code/genwin_functions.R", print.eval=TRUE)

## ## ## ## ## ## 
## arguments:  ##
## ## ## ## ## ##

# folder where results are stored
results_folder = "3-Identifying_Candidate_Loci/tables/"
# variables whose results we want to analyze with genwin
args = commandArgs(trailingOnly=TRUE)
var = args[1]

## ## ## ## ## ## 
## load data:  ##
## ## ## ## ## ##

# create a data frame with the results of the chosen variable
results_df = make_baypass_results_df_of_var(results_folder = results_folder, var = var)

# list of chromosomes
chr_list = unique(results_df$scaffold)

## ## ## ## ## ## ##
##   run GenWin   ##
## ## ## ## ## ## ##

# empty data frame to store results from each chromosome
all_spline = data.frame()

# loop through chr_list adding results to data frame
for (chr in chr_list){
  print(paste0("analyzing chromosome : ", chr))
  # analyze chromosome with splineAnalyze
  spline.data = run_splineanalyze(baypass_results_df = results_df, chr = chr)
  # store results in data frame with all spline results
  all_spline = rbind(all_spline, spline.data)
}

# add outliers
all_spline_outliers = add_wstat_outliers_column(all_spline)

# save table
write.table(x = all_spline_outliers,
            file = paste0(results_folder, var, "_genwin_windows.tsv"),
            quote = F,  col.names = T, row.names = F, sep= "\t")

# create a bed file of outliers
outliers_bed = all_spline_outliers %>% 
  filter(outlier == "yes") %>%
  select(scaffold, WindowStart, WindowStop)

# save outliers bed
write.table(x = outliers_bed,
            file = paste0(results_folder, var, "_genwin_windows_outliers.bed"),
            quote = F,  col.names = F, row.names = F, sep= "\t")
