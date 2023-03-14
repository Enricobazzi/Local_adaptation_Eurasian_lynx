############################################

library(tidyverse)
library(GenWin)

############################################

# create a data.frame with BayPass results for selected variable
make_baypass_results_df_of_var <- function(results_folder, var){
  
  # create an empty data.frame
  snps.table <- data.frame()
  
  # iterate through the 50 result tables for the variable
  for (n in 1:50){
    # upload the dataset table:
    var.snp <- read.table(paste0(results_folder, "AUX_", var, "_", n, "_summary_betai.out"),
                         h = T)
    # calculate the SNP database number sequence (1 every 50) 
    # NOTE: 2100553 is the total number of SNPs 
    # and 50 is the number of data sets I divided them into
    var.snp <- data.frame(var.snp, SNPnum = seq(n, 2100553, by = 50))
    # add the data set rows to the data frame
    snps.table <- rbind(snps.table, var.snp)
  }
  
  # order the table based on the SNP number
  snps.table <- snps.table %>% arrange(SNPnum)
  
  # Add SNP ID information from SNPIDs table
  SNPIDs <- read_tsv(paste0(results_folder, "finalset.maf5pc.SNPIDs"),
                     col_names = F, show_col_types = F)[,2-3] %>%
    rename("X2" = "scaffold", "X3" = "position")
  
  # make data frame with scaffold, position, snp number, and bayes factor score
  baypass_results_df <- data.frame(scaffold = SNPIDs$scaffold,
                                   position = SNPIDs$position,
                                   SNPnum = snps.table$SNPnum,
                                   BF.dB. = snps.table$BF.dB.,
                                   variable = var)
  
  return(baypass_results_df)
}

# analyze data of one chromosome with genwin splineAnalyze
run_splineanalyze <- function(baypass_results_df, chr, win_size = 10000, m = 3){
  # get data for chromosome
  data_chr <- baypass_results_df %>% filter(scaffold==chr)
  # get spline
  spline <- splineAnalyze(Y = data_chr$BF.dB., 
                          map = data_chr$position,
                          smoothness = win_size, 
                          mean = mean(baypass_results_df$BF.dB.),
                          s2 = var(baypass_results_df$BF.dB.),
                          method = m)
  # create a data frame from object created by splineAnalyze and additional info
  spline.data <- data.frame(scaffold = chr,
                            WindowStart = spline$windowData$WindowStart,
                            WindowStop = spline$windowData$WindowStop,
                            WindowLength = spline$windowData$WindowStop - spline$windowData$WindowStart,
                            variable = unique(baypass_results_df$variable),
                            SNPcount = spline$windowData$SNPcount,
                            MeanY = spline$windowData$MeanY,
                            Wstat = spline$windowData$Wstat)
  return(spline.data)
}

# add a column indicating if a window is an outlier or not
# outlier is defined by function quantile = produces sample quantiles corresponding to the given probabilities (smallest observation a prob of 0 and largest a prob of 1).
add_wstat_outliers_column <- function(all_spline, quant = 0.99){
  threshold <- as.numeric(quantile(all_spline$Wstat, probs=quant, na.rm=T))
  outlier_column <- ifelse(all_spline$Wstat > threshold, "yes", "no")
  return(data.frame(all_spline, outlier = outlier_column))
}
