# in this script I will prepare the table containing the values of the
# environmental variables for each population in a format useful for BayPass

library(tidyverse)

# read the persample table
persample_table <- read.delim(file = "../tables/allvars_persample_table.tsv",
                              header = T, sep = "\t") %>% column_to_rownames("sample")
# define populations in alphabetical order for BayPass
populations <- c("ca", "ki", "la", "mo", "tu", "ur", "vl", "ya")

# create empty df to fill with population averages
baypass_perpop_table <- data.frame(vars=colnames(persample_table))

# loop through populations
for (pop in populations) {
  
  # get a table with only the population's samples
  # mongolia has multiple samples with different pop names
  if (pop == "mo") {
    pop_table <- persample_table %>% filter(grepl("ka|og|to", rownames(persample_table)) == TRUE)
  } else {
    pop_table <- persample_table %>% filter(grepl(pop, rownames(persample_table)) == TRUE)
  }
  
  # get the population's average of each variable 
  avgs <- c()
  for (n in 1:(length(colnames(pop_table)))){
    avgs <- c(avgs, mean(pop_table[,n]))
  }
  baypass_perpop_table <- cbind(baypass_perpop_table, data.frame(x=avgs))
}

# adjust format
baypass_perpop_table <- baypass_perpop_table %>% column_to_rownames("vars") %>% setNames(populations)

# save table
write.table(x = baypass_perpop_table, file = "../tables/allvars_baypass_perpop_table.tsv", quote=FALSE,
            col.names = T, row.names = T, sep= "\t")
