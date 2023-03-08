library(readr)
library(dplyr)
library(ggplot2)

# Define the input Directory
wd_input <- "/Users/enricobazzicalupo/Documents/Selection_Eurasian_Lynx/SamTools_Depth/"

# Create a list of the sample files' names (the SAMTOOLS depth output files)
sample_files <- list.files(wd_input, pattern="*.depth$")

# Create an Empty dataframe to save values for each dataset for the final table
depth_per_sample <- data.frame()

# For every Sample file:
for (i in 1:length(sample_files)){
  
  # Import the sample file table
  input.depth <- read_delim(paste0(wd_input,sample_files[[i]]), col_names = F, delim = '\t')
  
  # Add a column (Total) which is the sum of all the depth columns (from the third to the last)
  input.depth$Total <- rowSums(input.depth[,3:ncol(input.depth)])
  
  # Create a frequency table of the values of the Total column
  freq_table_DF <- as.data.frame(table(input.depth$Total))
  
  # Define Dataset Name:
  population=unlist(strsplit(basename(sample_files[[i]]),"[.]"))[1]
  
  # Define the functions for mean and standard deviation,
  # and define maximum and minimum depth based on population:
  # minimum of 5x per individual is needed for "medium depth" populations (ca,ba)
  # and of 3x per individual for "low depth" populations (the rest)
  
  if (population=="ba" || population=="ca") {
    
    mymean <- mean(input.depth$Total)
    mysd <- sd(input.depth$Total)
    maxDepth_sd = mymean + (mysd * 1.5)
    minDepth <- 5 * (ncol(input.depth) - 3)
    
  } else {
    
    mymean <- mean(input.depth$Total)
    mysd <- sd(input.depth$Total)
    maxDepth_sd = mymean + (mysd * 1.5)
    minDepth <- 3 * (ncol(input.depth) - 3)
    
  }
  
  
  basename(sample_files[[i]])
  # Add dataset information to Dataframe
  depth_per_sample <- rbind(depth_per_sample,
                            data.frame(pop = population, nindividuals = (ncol(input.depth) - 3),
                                       mean = mymean, sd = mysd,
                                       maxDepthSD = maxDepth_sd, minDepth3x = minDepth))
  
  # Draw and save a graph of the distribution of depth values, with upper and lower depth limits
  ggplot(freq_table_DF, aes(x = as.numeric(Var1), y = Freq)) +
    geom_bar(stat = "identity", color = "black") +
    scale_x_continuous(breaks = 0:250*10, limits = c(0, maxDepth_sd*1.5)) +
    # scale_x_discrete(limits = c(0, maxDepth_DF*1.5)) +
    scale_y_continuous(expand=c(0,0)) +
    geom_vline(xintercept=maxDepth_sd,linetype="dashed", size=0.5) +
    geom_vline(xintercept=minDepth,linetype="dashed", size=0.5) +
    #geom_vline(xintercept=minDepth_DF_meanfolds, colour ="grey", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=maxDepth_DF_meanfolds, colour ="grey", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=minDepth_5x, colour ="red", linetype="dashed", size=0.5) +
    #geom_vline(xintercept=maxDepth_MEAN, colour ="red", linetype="dashed", size=0.5) +
    theme_classic() +
    theme(text = element_text(size=10))
  ggsave (filename = (paste0("graph_",basename(sample_files[[i]]),".pdf")), path = wd_input)
}

# Print the table to a file
write.table(x = depth_per_sample,file = paste0(wd_input,"depth_per_sample.csv"),
            quote=FALSE, col.names = T, row.names = FALSE, sep= ",")