# necessary libraries:
library(tidyverse)
library(ade4)
library(adegenet)
library(made4)
library(RColorBrewer)
library(viridis)
library(factoextra)
library(FactoMineR)

## general crap ####

get_gt_df_from_raw <- function(raw_table){
  # read plink raw data into a data.frame of the genotypes
  gl <- read.PLINK(raw_table)
  raw_df <- data.frame(as.matrix(gl))
  return(raw_df)
}

get_loc_from_samples <- function(sample_list){
  # get a vector of populations in the order of samples
  loc <- rep(NA, length(sample_list))
  loc[grep("ba", sample_list)] <- "Balkans"
  loc[grep("ca", sample_list)] <- "Caucasus"
  loc[grep("cr", sample_list)] <- "Carpathians"
  loc[grep("ka", sample_list)] <- "Mongolia"
  loc[grep("ki", sample_list)] <- "Kirov"
  loc[grep("la", sample_list)] <- "Latvia"
  loc[grep("no", sample_list)] <- "Norway"
  loc[grep("po", sample_list)] <- "NE-Poland"
  loc[grep("og", sample_list)] <- "Mongolia"
  loc[grep("to", sample_list)] <- "Mongolia"
  loc[grep("tu", sample_list)] <- "Tuva"
  loc[grep("ur", sample_list)] <- "Urals"
  loc[grep("vl", sample_list)] <- "Vladivostok"
  loc[grep("ya", sample_list)] <- "Yakutia"
  return(loc)
}

get_color_from_samples <- function(sample_list){
  # get a vector of colors in the order of samples
  col <- rep(NA, length(sample_list))
  col[grep("ba", sample_list)] <- "#A035AF"
  col[grep("ca", sample_list)] <- "#B8860b"
  col[grep("cr", sample_list)] <- "#CAB2D6"
  col[grep("ka", sample_list)] <- "#FDBF6F"
  col[grep("ki", sample_list)] <- "#440154FF"
  col[grep("la", sample_list)] <- "#B2DF8A"
  col[grep("no", sample_list)] <- "#3B528BFF"
  col[grep("po", sample_list)] <- "#21908CFF"
  col[grep("og", sample_list)] <- "#FDBF6F"
  col[grep("to", sample_list)] <- "#FDBF6F"
  col[grep("tu", sample_list)] <- "#FF7F00"
  col[grep("ur", sample_list)] <- "#0F4909"
  col[grep("vl", sample_list)] <- "#FB9A99"
  col[grep("ya", sample_list)] <- "#E31A1C"
  return(col)
}

## pca and coinertia ####

build_df_from_pca <- function(pca_obj){
  # build a data.frame of samples (PCs, populations, colors)
  # from a PCA object
  pca_df <- as.data.frame(pca_obj$ind$coord)
  pca_df$population <- get_loc_from_samples(row.names(pca_df))
  pca_df$color <- get_color_from_samples(row.names(pca_df))
  return(pca_df)
}

plot_pca <- function(pca_df, variance_percents, x = 1, y = 2){
  # produce a pca biplot from pca data.frame and variance percentages
  pca_plot <- ggplot() +
    geom_point(data = pca_df, aes(x = pca_df[,x], y = pca_df[,y]),
               fill = pca_df[,"color"], shape = 21) +
    theme_bw(base_size = 14, base_family = "Times") +
    theme(panel.background = element_blank(),
          # panel.grid = element_blank(),
          plot.background = element_blank(),
          #axis.text.x = element_blank(),
          #axis.ticks.x=element_blank()
          ) +
    # xlim(min = min_lim, max = max_lim) +
    # ylim(min = min_lim, max = max_lim) +
    xlab(paste0("PC", x, " - ", round(variance_percents[x], 2), "%")) +
    ylab(paste0("PC", y, " - ", round(variance_percents[y], 2), "%"))
  
  return(pca_plot)
}

build_coi_df_from_gt_dfs <- function(gt_df1, gt_df2){
  # gt_df1 = sd_candidate_snps
  # gt_df2 = neutral_snps
  cia1 <- cia(t(gt_df2), t(gt_df1),
      cia.nf=3, cia.scan=FALSE, nsc=TRUE)
  print(paste("co-inertia RV is", cia1$coinertia$RV))
  cia_df <- data.frame(sample=(rownames(cia1$coinertia$mX)),
                       start_x=(cia1$coinertia$mX$NorS1),
                       start_y=(cia1$coinertia$mX$NorS2),
                       start_z=(cia1$coinertia$mX$NorS3),
                       end_x=(cia1$coinertia$mY$NorS1),
                       end_y=(cia1$coinertia$mY$NorS2),
                       end_z=(cia1$coinertia$mY$NorS3))
  cia_df$population <- get_loc_from_samples(cia_df$sample)
  cia_df$color <- get_color_from_samples(cia_df$sample)
  
}


ggplot() + 
  geom_point(data=cia_df, mapping=aes(x=start_x, y=start_y), size=0.5, fill=cia_df$color) +
  geom_segment(data=cia_df, mapping=aes(x=start_x, y=start_y, xend=end_x, yend=end_y), 
               arrow = arrow(length = unit(0.25,"cm")), size=0.4, color="darkgrey", alpha=0.8) + 
  geom_point(data=cia_df, mapping=aes(x=end_x, y=end_y), size=3, shape=21, fill=cia_df$color, alpha=0.6) +
  # scale_y_continuous(limits=c(-1.25,0.65))+
  theme_minimal()
