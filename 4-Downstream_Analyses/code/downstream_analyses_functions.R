# necessary libraries:
library(tidyverse)
library(ade4)
library(adegenet)
library(made4)
library(RColorBrewer)
library(viridis)
library(factoextra)
library(FactoMineR)
library(raster)

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

## ## ## ## ## ## ## 
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
  loc[grep("vl", sample_list)] <- "Primorsky Krai"
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

get_number_from_samples <- function(sample_list){
  num <- rep(NA, length(sample_list))
  num[grep("ba", sample_list)] <- 1
  num[grep("ca", sample_list)] <- 3
  num[grep("cr", sample_list)] <- 2
  num[grep("ka", sample_list)] <- 6
  num[grep("ki", sample_list)] <- 4
  num[grep("la", sample_list)] <- 5
  num[grep("no", sample_list)] <- 8
  num[grep("po", sample_list)] <- 7
  num[grep("og", sample_list)] <- 6
  num[grep("to", sample_list)] <- 6
  num[grep("tu", sample_list)] <- 9
  num[grep("ur", sample_list)] <- 10
  num[grep("vl", sample_list)] <- 11
  num[grep("ya", sample_list)] <- 12
  return(num)
}

calculate_variance_percents_from_eigenvals <- function(eigenvals){
  variance_percents = round(x = (eigenvals / sum(eigenvals) * 100), 2)
}

translate_var_from_dict <- function(variables){
  return(c(as.character(var_dict[variables])))
}

## ## ## ## ## ## ## ## #
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
               fill = pca_df[,"color"], shape = 21, size = 2.5) +
    theme_bw() +
    theme(panel.background = element_blank(),
          # panel.grid = element_blank(),
          plot.background = element_blank(),
          #axis.text.x = element_blank(),
          #axis.ticks.x=element_blank()
          ) +
    # xlim(min = min_lim, max = max_lim) +
    # ylim(min = min_lim, max = max_lim) +
    
    xlab(paste0("PC", x, " - ", round(variance_percents[x], 2), "%")) +
    ylab(paste0("PC", y, " - ", round(variance_percents[y], 2), "%")) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9))
  
  return(pca_plot)
}

# run coinertia from gt data.frames 
# - default pca with 4 axes per df and coin
run_coin_from_gt_dfs <- function(gt_df1, gt_df2,
                                 method = "pca",
                                 nf_1 = 4, nf_2 = 4,
                                 coin_nf = 4){
  if (method == "pca"){
    dudi_1 = dudi.pca(gt_df1, scannf = FALSE, nf = nf_1)
    dudi_2 = dudi.pca(gt_df2, scannf = FALSE, nf = nf_2)
  }
  coin = coinertia(dudi_1, dudi_2, scannf = FALSE, nf = coin_nf)
  return(coin)
}

build_coin_df_from_coin <- function(coin){
  coin_df = data.frame(
  neutral_1 = coin$mX$NorS1,
  neutral_2 = coin$mX$NorS2,
  neutral_3 = coin$mX$NorS3,
  neutral_4 = coin$mX$NorS4,
  candidate_1 = coin$mY$NorS1,
  candidate_2 = coin$mY$NorS2,
  candidate_3 = coin$mY$NorS3,
  candidate_4 = coin$mY$NorS4,
  color = get_color_from_samples(rownames(neutral_pca_df)),
  population = get_loc_from_samples(rownames(neutral_pca_df))
)
}

plot_coinertia <- function(coin_df, x = 1, y = 2,
                           variance_percents){
  df1_ax1 = x
  df1_ax2 = y
  df2_ax1 = x + 4
  df2_ax2 = y + 4
  
  coin_plot <- ggplot() + 
    geom_point(data=coin_df, mapping=aes(x=coin_df[,df1_ax1], y=coin_df[,df1_ax2]), size=0.5, fill=coin_df$color) +
    geom_segment(data=coin_df, mapping=aes(x=coin_df[,df1_ax1], y=coin_df[,df1_ax2], xend=coin_df[,df2_ax1], yend=coin_df[,df2_ax2]), 
                 arrow = arrow(length = unit(0.25,"cm")), size=0.4, color="darkgrey", alpha=0.8) + 
    geom_point(data=coin_df, mapping=aes(x=coin_df[,df2_ax1], y=coin_df[,df2_ax2]), size=3, shape=21, fill=coin_df$color, alpha=0.6) +
    theme_bw(base_size = 14, base_family = "Times") +
    theme(panel.background = element_blank(),
          # panel.grid = element_blank(),
          plot.background = element_blank(),
          #axis.text.x = element_blank(),
          #axis.ticks.x=element_blank()
    ) +
    # xlim(min = min_lim, max = max_lim) +
    # ylim(min = min_lim, max = max_lim) +
    xlab(paste0("co-inertia axis ", x, " - ", variance_percents[x], "%")) +
    ylab(paste0("co-inertia axis ", y, " - ", round(variance_percents[y], 2), "%"))
  
  return(coin_plot)
}

get_distance_df_from_snps <- function(gt_df){
  # distance matrix
  D <- as.matrix(dist(gt_df))
  # empty distance_dfs
  distance_df <- data.frame()
  pop_distance_df <- data.frame()
  # for every population
  for (pop in c("ba", "ca", "cr", "ki", "la", "no", "po", "tu", "ur", "vl", "ya")){
    # get idate part of table
    pairwise_dist = c(D[rownames(D)[grep(pop, rownames(D))],])
    names_labels = labels(D[rownames(D)[grep(pop, rownames(D))],])
    combs = expand.grid(names_labels[[1]], names_labels[[2]])
    color = get_color_from_samples(as.character(combs$Var2))
    location = get_loc_from_samples(as.character(combs$Var2))
    pairwise_dist <- data.frame(p1 = get_loc_from_samples(pop),
                                p2 = location,
                                pairwise_dist = pairwise_dist,
                                color = color)
    distance_df <- rbind(distance_df, pairwise_dist)
    # add mean population distances
    # (same value for each sample in population pairs)
    for (loc in unique(pairwise_dist$p2)){
      pop_dist <- pairwise_dist %>% filter(p2 == loc & pairwise_dist > 0)
      pop_mean_dist <- mean(pop_dist$pairwise_dist)
      pop_dist$pop_mean_dist <- pop_mean_dist
      pop_dist <- pop_dist[which((!duplicated(pop_dist$p2)) == "TRUE"), ]
      pop_distance_df <- rbind(pop_dist, pop_distance_df)
    }
  }
  return(list(distance_df, pop_distance_df))
}
  

 ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## generalized dissimilarity modeling ####

get_snp_dist <- function(gt_df){
  # generate Euclidean distance matrix
  snp.dist <- dist(gt_df, diag = T, upper = T)
  #rescale by dividing by max value
  snp.dist.1 <- snp.dist/(max(snp.dist))
  snp.dist.1 <- as.matrix(snp.dist.1)
  # make the row names an actual column named "ID"
  snp.dist.1 <- cbind(ID = rownames(snp.dist.1), snp.dist.1)
  #remove prior row names
  rownames(snp.dist.1) <- NULL  
  # return distance object
  return(snp.dist.1)
}

get_clim_points <- function(sample_list){
  # get environmental data
  clim.data <- read_delim("2-Preparing_Environmental_Data/tables/allvars_persample_table.tsv",
                          col_names = T, delim = "\t") %>%
    filter(sample %in% sample_list == T) %>%
    column_to_rownames(var = "sample")
  
  # get coordinate data
  coord_table <- read_delim("2-Preparing_Environmental_Data/tables/wholeset_coordinates.csv",
                            col_names = T, delim = ',') %>%
    filter(id %in% sample_list == T)
  
  # combine environment and coordinates
  clim.points <- data.frame(clim.data,
                            x=as.numeric(coord_table$longitude),
                            y=as.numeric(coord_table$latitude)) %>% 
    rownames_to_column(var="ID")
  
  # return combined dataset
  return(clim.points)
}

# prepare gdm input
prepare_gdm_input <- function(snp.dist.1, clim.points){
  gdm.input <- gdm::formatsitepair(bioData = snp.dist.1, bioFormat = 3, 
                                           predData = clim.points, siteColumn = "ID", 
                                           XColumn = "x", YColumn = "y")
  return(gdm.input)
}

# plot importance
plot_importance <- function(gdm_obj){
  # plot a sorted horizontal bar plot of variable importance
  ssplines <- isplineExtract(gdm_obj)
  importance <- data.frame(ssplines$y) %>% summarise_if(is.numeric, max)
  importance_df <- data.frame(importance=t(importance)) %>%
    rownames_to_column("variable")
  importance_df$variable <- translate_var_from_dict(importance_df$variable)
  importance_df <- importance_df[order(-importance_df$importance),]
  importance_df$importance <- round(
    x = importance_df$importance / sum(importance_df$importance) * 100,
    digits = 2)
  p <- ggplot(data = filter(importance_df, importance > 0.01),
              aes(x=variable, y=importance)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=c(rev(importance_df[,1]))) +
    #scale_y_continuous(limits = c(0,1))+
    theme_minimal(base_size = 24) +
    theme(axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 0),
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9)) +
    ylab("Percentage of importance") +
    coord_flip()
  return(p)
}

# get raw climate layer raster
get_clim_layer <- function(var_list){
  # might add a var_list (e.g. c("bio2", "bio9", "bio16")) as argument
  # if I need it to be more reproducible ( -> stack(worldclim[[var_list]] )
  worldclim <- getData("worldclim", path = "2-Preparing_Environmental_Data/tables/",
                       var = "bio", res = 10)
  if ("jan_mean_depth" %in% var_list & "mean_snow_days" %in% var_list){
    jan_mean_depth <- raster("2-Preparing_Environmental_Data/tables/jan_mean_depth.tif")
    mean_snow_days <- raster("2-Preparing_Environmental_Data/tables/mean_snow_days.tif")
    clim.layer <- stack(worldclim[[var_list]],
                        projectRaster(jan_mean_depth, worldclim),
                        projectRaster(mean_snow_days, worldclim))
  } else if ("jan_mean_depth" %in% var_list){
    jan_mean_depth <- raster("2-Preparing_Environmental_Data/tables/jan_mean_depth.tif")
    clim.layer <- stack(worldclim[[var_list]],
                        projectRaster(jan_mean_depth, worldclim))
  } else if ("mean_snow_days" %in% var_list){
    mean_snow_days <- raster("2-Preparing_Environmental_Data/tables/mean_snow_days.tif")
    clim.layer <- stack(worldclim[[var_list]],
                        projectRaster(mean_snow_days, worldclim))
  }
  return(clim.layer)
}

# raster pca
run_pca_raster <- function(raster_obj){
  # extract values from raster and run pca on them
  rast_values <- na.omit(getValues(raster_obj))
  rast_pca <- prcomp(rast_values, center = TRUE, scale. = FALSE)
  return(rast_pca)
}

# write pca table
write_pca_table <- function(raster_pca_obj, file_name){
  # take a pca object obtained from a raster and save the table of its rotations
  pca_table <- data.frame(raster_pca_obj$rotation) %>%
    mutate_if(is.numeric, round, digits=3) %>% 
    rownames_to_column(var="predictor")
  write.table(x = pca_table,
              file = file_name,
              quote=FALSE,  col.names = T,
              row.names = FALSE, sep= "\t")
}

# build rbg pca raster
create_rgb_pca_raster <- function(trans_clim_layers, raster_pca,
                                  extent = c(-10, 180, 20, 80)){
  # create the pca raster
  pca.rast <- predict(trans_clim_layers, raster_pca,
                      index=1:ncol(raster_pca$rotation))
  
  # transform values to rgb
  for (n in 1:ncol(raster_pca$rotation)){
    pca.rast[[n]] <- (pca.rast[[n]]-pca.rast[[n]]@data@min) / (pca.rast[[n]]@data@max-pca.rast[[n]]@data@min)*255
  }

  # crop to desired extent
  pca.rast.crop <- crop(pca.rast, extent)
  
  return(pca.rast.crop)
}

# make rgb map
plot_rgb_map <- function(sample_list, pca.rast.crop,
                         r = 1, g = 2, b = 3){
  extent <- c(-10, 180, 20, 80)
  distr.map <- rgdal::readOGR("2-Preparing_Environmental_Data/tables/redlist_species_data/data_0.shp")
  loc <- get_loc_from_samples(sample_list)
  num <- get_number_from_samples(sample_list)
  ola <- data.frame(sample = sample_list, pop = loc, n = num)
  eco <- levels(ola$pop)
  bg <- get_color_from_samples(sample_list)
  #bg <- cols <- c("#A035AF",
  #                "#CAB2D6",
  #                "#B8860b",
  #                "#440154FF",
  #                "#B2DF8A",
  #                "#FDBF6F",
  #                "#21908CFF",
  #                "#3B528BFF",
  #                "#FF7F00",
  #                "#0F4909",
  #                "#FB9A99",
  #                "#E31A1C")
  r2 <- crop(pca.rast.crop, extent(distr.map))
  r3 <- mask(r2, distr.map)
  
  # Empty plot to ???
  empty <- raster(crs=pca.rast.crop)
  values(empty) <- 0
  empty.crop <- crop(empty, extent)
  # line around land
  ras <- pca.rast.crop[[1]] > -Inf
  pp <- rasterToPolygons(ras, dissolve=TRUE)
  
  
  # plot 
  # plotRGB(pca.rast.crop, r = r, g = g, b = b, bgalpha=0)
  plot(empty.crop, col=c('lightcyan2'), legend=F, asp=NA)
  plotRGB(pca.rast.crop, r = r, g = g, b = b, bgalpha=0, add=T)
  # plotRGB(r3, r = r, g = g, b = b, bgalpha=0, add=T)
  plot(pp, lwd=1.2, border='black', add=TRUE)
  plot(distr.map, lwd=1.2, border='darkgreen', add=T)
  points(clim.points$x,clim.points$y,
         pch=1, lwd=1, cex=1, 
         bg=bg[match(clim.points$ID, sample_list)]) # the lynxes
#  points(clim.points$x,clim.points$y,
#         pch=4, lwd=0.7, cex=0.8, col = rgb(red = 0, green = 0,
#                                            blue = 0, alpha = 0.5),
#         bg=bg[match(clim.points$ID, sample_list)]) # the lynxes
}

# difference function - included in plot_procrust_resid_map
# maybe should be kept separate?
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[c(mapCells)] <- residMap
  return(list(max(residMap), rast))
}

# plot difference function
plot_procrust_resid_map <- function(neutral_pca, candidate_pca, rast,
                                    samples, clim.points){
  # calculate procrustes residuals between two pcas
  diffProcrust <- procrustes(neutral_pca, candidate_pca,
                             scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  
  # extract cells in raster to be filled with residual values
  mapCells_df <- data.frame(getValues(rast)) %>%
    rownames_to_column("cells") %>% drop_na()
  mapCells <- as.numeric(mapCells_df$cells)
  
  # substitute raster values with residual values
  rast[c(mapCells)] <- residMap
  
  # prepare plot - extent
  extent <- c(-10, 180, 20, 80)
  # prepare plot - distribution outline
  distr.map <- rgdal::readOGR("2-Preparing_Environmental_Data/tables/redlist_species_data/data_0.shp")
  # prepare plot - crop raster to extent
  diff_map.crop <- crop(rast, extent)
  # diff_map_r2 <- crop(diff_map.crop, extent(distr.map))
  # diff_map_r3 <- mask(diff_map_r2, distr.map)
  # prepare plot - colors
  dif_cols <- get_color_from_samples(samples)
  
  # create plot
  # plot(diff_map.crop, alpha = 0.7, col = rev(brewer.pal(6, "RdYlBu")),
  #      legend = F, box = F)
  # plot(diff_map_r3, col = rev(brewer.pal(6, "RdYlBu")),
  #      legend = T, margins = F, add = T)
  plot(diff_map.crop, alpha = 0.7, col = rev(brewer.pal(10, "RdYlBu")),
       legend = T, box = F, asp=NA)
  plot(distr.map, lwd = 0.8, add = T)
  points(clim.points$x, clim.points$y,
         pch = 1, lwd = 0.7, cex = 0.8,
         bg = dif_cols[match(clim.points$ID, samples)]) # the lynxes
}

# biplot of the gdm pca
prepare_gdm_biplot_list <- function(trans_clim_layers, raster_pca){
  plotcoords <- data.frame(raster_pca$x)
  
  envloading <- data.frame(raster_pca$rotation)
  row.names(envloading) <- translate_var_from_dict(row.names(envloading))
  
  pca.rast <- predict(trans_clim_layers, raster_pca,
                      index=1:ncol(raster_pca$rotation))
  
  # transform values to rgb
  for (n in 1:ncol(raster_pca$rotation)){
    pca.rast[[n]] <- (pca.rast[[n]]-pca.rast[[n]]@data@min) / (pca.rast[[n]]@data@max-pca.rast[[n]]@data@min) * 255
  }
  
  legend.colors <- na.omit(getValues(pca.rast[[c(1,2,3)]]))
  
  colors <- c()
  
  for (r in 1:NROW(legend.colors)){
    message('\r', round((r / NROW(legend.colors) * 100), 2), appendLF = FALSE)
    red <- legend.colors[r, 1] / 255
    green <- legend.colors[r, 2] / 255
    blue <- legend.colors[r, 3] / 255
    colo <- rgb(red, green, blue, alpha = 1)
    colors <- c(colors,colo)
  }
  return(list(plotcoords, envloading, colors))
}

plot_gdm_pca_biplot <- function(gdm_biplot_list){
  
  plotcoords <- gdm_biplot_list[[1]]
  envloading <- gdm_biplot_list[[2]]
  
  envloading_x <- envloading[which(abs(envloading[,1]) > 0.3), c(1,2)]
  envloading_y <- envloading[which(abs(envloading[,2]) > 0.3), c(1,2)]
  envloading <- rbind(envloading_x, envloading_y)

  colors <- gdm_biplot_list[[3]]
  
  PCA <- ggplot() +
    geom_point(aes(x = plotcoords$PC1, y = -plotcoords$PC2),
               size = 6, shape = 20, color = colors) +
    # geom_point(aes(x = plotcoords$PC1, y = plotcoords$PC2),
    #            size = 4, shape = 20) +
    # scale_x_continuous(limits = c(-0.25,0.7)) +
    geom_segment(aes(x = rep(0, nrow(envloading)),
                     y = rep(0, nrow(envloading)),
                     xend = c(envloading[,1]/3),
                     yend = -c(envloading[,2]/3)),
                 color = "black",
                 arrow = arrow(angle = 20, length = unit(0.4,"cm"),
                               ends = "last", type = "open"),
                 size = 1.2) +
    theme_bw() +
    theme_void() + # remove background, grid, numeric labels
    theme(panel.background = element_blank(),
          # panel.grid = element_blank(),
          plot.background = element_blank())
  return(PCA)
}
