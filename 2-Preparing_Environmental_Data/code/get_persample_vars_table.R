# in this script I will prepare the table containing the values of the
# environmental variables for:
#  - each sample (persample_table) 

# prepare R:
library(raster)
library(tidyverse)
# library(grDevices)
library(rgeos)
library(geobuffer)

## CREATE ENVIRONMENTAL DATA RASTER
# load WorldClim data:
worldclim <- getData("worldclim", path = "../tables/", var = "bio", res = 10)
# load january depth and reproject on worldclim cs
jan_mean_depth <- raster("../tables/jan_mean_depth.tif")
jan_mean_depth_new <- projectRaster(jan_mean_depth, worldclim)
# load mean snow days and reproject on worldclim cs
mean_snow_days <- raster("../tables/mean_snow_days.tif")
mean_snow_new <- projectRaster(mean_snow_days, worldclim)
# combine climate layers
clim.layer <- stack(worldclim,jan_mean_depth_new,mean_snow_new)

## LOAD SAMPLES AND COORDINATES
# load coordinate data:
coord_table <- read_delim("../tables/wholeset_coordinates.csv", col_names = T, delim = ',') 
# define samples
samples <- coord_table$id %>% unique
# remove bad samples and ones without coordinates
elements_2_remove <- c("c_ll_cr_0211", "c_ll_ba_0233","c_ll_ba_0216", 
                       "h_ll_ba_0214", "h_ll_ba_0215")
samples = samples[!(samples %in% elements_2_remove)]

# EXTRACT VALUES FOR TABLES
# create a table with all variable values in a buffer of 100km around each sample
# see https://www.iucnredlist.org/species/12519/121707666#habitat-ecology for home range

# create a table to bind to the rest of the calculated data
table <- data.frame()

# now lets loop though all samples:
for (sample in samples){
  # get sample coordinates
  coord_sample <- coord_table %>% filter(coord_table$id == sample)
  # in a dataframe with columns x and y
  coords <- data.frame(x=as.numeric(coord_sample$longitude),
                       y=as.numeric(coord_sample$latitude))
  # get buffer around coordinates
  buff <- geobuffer_pts(xy = coords, dist_m = 100000)

  # create df for sample's row in general table
  row <- data.frame(sample=sample)
  # extract each variable's mean values in buffer
  for (k in 1:nlayers(clim.layer)){
    # subset for just one biovariable
    DEM <- clim.layer[[k]]
    cat(paste(sample, ":", names(DEM), "\r")) # print the name of the biovariable
    flush.console()
    # Get mean of biovariable value of all points within the polygon
    means <- raster::extract(
      DEM,                  # raster layer
      buff,               # SPDF with centroids for buffer
      fun=mean,             # get the mean
      na.rm = TRUE,         # remove NAs
      df=TRUE)              # return a dataframe
    
    # add the mean to the population table
    row <- cbind(row, data.frame(means[2]))
  }
  table <- rbind(table, row)
}

# write persample table
write.table(x = table, file = "../tables/allvars_persample_table.tsv", quote=FALSE, col.names = T, row.names = FALSE, sep= "\t")
