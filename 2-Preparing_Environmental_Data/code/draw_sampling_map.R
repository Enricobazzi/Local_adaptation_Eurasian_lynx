library(tidyverse)
library(raster)
library(rgdal)

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

# extent
extent <- c(-10, 180, 20, 80)

# distribution map
distr.map <- rgdal::readOGR("2-Preparing_Environmental_Data/tables/redlist_species_data/data_0.shp")

# land by turning bio1 into 1 as it has water as NA
worldclim <- getData("worldclim", path = "2-Preparing_Environmental_Data/tables/",
                     var = "bio", res = 10)
world <- worldclim[[1]] > -Inf
world.crop <- crop(world, extent)

# green color palette for land
pal <- colorRampPalette(c("darkolivegreen4","darkolivegreen4"))

# land borders
pp <- rasterToPolygons(world, dissolve=TRUE)

# empty raster for sea color
empty <- raster(crs=world)
values(empty) <- 0
empty.crop <- crop(empty, extent)

# land map with mountains from Natural Earth
world_map <- raster("2-Preparing_Environmental_Data/tables/HYP_HR_SR_W/HYP_HR_SR_W.tif")
world_map.crop <- crop(world_map, extent)

# samples
coord_table <- read_delim("2-Preparing_Environmental_Data/tables/wholeset_coordinates.csv",
                          col_names = T, delim = ',')

# Plot the grayscale world map
pdf("2-Preparing_Environmental_Data/plots/distribution_map.pdf", width = 12.5)

plot(empty.crop, col=c('lightcyan2'), asp=NA, legend = F)

plot(world.crop, col=pal(3), legend = F, add = T)

plot(world_map.crop, col=gray(seq(0, 1, length=100), alpha = 0.4),
    legend = F, add = T)

plot(crop(pp, extent), lwd=1.2, border='black', add=TRUE)

plot(crop(distr.map, extent), lwd=1.2, border='darkgreen', add=T)

points(coord_table$longitude, coord_table$latitude,
       pch=21, lwd=0.7, cex=1, 
       bg=get_color_from_samples(coord_table$id))

dev.off()

