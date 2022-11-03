library(raster)
library(tidyverse)

# create empty raster stack for yearly mean of january snow depth:
jan_mean_depth_stack <- raster()
# create empty raster stack for yearly mean days with snow cover:
snow_days_stack <- raster()

for (y in 1999:2018){
  print(y)
  # load year stack
  data <- stack(paste0("../tables/snow_data_daily/cmc_sdepth_dly_", y, "_v01.2.tif"))
  
  ## JANUARY DEPTH
  # get data of january only
  print("calculating january depth")
  data_jan <- data[[c(1:31)]]
  # create a raster that is the mean of january depth values for year y
  mean_jan <- calc(data_jan, fun = mean)
  # add the year's mean january depth stack to raster stack of yearly values
  jan_mean_depth_stack <- stack(jan_mean_depth_stack, mean_jan)
  
  ## MEAN DAYS WITH SNOW COVER
  # binarize data (day has snow = 1 or not = 0)
  print("calculating snow cover days")
  binary_data <- data > 0
  # create a raster that is the number of days with snow cover for year y
  sum_snow <- calc(binary_data, fun = sum)
  # add the year's number of days with snow cover stack to raster stack of yearly values
  snow_days_stack <- stack(snow_days_stack, sum_snow)
}
# create a raster that is the 20 year average of yearly january mean depth values 
jan_mean_depth <- calc(jan_mean_depth_stack, fun = mean)
# plot(jan_mean_depth)
writeRaster(jan_mean_depth, filename = "../tables/jan_mean_depth.tif", overwrite=TRUE)

# create a raster that is the 20 year average of yearly number of days with snow cover
mean_snow_days <- calc(snow_days_stack, fun = mean)
# plot(mean_snow_days)
writeRaster(mean_snow_days, filename = "../tables/mean_snow_days.tif", overwrite=TRUE)
