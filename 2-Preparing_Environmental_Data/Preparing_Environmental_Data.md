---
title: "Preparing_Environmental_Data"
author: "Enrico"
date: "2022-10-17"
output: html_document
---

### Extracting and analyzing environmental data

In this MD I will describe the way that I extracted and analyzed the environmental data for our study.

The WorldClim data is easily downloaded using the raster package in R, but I needed to prepare the data for the two snow variables we were interested in, the average snow depth for the month of January and the yearly average number of days with snow cover. To do so, I downloaded snow depth data for from the Northern Hemisphere subset of the Canadian Meteorological Centre operational global daily snow depth analysis ([Brown, Ross & Brasnett, Bruce, 2010](https://nsidc.org/data/nsidc-0447/versions/1)) for the years bettween 1999 and 2018. Using a custom made [R script](code/create_snow_rasters.R) I transformed the daily data into rasters of the desired variables with the same coordinate system as the WorldClim data.

With the data rasters ready I used another custom [R script](code/get_perpop_vars_table_for_baypass.R) to generate a table with each sample's value for each of the 21 variables (19 worldclim + 2 snow). Each sample's value was calculated as the average of the values in a circle of the area of 100km2 around the sample's coordinate point (see the [IUCN page](https://www.iucnredlist.org/species/12519/121707666#habitat-ecology) for references regarding home range sizes).

As BayPass is a population based test (compared to the individual approach we adopt in the RDA), we need population averages for each of the variables, which were calculated and formatted as needed by BayPass using a custom [R script](code/get_perpop_vars_table_for_baypass.R).

```
cd /home/ebazzicalupo/Local_adaptation_Eurasian_lynx/2-Preparing_Environmental_Data/tables

varLIST=($(cut -f1 allvars_baypass_perpop_table.tsv | grep -v "vars"))

for var in ${varLIST[@]}
 do
  echo "creating ${var} table"
  grep -w "${var}" allvars_baypass_perpop_table.tsv |
  cut -f2-  | tr '\t' ' ' \
  > ${var}_data.txt
done
```

### BONUS: distribution map

The script with which I drew the distribution map is [here](code/draw_sampling_map.R)!!
