source("4-Downstream_Analyses/code/downstream_analyses_functions.R", print.eval=TRUE)
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fgeb.13459&file=geb13459-sup-0001-AppendixS1.pdf
library(gdm)
library(vegan)

# parse arguments - formula file line
args = commandArgs(trailingOnly=TRUE)
formula_file_line = args[1]

# read formula table
formula_file = "3-Identifying_Candidate_Loci/tables/rda_exploration_formulas.txt"
form = read.table(formula_file, header = T, sep = "\t")
v_form = c(form[formula_file_line,1], form[formula_file_line,2])
v_form_complete = paste0("vars_", v_form[1], ".cond_", v_form[2])

## ## ## ## ## ## ## 
## prepare data ####

# raw files
sd_candidate_snps_raw = paste0("4-Downstream_Analyses/tables/",
                               v_form_complete, ".sd_candidate_snps.raw")
neutral_snps_raw = paste0("4-Downstream_Analyses/tables/",
                          v_form_complete, ".neutral_snps.raw")

# get gt data frames
sd_candidate_snps = get_gt_df_from_raw(raw_table = sd_candidate_snps_raw)
neutral_snps = get_gt_df_from_raw(raw_table = neutral_snps_raw)
# neutral_snps = neutral_snps[, sample(ncol(neutral_snps), 1000)]
# calculate euclidian distance between samples based on SNP data (dissimilarity)
neutral_snp.dist.1 = get_snp_dist(gt_df = neutral_snps)
candidate_snp.dist.1 = get_snp_dist(gt_df = sd_candidate_snps)

# samples
samples = row.names(sd_candidate_snps)

# climate points - df with all variables and x,y coordinates
clim.points = get_clim_points(sample_list = samples)

# variables to be analyzed
# variables = strsplit(x = v_form[1], split = "\\+")[[1]]
variables = colnames(clim.points)[c(2:22)]
# prepare specific clim.points with only variables to be analyzed
# clim.points = clim.points[,c("ID", variables, "x", "y")]

# prepare gdm input
neutral_gdm.input = prepare_gdm_input(snp.dist.1 = neutral_snp.dist.1,
                                      clim.points = clim.points)
candidate_gdm.input = prepare_gdm_input(snp.dist.1 = candidate_snp.dist.1,
                                        clim.points = clim.points)


## ## ## ## ## ## ## 
## running gdm  ####

# run gdm
neutral_gdm <- gdm(neutral_gdm.input, geo = T, splines = NULL, knots = NULL)
candidate_gdm <- gdm(candidate_gdm.input, geo = F, splines = NULL, knots = NULL)
geo_candidate_gdm <- gdm(candidate_gdm.input, geo = T, splines = NULL, knots = NULL)

# print summaries
# print(summary(neutral_gdm))
# print(summary(candidate_gdm))
# print(summary(geo_candidate_gdm))

# print power of GDM
print(paste("overall power of neutral gdm:", round(neutral_gdm$explained, 2)))
print(paste("overall power of candidate gdm:", round(candidate_gdm$explained, 2)))
print(paste("overall power of candidate gdm:", round(geo_candidate_gdm$explained, 2)))

## ## ## ## ## ## ## ## ## ##
##  process gdm results  ####

# load raw climate data rasters
clim.layer = get_clim_layer(var_list = variables)
# xm <- matrix(xFromCell(s3, c(1:1944000)), nrow = 900, byrow = TRUE)
# x <- raster(xm, xmn = -180, xmx = 180, ymn = -60, ymx = 90)
# names(x) <- "xCoord"
# projection(x) <- "+proj=longlat +datum=WGS84 +no_defs"
# # plot(x)
# ym <- matrix(yFromCell(s3, c(1:1944000)), nrow = 900, byrow = TRUE)
# y <- raster(ym, xmn = -180, xmx = 180, ymn = -60, ymx = 90)
# projection(y) <- "+proj=longlat +datum=WGS84 +no_defs"
# names(y) <- "yCoord"
# # plot(y)
# xy_layer <- stack(x, y, clim.layer)
# s3 <- mask(xy_layer, calc(xy_layer, fun = sum))
s3 <- mask(clim.layer, calc(clim.layer, fun = sum))

# worldclim <- getData("worldclim", path = "2-Preparing_Environmental_Data/tables/",
#                      var = "bio", res = 10)
# jan_mean_depth <- raster("2-Preparing_Environmental_Data/tables/jan_mean_depth.tif")
# mean_snow_days <- raster("2-Preparing_Environmental_Data/tables/mean_snow_days.tif")
# clim.layer <- stack(worldclim,
#                     projectRaster(jan_mean_depth, worldclim),
#                     projectRaster(mean_snow_days, worldclim))

# transform climate rasters based on gdm
neutral_clim.trans = gdm.transform(neutral_gdm, s3)
candidate_clim.trans = gdm.transform(candidate_gdm, s3)
# geo_candidate_clim.trans = gdm.transform(geo_candidate_gdm, s3)
# geo_candidate_clim.trans = outputRasts

# run PCA on transformed raster values
neutral_pca = run_pca_raster(raster_obj = neutral_clim.trans)
candidate_pca = run_pca_raster(raster_obj = candidate_clim.trans)
# geo_candidate_pca = run_pca_raster(raster_obj = geo_candidate_clim.trans)

# write predictor tables of PCAs
write_pca_table(raster_pca_obj = neutral_pca,
                file_name = paste0("4-Downstream_Analyses/tables/",
                                   v_form_complete, ".gdm_neutral_pca_table.tsv"))
write_pca_table(raster_pca_obj = candidate_pca,
                file_name = paste0("4-Downstream_Analyses/tables/",
                                   v_form_complete, ".gdm_candidate_pca_table.tsv"))

# create a raster with each layer representing the pca transformation of 
# the original raster with each layer representing one principal component
neutral.pca.rast.crop = create_rgb_pca_raster(trans_clim_layers = neutral_clim.trans,
                      raster_pca = neutral_pca)

candidate.pca.rast.crop = create_rgb_pca_raster(trans_clim_layers = candidate_clim.trans,
                                              raster_pca = candidate_pca)

# additional_raster_layer = candidate.pca.rast.crop$layer.1 
# additional_raster_layer[additional_raster_layer > -9999] <- 0
# candidate.pca.rast.crop = stack(candidate.pca.rast.crop, additional_raster_layer)

## ## ## ## ## ## ## ## ## ##
##  plotting gdm results ####

# plot variable importance in neutral model
neutral_importance_plot = plot_importance(neutral_gdm)
ggsave(filename = "4-Downstream_Analyses/plots/gdm_line1_allvars_neutral_importance.pdf",
       plot = neutral_importance_plot,
       width = 60, height = 80, units = "mm")

# plot variable importance in candidate model
candidate_importance_plot = plot_importance(geo_candidate_gdm)
ggsave(filename = "4-Downstream_Analyses/plots/gdm_line1_allvars_candidate_importance.pdf",
       plot = candidate_importance_plot,
       width = 60, height = 80, units = "mm")

# plot candidate rgb map
pdf("4-Downstream_Analyses/plots/gdm_line1_allvars_candidate_rgb_123.pdf",
    width = 9, height = 6)
plot_rgb_map(sample_list = samples, pca.rast.crop = candidate.pca.rast.crop,
             r = 1, g = 2, b = 3)
dev.off()

# plot neutral rgb map
pdf("4-Downstream_Analyses/plots/gdm_line1_allvars_neutral_rgb_123.pdf",
    width = 9, height = 6)
plot_rgb_map(sample_list = samples, pca.rast.crop = neutral.pca.rast.crop,
             r = 1, g = 2, b = 3)
dev.off()

# plot neutral vs candidate residuals raster
pdf("4-Downstream_Analyses/plots/gdm_line1_allvars_residuals.pdf",
    width = 9, height = 6)
plot_procrust_resid_map(neutral_pca = neutral_pca, candidate_pca = candidate_pca,
                        rast = neutral_clim.trans[[1]],
                        samples = samples, clim.points = clim.points)
dev.off()


## LATVIA ZOOM
latvia_poland_extent <- c(17, 31, 50, 62)

empty.crop.latvia <- crop(empty.crop, latvia_poland_extent)
candidate.pca.rast.crop.latvia <- crop(candidate.pca.rast.crop, latvia_poland_extent)
neutral.pca.rast.crop.latvia <- crop(neutral.pca.rast.crop, latvia_poland_extent)
pp.latvia <- crop(pp, latvia_poland_extent)
distr.map.latvia <- crop(distr.map, latvia_poland_extent)
clim.points.latvia <- clim.points[grep("la|po", clim.points$ID), ]


pdf("4-Downstream_Analyses/plots/gdm_LATVIA_line1_allvars_neutral_rgb_123.pdf",
    width = 4, height = 4)
plot(empty.crop.latvia, col=c('lightcyan2'), legend=F, asp=NA)
plotRGB(neutral.pca.rast.crop.latvia, r = 1, g = 2, b = 3, bgalpha=0, add = T)
# plot(world_map.crop, col=gray(seq(0, 1, length=100), alpha = 0.4),
#     legend = F, add=T)
plot(pp.latvia, lwd=1.2, border='black', add=TRUE)
plot(distr.map.latvia, lwd=1.2, border='darkgreen', add=T)
points(clim.points.latvia$x, clim.points.latvia$y,
       pch=1, lwd=1.3, cex=1.3) # the lynxes
# text(clim.points.latvia$x-1, clim.points.latvia$y,
#      labels=clim.points.latvia$ID, cex = 0.3)
dev.off()

#
pdf("4-Downstream_Analyses/plots/gdm_LATVIA_line1_allvars_candidate_rgb_123.pdf",
    width = 4, height = 4)
plot(empty.crop.latvia, col=c('lightcyan2'), legend=F, asp=NA)
plotRGB(candidate.pca.rast.crop.latvia, r = 1, g = 2, b = 3, bgalpha=0, add = T)
# plot(world_map.crop, col=gray(seq(0, 1, length=100), alpha = 0.4),
#     legend = F, add=T)
plot(pp.latvia, lwd=1.2, border='black', add=TRUE)
plot(distr.map.latvia, lwd=1.2, border='darkgreen', add=T)
points(clim.points.latvia$x, clim.points.latvia$y,
       pch=1, lwd=1.3, cex=1.3) # the lynxes
# text(clim.points.latvia$x-1, clim.points.latvia$y,
#      labels=clim.points.latvia$ID, cex = 0.3)
dev.off()
## which LATVIA are in Courland Peninsula
# clim.points.latvia$ID[which(clim.points.latvia$x < 24 & clim.points.latvia$y > 55)]
# "c_ll_la_0044" "c_ll_la_0048" "c_ll_la_0052" "c_ll_la_0053"

# biplots ####
neutral_biplot_list = prepare_gdm_biplot_list(trans_clim_layers = neutral_clim.trans,
                                     raster_pca = neutral_pca)

envloading_x <- neutral_biplot_list[[2]][which(abs(neutral_biplot_list[[2]][,1]) > 0.3), c(1,2)]
envloading_y <- neutral_biplot_list[[2]][which(abs(neutral_biplot_list[[2]][,2]) > 0.3), c(1,2)]
envloading <- rbind(envloading_x, envloading_y)

neutral_biplot = plot_gdm_pca_biplot(gdm_biplot_list = neutral_biplot_list) +
  annotate("text",
           x = c(envloading[1,1]/1.7,  # xCoord
                 envloading[2,1]+0.1,  # yCoord
                 envloading[3,1]/2.7), # bio3
           y = -c(envloading[1,2]/2.7,  # xCoord
                 envloading[2,2]+0.3,  # yCoord
                 envloading[3,2]/2.9), # bio3
           label = rownames(envloading),
           color = "black", size = 10)

ggsave(filename = "4-Downstream_Analyses/plots/gdm_line1_allvars_neutral_pca_biplot.pdf",
       plot = neutral_biplot,
       width = 8, height = 6, units = "in")


candidate_biplot_list = prepare_gdm_biplot_list(trans_clim_layers = candidate_clim.trans,
                                              raster_pca = candidate_pca)

envloading_x <- candidate_biplot_list[[2]][which(abs(candidate_biplot_list[[2]][,1]) > 0.3), c(1,2)]
envloading_y <- candidate_biplot_list[[2]][which(abs(candidate_biplot_list[[2]][,2]) > 0.3), c(1,2)]
envloading <- envloading_x # both x and y are the same two variables

candidate_biplot = plot_gdm_pca_biplot(gdm_biplot_list = candidate_biplot_list) +
  annotate("text",
           x = c(envloading[1,1]/5.5,  # T_range_day
                 envloading[2,1]/2.7), # Mean_snow_days
           y = c(envloading[1,2]/2.3,  # T_range_day
                 envloading[2,2]/2.7), # Mean_snow_days
           label = rownames(envloading),
           color = "black", size = 10)


ggsave(filename = "4-Downstream_Analyses/plots/gdm_line1_allvars_candidate_pca_biplot.pdf",
       plot = candidate_biplot,
       width = 8, height = 6, units = "in")
