library(tidyverse)
library(psych)
library(vegan)
library(adegenet)
library(viridis)
library(RColorBrewer)

# load coordinate data
coord_table <- read_delim("2-Preparing_Environmental_Data/tables/wholeset_coordinates.csv",
                          col_names = T, delim = ',') 

# load genetic data
gt_data <- read.PLINK("1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))

# load matrix of environmental data
env.predictors <- read_delim("2-Preparing_Environmental_Data/tables/allvars_persample_table.tsv",
                             col_names = T, delim = "\t") %>%
  filter(sample %in% rownames(gt_data_tsv) == T) %>%
  column_to_rownames(var = "sample")

# run RDA

# reviewer 2 suggestion:
# stepwise selection approach needs the RDA object to evaluate variables in the model

# I start with running one RDA of a model with only an intercept and one RDA with all the variables
# (except for snow! - it would be eliminated if included but I want to keep it - see below how)
rda.0 <- rda (gt_data ~ 1, data = env.predictors[,-c(20,21)], scale=T) # model containing only gt matrix and intercept
rda.all <- rda (gt_data ~ ., data = env.predictors[,-c(20,21)], scale=T) # model including all variables from matrix chem1 (the dot after tilda (~) means ALL!)

# ordi2step as in tutorial by Capblancq
rda.sel <- ordiR2step(rda.0, rda.all, Pin = 0.01, R2permutations = 1000,
                      R2scope = T, trace = 10)

# rda.sel
# Call: rda(formula = gt_data ~ bio9 + bio2 + bio15 + bio16 + bio13, data
#           = env.predictors[, -c(20, 21)], scale = T)
# 
# Inertia Proportion Rank
# Total         9.660e+05  1.000e+00     
# Constrained   1.683e+05  1.743e-01    5
# Unconstrained 7.977e+05  8.257e-01   63
# Inertia is correlations 
# 
# Eigenvalues for constrained axes:
#   RDA1  RDA2  RDA3  RDA4  RDA5 
# 76330 43198 21026 15798 11997 
# 
# Eigenvalues for unconstrained axes:
#   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
# 58572 39751 27792 22878 18975 18061 16258 15482 
# (Showing 8 of 63 unconstrained eigenvalues)

# rda.sel$anova
#                   R2.adj Df    AIC      F Pr(>F)   
# + bio9          0.027322  1 950.94 2.9101  0.002 **
# + bio2          0.059092  1 949.62 3.2622  0.002 **
# + bio15         0.088930  1 948.34 3.1616  0.002 **
# + bio16         0.098198  1 948.56 1.6680  0.004 **
# + bio13         0.108735  1 948.67 1.7566  0.002 **
# <All variables> 0.154015                           
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## VIF calculations

usdm::vifcor(env.predictors[,c(9, 2, 15, 16, 13)], th=10)
# ---------- VIFs of the remained variables -------- 
#   Variables       VIF
# 1      bio9  3.040581
# 2      bio2  4.481904
# 3     bio15  4.639044
# 4     bio16 57.069018
# 5     bio13 59.442046

# remove bio13! 
usdm::vifcor(env.predictors[,c(9, 2, 15, 16)], th=10)
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1      bio9 2.456010
# 2      bio2 3.561682
# 3     bio15 4.362175
# 4     bio16 1.266282

# remove bio15! 
usdm::vifcor(env.predictors[,c(9, 2, 16)], th=10)
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1      bio9 1.951767
# 2      bio2 1.839526
# 3     bio16 1.084348

## adding snow variable

# add mean_snow_days and jan_depth separately and calculate model's adj.R2 for each:
rda.sel.ns <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16, data = env.predictors, scale = T)
# RsquareAdj(rda.sel.ns)
# $r.squared
# [1] 0.1512455
# 
# $adj.r.squared
# [1] 0.09819832

rda.sel.sd <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days,
                  data = env.predictors, scale = T)
# RsquareAdj(rda.sel.sd)
# $r.squared
# [1] 0.1660375
# 
# $adj.r.squared
# [1] 0.09984998



# vif with the addition of mean snow days
usdm::vifcor(env.predictors[,c(9, 2, 15, 16, 21)], th=10)
# max correlation ( mean_snow_days ~ bio9 ):  -0.8321142 
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1           bio9 6.049206
# 2           bio2 4.303837
# 3          bio15 5.410557
# 4          bio16 1.318614
# 5 mean_snow_days 4.119881

rda.sel.jd <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16 + jan_mean_depth,
                  data = env.predictors, scale = T)
# RsquareAdj(rda.sel.jd)
# $r.squared
# [1] 0.1653665
# 
# $adj.r.squared
# [1] 0.09912574


# vif with the addition of january depth
usdm::vifcor(env.predictors[,c(9, 2, 15, 16, 20)], th=10)
# max correlation ( bio2 ~ bio9 ):  -0.6729578 
# 
# ---------- VIFs of the remained variables -------- 
#   Variables      VIF
# 1           bio9 3.015221
# 2           bio2 3.837336
# 3          bio15 6.408592
# 4          bio16 1.346865
# 5 jan_mean_depth 1.546143

## anovas

# anova(rda.sel.sd)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days, data = env.predictors, scale = T)
# Df Variance      F Pr(>F)
# Model     5   160397 2.5086  0.001 ***
# Residual 63   805632
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# anova(rda.sel.ns)

# # correlograms
pdf("2-Preparing_Environmental_Data/plots/model_vars_pairs_panels.pdf",
    width = 13, height = 13)
pairs.panels(env.predictors[,c(9, 2, 13, 15, 16, 20, 21)], scale=T, cex.cor = 1.5)
dev.off()


# # clustering dendrogram
# cor_matrix <- abs(cor(env.predictors))
# var_clusters <- hclust(dist(cor_matrix), method = "ward.D")
# plot(var_clusters)


