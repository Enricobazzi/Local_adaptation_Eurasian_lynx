library(tidyverse)
library(psych)
library(vegan)
library(adegenet)
library(viridis)
library(RColorBrewer)
library(corrplot)

# load neutral genetic data
gt_data_neutral <- read.PLINK("1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.nogenes.raw")
gt_data_neutral_tsv <- data.frame(as.matrix(gt_data_neutral))

# load whole genome genetic data
gt_data <- read.PLINK("1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))

# load matrix of environmental data
env.predictors <- read_delim("2-Preparing_Environmental_Data/tables/allvars_persample_table.tsv",
                             col_names = T, delim = "\t") %>%
  filter(sample %in% rownames(gt_data_tsv) == T) %>%
  column_to_rownames(var = "sample")

# load coordinate data
coord_table <- read_delim("2-Preparing_Environmental_Data/tables/wholeset_coordinates.csv",
                          col_names = T, delim = ',') %>%
  filter(id %in% rownames(gt_data_tsv) == T) %>%
  column_to_rownames(var = "id")

### inferring population structure ###

# run PCA on genotype data (rda with no variables == PCA)
pca <- rda(gt_data_neutral, scale=T)

# eigenvalues for unconstrained axes:
#   PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8 
# 60297 34985 16033 12305 11991  9708  9116  8678 

# first three PCs seems good:
PCs <- scores(pca, choices=c(1:3), display="sites", scaling=0)
# have PC table with same row order as env.predictors for cbind
PCs <- PCs[match(rownames(env.predictors), rownames(PCs)),]

### building the matrix for the model ###

# build matrix with all components I want to include in full model => selected variables + genetic structure + geography

# first a data.matrix to get all values as numeric
variables <- data.matrix(cbind(data.frame(PCs), data.frame(env.predictors[,c(9, 2, 15, 16, 21)]), data.frame(coord_table[,c(1,2)])))

# plot correlations:
M <- cor(variables)
pdf("2-Preparing_Environmental_Data/plots/model_vars_correlogram.pdf", width = 8, height = 8)
corrplot(M, type="upper", col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

# convert variables to data frame again
variables <- data.frame(variables)

### run full model ###
pRDAfull <- rda(gt_data ~ PC1 + PC2 + PC3 + longitude + latitude + bio9 + bio2 + bio15 + bio16 + mean_snow_days, 
                data = variables, scale = T)

# RsquareAdj(pRDAfull)
# $r.squared
# [1] 0.2615709
# 
# $adj.r.squared
# [1] 0.1342556

# anova(pRDAfull)
#
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gt_data ~ PC1 + PC2 + PC3 + longitude + latitude + bio9 + bio2 + bio15 + bio16 + mean_snow_days, data = variables, scale = T)
#          Df Variance      F Pr(>F)
# Model    10   252685 2.0545  0.001 ***
# Residual 58   713344
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1


## run pure climate model ##
pRDAclim <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days + 
                  Condition(PC1 + PC2 + PC3 + longitude + latitude),
                data = variables, scale = T)

# RsquareAdj(pRDAclim)
# $r.squared
# [1] 0.08622409
# $adj.r.squared
# [1] 0.0243574

# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days + Condition(PC1 + PC2 + PC3 + longitude + latitude), data = variables, scale = T)
#          Df Variance      F Pr(>F)
# Model     5    83295 1.3545  0.001 ***
# Residual 58   713344
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## pure neutral population structure model

pRDAstruct <- rda(gt_data ~ PC1 + PC2 + PC3 + 
                  Condition(longitude + latitude + bio9 + bio2 + bio15 + bio16 + mean_snow_days),
                data = variables, scale = T)

# RsquareAdj(pRDAstruct)
# $r.squared
# [1] 0.05713851
# $adj.r.squared
# [1] 0.02111779

# anova(pRDAstruct)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gt_data ~ PC1 + PC2 + PC3 + Condition(longitude + latitude + bio9 + bio2 + bio15 + bio16 + mean_snow_days), data = variables, scale = T)
#          Df Variance     F Pr(>F)
# Model     3    55197 1.496  0.001 ***
# Residual 58   713344
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pRDAgeog <- rda(gt_data ~ longitude + latitude +
                  Condition(PC1 + PC2 + PC3 + bio9 + bio2 + bio15 + bio16 + mean_snow_days),
                data = variables, scale = T)

# RsquareAdj(pRDAgeog)
# $r.squared
# [1] 0.04364243
# $adj.r.squared
# [1] 0.02060328

# anova(pRDAgeog)
# Permutation test for rda under reduced model
# Permutation: free
# Number of permutations: 999
# 
# Model: rda(formula = gt_data ~ longitude + latitude + Condition(PC1 + PC2 + PC3 + bio9 + bio2 + bio15 + bio16 + mean_snow_days), data = variables, scale = T)
#          Df Variance     F Pr(>F)
# Model     2    42160 1.714  0.001 ***
# Residual 58   713344
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



