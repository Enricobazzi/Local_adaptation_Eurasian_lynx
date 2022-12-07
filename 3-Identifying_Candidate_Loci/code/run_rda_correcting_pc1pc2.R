library(tidyverse)
library(vegan)
library(adegenet)
library(viridis)
library(RColorBrewer)
library(robust)
library(qvalue)

# outliers based on qvalue
rdadapt<-function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# outliers based on loading deviation
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]              # locus names in these tails
}

#######################
## prepare data

## load genetic data and per sample environmental data table
# load whole genome genetic data
gt_data <- read.PLINK("1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.raw")
gt_data_tsv <- data.frame(as.matrix(gt_data))

# load matrix of environmental data
env.predictors <- read_delim("2-Preparing_Environmental_Data/tables/allvars_persample_table.tsv",
                             col_names = T, delim = "\t") %>%
  filter(sample %in% rownames(gt_data_tsv) == T) %>%
  column_to_rownames(var = "sample")

## to add PCs to environmental variables (see variance_partitioning_in_rda_model.R):
# load neutral genetic data
gt_data_neutral <- read.PLINK("1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.nogenes.raw")
gt_data_neutral_tsv <- data.frame(as.matrix(gt_data_neutral))

# run PCA on genotype data (rda with no variables == PCA)
pca <- rda(gt_data_neutral, scale=T)

# as I will correct for PC1 and PC2 I keep these two PCs only
PCs <- scores(pca, choices=c(1:2), display="sites", scaling=0)
# have PC table with same row order as env.predictors for cbind
PCs <- PCs[match(rownames(env.predictors), rownames(PCs)),]

# prepare variables the data frame
variables <- data.matrix(cbind(data.frame(PCs), data.frame(env.predictors[,c(9, 2, 15, 16, 21)])))
variables <- data.frame(variables)

#######################
## Run RDA

RDA <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days + 
                  Condition(PC1 + PC2), data = variables, scale = T)

# RDA <- rda(gt_data ~ bio9 + bio2 + bio15 + bio16 + mean_snow_days + 
#              Condition(PC2), data = variables, scale = T)
# RDA <- rda(gt_data ~ bio9 + bio2 + bio16 + mean_snow_days + 
#             Condition(PC1 + PC2), data = variables, scale = T)

RDA
#######################
## outlier calculations

# Chosing the best number of axes for the analysis
pdf("3-Identifying_Candidate_Loci/plots/rda_inertia_per_axis.pdf", width = 8, height = 8)
ggplot() +
  geom_line(aes(x=c(1:length(RDA$CCA$eig)), y=as.vector(RDA$CCA$eig)), linetype="dotted",
            size = 1.5, color="darkgrey") + geom_point(aes(x=c(1:length(RDA$CCA$eig)), y=as.vector(RDA$CCA$eig)), size = 3,
                                                       color="darkgrey") +
  scale_x_discrete(name = "Ordination axes", limits=c(1:9)) + ylab("Inertia") +
  theme_bw()
dev.off()

# from inspection of inertia per RDA axis 2/3 axes were chosen

################################################################################
# Using rdadapt function with K=3 gives no outliers!
# Using rdadapt function with K=2
res_rdadapt <- rdadapt(RDA, 2)
# outliers loci (q.values < 0.1)
length(which(res_rdadapt[,2] < 0.1))
hist(res_rdadapt[c(which(res_rdadapt[,2] < 0.1)),2])
# gives very few outliers! maybe going with SD is better?
# let's see what happens with these few outliers

nsig.axes <- 2
load.rda <- scores(RDA, choices=c(1:nsig.axes), display="species")
cand <- data.frame(load.rda[c(which(res_rdadapt[,2] < 0.1)),])
cand <- (cbind(cand, res_rdadapt[c(which(res_rdadapt[,2] < 0.1)),])) %>% 
  rownames_to_column("snp")

ncand <- NROW(cand)

# add correlation of candidate SNPs with all environmental variables
foo <- matrix(nrow=(ncand), ncol=NCOL(variables[,-c(1,2)]))
colnames(foo) <- colnames(variables[,-c(1,2)])
for (i in 1:length(cand$snp)) {
  nam <- gsub(":", ".", cand[i,1])
  snp.gen <- gt_data_tsv[,nam]
  foo[i,] <- apply(variables[,-c(1,2)], 2, function(x) cor(x, snp.gen))
}
cand <- cbind.data.frame(cand,foo)
cand$predictor <- NA
cand$correlation <- NA
# assign SNPs to predictor based on highest correlation value
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,NCOL(bar)-1] <- names(which.max(abs(bar[,c((4+nsig.axes):(NCOL(bar)-2))]))) # gives the variable
  cand[i,NCOL(bar)] <- max(abs(bar[,c((4+nsig.axes):(NCOL(bar)-2))]))              # gives the correlation
}

################################################################################
#
## find outliers as extreme loading values for significant axes
nsig.axes <- 2
# extract significant axes of RDA
load.rda <- scores(RDA, choices=c(1:nsig.axes), display="species")
# create a new data.frame with candidate snps and loading values
cand <- data.frame()
for(i in 1:nsig.axes){
  candN <- outliers(load.rda[,i],3)
  candN <- cbind.data.frame(rep(i,times=length(candN)), names(candN), unname(candN))
  colnames(candN) <- c("axis","snp","loading")
  cand <- rbind(cand, candN)
}
cand$snp <- as.character(cand$snp)

# check how many SNPs are duplicated (snps that are outliers for more than one axis)
NROW(cand[duplicated(cand$snp),])
# only one! -> Super_Scaffold_4:69703703_G

# remove duplicate entries
cand <- cand[!duplicated(cand$snp),]
# total number of unique candidate SNPs = 8413
ncand <- NROW(cand)

# THERE ARE ONLY 7 OUTLIER SNPS FROM RDA AXIS 1 - DOESN'T SEEM VERY GOOD

# add correlation of candidate SNPs with all environmental variables
foo <- matrix(nrow=(ncand), ncol=NCOL(variables[,-c(1,2)]))
colnames(foo) <- colnames(variables[,-c(1,2)])
for (i in 1:length(cand$snp)) {
  nam <- gsub(":", ".", cand[i,2])
  snp.gen <- gt_data_tsv[,nam]
  foo[i,] <- apply(variables[,-c(1,2)], 2, function(x) cor(x, snp.gen))
}
cand <- cbind.data.frame(cand,foo)  

# assign SNPs to predictor based on highest correlation value
cand$predictor <- NA
cand$correlation <- NA
for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,(4+(NCOL(variables[,-c(1,2)])))] <- names(which.max(abs(bar[4:(3+(NCOL(variables[,-c(1,2)])))]))) # gives the variable
  cand[i,(5+(NCOL(variables[,-c(1,2)])))] <- max(abs(bar[4:(3+(NCOL(variables[,-c(1,2)])))]))              # gives the correlation
}

######################################
## plot RDA results
sel <- cand$snp
env <- cand$predictor
#
# choose and add colors for each predictor
env[env=="bio2"] <- 'purple'
env[env=="bio9"] <- 'chartreuse3'
env[env=="bio16"] <- 'coral3'
env[env=="mean_snow_days"] <- 'dodgerblue4'
env[env=="bio15"] <- 'gold'

# pull all the SNP names
col.pred <- rownames(RDA$CCA$v)

# replace snp name with its color code (only candidates here)
for (i in 1:length(sel)) {
  foo <- match(sel[i],col.pred)
  col.pred[foo] <- env[i]
}
# replace non-candidates with grey
col.pred[grep("affold",col.pred)] <- '#f1eef6'

# create an empty (transparent) list for when drawing only candidates
empty <- col.pred
empty[grep("#f1eef6",empty)] <- rgb(0,1,0, alpha=0) # transparent
# and same but with grey for candidates
empty.outline <- ifelse(empty=="#00FF0000","#00FF0000","gray32")
# list all colors
bg <- c('purple','chartreuse3','coral3','dodgerblue4','gold')

pdf(file = paste0("try_rda1_rda2_snps.pdf"),
    width = 8,
    height = 8)

plot(RDA, type="n", 
     xlim=c(-(abs(load.rda[which.min(load.rda[,1]),1])*1.5),(abs(load.rda[which.max(load.rda[,1]),1])*1.5)), 
     ylim=c(-(abs(load.rda[which.min(load.rda[,2]),1])*1.5),(abs(load.rda[which.max(load.rda[,2]),1])*1.5)))
points(RDA, display="species", pch=4, cex=0.7, col="gray32", bg=col.pred)
points(RDA, display="species", pch=21, cex=1, col=empty.outline, bg=empty)
text(RDA, scaling=3, display="bp", col="#0868ac", cex=1)
legend("bottomright", legend=c("bio2", "bio9", "bio16", "mean_snow_days", "bio15"), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)

dev.off()


# Manhattan plot (outliers : q.value < 0.1 are colored in orange)

# Projection of loci in the RDA space

# transform our candidate table into a BED for downstream analyses
