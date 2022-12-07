############################################
library(tidyverse)
library(vegan)
library(adegenet)
library(viridis)
library(RColorBrewer)
library(robust)
library(qvalue)

############################################

# load whole genome genetic data
read_raw_file <- function(raw_file) {
  gt_data <- read.PLINK(raw_file)
  return(gt_data)
}

get_tsv_from_gt_data <- function(gt_data){
  return(data.frame(as.matrix(gt_data)))
}

# load matrix of environmental data
read_env_data <- function(env_table, sample_list, sample_colname = "sample"){
  env.predictors <- read_delim(env_table, col_names = T,
                               delim = "\t", show_col_types = FALSE) %>%
  filter(sample %in% sample_list == T) %>%
  column_to_rownames(var = sample_colname)
  return(env.predictors)
}

## to add PCs to environmental variables (see variance_partitioning_in_rda_model.R):
get_pcs_from_ade <- function(gt_data, npcs = c(1:2)){
  # run PCA on genotype data (rda with no variables == PCA)
  PCA <- rda(gt_data, scale=T)
  # as I will correct for PC1 and PC2 I keep these two PCs only
  PCs <- scores(PCA, choices=npcs, display="sites", scaling=0)
  # have PC table with same row order as env.predictors for cbind
  PCs <- PCs[match(rownames(env.predictors), rownames(PCs)),]
  return(data.frame(PCs))
}

# prepare variables data frame
get_vars_df <- function(PCs, env.predictors, vars_list){
  variables <- data.matrix(cbind(data.frame(PCs), data.frame(env.predictors[,vars_list])))
  return(data.frame(variables))
}

# read formula from input file
read_formula_from_file <- function(formula_file, N){
  form <- read.table(formula_file, header = T, sep = "\t")
  v_form <- c(form[N,1], form[N,2])
  return(v_form)
}

# run rda
run_rda <- function(gt_data, variables, vars_list, cond_list) {
  
  Formula <- formula(paste0(deparse(substitute(gt_data)), " ~ ", paste(vars_list, collapse = " + "), 
                            " + Condition(", paste(cond_list, collapse = " + "), ")"))
  RDA <- rda(Formula, data = variables, scale = T)
  
  return(RDA)
}

# plot amount of inertia per ordination axis
plot_ord_x_inert_y <- function(plot_name, RDA){
  pdf(plot_name, width = 8, height = 8)
  plot = ggplot() +
    geom_line(aes(x=c(1:length(RDA$CCA$eig)), y=as.vector(RDA$CCA$eig)), linetype="dotted",
              size = 1.5, color="darkgrey") + geom_point(aes(x=c(1:length(RDA$CCA$eig)), y=as.vector(RDA$CCA$eig)), size = 3,
                                                         color="darkgrey") +
    scale_x_discrete(name = "Ordination axes", limits=c(1:9)) + ylab("Inertia") +
    theme_bw()
  print(plot)
  dev.off()
}

# outliers based on qvalue
rdadapt <- function(RDA,K) {
  zscores<-RDA$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  q.values_rdadapt<-qval$qvalues
  return(data.frame(p.values=reschi2test, q.values=q.values_rdadapt))
}

# get candidate (outlier SNPs) dataframe
get_cand_from_rdadapt <- function(RDA, p_q_vals, K, limit = 0.1){
  nsig.axes <- K
  load.rda <- scores(RDA, choices=c(1:nsig.axes), display="species")
  cand <- data.frame(load.rda[c(which(p_q_vals[,2] < limit)),])
  cand <- (cbind(cand, p_q_vals[c(which(p_q_vals[,2] < limit)),])) %>% 
    rownames_to_column("snp")
  ncand <- NROW(cand)
  print(paste("There are", ncand, "candidate SNPs!"))
  
  return(cand)
}

# outliers based on loading deviation
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]]              # locus names in these tails
}

# get candidate (outlier SNPs) dataframe from loadings exceeding X*SD
get_cand_from_sd <- function(RDA, K, limit = 3){
  ## find outliers as extreme loading values for significant axes
  nsig.axes <- K
  # extract significant axes of RDA
  load.rda <- scores(RDA, choices=c(1:nsig.axes), display="species")
  # create a new data.frame with candidate snps and loading values
  cand <- data.frame()
  for(i in 1:nsig.axes){
    candN <- outliers(load.rda[,i], limit)
    candN <- cbind.data.frame(rep(i,times=length(candN)), names(candN), unname(candN))
    colnames(candN) <- c("axis","snp","loading")
    cand <- rbind(cand, candN)
  }
  cand$snp <- as.character(cand$snp)
  
  # check how many SNPs are duplicated (snps that are outliers for more than one axis)
  
  print(paste("There are", NROW(cand[duplicated(cand$snp),]), "duplicate SNPs!"))
  # only one! -> Super_Scaffold_4:69703703_G
  
  # remove duplicate entries
  cand <- cand[!duplicated(cand$snp),]
  # total number of unique candidate SNPs = 8413
  ncand <- NROW(cand)
  print(paste("There are", ncand, "candidate SNPs!"))
  return(cand)
}


# add correlation of candidate SNPs with all environmental variables
get_corr_snps_vars <- function(cand, variables, vars_list){
  ncand <- NROW(cand)
  foo <- matrix(nrow=(ncand), ncol=length(vars_list)+2)
  colnames(foo) <- c(vars_list, "predictor", "correlation")
  for (i in 1:length(cand$snp)) {
    nam <- gsub(":", ".", cand[i,"snp"])
    snp.gen <- gt_data_tsv[,nam]
    foo[i,1:length(vars_list)] <- apply(variables[,which(colnames(variables) %in% vars_list)],
                                        2, function(x) cor(x, snp.gen))
    foo[i,length(vars_list)+1] <- vars_list[which.max(abs(as.numeric(foo[i,1:length(vars_list)])))]
    foo[i,length(vars_list)+2] <- max(abs(as.numeric(foo[i,1:length(vars_list)])))
  }
  return(foo)
}

# cand <- cbind.data.frame(cand,foo)


## plot RDA results

plot_cand_data <- function(cand, RDA, plot_name, x = 1, y = 2){
  sel <- cand$snp
  env <- cand$predictor
  vars <- unique(env)
  
  for (i in 1:length(vars)){
    var = vars[i]
    col = brewer.pal(9, "Set1")[i] 
    env[env==var] = col
  }
  
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
  bg <- unique(env)
  
  pdf(file = plot_name,
      width = 8,
      height = 8)
  
  plot(RDA, type="n", choices=c(x,y), 
       xlim=c(-0.1,0.1), 
       ylim=c(-0.1,0.1))
  points(RDA, display="species", pch=4, cex=0.7, col="gray32", bg=col.pred, choices=c(x,y))
  points(RDA, display="species", pch=21, cex=1, col=empty.outline, bg=empty, choices=c(x,y))
  text(RDA, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(x,y))
  legend("bottomright", legend=vars, bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
  
  dev.off()
}

# plot samples instead of SNPs
plot_samples_rda <- function(plot_name, sample_list, RDA, x = 1, y = 2){
  # Add populations and colors
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
  
  ola <- data.frame(sample = sample_list, pop = loc, n = num)
  eco <- levels(ola$pop)
  bg <- cols <- c("#A035AF",
                  brewer.pal(12,"Paired")[9],
                  "#B8860b",
                  viridis_pal()(5)[1],
                  brewer.pal(12,"Paired")[3],
                  brewer.pal(12,"Paired")[7],
                  viridis_pal()(5)[3],
                  viridis_pal()(5)[2],
                  brewer.pal(12,"Paired")[8],
                  "#0F4909",
                  brewer.pal(12,"Paired")[5],
                  brewer.pal(12,"Paired")[6])
 
  pdf(file = plot_name,
      width = 8,
      height = 8)
  # plot
  plot(RDA, type="n", scaling=3, choices=c(x,y))
  # points(rda, display="species", pch=4, cex=0.7, col="gray32", scaling=3)
  points(RDA, display="sites", pch=21, cex=1.5, col="gray32", scaling=3, bg=bg[num], choices=c(x,y))
  text(RDA, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(x,y))
  # legend("topright", legend=eco, bty="n", col="gray32", pch=21, cex=0.7, pt.bg=bg)
  dev.off()
}


######################################
# plot RDA samples + SNPs in same plot
plot_rda_snps_and_samples <- function(RDA, cand, x = 1, y = 2, plot_name){
  # data frame for plotting samples
  sample_df <- data.frame(scores(RDA, choices=c(1:3), display="sites", scaling="none"))
  
  # add populations
  loc <- rep(NA, length(rownames(sample_df)))
  loc[grep("ba", rownames(sample_df))] <- "Balkans"
  loc[grep("ca", rownames(sample_df))] <- "Caucasus"
  loc[grep("cr", rownames(sample_df))] <- "Carpathians"
  loc[grep("ka", rownames(sample_df))] <- "Mongolia"
  loc[grep("ki", rownames(sample_df))] <- "Kirov"
  loc[grep("la", rownames(sample_df))] <- "Latvia"
  loc[grep("no", rownames(sample_df))] <- "Norway"
  loc[grep("po", rownames(sample_df))] <- "NE-Poland"
  loc[grep("og", rownames(sample_df))] <- "Mongolia"
  loc[grep("to", rownames(sample_df))] <- "Mongolia"
  loc[grep("tu", rownames(sample_df))] <- "Tuva"
  loc[grep("ur", rownames(sample_df))] <- "Urals"
  loc[grep("vl", rownames(sample_df))] <- "Vladivostok"
  loc[grep("ya", rownames(sample_df))] <- "Yakutia"
  
  sample_df$population <- loc
  
  # add colors
  col <- rep(NA, length(rownames(sample_df)))
  col[grep("ba", rownames(sample_df))] <- "#A035AF"
  col[grep("ca", rownames(sample_df))] <- "#B8860b"
  col[grep("cr", rownames(sample_df))] <- "#CAB2D6"
  col[grep("ka", rownames(sample_df))] <- "#FDBF6F"
  col[grep("ki", rownames(sample_df))] <- "#440154FF"
  col[grep("la", rownames(sample_df))] <- "#B2DF8A"
  col[grep("no", rownames(sample_df))] <- "#3B528BFF"
  col[grep("po", rownames(sample_df))] <- "#21908CFF"
  col[grep("og", rownames(sample_df))] <- "#FDBF6F"
  col[grep("to", rownames(sample_df))] <- "#FDBF6F"
  col[grep("tu", rownames(sample_df))] <- "#FF7F00"
  col[grep("ur", rownames(sample_df))] <- "#0F4909"
  col[grep("vl", rownames(sample_df))] <- "#FB9A99"
  col[grep("ya", rownames(sample_df))] <- "#E31A1C"
  
  sample_df$color <- col
  
  # data frame for plotting SNPs
  snps_df <- data.frame(scores(RDA, choices=c(1:3), display="species", scaling="none"))
  snps_df$names <- row.names(snps_df)
  snps_df$type <- "Neutral"
  snps_df$type[snps_df$names %in% cand$snp] <- "Candidate"
  
  # add color
  sel <- cand$snp
  env <- cand$predictor
  vars <- unique(env)
  # list of all predictors
  vars <- unique(env)
  # assign a color to each predictor
  for (i in 1:length(vars)){
    var = vars[i]
    col = brewer.pal(9, "Set1")[i] 
    env[env==var] = col
  }
  # pull all the SNP names
  col.pred <- rownames(RDA$CCA$v)
  # replace snp name with its color code (only candidates here)
  for (i in 1:length(sel)) {
    foo <- match(sel[i],col.pred)
    col.pred[foo] <- env[i]
  }
  # replace non-candidates with grey
  col.pred[grep("affold",col.pred)] <- '#f1eef6'
  # add color column
  snps_df$color <- col.pred
  
  # divide SNPs into 2 data.frames
  neutral_snps_df <- snps_df %>% filter(type == "Neutral")
  candidate_snps_df <- snps_df %>% filter(type == "Candidate")
  
  # variables data frame
  TAB_var <- data.frame(scores(RDA, choices=c(1:3), display="bp"))
  
  # plot limits
  min_lim = -max(c(abs(sample_df[,x]), abs(sample_df[,y])))*1.1
  max_lim = max(c(abs(sample_df[,x]), abs(sample_df[,y])))*1.1
  
  # plot the RDA
  rda_plot <- ggplot() +
    geom_hline(yintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
    geom_vline(xintercept=0, linetype="dashed", color = gray(.80), size=0.3) +
    
    geom_point(aes(x = neutral_snps_df[,x] * 25, y = neutral_snps_df[,y] * 25),
               shape = 4, cex = 1, fill = neutral_snps_df$color, alpha = 0.5) +
    
    geom_segment(aes(xend=TAB_var[,x]/2, yend=TAB_var[,y]/2, x=0, y=0),
                 colour="black", size=0.15, linetype=1,
                 arrow=arrow(length = unit(0.01, "npc"))) +
    geom_text(aes(x=1.05*(TAB_var[,x]/2), y=1.05*(TAB_var[,y]/2), label = row.names(TAB_var)),
              size = 2.5, family = "Times") +
    
    geom_point(aes(x = candidate_snps_df[,x] * 25, y = candidate_snps_df[,y] * 25),
               shape = 21, cex = 1.5, fill = candidate_snps_df$color, alpha = 0.8) +
    
    geom_point(aes(x = sample_df[,x], y = sample_df[,y]),
               shape = 21, cex = 3, fill = sample_df$color, alpha = 0.8) +
    
    theme_bw(base_size = 14, base_family = "Times") +
    theme(panel.background = element_blank(), panel.grid = element_blank(),
          plot.background = element_blank()) +
    # xlim(min = min_lim, max = max_lim) +
    # ylim(min = min_lim, max = max_lim) +
    xlab(paste("RDA", x)) + ylab(paste("RDA", y))
  
  ggsave(filename = plot_name, plot = rda_plot, width = 8, height = 8, units = "in")
  
}

######################################
# manhattan plots

prepare_snps_dataframe_for_plotting <- function(RDA, K, cand_qvals, cand_sd){

  # data frame for plotting SNPs
  snps_df <- data.frame(scores(RDA, choices=c(1:K), display="species", scaling="none"))
  snps_df$names <- row.names(snps_df)
  snps_df$type <- "Neutral"
  snps_df$type[snps_df$names %in% cand_sd$snp] <- "Candidate_SD"
  snps_df$type[snps_df$names %in% cand_qvals$snp] <- "Candidate_Qval"
  
  # add color - use cand_sd as they are always a superset of cand_qvals
  sel <- cand_sd$snp
  correl <- cand_sd$predictor
  env <- cand_sd$predictor
  vars <- unique(env)
  # list of all predictors
  vars <- unique(env)
  all_vars <- c("bio9", "bio2", "bio15", "bio16", "mean_snow_days", "jan_mean_depth")
  # assign a color to each predictor
  for (i in 1:length(all_vars)){
    if (all_vars[i] %in% vars){
      var = all_vars[i]
      col = brewer.pal(9, "Set1")[i] 
    }
    env[env==var] = col
  }
  # pull all the SNP names
  var.pred <- rownames(RDA$CCA$v)
  col.pred <- rownames(RDA$CCA$v)
  # replace snp name with its color code (only candidates here)
  for (i in 1:length(sel)) {
    foo <- match(sel[i],col.pred)
    col.pred[foo] <- env[i]
    var.pred[foo] <- correl[i]
  }
  # replace non-candidates with grey
  col.pred[grep("affold",col.pred)] <- '#f1eef6'
  var.pred[grep("affold",var.pred)] <- 'none'
  
  # add color, correlation and snp number columns
  snps_df$color <- col.pred
  snps_df$correlation <- var.pred
  snps_df$snpnum <- 1:nrow(snps_df)
  
  # add chromosome and position
  names_list = unlist(strsplit(snps_df$names, ":"))
  chr_list = names_list[seq(1, length(names_list), 2)]
  pos_list = unlist(strsplit(names_list[seq(2, length(names_list), 2)], "_"))
  pos_list = pos_list[seq(1, length(pos_list), 2)]
  
  snps_df$chromosome <- chr_list
  snps_df$position <- pos_list
  
  return(snps_df)
}

plot_rda_manhattan <- function(snps_df, K, plots_prefix){
  
  # divide SNPs into 3 dfs
  no_sel_df <- snps_df %>% filter(type == "Neutral")
  cand_sd_df <- snps_df %>% filter(type == "Candidate_SD")
  cand_qval_df <- snps_df %>% filter(type == "Candidate_Qval")
  
  # 
  for (k in 1:K){
    plot_name = paste0(plots_prefix, "K", K, ".rda", k, ".manhattan_plot.pdf")
    
    man_plot <- ggplot() +
      geom_point(data = no_sel_df, aes(position, y = no_sel_df[,k]),
                 shape = 4, cex = 1, fill = no_sel_df$color, alpha = 0.5) +
      geom_point(data = cand_sd_df, aes(x = position, y = cand_sd_df[,k]),
                 shape = 21, cex = 1.5, fill = cand_sd_df$color, alpha = 0.8) +
      geom_point(data = cand_qval_df, aes(x = position, y = cand_qval_df[,k]),
                 shape = 22, cex = 1.8, fill = cand_qval_df$color, alpha = 0.8) +
      theme_bw(base_size = 14, base_family = "Times") +
      theme(panel.background = element_blank(), panel.grid = element_blank(),
            plot.background = element_blank(), axis.text.x = element_blank(),
            axis.ticks.x=element_blank()) +
      # xlim(min = min_lim, max = max_lim) +
      # ylim(min = min_lim, max = max_lim) +
      xlab("SNPs") + ylab(paste(colnames(no_sel_df)[k], "loadings")) +
      facet_wrap(~ chromosome, scales = "free_x")

      ggsave(filename = plot_name, plot = man_plot, width = 16, height = 8, units = "in")
    
  }
  
}

# compare_rda_results set names
get_name_set_abbreviation <- function(name_set){
  if (name_set == "(bio9+bio2+bio15+bio16+mean_snow_days)-(PC1+PC2)") {
    name_set_abb = "sd_bio15_pc12"
  } else if (name_set == "(bio9+bio2+bio15+bio16+mean_snow_days)-(PC2)") {
    name_set_abb = "sd_bio15_pc2"
  } else if (name_set == "(bio9+bio2+bio16+mean_snow_days)-(PC1+PC2)") {
    name_set_abb = "sd_pc12"
  } else if (name_set == "(bio9+bio2+bio16+mean_snow_days)-(PC2)") {
    name_set_abb = "sd_pc2"
  } else if (name_set == "(bio9+bio2+bio16+jan_mean_depth)-(PC1+PC2)") {
    name_set_abb = "jd_pc12"
  } else if (name_set == "(bio9+bio2+bio16+jan_mean_depth)-(PC2)") {
    name_set_abb = "jd_pc2"
  }
  return(name_set_abb)
}
