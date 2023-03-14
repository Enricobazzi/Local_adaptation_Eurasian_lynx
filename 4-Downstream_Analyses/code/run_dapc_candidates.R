source("4-Downstream_Analyses/code/downstream_analyses_functions.R", print.eval=TRUE)

## DAPC of adaptive loci

# parse arguments - formula file line
args = commandArgs(trailingOnly=TRUE)
formula_file_line = args[1]

# read formula table
formula_file = "3-Identifying_Candidate_Loci/tables/rda_exploration_formulas.txt"
form = read.table(formula_file, header = T, sep = "\t")
v_form = c(form[formula_file_line,1], form[formula_file_line,2])
v_form_complete = paste0("vars_", v_form[1], ".cond_", v_form[2])

# raw files
sd_candidate_snps_raw = paste0("4-Downstream_Analyses/tables/", 
                               v_form_complete, ".sd_candidate_snps.raw")
# get gt data frames
sd_candidate_snps = get_gt_df_from_raw(raw_table = sd_candidate_snps_raw)

# interactive DAPC
grp <- find.clusters(sd_candidate_snps)

# DAPC based on criterion
grp <- find.clusters(sd_candidate_snps, n.pca = 103,
                     choose.n.clust = F, criterion = "goodfit")
levels(grp$grp)
# preditermined DAPC
grp <- find.clusters(sd_candidate_snps, n.pca = 103,
                     n.clust=4, max.n.clust = 20)

# groups
candidate_grps <- grp$grp

xval <- xvalDapc(sd_candidate_snps, candidate_grps, n.pca.max = 200,
                 result = "overall",
                 n.pca = NULL, n.rep = 100, xval.plot = TRUE)

candidate_DAPC <- dapc(sd_candidate_snps, candidate_grps,
                       n.pca = xval$DAPC$n.pca,
                       n.da = xval$DAPC$n.da)

candidate_dapc_matrix <- data.frame(candidate_DAPC$ind.coord)
discr_func_var_explained <- round((candidate_DAPC$eig/sum(candidate_DAPC$eig) * 100), 2)


####
candidate_DAPC$posterior

myCol <- c(brewer.pal(n = 8, name = "Dark2"))
scatter(candidate_DAPC, col=myCol,
        fill=get_color_from_samples(rownames(candidate_dapc_matrix)),
        scree.da=FALSE, xax = 1, yax = 2,
        cell=1.5, cex=2, bg="white", cstar=1)

ggplot() +
  geom_point(aes(x = candidate_dapc_matrix$LD1, y = candidate_dapc_matrix$LD2),
             shape = 21,
             fill = get_color_from_samples(rownames(candidate_dapc_matrix)))
ggplot() +
  geom_point(aes(x = candidate_dapc_matrix$LD1, y = candidate_dapc_matrix$LD3),
             shape = 21,
             fill = get_color_from_samples(rownames(candidate_dapc_matrix)))

library(scatterplot3d)

# Use scatterplot3d to create a 3D scatter plot
scatterplot3d(candidate_dapc_matrix$LD1,
              candidate_dapc_matrix$LD2,
              candidate_dapc_matrix$LD3,
              type = "h", pch = 21, box = F,
              bg = get_color_from_samples(rownames(candidate_dapc_matrix)))


run_dapc_fast <- function(n.clust){
  # preditermined DAPC
  grp <- find.clusters(sd_candidate_snps, n.pca = 103,
                       n.clust=n.clust, max.n.clust = 20)
  candidate_grps <- grp$grp
  
  # xval <- xvalDapc(sd_candidate_snps, candidate_grps, n.pca.max = 200,
  #                  result = "overall",
  #                  n.pca = NULL, n.rep = 100, xval.plot = FALSE)
  
  candidate_DAPC <- dapc(sd_candidate_snps, candidate_grps,
                         n.pca = 10,
                         n.da = 3)
  
  candidate_dapc_matrix <- data.frame(candidate_DAPC$ind.coord)
  
  scatterplot3d(candidate_dapc_matrix$LD1,
                candidate_dapc_matrix$LD2,
                candidate_dapc_matrix$LD3,
                type = "h", pch = 21, box = F,
                bg = get_color_from_samples(rownames(candidate_dapc_matrix)))
}
for (i in 1:30){
  run_dapc_fast(6)
}

candidate_DAPC
grp$size
