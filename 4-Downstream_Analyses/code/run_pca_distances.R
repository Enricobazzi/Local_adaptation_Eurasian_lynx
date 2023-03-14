source("4-Downstream_Analyses/code/downstream_analyses_functions.R", print.eval=TRUE)

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
# qvals_candidate_snps_raw = "4-Downstream_Analyses/tables/qvals_candidate_snps.raw"
neutral_snps_raw = paste0("4-Downstream_Analyses/tables/",
                          v_form_complete, ".neutral_snps.raw")

# get gt data frames
sd_candidate_snps = get_gt_df_from_raw(raw_table = sd_candidate_snps_raw)
# qvals_candidate_snps = get_gt_df_from_raw(raw_table = qvals_candidate_snps_raw)
neutral_snps = get_gt_df_from_raw(raw_table = neutral_snps_raw)
# neutral_snps = neutral_snps[, sample(ncol(neutral_snps), 1000)]

# run PCA 
sd_candidate_pca <- FactoMineR::PCA(sd_candidate_snps, graph = F)
# qvals_candidate_pca <- FactoMineR::PCA(qvals_candidate_snps, graph = F)
neutral_pca <- FactoMineR::PCA(neutral_snps, graph = F)

# build dataframes of PCA results
sd_candidate_pca_df = build_df_from_pca(pca_obj = sd_candidate_pca)
# qvals_candidate_pca_df = build_df_from_pca(pca_obj = qvals_candidate_pca)
neutral_pca_df = build_df_from_pca(pca_obj = neutral_pca)


# calculate percents of variance explained by pcs
sd_candidate_variance_percents <- sd_candidate_pca$eig[,2]
# qvals_candidate_variance_percents <- qvals_candidate_pca$eig[,2]
neutral_variance_percents <- neutral_pca$eig[,2]

# plot PCA results
sd_pca_plot = plot_pca(pca_df = sd_candidate_pca_df,
                       variance_percents = sd_candidate_variance_percents)
sd_pca_plot

ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
                         v_form_complete, ".pca_pc12_candidates.pdf"),
       plot = sd_pca_plot, width = 85, height = 80, units = "mm")

sd_pca_plot_3_4 = plot_pca(pca_df = sd_candidate_pca_df,
                       variance_percents = sd_candidate_variance_percents,
                       x = 3, y = 4)
sd_pca_plot_3_4

ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
                         v_form_complete, ".pca_pc34_candidates.pdf"),
       plot = sd_pca_plot_3_4, width = 85, height = 80, units = "mm")

sd_pca_plot_1_3 = plot_pca(pca_df = sd_candidate_pca_df,
                           variance_percents = sd_candidate_variance_percents,
                           y = 3)
sd_pca_plot_1_3

sd_pca_plot_1_4 = plot_pca(pca_df = sd_candidate_pca_df,
                           variance_percents = sd_candidate_variance_percents,
                           y = 4)
sd_pca_plot_1_4


# qvals_pca_plot = plot_pca(pca_df = qvals_candidate_pca_df,
#                        variance_percents = qvals_candidate_variance_percents)
# qvals_pca_plot
# 
# qvals_pca_plot_1_3 = plot_pca(pca_df = qvals_candidate_pca_df,
#                           variance_percents = qvals_candidate_variance_percents,
#                           y = 3)
# qvals_pca_plot_1_3

neutral_pca_plot = plot_pca(pca_df = neutral_pca_df,
                          variance_percents = neutral_variance_percents)
neutral_pca_plot

ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
                         v_form_complete, ".pca_pc12_neutral.pdf"),
       plot = neutral_pca_plot, width = 85, height = 80, units = "mm")


neutral_pca_plot_1_3 = plot_pca(pca_df = neutral_pca_df,
                            variance_percents = neutral_variance_percents,
                            x = 1, y = 3)
neutral_pca_plot_1_3

ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
                         v_form_complete, ".pca_pc13_neutral.pdf"),
       plot = neutral_pca_plot_1_3, width = 85, height = 80, units = "mm")

# coinertia between neutral and candidates 
# coin = run_coin_from_gt_dfs(gt_df1 = neutral_snps, gt_df2 = sd_candidate_snps)
# print(paste("co-inertia RV is", round(coin$RV, 3)))
# variance_percents = calculate_variance_percents_from_eigenvals(coin$eig)
# 
# 
# # build coinertia data.frame and plot
# coin_df = build_coin_df_from_coin(coin)
# coin_var_percent = calculate_variance_percents_from_eigenvals(eigenvals = coin$eig)
# coin_plot = plot_coinertia(coin_df = coin_df,
#                            variance_percents = coin_var_percent)
# coin_plot
# ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
#                          v_form_complete, ".coinertia_plot.pdf"),
#        plot = coin_plot, width = 8, height = 8, units = "in")

## crap
# coin$aX
# coin$aY
# 
# df1_ax1 = 1
# df1_ax2 = 2
# df2_ax1 = 5
# df2_ax2 = 6
# 
# coin_plot <- ggplot() + 
#   geom_point(data=coin_df, mapping=aes(x=coin_df[,df1_ax1], y=coin_df[,df1_ax2]), size=0.5, fill=coin_df$color) +
#   geom_segment(data=coin_df, mapping=aes(x=coin_df[,df1_ax1], y=coin_df[,df1_ax2], xend=coin_df[,df2_ax1], yend=coin_df[,df2_ax2]), 
#                arrow = arrow(length = unit(0.25,"cm")), size=0.4, color="darkgrey", alpha=0.8) + 
#   geom_point(data=coin_df, mapping=aes(x=coin_df[,df2_ax1], y=coin_df[,df2_ax2]), size=3, shape=21, fill=coin_df$color, alpha=0.6) +
#   theme_bw(base_size = 14, base_family = "Times") +
#   theme(panel.background = element_blank(),
#         # panel.grid = element_blank(),
#         plot.background = element_blank(),
#         #axis.text.x = element_blank(),
#         #axis.ticks.x=element_blank()
#   ) +
#   # xlim(min = min_lim, max = max_lim) +
#   # ylim(min = min_lim, max = max_lim) +
#   xlab(paste0("co-inertia axis ", x, " - ", variance_percents[x], "%")) +
#   ylab(paste0("co-inertia axis ", y, " - ", round(variance_percents[y], 2), "%"))
# 
# coin_plot


# get gt data frames
# sd_candidate_snps = get_gt_df_from_raw(raw_table = sd_candidate_snps_raw)
# # qvals_candidate_snps = get_gt_df_from_raw(raw_table = qvals_candidate_snps_raw)
# neutral_snps = get_gt_df_from_raw(raw_table = neutral_snps_raw)
neutral_snps = neutral_snps[, sample(ncol(neutral_snps), 363)]

# get distance data frames
sd_candidate_dist_df = get_distance_df_from_snps(gt_df = sd_candidate_snps)
sd_candidate_dist_df[[1]]$type = "Candidate"
neutral_dist_df = get_distance_df_from_snps(gt_df = neutral_snps)
neutral_dist_df[[1]]$type = "Neutral"
# pair_dist_df = rbind(sd_candidate_dist_df[[1]], neutral_dist_df[[1]])

# upperline = mean(neutral_dist_df[[1]]$pairwise_dist) + (2 * sd(neutral_dist_df[[1]]$pairwise_dist))
# lowerline = mean(neutral_dist_df[[1]]$pairwise_dist) - (2 * sd(neutral_dist_df[[1]]$pairwise_dist))

# plot
# dist_plot <- ggplot(data = pair_dist_df) +
#   geom_jitter(aes(x = p1, y = pairwise_dist),
#               shape=16, position=position_jitter(0.2),
#               color = pair_dist_df$color) +
#  # geom_hline(yintercept = upperline, linetype = "dashed") + 
#  # geom_hline(yintercept = lowerline, linetype = "dashed") +
#   labs(y = "Pairwise Distances") +
#   scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
#   theme_bw() +
#   theme(panel.background = element_blank(),
#         plot.background = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_text(size = 12),
#         axis.text.x = element_text(size = 9)) +
#         #axis.text.x = element_text(size = 9, angle = 45, vjust = 0.5, hjust=1)) +
#   facet_wrap(~ type, scales = "free_x")
#   
# dist_plot
# 
# ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
#                          "dist_plot.pdf"),
#        plot = dist_plot, width = 180, height = 70, units = "mm")

## which LATVIA are most distant in PCA
# sd_candidate_pca_df.latvia = sd_candidate_pca_df[grep("la", rownames(sd_candidate_pca_df)),]
# rownames(sd_candidate_pca_df.latvia)[which(sd_candidate_pca_df.latvia$Dim.1 > 32)]
# "c_ll_la_0048" "c_ll_la_0052" "c_ll_la_0053"

df <- sd_candidate_dist_df[[1]]
overall <- data.frame(p1 = "Overall", p2 = df$p2,
                      pairwise_dist = df$pairwise_dist,
                      color = df$color, type = df$type)
df1 <- rbind(df, overall)
df1$fact = factor(df1$p1, levels=c("Balkans", "Caucasus", "Carpathians", "Kirov",
                                   "Latvia", "Norway", "NE-Poland", "Tuva",
                                   "Urals", "Primorsky Krai", "Yakutia", "Overall"))

p <- ggplot(data = df1) +
  geom_density(aes(x = pairwise_dist, fill = p2, alpha = 0.3)) +
  scale_fill_manual(values = unique(df1$color)[c(1,3,2,5,6,4,8,7,9,10,11,12)]) +
  theme_bw() +
  theme(panel.background = element_blank(),
        plot.background = element_blank()) +
  facet_wrap(~ fact, scales = "free_y")
ggsave(filename = paste0("4-Downstream_Analyses/plots/",
                         "adaptive_dista_distri_panels.pdf"),
       plot = p, width = 205, height = 100, units = "mm")

################################################################################



df2 <- neutral_dist_df[[1]]
ggplot(data = df2) +
  geom_density(aes(x = df2$pairwise_dist, fill = df2$p2, alpha = 0.3)) +
  scale_fill_manual(values = unique(df2$color)[c(1,3,2,5,6,4,8,7,9,10,11,12)]) +
  facet_wrap(~ p1, scales = "free_y")

df3 <- data.frame(p1 = df$p1, p2 = df$p2,
                  delta_dist = df$pairwise_dist - df2$pairwise_dist,
                  color = df$color) 
overall <- data.frame(p1 = "Overall", p2 = df$p2,
                      delta_dist = df$pairwise_dist - df2$pairwise_dist,
                      color = df$color)
df3 <- rbind(df3, overall)

df3$fact = factor(df3$p1, levels=c("Balkans", "Caucasus", "Carpathians", "Kirov",
                                   "Latvia", "Norway", "NE-Poland", "Tuva",
                                   "Urals", "Primorsky Krai", "Yakutia", "Overall"))

dist_plot <- ggplot(data = df3) +
  geom_density(aes(x = delta_dist, fill = p2, alpha = 0.3)) +
  scale_fill_manual(values = unique(df2$color)[c(1,3,2,5,6,4,8,7,9,10,11,12)]) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(panel.background = element_blank(),
        plot.background = element_blank()) +
  facet_wrap(~ fact, scales = "free_y")
dist_plot
ggsave(filename = paste0("4-Downstream_Analyses/plots/", 
                         "dist_plot.pdf"),
       plot = dist_plot, width = 205, height = 100, units = "mm")

## if i want to add some sort of population mean differentiation and delta
## between populations I have some work done in the [[2]] dist_dfs

