source("4-Downstream_Analyses/code/downstream_analyses_functions.R", print.eval=TRUE)

sd_candidate_snps_raw = "4-Downstream_Analyses/tables/sd_candidate_snps.raw"
qvals_candidate_snps_raw = "4-Downstream_Analyses/tables/qvals_candidate_snps.raw"
neutral_snps_raw = "4-Downstream_Analyses/tables/neutral_snps.raw"

sd_candidate_snps = get_gt_df_from_raw(raw_table = sd_candidate_snps_raw)
qvals_candidate_snps = get_gt_df_from_raw(raw_table = qvals_candidate_snps_raw)
neutral_snps = get_gt_df_from_raw(raw_table = neutral_snps_raw)

sd_candidate_pca <- FactoMineR::PCA(sd_candidate_snps, graph = F)
qvals_candidate_pca <- FactoMineR::PCA(qvals_candidate_snps, graph = F)
neutral_pca <- FactoMineR::PCA(neutral_snps, graph = F)

sd_candidate_pca_df = build_df_from_pca(pca_obj = sd_candidate_pca)
qvals_candidate_pca_df = build_df_from_pca(pca_obj = qvals_candidate_pca)
neutral_pca_df = build_df_from_pca(pca_obj = neutral_pca)

sd_candidate_variance_percents <- sd_candidate_pca$eig[,2]
qvals_candidate_variance_percents <- qvals_candidate_pca$eig[,2]
neutral_variance_percents <- neutral_pca$eig[,2]

sd_pca_plot = plot_pca(pca_df = sd_candidate_pca_df,
                       variance_percents = sd_candidate_variance_percents)
sd_pca_plot

sd_pca_plot_1_3 = plot_pca(pca_df = sd_candidate_pca_df,
                       variance_percents = sd_candidate_variance_percents,
                       y = 3)
sd_pca_plot_1_3

qvals_pca_plot = plot_pca(pca_df = qvals_candidate_pca_df,
                       variance_percents = qvals_candidate_variance_percents)
qvals_pca_plot

qvals_pca_plot_1_3 = plot_pca(pca_df = qvals_candidate_pca_df,
                          variance_percents = qvals_candidate_variance_percents,
                          y = 3)
qvals_pca_plot_1_3


neutral_pca_plot = plot_pca(pca_df = neutral_pca_df,
                          variance_percents = neutral_variance_percents)
neutral_pca_plot

neutral_pca_plot_1_3 = plot_pca(pca_df = neutral_pca_df,
                            variance_percents = neutral_variance_percents,
                            y = 3)
neutral_pca_plot_1_3

