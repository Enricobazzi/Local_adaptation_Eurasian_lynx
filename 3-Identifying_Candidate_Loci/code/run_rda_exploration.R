source("3-Identifying_Candidate_Loci/code/rda_functions.R", print.eval=TRUE)

## ## ## ## ## ## 
## arguments:  ##
## ## ## ## ## ##
args = commandArgs(trailingOnly=TRUE)
formula_file_line = args[1]
raw_file = "1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.raw"
env_table = "2-Preparing_Environmental_Data/tables/allvars_persample_table.tsv"
neutral_raw_file = "1-Preparing_Genetic_Data/tables/finalset.maf5pc.pruned.nogenes.raw"
formula_file = "3-Identifying_Candidate_Loci/tables/rda_exploration_formulas.txt"





# variables and conditions to build RDA formula
v_form = read_formula_from_file(formula_file = formula_file, N = formula_file_line)
vars_list = strsplit(v_form[1], split = "+", fixed = T)[[1]]
cond_list = strsplit(v_form[2], split = "+", fixed = T)[[1]]

# prefix of plots
plots_prefix = paste0("3-Identifying_Candidate_Loci/plots/",
                      "vars_", str_replace_all(paste(vars_list, collapse = "-"),
                                           "_", ""), ".",
                      "cond_", paste(cond_list, collapse = "-"), ".")

# prefix of tables
table_prefix = paste0("3-Identifying_Candidate_Loci/tables/",
                      "vars_", str_replace_all(paste(vars_list, collapse = "-"),
                                               "_", ""), ".",
                      "cond_", paste(cond_list, collapse = "-"), ".")

## ## ## ## ## ## 
## load data:  ##
## ## ## ## ## ## 

# load genetic data
gt_data = read_raw_file(raw_file = raw_file)
gt_data_tsv = get_tsv_from_gt_data(gt_data = gt_data)
# sample list
sample_list = rownames(gt_data_tsv)
# environmental variables
env.predictors = read_env_data(env_table = env_table, sample_list = sample_list)
# load neutral genetic data
neutral_gt_data = read_raw_file(raw_file = neutral_raw_file)
# PCs of neutral genetic structure
PCs = get_pcs_from_ade(gt_data = neutral_gt_data)
# variables dataframe = PCs + variables I want to include in model
variables = get_vars_df(PCs, env.predictors, vars_list = vars_list)

## ## ## ## ## ## ##
## run analysis:  ##
## ## ## ## ## ## ##

# run RDA
RDA = run_rda(gt_data = gt_data, variables = variables,
              vars_list = vars_list, cond_list = cond_list) #######

# plot inertia per axis
plot_ord_x_inert_y(plot_name = paste0(plots_prefix, "inertia_per_axis.pdf"),
                   RDA = RDA)

## ## ## ## ## ## ## ## ## ## ## ##
## get results for different Ks: ##
## ## ## ## ## ## ## ## ## ## ## ##

# axes to analyze for candidates
for (K in 2:3){
  
  ## candidate table using using q values: ##
  
  # calculate p and q values
  p_q_vals = rdadapt(RDA = RDA, K = K)
  # candidates from qvalues
  cand_qvals = get_cand_from_rdadapt(RDA = RDA, p_q_vals = p_q_vals, K = K, limit = 0.1)
  # stop if no outliers are found
  if (NROW(cand_qvals) == 0){
    cat("no outliers found",
        file = paste0(table_prefix,
                      deparse(substitute(cand_qvals)), ".",
                      "K", K, ".", "candidates.tsv"),
        sep="\n")
  } else {
    # correlation of candidates with env variables
    corrs = get_corr_snps_vars(cand = cand_qvals,
                               variables = variables, vars_list = vars_list)
    # paste together and save table
    cand_qvals <- cbind.data.frame(cand_qvals, corrs[,c(ncol(corrs)-1,ncol(corrs))])
    write.table(x = cand_qvals, file = paste0(table_prefix,
                                              deparse(substitute(cand_qvals)), ".",
                                              "K", K, ".", "candidates.tsv"),
                row.names = F, sep = "\t", quote = F)
  }  
  
  ## candidate table using using SD outliers: ##

  # candidates from sd outliers
  cand_sd = get_cand_from_sd(RDA = RDA, K = K, limit = 3)
  # correlation of candidates with env variables
  corrs_sd = get_corr_snps_vars(cand = cand_sd,
                                variables = variables, vars_list = vars_list)
  # paste together and save table
  cand_sd <- cbind.data.frame(cand_sd, corrs_sd[,c(ncol(corrs_sd)-1,ncol(corrs_sd))])
  write.table(x = cand_sd, file = paste0(table_prefix,
                                         deparse(substitute(cand_sd)), ".",
                                         "K", K, ".","candidates.tsv"),
              row.names = F, sep = "\t")
  
  ## plotting RDA loadings as manhattan plots
  snps_df = prepare_snps_dataframe_for_plotting(RDA = RDA, K = K, 
                                      cand_qvals = cand_qvals,
                                      cand_sd = cand_sd)
  # generate manhattan plots
  plot_rda_manhattan(snps_df = snps_df, K = K,
                     plots_prefix = plots_prefix)

  ## plotting RDAs
  if (K == 2) {
    if (NROW(cand_qvals) > 0){
      # qval plot
      plot_rda_snps_and_samples(RDA = RDA, cand = cand_qvals, 
                                plot_name = paste0(plots_prefix, 
                                                   deparse(substitute(cand_qvals)),
                                                   ".K2.RDA1-RDA2.rda.pdf"))
    }
    # SD plot
    plot_rda_snps_and_samples(RDA = RDA, cand = cand_sd,
                              plot_name = paste0(plots_prefix, 
                                                 deparse(substitute(cand_sd)),
                                                 ".K2.RDA1-RDA2.rda.pdf"))
  } else if (K == 3){
    if (NROW(cand_qvals) > 0){
      # qval plot RDA1 RDA2
      plot_rda_snps_and_samples(RDA = RDA, cand = cand_qvals, 
                                plot_name = paste0(plots_prefix, 
                                                   deparse(substitute(cand_qvals)),
                                                   ".K3.RDA1-RDA2.rda.pdf"))
      # qval plot RDA1 RDA3
      plot_rda_snps_and_samples(RDA = RDA, cand = cand_qvals, x = 1, y = 3,
                                plot_name = paste0(plots_prefix, 
                                                   deparse(substitute(cand_qvals)),
                                                   ".K3.RDA1-RDA3.rda.pdf"))
    }  
    # SD plot RDA1 RDA2
    plot_rda_snps_and_samples(RDA = RDA, cand = cand_sd,
                              plot_name = paste0(plots_prefix, 
                                                 deparse(substitute(cand_sd)),
                                                 ".K3.RDA1-RDA2.rda.pdf"))
    # SD plot RDA1 RDA3
    plot_rda_snps_and_samples(RDA = RDA, cand = cand_sd, x = 1, y = 3,
                              plot_name = paste0(plots_prefix, 
                                                 deparse(substitute(cand_sd)),
                                                 ".K3.RDA1-RDA3.rda.pdf"))
  }
}
