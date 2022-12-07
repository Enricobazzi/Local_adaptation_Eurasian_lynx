# compare the SNPs in the different tables of candidates
source("3-Identifying_Candidate_Loci/code/rda_functions.R", print.eval=TRUE)
library(venn)

args = commandArgs(trailingOnly=TRUE)
formula_file_line_1 = args[1]
formula_file_line_2 = args[2]
formula_file_line_3 = args[3]

formula_file = "3-Identifying_Candidate_Loci/tables/rda_exploration_formulas.txt"


# variables and conditions to build RDA formula
v_form_1 = read_formula_from_file(formula_file = formula_file, N = formula_file_line_1)
vars_list_1 = strsplit(v_form_1[1], split = "+", fixed = T)[[1]]
cond_list_1 = strsplit(v_form_1[2], split = "+", fixed = T)[[1]]

v_form_2 = read_formula_from_file(formula_file = formula_file, N = formula_file_line_2)
vars_list_2 = strsplit(v_form_2[1], split = "+", fixed = T)[[1]]
cond_list_2 = strsplit(v_form_2[2], split = "+", fixed = T)[[1]]

v_form_3 = read_formula_from_file(formula_file = formula_file, N = formula_file_line_3)
vars_list_3 = strsplit(v_form_3[1], split = "+", fixed = T)[[1]]
cond_list_3 = strsplit(v_form_3[2], split = "+", fixed = T)[[1]]

# prefix of tables
table_prefix_1 = paste0("3-Identifying_Candidate_Loci/tables/",
                      "vars_", str_replace_all(paste(vars_list_1, collapse = "-"),
                                               "_", ""), ".",
                      "cond_", paste(cond_list_1, collapse = "-"), ".")

table_prefix_2 = paste0("3-Identifying_Candidate_Loci/tables/",
                      "vars_", str_replace_all(paste(vars_list_2, collapse = "-"),
                                               "_", ""), ".",
                      "cond_", paste(cond_list_2, collapse = "-"), ".")

table_prefix_3 = paste0("3-Identifying_Candidate_Loci/tables/",
                        "vars_", str_replace_all(paste(vars_list_3, collapse = "-"),
                                                 "_", ""), ".",
                        "cond_", paste(cond_list_3, collapse = "-"), ".")

# candidate table 1
if ("PC1" %in% cond_list_1){
  K = 2
} else {
  K = 3
}
table_1_qvals = read.table(paste0(table_prefix_1,
                  deparse(substitute(cand_qvals)), ".",
                  "K", K, ".", "candidates.tsv"), header = T)

table_1_sd = read.table(paste0(table_prefix_1,
                                  deparse(substitute(cand_sd)), ".",
                                  "K", K, ".", "candidates.tsv"), header = T)

# candidate table 2
if ("PC1" %in% cond_list_2){
  K = 2
} else {
  K = 3
}
table_2_qvals = read.table(paste0(table_prefix_2,
                                  deparse(substitute(cand_qvals)), ".",
                                  "K", K, ".", "candidates.tsv"), header = T)

table_2_sd = read.table(paste0(table_prefix_2,
                               deparse(substitute(cand_sd)), ".",
                               "K", K, ".", "candidates.tsv"), header = T)

# candidate table 3
if ("PC1" %in% cond_list_3){
  K = 2
} else {
  K = 3
}
table_3_qvals = read.table(paste0(table_prefix_3,
                                 deparse(substitute(cand_qvals)), ".",
                                 "K", K, ".", "candidates.tsv"), header = T)

table_3_sd = read.table(paste0(table_prefix_3,
                               deparse(substitute(cand_sd)), ".",
                               "K", K, ".", "candidates.tsv"), header = T)

name_set_1 = paste0("(", paste(v_form_1, collapse = ")-("), ")")
name_set_1 = get_name_set_abbreviation(name_set_1)
name_set_2 = paste0("(", paste(v_form_2, collapse = ")-("), ")")
name_set_2 = get_name_set_abbreviation(name_set_2)
name_set_3 = paste0("(", paste(v_form_3, collapse = ")-("), ")")
name_set_3 = get_name_set_abbreviation(name_set_3)

# venn diagram of sd candidates
pdf(file = paste0("3-Identifying_Candidate_Loci/plots/",
                 "venn_diagrams.",
                 "line_", formula_file_line_1, ".",
                 "line_", formula_file_line_2, ".",
                 "line_", formula_file_line_3, ".",
                 "sd_candidates.pdf"),
    width = 8,
    height = 8)

venn(list(table_1_sd$snp, table_2_sd$snp, table_3_sd$snp),
     snames = c(name_set_1, name_set_2, name_set_3), 
     zcolor = 'style',
     ilcs = 1, sncs = 0.8,
     box = F)

dev.off()

# venn diagram of qval candidates
pdf(file = paste0("3-Identifying_Candidate_Loci/plots/",
                  "venn_diagrams.",
                  "line_", formula_file_line_1, ".",
                  "line_", formula_file_line_2, ".",
                  "line_", formula_file_line_3, ".",
                  "qval_candidates.pdf"),
    width = 8,
    height = 8)

venn(list(table_1_qvals$snp, table_2_qvals$snp, table_3_qvals$snp),
     snames = c(name_set_1, name_set_2, name_set_3), 
     zcolor = 'style',
     ilcs = 1, sncs = 0.8,
     box = F)

dev.off()

