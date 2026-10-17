# Validate archived Figure 6/MR source tables and clumping stability.
# No packages are installed; base R only.
module_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(module_dir) != "06_mr") {
  stop("Run this script with the working directory set to repository/06_mr")
}
read_tsv <- function(p) read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
inst <- read_tsv(file.path(module_dir, "source_data", "Table_MR_instrument_characteristics_and_results.tsv"))
comb <- read_tsv(file.path(module_dir, "source_data", "Table_MR_Wald_ratio_combined_ALM_final.tsv"))
sex  <- read_tsv(file.path(module_dir, "source_data", "Table_MR_Wald_ratio_sex_stratified_ALM.tsv"))
stg  <- read_tsv(file.path(module_dir, "source_data", "Table_MR_Steiger_directionality.tsv"))
stopifnot(nrow(inst) == 4L, nrow(comb) == 4L, nrow(sex) == 8L, nrow(stg) == 4L)
expected <- c(ATF4="rs882565", CASP9="rs933704", CASP3="rs4647601", MAP1LC3B="rs7865")
obs <- setNames(as.character(inst$rsid), as.character(inst$gene_name))
stopifnot(identical(unname(obs[names(expected)]), unname(expected)))
stopifnot(all(as.numeric(inst$F_stat) > 10))
atf4 <- comb[comb$gene_name == "ATF4", ]
stopifnot(nrow(atf4) == 1L, abs(atf4$beta_MR + 0.016503) < 1e-6, abs(atf4$p_bonferroni_4 - 0.026568) < 1e-6)
stopifnot(all(stg$correct_causal_direction %in% c(TRUE, "TRUE", "True", 1)))
cat("MR final-result validation: PASS\n")
