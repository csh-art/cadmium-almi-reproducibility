# Reconstructed M11 frozen-chain validator
# This is NOT recovered original upstream analysis code.
# Run from the repository root:
#   source("06_mr/scripts/02_validate_frozen_M11_chain.R", encoding="UTF-8")

options(stringsAsFactors = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(root, "06_mr"))) {
  stop("Run this script from the repository root.", call. = FALSE)
}

rd <- function(path) {
  read.delim(file.path(root, path), stringsAsFactors = FALSE, check.names = FALSE)
}

mr05 <- rd("06_mr/provenance/harmonization/MR_05_strict_candidate_instruments.tsv")
mr06 <- rd("06_mr/provenance/harmonization/MR_06_combined_ALM_exact_matches.tsv")
mr06b <- rd("06_mr/provenance/harmonization/MR_06b_combined_ALM_missing_variants.tsv")
mr06c <- rd("06_mr/provenance/harmonization/MR_06c_combined_ALM_position_recovered.tsv")
mr07 <- rd("06_mr/provenance/harmonization/MR_07_combined_ALM_harmonised_candidates.tsv")
mr09 <- rd("06_mr/provenance/clumping/MR_09_PLINK_clump_mapping.tsv")
mr10 <- rd("06_mr/source_data/MR_10_final_independent_instruments_combined_ALM.tsv")
combined <- rd("06_mr/source_data/Table_MR_Wald_ratio_combined_ALM_final.tsv")
sex <- rd("06_mr/source_data/Table_MR_Wald_ratio_sex_stratified_ALM.tsv")
heterogeneity <- rd("06_mr/source_data/Table_MR_female_male_effect_difference_final.tsv")
ld_sens <- rd("06_mr/source_data/Table_MR_LD_clumping_sensitivity.tsv")
steiger <- rd("06_mr/source_data/Table_MR_Steiger_directionality.tsv")

expected_genes <- c("ATF4","CASP9","CASP3","MAP1LC3B")
expected_rsids <- c("rs882565","rs933704","rs4647601","rs7865")

stopifnot(
  nrow(mr05) == 342L,
  nrow(mr06) == 327L,
  nrow(mr06b) == 15L,
  nrow(mr06c) == 6L,
  nrow(mr07) == 333L,
  nrow(mr09) == 309L,
  nrow(mr10) == 4L,
  setequal(mr10$gene_name, expected_genes),
  setequal(mr10$rsid, expected_rsids)
)

wald <- transform(
  mr10,
  beta_MR_recalc = beta_outcome / beta_exposure,
  se_MR_recalc = se_outcome / abs(beta_exposure),
  F_recalc = (beta_exposure / se_exposure)^2
)
wald <- merge(
  wald,
  combined[, c("gene_name","rsid","beta_MR","se_MR","F_stat")],
  by = c("gene_name","rsid"),
  suffixes = c("_instrument","_final"),
  sort = FALSE
)

stopifnot(
  max(abs(wald$beta_MR_recalc - wald$beta_MR)) < 1e-6,
  max(abs(wald$se_MR_recalc - wald$se_MR)) < 1e-6,
  max(abs(wald$F_recalc - wald$F_stat_final)) < 5e-3
)

by_thr <- split(ld_sens, ld_sens$r2_threshold)
stopifnot(
  all(c("0.001","0.01") %in% names(by_thr)),
  setequal(by_thr[["0.001"]]$rsid, expected_rsids),
  setequal(by_thr[["0.01"]]$rsid, expected_rsids)
)

stopifnot(
  nrow(steiger) == 4L,
  all(as.logical(steiger$correct_causal_direction))
)

cat("M11 frozen-chain validation: PASS\n")
cat("342 strict candidates -> 333 harmonized -> 309 EUR-present -> 4 retained instruments\n")
cat("Primary and r2=0.01 sensitivity clumping retain the same four instruments.\n")
