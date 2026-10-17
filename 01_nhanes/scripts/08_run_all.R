scripts <- c("01_build_core_dataset.R","02_build_covariates.R","03_finalize_analysis_dataset.R","04_weighted_quartiles.R","05_primary_regression.R","06_sensitivity_analyses.R","07_validate_expected_results.R")
for (s in scripts) {
  cat("\n=== Running ",s," ===\n",sep="")
  source(file.path("01_nhanes","scripts",s),encoding="UTF-8")
}
cat("\nNHANES reconstructed pipeline completed.\n")
