cat("R script started\n")
flush.console()

library(arrow)

file <- "E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M11_MR/02_analysis_files/Muscle_Skeletal.v11.eQTLs.signif_pairs.parquet"

cat("File exists:", file.exists(file), "\n")
cat("File size MB:", round(file.info(file)$size / 1024^2, 2), "\n")
flush.console()

ds <- open_dataset(file, format = "parquet")

cat("\nDATASET INFORMATION:\n")
print(ds)

cat("\nCOLUMN NAMES:\n")
print(names(ds))

cat("\nScript finished\n")
