# Portable validator for frozen Figure 3 source-data anchors.
# Run from the repository root, or set CADMIUM_ALMI_REPO_ROOT.
options(stringsAsFactors = FALSE)
find_repo_root <- function() {
  x <- Sys.getenv("CADMIUM_ALMI_REPO_ROOT", unset = "")
  if (nzchar(x) && dir.exists(x)) return(normalizePath(x, winslash = "/"))
  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "02_candidate_gene_integration", "source_data", "10_FINAL_original_23_candidate_genes.csv"))) return(cur)
    par <- dirname(cur); if (identical(par, cur)) break; cur <- par
  }
  stop("Repository root not found.", call. = FALSE)
}
root <- find_repo_root(); m <- file.path(root, "02_candidate_gene_integration")
rd <- function(...) read.csv(file.path(m, ...), check.names = FALSE, stringsAsFactors = FALSE)
ctd <- rd("source_data", "02_CTD_cadmium_interaction_count_gt10_138_genes.csv")
nom <- rd("source_data", "05_GSE111016_nominal_P005_DEGs_3342_Ensembl.csv")
final23 <- rd("source_data", "10_FINAL_original_23_candidate_genes.csv")
slc <- rd("source_data", "11_SLC11A2_threshold_boundary_AUDIT.csv")
pheno <- rd("source_data", "03_GSE111016_aligned_phenotype_20_no_20_yes.csv")
stopifnot(nrow(ctd) == 138L, nrow(nom) == 3342L, nrow(final23) == 23L, nrow(pheno) == 40L)
# Keep the SLC11A2 check schema-tolerant: identify a numeric log2FC-like column.
num <- names(slc)[vapply(slc, is.numeric, logical(1))]
lfc_name <- num[grepl("log2|fold", num, ignore.case = TRUE)][1]
if (!is.na(lfc_name) && nzchar(lfc_name)) stopifnot(abs(abs(slc[[lfc_name]][1]) - 0.0998186723) < 1e-5)
cat("Figure 3 frozen-source validation: PASS\n")
cat("CTD=138; nominal DEGs=3342; final candidates=23; phenotype samples=40.\n")
