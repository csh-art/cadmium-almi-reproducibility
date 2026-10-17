# Portable validator for frozen Figure 5 human-readable source data.
# Run from the repository root, or set CADMIUM_ALMI_REPO_ROOT.
options(stringsAsFactors = FALSE)
find_repo_root <- function() {
  x <- Sys.getenv("CADMIUM_ALMI_REPO_ROOT", unset = "")
  if (nzchar(x) && dir.exists(x)) return(normalizePath(x, winslash = "/"))
  cur <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "05_mitochondrial_pathways", "source_data", "04_Figure5C_correlations.csv"))) return(cur)
    par <- dirname(cur); if (identical(par, cur)) break; cur <- par
  }
  stop("Repository root not found.", call. = FALSE)
}
root <- find_repo_root(); s <- file.path(root, "05_mitochondrial_pathways", "source_data")
rd <- function(x) read.csv(file.path(s, x), check.names = FALSE, stringsAsFactors = FALSE)
group <- rd("01_group_df.csv")
cor <- rd("04_Figure5C_correlations.csv")
gsea4 <- rd("05_Figure5D_target_GSEA_statistics.csv")
full <- rd("14_Figure5A_full_GO_BP_GSEA_result.csv")
hubs <- rd("09_hub_genes.csv")
paths <- rd("11_Figure5C_pathways_after_hub_removal.csv")
stopifnot(nrow(group) == 40L, nrow(cor) == 44L, nrow(gsea4) == 4L, nrow(full) == 5880L, nrow(hubs) == 11L)
# Schema-tolerant hub-removal check: compare every character cell in path membership to hub symbols.
hub_symbols <- unique(unlist(hubs, use.names = FALSE)); hub_symbols <- hub_symbols[grepl("^[A-Z0-9-]+$", hub_symbols)]
path_chars <- unique(unlist(paths, use.names = FALSE)); path_chars <- as.character(path_chars)
stopifnot(length(intersect(hub_symbols, path_chars)) == 0L)
cat("Figure 5 frozen-source validation: PASS\n")
cat("40 samples; 44 correlations; 4 target GSEA pathways; 5880 GO-BP rows; 11 hubs removed from pathway memberships.\n")
