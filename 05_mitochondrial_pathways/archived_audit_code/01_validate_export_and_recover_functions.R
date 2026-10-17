# ============================================================
# M05 repository audit and source-data export
# Read-only with respect to the original analysis objects.
# Outputs are written only to:
#   <M05 project>/02_analysis_files/00_repository_audit/
#
# Purpose:
# 1) Validate the frozen Figure 5 canonical objects.
# 2) Export human-readable source-data CSVs from the frozen core.
# 3) Verify that the no-hub pathway sets exclude all 11 hub genes.
# 4) Recover relevant saved functions from project_workspace.RData
#    without modifying the original workspace.
#
# No packages are installed by this script.
# ============================================================

options(stringsAsFactors = FALSE)

m05_root <- "E:/1Data/Cd_Sarcopenia_NHANES/new_workspace/03_results/M05_PathwayValidation_GSEA_ssGSEA_Correlation"
analysis_dir <- file.path(m05_root, "02_analysis_files")
backup_dir <- file.path(m05_root, "backup")
out_dir <- file.path(analysis_dir, "00_repository_audit")
source_dir <- file.path(out_dir, "source_data")
function_dir <- file.path(out_dir, "recovered_functions")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(function_dir, recursive = TRUE, showWarnings = FALSE)

core_path <- file.path(analysis_dir, "canonical_core", "M05_Figure5_CANONICAL_CORE_FINAL.RData")
bc_path <- file.path(analysis_dir, "canonical_inputs", "M05_BC_GSE111016_log2CPM_group_canonical_input.rds")
gsea_path <- file.path(analysis_dir, "canonical_inputs", "M05A_GO_BP_GSEA_clusterProfiler_object_v2.rds")
workspace_path <- file.path(backup_dir, "project_workspace.RData")

required_files <- c(core_path, bc_path, gsea_path)
if (any(!file.exists(required_files))) {
  stop("Missing required M05 file(s):\n", paste(required_files[!file.exists(required_files)], collapse="\n"), call.=FALSE)
}

# ---------- helpers ----------
obj_summary <- function(name, x) {
  dims <- dim(x)
  data.frame(
    Object = name,
    Class = paste(class(x), collapse = " | "),
    Length = length(x),
    Nrow = if (length(dims) >= 1) dims[1] else NA_integer_,
    Ncol = if (length(dims) >= 2) dims[2] else NA_integer_,
    Size_MB = round(as.numeric(object.size(x))/1024^2, 4),
    stringsAsFactors = FALSE
  )
}

write_matrix_csv <- function(x, path, row_name_col="Row") {
  x <- as.matrix(x)
  out <- data.frame(stringsAsFactors = FALSE)
  out[[row_name_col]] <- rownames(x)
  out <- cbind(out, as.data.frame(x, check.names = FALSE, stringsAsFactors = FALSE))
  utils::write.csv(out, path, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

list_to_long <- function(x, list_name="Pathway", value_name="Gene") {
  if (!is.list(x)) return(data.frame())
  pieces <- lapply(names(x), function(nm) {
    vals <- as.character(x[[nm]])
    data.frame(Pathway = nm, Gene = vals, stringsAsFactors = FALSE)
  })
  if (!length(pieces)) return(data.frame())
  do.call(rbind, pieces)
}

safe_deparse_function <- function(fn, name, dir_path) {
  code <- deparse(fn, control = "all", width.cutoff = 500L)
  if (!length(code)) return(NA_character_)
  code[1] <- paste0(name, " <- ", code[1])
  path <- file.path(dir_path, paste0(name, ".R"))
  writeLines(code, path, useBytes = TRUE)
  path
}

# ---------- load frozen core ----------
core_env <- new.env(parent = globalenv())
loaded_objects <- load(core_path, envir = core_env)
required_core <- c(
  "group_df", "hub_genes", "strict_hub_genes", "ssgsea_scores",
  "figure5b_statistics", "figure5c_correlations", "target_gsea_statistics",
  "figure5d_es", "figure5d_hits", "rank_vec"
)
missing_core <- required_core[!vapply(required_core, exists, logical(1), envir=core_env, inherits=FALSE)]

inventory <- do.call(rbind, lapply(loaded_objects, function(nm) obj_summary(nm, get(nm, envir=core_env, inherits=FALSE))))
utils::write.csv(inventory, file.path(out_dir, "00_M05_core_object_inventory.csv"), row.names=FALSE, fileEncoding="UTF-8")

if (length(missing_core)) {
  stop("Frozen core lacks required object(s): ", paste(missing_core, collapse=", "), call.=FALSE)
}

# ---------- validate locked anchors ----------
group_df <- as.data.frame(get("group_df", envir=core_env), stringsAsFactors=FALSE)
hub_genes <- as.character(get("hub_genes", envir=core_env))
strict_hub_genes <- as.character(get("strict_hub_genes", envir=core_env))
ssgsea_scores <- as.matrix(get("ssgsea_scores", envir=core_env))
figure5b_statistics <- as.data.frame(get("figure5b_statistics", envir=core_env), stringsAsFactors=FALSE)
figure5c_correlations <- as.data.frame(get("figure5c_correlations", envir=core_env), stringsAsFactors=FALSE)
target_gsea_statistics <- as.data.frame(get("target_gsea_statistics", envir=core_env), stringsAsFactors=FALSE)
figure5d_es <- as.data.frame(get("figure5d_es", envir=core_env), stringsAsFactors=FALSE)
figure5d_hits <- as.data.frame(get("figure5d_hits", envir=core_env), stringsAsFactors=FALSE)
rank_vec <- get("rank_vec", envir=core_env)

expected_pathways <- c(
  "Oxidative phosphorylation",
  "Respiratory electron transport chain",
  "Mitochondrial ATP synthesis",
  "Mitochondrial translation"
)

checks <- data.frame(
  Check = c(
    "group_df has 40 samples",
    "Control n = 20",
    "Sarcopenia n = 20",
    "ssGSEA has 40 sample columns",
    "ssGSEA contains four locked mitochondrial pathways",
    "hub_genes length = 11",
    "strict_hub_genes length = 2",
    "strict hubs are ATF4 and CASP9",
    "Figure 5C has 44 gene-pathway correlations",
    "target GSEA statistics contain four pathways",
    "pathways_nohub object exists",
    "ssgsea_scores_nohub object exists"
  ),
  Pass = c(
    nrow(group_df) == 40L,
    sum(as.character(group_df$Group) == "Control") == 20L,
    sum(as.character(group_df$Group) == "Sarcopenia") == 20L,
    ncol(ssgsea_scores) == 40L,
    all(expected_pathways %in% rownames(ssgsea_scores)),
    length(hub_genes) == 11L,
    length(strict_hub_genes) == 2L,
    setequal(strict_hub_genes, c("ATF4", "CASP9")),
    nrow(figure5c_correlations) == 44L,
    nrow(target_gsea_statistics) == 4L,
    exists("pathways_nohub", envir=core_env, inherits=FALSE),
    exists("ssgsea_scores_nohub", envir=core_env, inherits=FALSE)
  ),
  stringsAsFactors = FALSE
)

# ---------- export source data ----------
utils::write.csv(group_df, file.path(source_dir, "01_group_df.csv"), row.names=FALSE, fileEncoding="UTF-8")
write_matrix_csv(ssgsea_scores, file.path(source_dir, "02_Figure5B_ssGSEA_scores.csv"), "Pathway")
utils::write.csv(figure5b_statistics, file.path(source_dir, "03_Figure5B_statistics.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(figure5c_correlations, file.path(source_dir, "04_Figure5C_correlations.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(target_gsea_statistics, file.path(source_dir, "05_Figure5D_target_GSEA_statistics.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(figure5d_es, file.path(source_dir, "06_Figure5D_running_ES.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(figure5d_hits, file.path(source_dir, "07_Figure5D_gene_set_hits.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(data.frame(Gene=names(rank_vec), Rank_metric=as.numeric(rank_vec), stringsAsFactors=FALSE), file.path(source_dir, "08_Figure5D_rank_vector.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
utils::write.csv(data.frame(Gene=hub_genes, Hub_class=ifelse(hub_genes %in% strict_hub_genes, "Strict", "Robust"), stringsAsFactors=FALSE), file.path(source_dir, "09_hub_genes.csv"), row.names=FALSE, fileEncoding="UTF-8")

optional_df_objects <- c("pathway_size_table", "figure5a_display")
for (nm in optional_df_objects) {
  if (exists(nm, envir=core_env, inherits=FALSE)) {
    obj <- as.data.frame(get(nm, envir=core_env), stringsAsFactors=FALSE)
    utils::write.csv(obj, file.path(source_dir, paste0("10_", nm, ".csv")), row.names=FALSE, na="", fileEncoding="UTF-8")
  }
}

hub_overlap_pass <- NA
hub_overlap_rows <- data.frame()
if (exists("pathways_nohub", envir=core_env, inherits=FALSE)) {
  pathways_nohub <- get("pathways_nohub", envir=core_env)
  pathways_nohub_long <- list_to_long(pathways_nohub)
  if (nrow(pathways_nohub_long)) {
    utils::write.csv(pathways_nohub_long, file.path(source_dir, "11_Figure5C_pathways_after_hub_removal.csv"), row.names=FALSE, fileEncoding="UTF-8")
    hub_overlap_rows <- pathways_nohub_long[pathways_nohub_long$Gene %in% hub_genes, , drop=FALSE]
    hub_overlap_pass <- nrow(hub_overlap_rows) == 0L
  }
}

if (exists("ssgsea_scores_nohub", envir=core_env, inherits=FALSE)) {
  write_matrix_csv(get("ssgsea_scores_nohub", envir=core_env), file.path(source_dir, "12_Figure5C_ssGSEA_scores_after_hub_removal.csv"), "Pathway")
}

if (exists("figure5a_display", envir=core_env, inherits=FALSE)) {
  utils::write.csv(as.data.frame(get("figure5a_display", envir=core_env), stringsAsFactors=FALSE), file.path(source_dir, "13_Figure5A_display_terms_from_core.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
}

# ---------- export full stored GO-BP GSEA result ----------
gsea_export_status <- "NOT_ATTEMPTED"
if (requireNamespace("clusterProfiler", quietly=TRUE)) {
  gsea_obj <- readRDS(gsea_path)
  if (methods::is(gsea_obj, "gseaResult") && "result" %in% methods::slotNames(gsea_obj)) {
    gsea_df <- as.data.frame(methods::slot(gsea_obj, "result"), stringsAsFactors=FALSE)
    utils::write.csv(gsea_df, file.path(source_dir, "14_Figure5A_full_GO_BP_GSEA_result.csv"), row.names=FALSE, na="", fileEncoding="UTF-8")
    gsea_export_status <- paste0("PASS: ", nrow(gsea_df), " rows exported")
  } else {
    gsea_export_status <- paste0("WARN: unexpected GSEA object class = ", paste(class(gsea_obj), collapse=" | "))
  }
} else {
  gsea_export_status <- "SKIPPED: clusterProfiler is not installed in the current R library"
}

# ---------- inspect canonical B/C input ----------
bc_input <- readRDS(bc_path)
bc_names <- names(bc_input)
bc_inventory <- do.call(rbind, lapply(bc_names, function(nm) obj_summary(nm, bc_input[[nm]])))
utils::write.csv(bc_inventory, file.path(out_dir, "01_M05_BC_input_inventory.csv"), row.names=FALSE, fileEncoding="UTF-8")

# ---------- recover relevant functions from workspace backup ----------
function_manifest <- data.frame(Function=character(), Size_MB=numeric(), Output_file=character(), stringsAsFactors=FALSE)
workspace_status <- "workspace backup not found"
if (file.exists(workspace_path)) {
  ws_env <- new.env(parent = globalenv())
  ws_names <- load(workspace_path, envir=ws_env)
  fn_names <- ws_names[vapply(ws_names, function(nm) is.function(get(nm, envir=ws_env, inherits=FALSE)), logical(1))]
  relevant_fn <- fn_names[grepl("ssgsea|figure5|gsea|pathway|canonical|core", fn_names, ignore.case=TRUE)]
  if (length(relevant_fn)) {
    rows <- lapply(relevant_fn, function(nm) {
      fn <- get(nm, envir=ws_env, inherits=FALSE)
      path <- safe_deparse_function(fn, nm, function_dir)
      data.frame(Function=nm, Size_MB=round(as.numeric(object.size(fn))/1024^2,4), Output_file=path, stringsAsFactors=FALSE)
    })
    function_manifest <- do.call(rbind, rows)
    workspace_status <- paste0("PASS: recovered ", nrow(function_manifest), " relevant function(s)")
  } else {
    workspace_status <- "WARN: no relevant functions matched the recovery pattern"
  }
}
utils::write.csv(function_manifest, file.path(out_dir, "02_M05_recovered_function_manifest.csv"), row.names=FALSE, fileEncoding="UTF-8")

# ---------- final report ----------
checks <- rbind(checks, data.frame(Check="No hub gene remains in pathways_nohub", Pass=isTRUE(hub_overlap_pass), stringsAsFactors=FALSE))
utils::write.csv(checks, file.path(out_dir, "03_M05_validation_checks.csv"), row.names=FALSE, fileEncoding="UTF-8")

historical_code <- file.path(m05_root, "03_outputs_figures_tables", "figure_re", "M05_BC_GSE111016_canonical_reconstruction_code.txt")
historical_code_size <- if (file.exists(historical_code)) file.info(historical_code)$size else NA_real_

report <- c(
  "============================================================",
  "M05 FIGURE 5 REPOSITORY AUDIT",
  "============================================================",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("Core file: ", core_path),
  paste0("B/C canonical input: ", bc_path),
  paste0("GO-BP GSEA object: ", gsea_path),
  "",
  paste0("Frozen core objects: ", length(loaded_objects)),
  paste0("Missing required core objects: ", ifelse(length(missing_core), paste(missing_core, collapse=", "), "none")),
  paste0("GSEA full-result export: ", gsea_export_status),
  paste0("Workspace function recovery: ", workspace_status),
  paste0("Historical canonical reconstruction-code file size: ", ifelse(is.na(historical_code_size), "not found", paste0(historical_code_size, " bytes"))),
  "",
  "Validation checks:",
  paste0(ifelse(checks$Pass, "PASS", "FAIL"), " | ", checks$Check),
  "",
  paste0("Export folder: ", out_dir),
  "",
  "Please upload these files/folders to ChatGPT after the run:",
  "- 00_M05_core_object_inventory.csv",
  "- 01_M05_BC_input_inventory.csv",
  "- 02_M05_recovered_function_manifest.csv",
  "- 03_M05_validation_checks.csv",
  "- 04_M05_validation_report.txt",
  "- recovered_functions/ (especially run_ssgsea.R if present)",
  "- source_data/ (ZIP the folder; these are small CSV exports)",
  "============================================================"
)
writeLines(report, file.path(out_dir, "04_M05_validation_report.txt"), useBytes=TRUE)

cat(paste(report, collapse="\n"), "\n")

invisible(list(
  core_inventory=inventory,
  bc_inventory=bc_inventory,
  checks=checks,
  recovered_functions=function_manifest,
  output_dir=out_dir
))
