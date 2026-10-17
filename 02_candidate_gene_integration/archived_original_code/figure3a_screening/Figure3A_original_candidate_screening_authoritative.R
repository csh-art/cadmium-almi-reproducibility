# Authoritative original Figure 3A candidate screening
# Complete standalone reproducible R code

rm(list = ls())
graphics.off()

run_original_figure3A_screening <- function () 
{
    required_packages <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db")
    missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing_packages) > 0) {
        stop("以下R包尚未安装：\n", paste(missing_packages, collapse = ", "), "\n\n请使用BiocManager安装缺失的Bioconductor包。")
    }
    project_dir <- paste0("E:/1Data/Cd_Sarcopenia_NHANES/", "new_workspace/03_results/", "M03_CandidateGenes_GO_PPI_GSE111016")
    raw_dir <- file.path(project_dir, "01_raw_data")
    analysis_root <- file.path(project_dir, "02_analysis_files")
    output_dir <- file.path(analysis_root, "Figure3A_original_screening_authoritative")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(raw_dir)) {
        stop("找不到原始数据目录：\n", raw_dir)
    }
    if (!dir.exists(analysis_root)) {
        stop("找不到分析目录：\n", analysis_root)
    }
    locate_input_file <- function(preferred_paths, target_basename, search_roots) {
        preferred_paths <- unique(preferred_paths[!is.na(preferred_paths) & preferred_paths != ""])
        preferred_existing <- preferred_paths[file.exists(preferred_paths)]
        if (length(preferred_existing) >= 1L) {
            return(normalizePath(preferred_existing[1], winslash = "/", mustWork = TRUE))
        }
        all_files <- unlist(lapply(search_roots, function(search_root) {
            if (!dir.exists(search_root)) {
                return(character(0))
            }
            list.files(search_root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
        }), use.names = FALSE)
        matches <- unique(all_files[basename(all_files) == target_basename])
        if (length(matches) == 0L) {
            stop("找不到输入文件：\n", target_basename, "\n\n搜索目录：\n", paste(search_roots, collapse = "\n"))
        }
        if (length(matches) > 1L) {
            legacy_matches <- matches[grepl("00_legacy_reference_inputs", matches, fixed = TRUE)]
            if (length(legacy_matches) == 1L) {
                matches <- legacy_matches
            }
            else {
                stop("找到多个同名输入文件，无法安全选择：\n", paste(matches, collapse = "\n"))
            }
        }
        normalizePath(matches[1], winslash = "/", mustWork = TRUE)
    }
    counts_file <- locate_input_file(preferred_paths = c(file.path(raw_dir, "GSE111016", "GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv.gz")), target_basename = "GSE111016_allSamplesCounts_htseqcov1_sss_forGEO.csv.gz", search_roots = c(raw_dir, analysis_root))
    phenotype_file <- locate_input_file(preferred_paths = c(file.path(analysis_root, "00_legacy_reference_inputs", "GSE111016_metadata", "GSE111016_pheno_metadata.rds"), file.path(analysis_root, "GSE111016_metadata", "GSE111016_pheno_metadata.rds")), target_basename = "GSE111016_pheno_metadata.rds", search_roots = c(analysis_root))
    ctd_counts_file <- locate_input_file(preferred_paths = c(file.path(analysis_root, "00_legacy_reference_inputs", "CTD_legacy", "cadmium_human_gene_interaction_counts.csv"), file.path(analysis_root, "CTD_legacy", "cadmium_human_gene_interaction_counts.csv")), target_basename = "cadmium_human_gene_interaction_counts.csv", search_roots = c(analysis_root))
    input_files <- c(CTD_interaction_counts = ctd_counts_file, GSE111016_counts = counts_file, GSE111016_phenotype = phenotype_file)
    input_information <- file.info(input_files)
    input_manifest <- data.frame(input_name = names(input_files), file_name = basename(input_files), full_path = unname(input_files), size_bytes = unname(input_information$size), modified_time = format(input_information$mtime, "%Y-%m-%d %H:%M:%S"), md5 = unname(tools::md5sum(input_files)), stringsAsFactors = FALSE)
    write.csv(input_manifest, file.path(output_dir, "01_input_file_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    screening_parameters <- data.frame(parameter = c("CTD_gene_rule", "CTD_interaction_count_cutoff", "DESeq2_low_expression_prefilter", "DESeq2_design", "DESeq2_reference_group", "DESeq2_comparison", "Figure3A_nominal_DEG_rule", "Candidate_priority_P_rule", "Candidate_priority_abs_log2FC_rule", "Intersection_identifier", "Expected_CTD_gene_count", "Expected_nominal_DEG_count", "Expected_final_candidate_count"), value = c("interaction_count > 10", "10", "None; all 65217 count-matrix rows enter DESeq2", 
        "~ sarcopenia", "no", "yes versus no", "pvalue < 0.05", "pvalue < 0.05", "abs(log2FoldChange) > 0.1", "Official gene symbol", "138", "3342", "23"), interpretation = c("Original CTD cadmium-gene screening rule", "Strictly greater than 10", "Do not apply rowSums/count filters before DESeq2", "One-factor sarcopenia model", "No-sarcopenia samples are the reference", "Positive log2FC means higher expression in sarcopenia", "Count displayed in the right circle of Figure 3A", "Nominal significance component of candidate prioritization", 
        "Original minimal effect-size component", "CTD and expression results are intersected after symbol mapping", "Required reproducibility check", "Required reproducibility check", "Required reproducibility check"), stringsAsFactors = FALSE)
    write.csv(screening_parameters, file.path(output_dir, "00_original_screening_parameters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    ctd_counts <- read.csv(ctd_counts_file, stringsAsFactors = FALSE, check.names = FALSE)
    required_ctd_columns <- c("GeneSymbol", "interaction_count")
    missing_ctd_columns <- setdiff(required_ctd_columns, names(ctd_counts))
    if (length(missing_ctd_columns) > 0) {
        stop("CTD文件缺少以下列：\n", paste(missing_ctd_columns, collapse = ", "))
    }
    ctd_counts$GeneSymbol <- trimws(as.character(ctd_counts$GeneSymbol))
    ctd_counts$interaction_count <- suppressWarnings(as.numeric(ctd_counts$interaction_count))
    ctd_counts <- ctd_counts[!is.na(ctd_counts$GeneSymbol) & ctd_counts$GeneSymbol != "" & !is.na(ctd_counts$interaction_count), , drop = FALSE]
    ctd_counts_aggregated <- stats::aggregate(interaction_count ~ GeneSymbol, data = ctd_counts, FUN = max)
    ctd_genes_gt10 <- ctd_counts_aggregated[ctd_counts_aggregated$interaction_count > 10, , drop = FALSE]
    ctd_genes_gt10 <- ctd_genes_gt10[order(-ctd_genes_gt10$interaction_count, ctd_genes_gt10$GeneSymbol), , drop = FALSE]
    rownames(ctd_genes_gt10) <- NULL
    if (nrow(ctd_genes_gt10) != 138L) {
        stop("CTD筛选没有复现138个基因。\n", "当前结果：", nrow(ctd_genes_gt10), "\n筛选规则：interaction_count > 10")
    }
    write.csv(ctd_genes_gt10, file.path(output_dir, "02_CTD_cadmium_interaction_count_gt10_138_genes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    writeLines(ctd_genes_gt10$GeneSymbol, file.path(output_dir, "02_CTD_cadmium_interaction_count_gt10_138_genes.txt"), useBytes = TRUE)
    phenotype_df <- as.data.frame(readRDS(phenotype_file), stringsAsFactors = FALSE)
    phenotype_names_clean <- trimws(gsub("﻿", "", names(phenotype_df)))
    names(phenotype_df) <- phenotype_names_clean
    required_phenotype_columns <- c("title", "geo_accession", "sarcopenia status:ch1")
    missing_phenotype_columns <- setdiff(required_phenotype_columns, names(phenotype_df))
    if (length(missing_phenotype_columns) > 0) {
        stop("表型文件缺少以下列：\n", paste(missing_phenotype_columns, collapse = ", "))
    }
    extract_sample_tag <- function(x) {
        x <- as.character(x)
        match_position <- regexpr("Sample\\s*[0-9]+", x, ignore.case = TRUE, perl = TRUE)
        match_length <- attr(match_position, "match.length")
        output <- rep(NA_character_, length(x))
        valid <- !is.na(x) & match_position > 0 & match_length > 0
        output[valid] <- substring(x[valid], first = match_position[valid], last = match_position[valid] + match_length[valid] - 1)
        output <- gsub("\\s+", " ", output)
        trimws(output)
    }
    counts_df <- read.csv(gzfile(counts_file), stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(counts_df) != 41L) {
        stop("count文件列数不是41。\n", "当前列数：", ncol(counts_df))
    }
    gene_id_raw <- as.character(counts_df[[1]])
    count_matrix <- as.matrix(counts_df[, -1, drop = FALSE])
    storage.mode(count_matrix) <- "numeric"
    if (nrow(count_matrix) != 65217L) {
        stop("count矩阵没有复现65,217个基因行。\n", "当前行数：", nrow(count_matrix))
    }
    if (ncol(count_matrix) != 40L) {
        stop("count矩阵样本数不是40。\n", "当前列数：", ncol(count_matrix))
    }
    if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
        stop("count矩阵包含NA或非有限值。")
    }
    if (any(count_matrix < 0)) {
        stop("count矩阵包含负值。")
    }
    if (any(count_matrix != round(count_matrix))) {
        stop("count矩阵包含非整数值。")
    }
    if (anyDuplicated(gene_id_raw)) {
        stop("count文件中存在重复的原始Ensembl ID。")
    }
    rownames(count_matrix) <- gene_id_raw
    storage.mode(count_matrix) <- "integer"
    count_sample_tags <- extract_sample_tag(colnames(count_matrix))
    phenotype_sample_tags <- extract_sample_tag(phenotype_df$title)
    if (anyNA(count_sample_tags) || anyNA(phenotype_sample_tags)) {
        stop("无法从count列名或表型title中提取Sample编号。")
    }
    phenotype_order <- match(count_sample_tags, phenotype_sample_tags)
    if (anyNA(phenotype_order)) {
        stop("部分count样本无法匹配到表型记录。")
    }
    phenotype_aligned <- phenotype_df[phenotype_order, , drop = FALSE]
    phenotype_aligned$sample_tag <- extract_sample_tag(phenotype_aligned$title)
    phenotype_aligned$sarcopenia <- factor(tolower(trimws(as.character(phenotype_aligned[["sarcopenia status:ch1"]]))), levels = c("no", "yes"))
    if (anyNA(phenotype_aligned$sarcopenia)) {
        stop("肌少症分组中存在no/yes以外的值。")
    }
    if (!identical(count_sample_tags, phenotype_aligned$sample_tag)) {
        stop("表型数据与count矩阵样本顺序不一致。")
    }
    group_counts <- table(phenotype_aligned$sarcopenia)
    if (as.integer(group_counts["no"]) != 20L || as.integer(group_counts["yes"]) != 20L) {
        stop("分组没有复现20个no和20个yes。\n", paste(names(group_counts), as.integer(group_counts), collapse = "; "))
    }
    phenotype_output <- data.frame(matrix_column = colnames(count_matrix), sample_tag = phenotype_aligned$sample_tag, geo_accession = phenotype_aligned$geo_accession, title = phenotype_aligned$title, sarcopenia = as.character(phenotype_aligned$sarcopenia), stringsAsFactors = FALSE)
    write.csv(phenotype_output, file.path(output_dir, "03_GSE111016_aligned_phenotype_20_no_20_yes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    col_data <- data.frame(sarcopenia = phenotype_aligned$sarcopenia, row.names = colnames(count_matrix))
    if (!identical(colnames(count_matrix), rownames(col_data))) {
        stop("DESeq2 countData列和colData行不一致。")
    }
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~sarcopenia)
    dds$sarcopenia <- stats::relevel(dds$sarcopenia, ref = "no")
    if (nrow(dds) != 65217L) {
        stop("进入DESeq2的基因数不是65,217。")
    }
    set.seed(1234)
    dds <- DESeq2::DESeq(dds, quiet = FALSE)
    deseq_result <- DESeq2::results(dds, contrast = c("sarcopenia", "yes", "no"), alpha = 0.05)
    deseq_df <- as.data.frame(deseq_result)
    deseq_df$EnsemblID_raw <- rownames(deseq_df)
    rownames(deseq_df) <- NULL
    deseq_df$EnsemblID <- sub("\\..*$", "", deseq_df$EnsemblID_raw)
    deseq_df <- deseq_df[, c("EnsemblID_raw", "EnsemblID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
    valid_result_df <- deseq_df[!is.na(deseq_df$pvalue) & is.finite(deseq_df$pvalue) & !is.na(deseq_df$log2FoldChange) & is.finite(deseq_df$log2FoldChange), , drop = FALSE]
    valid_pvalue_count <- nrow(valid_result_df)
    nominal_degs <- valid_result_df[valid_result_df$pvalue < 0.05, , drop = FALSE]
    nominal_degs <- nominal_degs[order(nominal_degs$pvalue, -abs(nominal_degs$log2FoldChange)), , drop = FALSE]
    rownames(nominal_degs) <- NULL
    nominal_deg_count <- nrow(nominal_degs)
    nominal_up_count <- sum(nominal_degs$log2FoldChange > 0)
    nominal_down_count <- sum(nominal_degs$log2FoldChange < 0)
    if (valid_pvalue_count != 53229L) {
        stop("有效P值基因数没有复现53,229。\n", "当前结果：", valid_pvalue_count)
    }
    if (nominal_deg_count != 3342L) {
        stop("nominal P < 0.05没有复现3,342个基因。\n", "当前结果：", nominal_deg_count, "\n\n请检查R包版本和输入文件。")
    }
    if (nominal_up_count != 1771L || nominal_down_count != 1571L) {
        stop("上下调数量没有复现1,771和1,571。\n", "当前上调：", nominal_up_count, "\n当前下调：", nominal_down_count)
    }
    prioritized_degs <- nominal_degs[abs(nominal_degs$log2FoldChange) > 0.1, , drop = FALSE]
    prioritized_degs <- prioritized_degs[order(prioritized_degs$pvalue, -abs(prioritized_degs$log2FoldChange)), , drop = FALSE]
    rownames(prioritized_degs) <- NULL
    symbol_map <- AnnotationDbi::mapIds(x = org.Hs.eg.db::org.Hs.eg.db, keys = unique(deseq_df$EnsemblID), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
    deseq_df$SYMBOL <- unname(symbol_map[deseq_df$EnsemblID])
    nominal_degs$SYMBOL <- unname(symbol_map[nominal_degs$EnsemblID])
    prioritized_degs$SYMBOL <- unname(symbol_map[prioritized_degs$EnsemblID])
    mapped_result_count <- sum(!is.na(deseq_df$SYMBOL) & deseq_df$SYMBOL != "")
    prioritized_symbol_all <- prioritized_degs[!is.na(prioritized_degs$SYMBOL) & prioritized_degs$SYMBOL != "", , drop = FALSE]
    prioritized_symbol_all$abs_log2FC <- abs(prioritized_symbol_all$log2FoldChange)
    prioritized_symbol_all <- prioritized_symbol_all[order(prioritized_symbol_all$SYMBOL, prioritized_symbol_all$pvalue, -prioritized_symbol_all$abs_log2FC), , drop = FALSE]
    prioritized_symbol_unique <- prioritized_symbol_all[!duplicated(prioritized_symbol_all$SYMBOL), , drop = FALSE]
    prioritized_symbol_unique <- prioritized_symbol_unique[order(prioritized_symbol_unique$pvalue, -prioritized_symbol_unique$abs_log2FC), , drop = FALSE]
    rownames(prioritized_symbol_unique) <- NULL
    nominal_symbol_all <- nominal_degs[!is.na(nominal_degs$SYMBOL) & nominal_degs$SYMBOL != "", , drop = FALSE]
    nominal_symbol_all$abs_log2FC <- abs(nominal_symbol_all$log2FoldChange)
    nominal_ctd_overlap_all <- nominal_symbol_all[nominal_symbol_all$SYMBOL %in% ctd_genes_gt10$GeneSymbol, , drop = FALSE]
    nominal_ctd_overlap_all <- nominal_ctd_overlap_all[order(nominal_ctd_overlap_all$SYMBOL, nominal_ctd_overlap_all$pvalue, -nominal_ctd_overlap_all$abs_log2FC), , drop = FALSE]
    nominal_ctd_overlap_unique <- nominal_ctd_overlap_all[!duplicated(nominal_ctd_overlap_all$SYMBOL), , drop = FALSE]
    nominal_ctd_overlap_unique <- nominal_ctd_overlap_unique[order(nominal_ctd_overlap_unique$SYMBOL), , drop = FALSE]
    rownames(nominal_ctd_overlap_unique) <- NULL
    if (nrow(nominal_ctd_overlap_unique) != 24L) {
        stop("CTD 138与nominal P < 0.05的直接交集没有复现24个基因。\n", "当前结果：", nrow(nominal_ctd_overlap_unique))
    }
    candidate_overlap_all <- prioritized_symbol_all[prioritized_symbol_all$SYMBOL %in% ctd_genes_gt10$GeneSymbol, , drop = FALSE]
    candidate_overlap_all <- candidate_overlap_all[order(candidate_overlap_all$SYMBOL, candidate_overlap_all$pvalue, -candidate_overlap_all$abs_log2FC), , drop = FALSE]
    candidate_overlap_unique <- candidate_overlap_all[!duplicated(candidate_overlap_all$SYMBOL), , drop = FALSE]
    candidate_overlap_unique$direction <- ifelse(candidate_overlap_unique$log2FoldChange > 0, "Up", "Down")
    expected_candidate_genes <- c("ABCG2", "ATF4", "CASP3", "CASP9", "CCL2", "CYCS", "CYP1A1", "DDIT3", "FOXO1", "FOXO3", "HIF1A", "HMGA2", "MAP1LC3B", "MFN2", "MMP9", "MTF1", "PHGDH", "PKM", "PPARGC1A", "SLC16A1", "SP1", "STAT3", "VEGFA")
    observed_candidate_genes <- sort(unique(candidate_overlap_unique$SYMBOL))
    if (length(observed_candidate_genes) != 23L) {
        stop("最终候选基因数量没有复现23。\n", "当前结果：", length(observed_candidate_genes))
    }
    if (!identical(observed_candidate_genes, sort(expected_candidate_genes))) {
        missing_expected <- setdiff(expected_candidate_genes, observed_candidate_genes)
        unexpected_observed <- setdiff(observed_candidate_genes, expected_candidate_genes)
        stop("最终23基因列表与原始列表不一致。\n\n", "缺少：", paste(missing_expected, collapse = ", "), "\n多出：", paste(unexpected_observed, collapse = ", "))
    }
    candidate_overlap_unique$gene_order <- match(candidate_overlap_unique$SYMBOL, expected_candidate_genes)
    candidate_final <- candidate_overlap_unique[order(candidate_overlap_unique$gene_order), , drop = FALSE]
    candidate_final <- candidate_final[, c("SYMBOL", "EnsemblID_raw", "EnsemblID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "abs_log2FC", "direction")]
    rownames(candidate_final) <- NULL
    slc11a2_rows <- nominal_symbol_all[nominal_symbol_all$SYMBOL == "SLC11A2", , drop = FALSE]
    slc11a2_rows <- slc11a2_rows[order(slc11a2_rows$pvalue, -slc11a2_rows$abs_log2FC), , drop = FALSE]
    if (nrow(slc11a2_rows) == 0L) {
        stop("在nominal P < 0.05结果中没有找到SLC11A2。")
    }
    slc11a2_rows$in_CTD_138 <- (slc11a2_rows$SYMBOL %in% ctd_genes_gt10$GeneSymbol)
    slc11a2_rows$passes_nominal_P005 <- (slc11a2_rows$pvalue < 0.05)
    slc11a2_rows$passes_abs_log2FC_gt_0.1 <- (slc11a2_rows$abs_log2FC > 0.1)
    slc11a2_rows$included_in_final_23 <- (slc11a2_rows$SYMBOL %in% candidate_final$SYMBOL)
    slc11a2_rows$exclusion_reason <- ifelse(!slc11a2_rows$passes_abs_log2FC_gt_0.1, paste0("Excluded by the original minimal effect-size rule: ", "|log2FoldChange| <= 0.1"), "")
    if (!all(slc11a2_rows$in_CTD_138) || !all(slc11a2_rows$passes_nominal_P005) || any(slc11a2_rows$passes_abs_log2FC_gt_0.1) || any(slc11a2_rows$included_in_final_23)) {
        stop("SLC11A2边界核查结果与原始筛选逻辑不一致。")
    }
    write.csv(deseq_df, file.path(output_dir, "04_GSE111016_DESeq2_all_results_65217_Ensembl.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(nominal_degs, file.path(output_dir, "05_GSE111016_nominal_P005_DEGs_3342_Ensembl.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(prioritized_degs, file.path(output_dir, "06_GSE111016_original_priority_P005_absLog2FC_gt0.1_Ensembl.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(prioritized_symbol_unique, file.path(output_dir, "07_GSE111016_original_priority_P005_absLog2FC_gt0.1_unique_Symbol.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(nominal_ctd_overlap_unique, file.path(output_dir, "08_CTD138_nominal_P005_direct_overlap_24_genes_AUDIT.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(candidate_overlap_all, file.path(output_dir, "09_CTD138_original_priority_overlap_all_Ensembl.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(candidate_final, file.path(output_dir, "10_FINAL_original_23_candidate_genes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    writeLines(expected_candidate_genes, file.path(output_dir, "10_FINAL_original_23_candidate_genes.txt"), useBytes = TRUE)
    write.csv(slc11a2_rows, file.path(output_dir, "11_SLC11A2_threshold_boundary_AUDIT.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    saveRDS(dds, file.path(output_dir, "14_GSE111016_DESeq2_dds_unfiltered_65217.rds"))
    saveRDS(deseq_result, file.path(output_dir, "14_GSE111016_DESeq2_results_object.rds"))
    screening_summary <- data.frame(item = c("Input count-matrix rows", "Input samples", "No-sarcopenia samples", "Sarcopenia samples", "CTD genes with interaction_count > 10", "DESeq2 rows with valid P values", "Nominal P < 0.05 results shown in Figure 3A", "Nominally upregulated genes", "Nominally downregulated genes", "Prioritized Ensembl rows: P < 0.05 and |log2FC| > 0.1", "Prioritized unique mapped symbols", "Direct CTD138 versus nominal-P overlap", "Original CTD138 versus prioritized-DEG overlap", 
        "Final candidate genes", "SLC11A2 in CTD 138", "SLC11A2 nominal P < 0.05", "SLC11A2 passes |log2FC| > 0.1", "Mapped DESeq2 rows"), value = c(nrow(count_matrix), ncol(count_matrix), as.integer(group_counts["no"]), as.integer(group_counts["yes"]), nrow(ctd_genes_gt10), valid_pvalue_count, nominal_deg_count, nominal_up_count, nominal_down_count, nrow(prioritized_degs), nrow(prioritized_symbol_unique), nrow(nominal_ctd_overlap_unique), nrow(candidate_final), length(expected_candidate_genes), 
        all(slc11a2_rows$in_CTD_138), all(slc11a2_rows$passes_nominal_P005), any(slc11a2_rows$passes_abs_log2FC_gt_0.1), mapped_result_count), expected_or_interpretation = c("65217", "40", "20", "20", "138", "53229", "3342; retained as the Figure 3A right-circle count", "1771", "1571", "Calculated and recorded at runtime", "Calculated and recorded at runtime", "24; audit result before applying the effect-size criterion", "23", "23", "TRUE", "TRUE", "FALSE", "Depends on the installed org.Hs.eg.db annotation version"), 
        stringsAsFactors = FALSE)
    write.csv(screening_summary, file.path(output_dir, "12_original_screening_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    screening_record <- list(record_name = "Authoritative original Figure 3A candidate screening", created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), project_dir = project_dir, output_dir = output_dir, input_manifest = input_manifest, parameters = screening_parameters, phenotype = phenotype_output, CTD_all_counts = ctd_counts_aggregated, CTD_138 = ctd_genes_gt10, DESeq2_results_all = deseq_df, nominal_DEGs_3342 = nominal_degs, original_prioritized_DEGs = prioritized_degs, original_prioritized_unique_symbols = prioritized_symbol_unique, 
        direct_nominal_overlap_24 = nominal_ctd_overlap_unique, final_candidate_genes_23 = candidate_final, SLC11A2_boundary_audit = slc11a2_rows, summary = screening_summary, final_candidate_gene_vector = expected_candidate_genes)
    saveRDS(screening_record, file.path(output_dir, "13_original_Figure3A_screening_record.rds"))
    session_information <- capture.output(sessionInfo())
    package_versions <- c(paste0("DESeq2: ", as.character(utils::packageVersion("DESeq2"))), paste0("AnnotationDbi: ", as.character(utils::packageVersion("AnnotationDbi"))), paste0("org.Hs.eg.db: ", as.character(utils::packageVersion("org.Hs.eg.db"))))
    writeLines(c("Original Figure 3A screening session information", paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "", package_versions, "", session_information), file.path(output_dir, "15_sessionInfo_and_package_versions.txt"), useBytes = TRUE)
    candidate_gene_text <- paste(expected_candidate_genes, collapse = ", ")
    readme_lines <- c("============================================================", "Figure 3A original candidate-gene screening record", "Figure 3A 原始候选基因筛选记录", "============================================================", "", paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "", "AUTHORITATIVE SCREENING LOGIC / 权威筛选逻辑", "", "1. CTD cadmium-related genes:", "   interaction_count > 10", "   Result: 138 genes.", "", "2. GSE111016 differential-expression analysis:", 
        "   DESeq2 design: ~ sarcopenia", "   Comparison: sarcopenia yes versus no", "   Reference group: no", "   All 65,217 count-matrix rows entered DESeq2.", "   No low-expression prefiltering was performed.", "", "3. Nominal DEG pool displayed in Figure 3A:", "   pvalue < 0.05", "   Result: 3,342 Ensembl-level DESeq2 results.", "   This is the value displayed in the right circle of Figure 3A.", "", "4. Original candidate-prioritization rule:", "   pvalue < 0.05 AND abs(log2FoldChange) > 0.1", 
        "", "5. Intersection:", "   CTD 138 genes were intersected with prioritized", "   GSE111016 genes after mapping to official gene symbols.", "   Result: 23 candidate genes.", "", "IMPORTANT INTERPRETATION / 重要解释", "", "The 23 genes are not the direct intersection obtained using", "pvalue < 0.05 alone. The direct nominal-P intersection contains", "24 genes, including SLC11A2.", "", "SLC11A2 satisfies nominal pvalue < 0.05 and is present in the", "CTD 138-gene set, but its absolute log2 fold change is below", 
        "0.1. It is therefore excluded automatically by the original", "minimal effect-size criterion; it was not manually deleted.", "", "Figure 3A should therefore be interpreted as an integration and", "candidate-prioritization schematic rather than a strict", "unfiltered two-set Venn diagram.", "", "中文说明：", "", "Figure 3A右侧的3342表示DESeq2分析中nominal P < 0.05的", "全部结果。最终23个候选基因在此基础上进一步满足", "|log2FC| > 0.1，并与CTD的138个基因在Gene Symbol层面", 
        "取交集。SLC11A2由于|log2FC|未超过0.1而自动排除，", "并非后期人工删除。", "", "FINAL 23 CANDIDATE GENES / 最终23个候选基因", "", candidate_gene_text, "", "FILES", "", "00_original_screening_parameters.csv", "    Complete thresholds and methodological settings.", "", "01_input_file_manifest.csv", "    Source paths, file sizes, modification times and MD5 values.", "", "02_CTD_cadmium_interaction_count_gt10_138_genes.csv/txt", "    Reproduced CTD 138-gene set.", 
        "", "03_GSE111016_aligned_phenotype_20_no_20_yes.csv", "    Verified sample alignment and group assignments.", "", "04_GSE111016_DESeq2_all_results_65217_Ensembl.csv", "    Complete unfiltered DESeq2 result table.", "", "05_GSE111016_nominal_P005_DEGs_3342_Ensembl.csv", "    Nominal DEG pool displayed as n = 3342 in Figure 3A.", "", "06_GSE111016_original_priority_P005_absLog2FC_gt0.1_Ensembl.csv", "    Original prioritized Ensembl-level DEG set.", "", "07_GSE111016_original_priority_P005_absLog2FC_gt0.1_unique_Symbol.csv", 
        "    Prioritized DEG set after symbol mapping and deduplication.", "", "08_CTD138_nominal_P005_direct_overlap_24_genes_AUDIT.csv", "    Audit-only direct overlap before the |log2FC| criterion.", "", "09_CTD138_original_priority_overlap_all_Ensembl.csv", "    CTD-prioritized overlap before final symbol deduplication.", "", "10_FINAL_original_23_candidate_genes.csv/txt", "    Authoritative final candidate-gene set.", "", "11_SLC11A2_threshold_boundary_AUDIT.csv", "    Explicit record explaining the SLC11A2 boundary case.", 
        "", "12_original_screening_summary.csv", "    Key counts and reproducibility checks.", "", "13_original_Figure3A_screening_record.rds", "    Consolidated R object containing all important tables.", "", "14_GSE111016_DESeq2_dds_unfiltered_65217.rds", "14_GSE111016_DESeq2_results_object.rds", "    Saved DESeq2 objects for future reuse.", "", "15_sessionInfo_and_package_versions.txt", "    R environment and package versions.", "", "Do not replace this original screening rule with a human-only", 
        "recount of the current CTD interaction table, and do not apply", "low-expression filtering before DESeq2 when reproducing Figure 3A.", "", "============================================================")
    readme_file <- file.path(output_dir, "00_README_original_Figure3A_screening.txt")
    writeLines(readme_lines, readme_file, useBytes = TRUE)
    generated_files <- list.files(output_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE)
    generated_information <- file.info(generated_files)
    output_manifest <- data.frame(file_name = basename(generated_files), size_bytes = unname(generated_information$size), size_MB = round(unname(generated_information$size)/1024^2, 4), modified_time = format(generated_information$mtime, "%Y-%m-%d %H:%M:%S"), md5 = unname(tools::md5sum(generated_files)), full_path = normalizePath(generated_files, winslash = "/", mustWork = FALSE), stringsAsFactors = FALSE)
    output_manifest <- output_manifest[order(output_manifest$file_name), , drop = FALSE]
    rownames(output_manifest) <- NULL
    write.csv(output_manifest, file.path(output_dir, "16_output_file_manifest.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    cat("\n============================================================\n", "Original Figure 3A screening reproduced successfully\n", "============================================================\n", "\nOutput directory:\n", output_dir, "\n\nCTD genes, interaction_count > 10:\n", nrow(ctd_genes_gt10), "\n\nDESeq2 rows entering analysis:\n", nrow(dds), "\n\nValid DESeq2 P values:\n", valid_pvalue_count, "\n\nNominal P < 0.05 genes shown in Figure 3A:\n", nominal_deg_count, "\n\nNominal upregulated genes:\n", 
        nominal_up_count, "\n\nNominal downregulated genes:\n", nominal_down_count, "\n\nDirect nominal-P CTD overlap:\n", nrow(nominal_ctd_overlap_unique), "\n\nFinal prioritized CTD overlap:\n", nrow(candidate_final), "\n\nFinal 23 genes:\n", paste(expected_candidate_genes, collapse = ", "), "\n\nREADME:\n", readme_file, "\n")
    list(project_dir = project_dir, output_dir = output_dir, readme_file = readme_file, input_manifest = input_manifest, screening_parameters = screening_parameters, screening_summary = screening_summary, CTD_138 = ctd_genes_gt10, nominal_DEGs_3342 = nominal_degs, original_prioritized_DEGs = prioritized_degs, direct_nominal_overlap_24 = nominal_ctd_overlap_unique, final_candidate_genes_23 = candidate_final, final_candidate_gene_vector = expected_candidate_genes, SLC11A2_boundary_audit = slc11a2_rows)
}

original_screening_result <- run_original_figure3A_screening()
