run_figure5_abcd_final <- function() {
    
    options(stringsAsFactors = FALSE)
    
    ############################################################
    ## 0. Project paths
    ############################################################
    
    project_root <- "E:/1Data/Cd_Sarcopenia_NHANES"
    
    m05_root <- file.path(
        project_root,
        "new_workspace",
        "03_results",
        "M05_PathwayValidation_GSEA_ssGSEA_Correlation"
    )
    
    analysis_dir <- file.path(
        m05_root,
        "02_analysis_files"
    )
    
    canonical_input_dir <- file.path(
        analysis_dir,
        "canonical_inputs"
    )
    
    figure5_analysis_dir <- file.path(
        analysis_dir,
        "figure5_rebuild_final"
    )
    
    canonical_core_dir <- file.path(
        analysis_dir,
        "canonical_core"
    )
    
    figure_dir <- file.path(
        m05_root,
        "03_outputs_figures_tables",
        "figure_re"
    )
    
    backup_dir <- file.path(
        m05_root,
        "backup"
    )
    
    dir.create(
        figure5_analysis_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    dir.create(
        canonical_core_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    dir.create(
        figure_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    cat("\n============================================================\n")
    cat("FIGURE 5A–D COMPLETE RECONSTRUCTION — FIXED\n")
    cat("============================================================\n")
    
    cat("\nAnalysis directory:\n")
    cat(figure5_analysis_dir, "\n")
    
    cat("\nFigure directory:\n")
    cat(figure_dir, "\n")
    
    cat("\nBackup will not be modified:\n")
    cat(backup_dir, "\n")
    
    ############################################################
    ## 1. Required packages
    ############################################################
    
    required_packages <- c(
        "GSVA",
        "AnnotationDbi",
        "org.Hs.eg.db",
        "ggplot2",
        "ggbeeswarm",
        "patchwork",
        "scales",
        "stringr"
    )
    
    missing_packages <- required_packages[
        !vapply(
            required_packages,
            requireNamespace,
            logical(1),
            quietly = TRUE
        )
    ]
    
    if (length(missing_packages) > 0L) {
        
        stop(
            "Required package(s) not installed: ",
            paste(
                missing_packages,
                collapse = ", "
            ),
            call. = FALSE
        )
    }
    
    ############################################################
    ## 2. Locate canonical input files
    ############################################################
    
    bc_input_path <- file.path(
        canonical_input_dir,
        "M05_BC_GSE111016_log2CPM_group_canonical_input.rds"
    )
    
    gsea_path_candidates <- c(
        file.path(
            canonical_input_dir,
            "M05A_GO_BP_GSEA_clusterProfiler_object_v2.rds"
        ),
        
        file.path(
            analysis_dir,
            "M05A_GO_BP_GSEA_clusterProfiler_object_v2.rds"
        ),
        
        file.path(
            analysis_dir,
            "legacy_source_materials",
            "M05A_GO_BP_GSEA_clusterProfiler_object_v2.rds"
        ),
        
        file.path(
            analysis_dir,
            "canonical_inputs",
            "M05A_offline_GO_BP_GSEA_clusterProfiler_object_v1.rds"
        ),
        
        file.path(
            analysis_dir,
            "M05A_offline_GO_BP_GSEA_clusterProfiler_object_v1.rds"
        )
    )
    
    gsea_path_candidates <- unique(
        gsea_path_candidates
    )
    
    existing_gsea_paths <- gsea_path_candidates[
        file.exists(
            gsea_path_candidates
        )
    ]
    
    if (!file.exists(bc_input_path)) {
        
        stop(
            "Canonical B/C input was not found:\n",
            bc_input_path,
            call. = FALSE
        )
    }
    
    if (length(existing_gsea_paths) == 0L) {
        
        stop(
            paste0(
                "No canonical GSEA RDS was found. Checked:\n",
                paste(
                    gsea_path_candidates,
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }
    
    gsea_path <- existing_gsea_paths[1]
    
    cat("\nCanonical B/C input:\n")
    cat(bc_input_path, "\n")
    
    cat("\nCanonical GSEA input:\n")
    cat(gsea_path, "\n")
    
    ############################################################
    ## 3. Load canonical B/C input
    ############################################################
    
    bc_input <- readRDS(
        bc_input_path
    )
    
    required_bc_elements <- c(
        "expr",
        "group_df",
        "hub_genes"
    )
    
    missing_bc_elements <- required_bc_elements[
        !required_bc_elements %in%
            names(bc_input)
    ]
    
    if (length(missing_bc_elements) > 0L) {
        
        stop(
            "Canonical B/C input lacks element(s): ",
            paste(
                missing_bc_elements,
                collapse = ", "
            ),
            call. = FALSE
        )
    }
    
    expr <- as.matrix(
        getElement(
            bc_input,
            "expr"
        )
    )
    
    group_df <- as.data.frame(
        getElement(
            bc_input,
            "group_df"
        ),
        stringsAsFactors = FALSE
    )
    
    rownames(expr) <- toupper(
        trimws(
            rownames(expr)
        )
    )
    
    colnames(expr) <- trimws(
        colnames(expr)
    )
    
    group_df$SampleID <- trimws(
        as.character(
            group_df$SampleID
        )
    )
    
    group_df$Group <- factor(
        as.character(
            group_df$Group
        ),
        levels = c(
            "Control",
            "Sarcopenia"
        )
    )
    
    group_df <- group_df[
        match(
            colnames(expr),
            group_df$SampleID
        ),
        ,
        drop = FALSE
    ]
    
    rownames(group_df) <- NULL
    
    hub_genes <- c(
        "ATF4",
        "CASP9",
        "FOXO1",
        "FOXO3",
        "MAP1LC3B",
        "PPARGC1A",
        "CASP3",
        "CYCS",
        "HIF1A",
        "STAT3",
        "CCL2"
    )
    
    strict_hub_genes <- c(
        "ATF4",
        "CASP9"
    )
    
    robust_hub_genes <- setdiff(
        hub_genes,
        strict_hub_genes
    )
    
    stopifnot(
        is.matrix(expr),
        is.numeric(expr),
        
        ncol(expr) == 40L,
        nrow(group_df) == 40L,
        
        identical(
            colnames(expr),
            group_df$SampleID
        ),
        
        sum(
            group_df$Group == "Control"
        ) == 20L,
        
        sum(
            group_df$Group == "Sarcopenia"
        ) == 20L,
        
        all(
            hub_genes %in%
                rownames(expr)
        )
    )
    
    cat("\nExpression dimensions:\n")
    print(dim(expr))
    
    cat("\nGroup distribution:\n")
    print(table(group_df$Group))
    
    ############################################################
    ## 4. Load canonical GSEA object
    ############################################################
    
    gsea_object <- readRDS(
        gsea_path
    )
    
    if (!methods::is(gsea_object, "gseaResult")) {
        
        stop(
            paste0(
                "Canonical GSEA input is not a gseaResult object. ",
                "Detected class: ",
                paste(
                    class(gsea_object),
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }
    
    required_gsea_slots <- c(
        "result",
        "geneList",
        "geneSets"
    )
    
    missing_gsea_slots <- required_gsea_slots[
        !required_gsea_slots %in%
            methods::slotNames(
                gsea_object
            )
    ]
    
    if (length(missing_gsea_slots) > 0L) {
        
        stop(
            "GSEA object lacks slot(s): ",
            paste(
                missing_gsea_slots,
                collapse = ", "
            ),
            call. = FALSE
        )
    }
    
    gsea_result <- as.data.frame(
        methods::slot(
            gsea_object,
            "result"
        ),
        stringsAsFactors = FALSE
    )
    
    rank_vec_original <- methods::slot(
        gsea_object,
        "geneList"
    )
    
    rank_vec <- as.numeric(
        rank_vec_original
    )
    
    names(rank_vec) <- names(
        rank_vec_original
    )
    
    gsea_gene_sets <- methods::slot(
        gsea_object,
        "geneSets"
    )
    
    required_gsea_columns <- c(
        "ID",
        "Description",
        "NES",
        "p.adjust"
    )
    
    missing_gsea_columns <- setdiff(
        required_gsea_columns,
        colnames(gsea_result)
    )
    
    if (length(missing_gsea_columns) > 0L) {
        
        stop(
            "GSEA result table lacks column(s): ",
            paste(
                missing_gsea_columns,
                collapse = ", "
            ),
            call. = FALSE
        )
    }
    
    valid_rank_entries <- (
        is.finite(rank_vec) &
            !is.na(names(rank_vec)) &
            nzchar(names(rank_vec))
    )
    
    rank_vec <- rank_vec[
        valid_rank_entries
    ]
    
    rank_vec <- rank_vec[
        !duplicated(
            names(rank_vec)
        )
    ]
    
    if (!all(diff(rank_vec) <= 0)) {
        
        rank_vec <- sort(
            rank_vec,
            decreasing = TRUE
        )
    }
    
    cat("\nGSEA result dimensions:\n")
    print(dim(gsea_result))
    
    cat("\nRanked gene-list length:\n")
    print(length(rank_vec))
    
    cat("\nStored gene-set count:\n")
    print(length(gsea_gene_sets))
    
    ############################################################
    ## 5. Define four mitochondrial pathways
    ############################################################
    
    target_pathways <- data.frame(
        Pathway = c(
            "Oxidative phosphorylation",
            "Respiratory electron transport chain",
            "Mitochondrial ATP synthesis",
            "Mitochondrial translation"
        ),
        
        Short = c(
            "Oxidative phosphorylation",
            "Respiratory ETC",
            "Mitochondrial ATP synthesis",
            "Mitochondrial translation"
        ),
        
        GO_ID = c(
            "GO:0006119",
            "GO:0022904",
            "GO:0042775",
            "GO:0032543"
        ),
        
        stringsAsFactors = FALSE
    )
    
    pathway_order <- target_pathways$Pathway
    pathway_short_order <- target_pathways$Short
    
    get_stored_gene_set <- function(
        go_id,
        pathway_description
    ) {
        
        available_names <- names(
            gsea_gene_sets
        )
        
        match_index <- which(
            available_names == go_id
        )
        
        if (length(match_index) == 0L) {
            
            match_index <- which(
                tolower(available_names) ==
                    tolower(pathway_description)
            )
        }
        
        if (length(match_index) == 0L) {
            
            match_index <- grep(
                go_id,
                available_names,
                fixed = TRUE
            )
        }
        
        if (length(match_index) == 0L) {
            
            stop(
                "Stored GSEA gene set not found for ",
                go_id,
                " / ",
                pathway_description,
                call. = FALSE
            )
        }
        
        unique(
            as.character(
                getElement(
                    gsea_gene_sets,
                    match_index[1]
                )
            )
        )
    }
    
    pathway_entrez <- vector(
        "list",
        nrow(target_pathways)
    )
    
    names(pathway_entrez) <- pathway_order
    
    for (
        pathway_index in seq_len(
            nrow(target_pathways)
        )
    ) {
        
        pathway_name <- target_pathways$Pathway[
            pathway_index
        ]
        
        pathway_go_id <- target_pathways$GO_ID[
            pathway_index
        ]
        
        pathway_entrez[
            pathway_name
        ] <- list(
            get_stored_gene_set(
                go_id = pathway_go_id,
                pathway_description = pathway_name
            )
        )
    }
    
    ############################################################
    ## 6. Map stored Entrez gene sets to symbols
    ############################################################
    
    annotation_database <- get(
        "org.Hs.eg.db",
        envir = asNamespace(
            "org.Hs.eg.db"
        )
    )
    
    pathway_symbols_raw <- vector(
        "list",
        length(pathway_entrez)
    )
    
    names(pathway_symbols_raw) <- names(
        pathway_entrez
    )
    
    pathway_mapping_records <- list()
    
    for (
        pathway_name in names(
            pathway_entrez
        )
    ) {
        
        entrez_ids <- unique(
            as.character(
                getElement(
                    pathway_entrez,
                    pathway_name
                )
            )
        )
        
        entrez_ids <- entrez_ids[
            !is.na(entrez_ids) &
                nzchar(entrez_ids)
        ]
        
        mapped_symbols <- AnnotationDbi::mapIds(
            annotation_database,
            keys = entrez_ids,
            keytype = "ENTREZID",
            column = "SYMBOL",
            multiVals = "first"
        )
        
        symbol_values <- unname(
            mapped_symbols[
                match(
                    entrez_ids,
                    names(mapped_symbols)
                )
            ]
        )
        
        symbol_values <- toupper(
            trimws(
                as.character(symbol_values)
            )
        )
        
        valid_symbols <- (
            !is.na(symbol_values) &
                nzchar(symbol_values)
        )
        
        pathway_symbols_raw[
            pathway_name
        ] <- list(
            unique(
                symbol_values[
                    valid_symbols
                ]
            )
        )
        
        pathway_mapping_records <- append(
            pathway_mapping_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    ENTREZID = entrez_ids,
                    SYMBOL = symbol_values,
                    Mapped = valid_symbols,
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    pathway_mapping_table <- do.call(
        rbind,
        pathway_mapping_records
    )
    
    pathways <- lapply(
        pathway_symbols_raw,
        function(symbol_values) {
            
            intersect(
                unique(symbol_values),
                rownames(expr)
            )
        }
    )
    
    pathways <- pathways[
        pathway_order
    ]
    
    pathway_size_table <- data.frame(
        Pathway = pathway_order,
        
        GO_ID = target_pathways$GO_ID,
        
        N_stored_Entrez = vapply(
            pathway_entrez[
                pathway_order
            ],
            length,
            integer(1)
        ),
        
        N_mapped_symbols = vapply(
            pathway_symbols_raw[
                pathway_order
            ],
            length,
            integer(1)
        ),
        
        N_detected_in_expression = vapply(
            pathways[
                pathway_order
            ],
            length,
            integer(1)
        ),
        
        stringsAsFactors = FALSE
    )
    
    if (
        any(
            pathway_size_table$N_detected_in_expression < 5L
        )
    ) {
        
        stop(
            "At least one pathway has fewer than five detected genes.",
            call. = FALSE
        )
    }
    
    write.csv(
        pathway_mapping_table,
        file.path(
            figure5_analysis_dir,
            "00_pathway_Entrez_to_symbol_mapping.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        pathway_size_table,
        file.path(
            figure5_analysis_dir,
            "01_four_pathway_gene_set_sizes.csv"
        ),
        row.names = FALSE
    )
    
    cat("\nCanonical pathway sizes:\n")
    
    print(
        pathway_size_table,
        row.names = FALSE
    )
    
    ############################################################
    ## 7. GSVA-compatible ssGSEA helper
    ############################################################
    
    run_ssgsea <- function(
        expression_matrix,
        gene_sets
    ) {
        
        new_api_available <- (
            "ssgseaParam" %in%
                getNamespaceExports(
                    "GSVA"
                )
        )
        
        new_api_error <- "Parameter-based API not available"
        
        if (new_api_available) {
            
            new_api_result <- tryCatch(
                {
                    parameter_object <- GSVA::ssgseaParam(
                        exprData = expression_matrix,
                        geneSets = gene_sets,
                        normalize = TRUE
                    )
                    
                    GSVA::gsva(
                        parameter_object,
                        verbose = FALSE
                    )
                },
                error = function(e) {
                    e
                }
            )
            
            if (!inherits(new_api_result, "error")) {
                
                return(
                    as.matrix(
                        new_api_result
                    )
                )
            }
            
            new_api_error <- conditionMessage(
                new_api_result
            )
        }
        
        legacy_api_result <- tryCatch(
            {
                GSVA::gsva(
                    expression_matrix,
                    gene_sets,
                    method = "ssgsea",
                    kcdf = "Gaussian",
                    ssgsea.norm = TRUE,
                    verbose = FALSE
                )
            },
            error = function(e) {
                e
            }
        )
        
        if (!inherits(legacy_api_result, "error")) {
            
            return(
                as.matrix(
                    legacy_api_result
                )
            )
        }
        
        stop(
            paste0(
                "Both GSVA ssGSEA interfaces failed.\n\n",
                "Parameter-based API:\n",
                new_api_error,
                "\n\nLegacy API:\n",
                conditionMessage(
                    legacy_api_result
                )
            ),
            call. = FALSE
        )
    }
    
    ############################################################
    ## 8. Figure 5A — GSEA mechanism dot plot
    ############################################################
    
    assign_mechanism <- function(description) {
        
        description_lower <- tolower(
            as.character(description)
        )
        
        if (
            grepl(
                paste0(
                    "mitochond|respiratory chain|electron transport|",
                    "oxidative phosphorylation|aerobic respiration|",
                    "cellular respiration|atp synthesis"
                ),
                description_lower
            )
        ) {
            return("Mitochondrial respiration")
        }
        
        if (grepl("apopt", description_lower)) {
            return("Apoptosis")
        }
        
        if (grepl("autophag", description_lower)) {
            return("Autophagy")
        }
        
        if (
            grepl(
                "endoplasmic reticulum|unfolded protein|er stress",
                description_lower
            )
        ) {
            return("ER stress")
        }
        
        if (
            grepl(
                paste0(
                    "inflamm|cytokine|interleukin|",
                    "tumor necrosis factor|nf-kappa"
                ),
                description_lower
            )
        ) {
            return("Inflammation")
        }
        
        if (
            grepl(
                "ubiquitin|proteasome|protein catabolic",
                description_lower
            )
        ) {
            return("Protein degradation")
        }
        
        if (
            grepl(
                "muscle|myoblast|myotube|atrophy",
                description_lower
            )
        ) {
            return("Muscle remodeling")
        }
        
        if (
            grepl(
                "oxidative|reactive oxygen",
                description_lower
            )
        ) {
            return("Oxidative stress")
        }
        
        NA_character_
    }
    
    gsea_result$ID <- as.character(
        gsea_result$ID
    )
    
    gsea_result$Description <- as.character(
        gsea_result$Description
    )
    
    gsea_result$NES <- as.numeric(
        gsea_result$NES
    )
    
    gsea_result$p.adjust <- as.numeric(
        gsea_result$p.adjust
    )
    
    gsea_result$Mechanism <- vapply(
        gsea_result$Description,
        assign_mechanism,
        character(1)
    )
    
    mechanism_df <- gsea_result[
        !is.na(
            gsea_result$Mechanism
        ) &
            is.finite(
                gsea_result$NES
            ) &
            is.finite(
                gsea_result$p.adjust
            ),
        ,
        drop = FALSE
    ]
    
    mechanism_df$FDR <- mechanism_df$p.adjust
    
    mechanism_df$neg_log10_FDR <- -log10(
        pmax(
            mechanism_df$FDR,
            .Machine$double.xmin
        )
    )
    
    mechanism_df$Direction <- ifelse(
        mechanism_df$NES >= 0,
        "Enriched in sarcopenia",
        "Enriched in control"
    )
    
    mechanism_df$rank_in_mechanism <- NA_integer_
    
    mechanism_groups <- unique(
        mechanism_df$Mechanism
    )
    
    for (mechanism_name in mechanism_groups) {
        
        mechanism_rows <- which(
            mechanism_df$Mechanism ==
                mechanism_name
        )
        
        mechanism_order <- order(
            mechanism_df$p.adjust[
                mechanism_rows
            ],
            -abs(
                mechanism_df$NES[
                    mechanism_rows
                ]
            )
        )
        
        mechanism_df$rank_in_mechanism[
            mechanism_rows[
                mechanism_order
            ]
        ] <- seq_along(
            mechanism_rows
        )
    }
    
    keep_mechanism_rows <- (
        (
            mechanism_df$Mechanism ==
                "Mitochondrial respiration" &
                mechanism_df$rank_in_mechanism <= 8L
        ) |
            (
                mechanism_df$Mechanism !=
                    "Mitochondrial respiration" &
                    mechanism_df$rank_in_mechanism <= 3L
            )
    )
    
    figure5a_display <- mechanism_df[
        keep_mechanism_rows,
        ,
        drop = FALSE
    ]
    
    figure5a_display <- figure5a_display[
        order(
            figure5a_display$p.adjust,
            -abs(
                figure5a_display$NES
            )
        ),
        ,
        drop = FALSE
    ]
    
    figure5a_display <- head(
        figure5a_display,
        18L
    )
    
    figure5a_display <- figure5a_display[
        order(
            figure5a_display$NES
        ),
        ,
        drop = FALSE
    ]
    
    figure5a_display$CleanDescription <- stringr::str_to_sentence(
        stringr::str_replace_all(
            figure5a_display$Description,
            "_",
            " "
        )
    )
    
    figure5a_display$CleanDescription <- factor(
        figure5a_display$CleanDescription,
        levels = figure5a_display$CleanDescription
    )
    
    mechanism_colors <- c(
        "Mitochondrial respiration" = "#2C7FB8",
        "Apoptosis" = "#D95F0E",
        "Autophagy" = "#31A354",
        "ER stress" = "#756BB1",
        "Inflammation" = "#E6550D",
        "Protein degradation" = "#636363",
        "Muscle remodeling" = "#E377C2",
        "Oxidative stress" = "#B15928"
    )
    
    p_A <- ggplot2::ggplot(
        figure5a_display,
        ggplot2::aes(
            x = NES,
            y = CleanDescription,
            size = neg_log10_FDR,
            color = Mechanism
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0,
            linewidth = 0.45,
            color = "grey40"
        ) +
        ggplot2::geom_point(
            alpha = 0.92
        ) +
        ggplot2::scale_color_manual(
            values = mechanism_colors,
            drop = FALSE
        ) +
        ggplot2::scale_size_continuous(
            name = expression(
                -log[10](FDR)
            ),
            range = c(
                2.8,
                7.2
            )
        ) +
        ggplot2::labs(
            title = paste0(
                "Pathway-level GSEA of sarcopenia-related ",
                "transcriptomic alterations"
            ),
            subtitle = paste0(
                "Negative NES indicates enrichment in controls and ",
                "relative suppression in sarcopenia"
            ),
            x = "Normalized enrichment score (NES)",
            y = NULL,
            color = "Mechanism"
        ) +
        ggplot2::theme_bw(
            base_size = 12
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = 14,
                hjust = 0
            ),
            plot.subtitle = ggplot2::element_text(
                size = 10.3,
                color = "grey35",
                hjust = 0,
                margin = ggplot2::margin(
                    b = 7
                )
            ),
            axis.text.y = ggplot2::element_text(
                size = 8.9,
                color = "black"
            ),
            axis.text.x = ggplot2::element_text(
                size = 10.2,
                color = "black"
            ),
            axis.title.x = ggplot2::element_text(
                size = 11.3,
                face = "bold"
            ),
            legend.position = "right",
            legend.title = ggplot2::element_text(
                face = "bold"
            ),
            legend.text = ggplot2::element_text(
                size = 8.8
            ),
            panel.grid.major.y = ggplot2::element_line(
                color = "grey90",
                linewidth = 0.25
            ),
            panel.grid.major.x = ggplot2::element_line(
                color = "grey92",
                linewidth = 0.25
            ),
            panel.grid.minor = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(
                8,
                10,
                8,
                8
            )
        )
    
    write.csv(
        mechanism_df,
        file.path(
            figure5_analysis_dir,
            "Figure5A_all_mechanism_related_GSEA_terms.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        figure5a_display,
        file.path(
            figure5_analysis_dir,
            "Figure5A_displayed_GSEA_terms.csv"
        ),
        row.names = FALSE
    )
    
    ############################################################
    ## 9. Figure 5B — sample-level ssGSEA
    ############################################################
    
    cat("\n============================================================\n")
    cat("RUNNING FIGURE 5B ssGSEA\n")
    cat("============================================================\n")
    
    ssgsea_scores <- run_ssgsea(
        expression_matrix = expr,
        gene_sets = pathways
    )
    
    ssgsea_scores <- ssgsea_scores[
        pathway_order,
        group_df$SampleID,
        drop = FALSE
    ]
    
    figure5b_records <- list()
    
    for (pathway_name in pathway_order) {
        
        pathway_scores <- as.numeric(
            ssgsea_scores[
                pathway_name,
                group_df$SampleID
            ]
        )
        
        pathway_sd <- stats::sd(
            pathway_scores,
            na.rm = TRUE
        )
        
        if (
            is.finite(pathway_sd) &&
            pathway_sd > 0
        ) {
            
            score_z <- as.numeric(
                scale(
                    pathway_scores
                )
            )
            
        } else {
            
            score_z <- rep(
                0,
                length(pathway_scores)
            )
        }
        
        figure5b_records <- append(
            figure5b_records,
            list(
                data.frame(
                    SampleID = group_df$SampleID,
                    Group = group_df$Group,
                    Pathway = pathway_name,
                    Score = pathway_scores,
                    Score_z = score_z,
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    figure5b_long <- do.call(
        rbind,
        figure5b_records
    )
    
    figure5b_long$Pathway <- factor(
        figure5b_long$Pathway,
        levels = pathway_order
    )
    
    figure5b_long$Group <- factor(
        figure5b_long$Group,
        levels = c(
            "Control",
            "Sarcopenia"
        )
    )
    
    figure5b_stat_records <- list()
    
    for (pathway_name in pathway_order) {
        
        pathway_data <- figure5b_long[
            as.character(
                figure5b_long$Pathway
            ) ==
                pathway_name,
            ,
            drop = FALSE
        ]
        
        control_scores <- pathway_data$Score[
            pathway_data$Group ==
                "Control"
        ]
        
        sarcopenia_scores <- pathway_data$Score[
            pathway_data$Group ==
                "Sarcopenia"
        ]
        
        wilcoxon_result <- stats::wilcox.test(
            sarcopenia_scores,
            control_scores,
            exact = FALSE
        )
        
        figure5b_stat_records <- append(
            figure5b_stat_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    n_Control = length(
                        control_scores
                    ),
                    n_Sarcopenia = length(
                        sarcopenia_scores
                    ),
                    median_Control = median(
                        control_scores,
                        na.rm = TRUE
                    ),
                    median_Sarcopenia = median(
                        sarcopenia_scores,
                        na.rm = TRUE
                    ),
                    delta_median = median(
                        sarcopenia_scores,
                        na.rm = TRUE
                    ) -
                        median(
                            control_scores,
                            na.rm = TRUE
                        ),
                    p_value = wilcoxon_result$p.value,
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    figure5b_statistics <- do.call(
        rbind,
        figure5b_stat_records
    )
    
    figure5b_statistics$FDR <- stats::p.adjust(
        figure5b_statistics$p_value,
        method = "BH"
    )
    
    figure5b_statistics$direction <- ifelse(
        figure5b_statistics$delta_median < 0,
        "Lower in sarcopenia",
        "Higher in sarcopenia"
    )
    
    figure5b_statistics$label <- ifelse(
        figure5b_statistics$FDR < 0.001,
        "FDR < 0.001",
        ifelse(
            figure5b_statistics$FDR < 0.01,
            "FDR < 0.01",
            ifelse(
                figure5b_statistics$FDR < 0.05,
                "FDR < 0.05",
                paste0(
                    "FDR = ",
                    signif(
                        figure5b_statistics$FDR,
                        2
                    )
                )
            )
        )
    )
    
    figure5b_annotation_records <- list()
    
    for (pathway_name in pathway_order) {
        
        pathway_plot_data <- figure5b_long[
            as.character(
                figure5b_long$Pathway
            ) ==
                pathway_name,
            ,
            drop = FALSE
        ]
        
        pathway_statistics <- figure5b_statistics[
            figure5b_statistics$Pathway ==
                pathway_name,
            ,
            drop = FALSE
        ]
        
        figure5b_annotation_records <- append(
            figure5b_annotation_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    y_position = max(
                        pathway_plot_data$Score_z,
                        na.rm = TRUE
                    ) + 0.55,
                    label = pathway_statistics$label[1],
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    figure5b_annotations <- do.call(
        rbind,
        figure5b_annotation_records
    )
    
    figure5b_annotations$Pathway <- factor(
        figure5b_annotations$Pathway,
        levels = pathway_order
    )
    
    group_colors <- c(
        "Control" = "#4E9BB8",
        "Sarcopenia" = "#D95F4B"
    )
    
    p_B <- ggplot2::ggplot(
        figure5b_long,
        ggplot2::aes(
            x = Group,
            y = Score_z
        )
    ) +
        ggplot2::geom_violin(
            ggplot2::aes(
                fill = Group
            ),
            width = 0.82,
            trim = FALSE,
            alpha = 0.28,
            color = NA
        ) +
        ggplot2::geom_boxplot(
            ggplot2::aes(
                fill = Group
            ),
            width = 0.23,
            outlier.shape = NA,
            alpha = 0.82,
            linewidth = 0.42,
            color = "black"
        ) +
        ggbeeswarm::geom_quasirandom(
            ggplot2::aes(
                color = Group
            ),
            width = 0.22,
            size = 1.8,
            alpha = 0.82,
            varwidth = FALSE
        ) +
        ggplot2::geom_text(
            data = figure5b_annotations,
            mapping = ggplot2::aes(
                x = 1.5,
                y = y_position,
                label = label
            ),
            inherit.aes = FALSE,
            size = 3.5
        ) +
        ggplot2::facet_wrap(
            ~Pathway,
            ncol = 2,
            scales = "free_y",
            labeller = ggplot2::labeller(
                Pathway = ggplot2::label_wrap_gen(
                    width = 28
                )
            )
        ) +
        ggplot2::scale_fill_manual(
            values = group_colors
        ) +
        ggplot2::scale_color_manual(
            values = group_colors
        ) +
        ggplot2::labs(
            title = "Sample-level mitochondrial pathway activity",
            subtitle = "GO-BP-based ssGSEA analysis in GSE111016",
            x = NULL,
            y = "Pathway activity score, z-scaled"
        ) +
        ggplot2::theme_classic(
            base_size = 12.2
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = 15,
                hjust = 0
            ),
            plot.subtitle = ggplot2::element_text(
                size = 10.5,
                color = "grey35",
                hjust = 0,
                margin = ggplot2::margin(
                    b = 8
                )
            ),
            legend.position = "top",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(
                size = 10.5
            ),
            strip.background = ggplot2::element_rect(
                fill = "grey94",
                color = "grey70",
                linewidth = 0.45
            ),
            strip.text = ggplot2::element_text(
                face = "bold",
                size = 10.2,
                margin = ggplot2::margin(
                    t = 5,
                    b = 5
                )
            ),
            axis.title.y = ggplot2::element_text(
                face = "bold",
                size = 11.5
            ),
            axis.text.y = ggplot2::element_text(
                color = "black",
                size = 9.8
            ),
            axis.text.x = ggplot2::element_text(
                color = "black",
                size = 9.8,
                angle = 18,
                hjust = 1
            ),
            panel.spacing = grid::unit(
                0.9,
                "lines"
            ),
            plot.margin = ggplot2::margin(
                8,
                10,
                8,
                8
            )
        ) +
        ggplot2::coord_cartesian(
            clip = "off"
        )
    
    write.csv(
        ssgsea_scores,
        file.path(
            figure5_analysis_dir,
            "Figure5B_four_pathway_ssGSEA_scores.csv"
        )
    )
    
    write.csv(
        figure5b_long,
        file.path(
            figure5_analysis_dir,
            "Figure5B_plot_data.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        figure5b_statistics,
        file.path(
            figure5_analysis_dir,
            "Figure5B_Wilcoxon_BH_statistics.csv"
        ),
        row.names = FALSE
    )
    
    ############################################################
    ## 10. Figure 5C — 11-gene correlation heatmap
    ############################################################
    
    cat("\n============================================================\n")
    cat("RUNNING FIGURE 5C 11-GENE CORRELATION ANALYSIS\n")
    cat("============================================================\n")
    
    pathways_nohub <- lapply(
        pathways,
        function(gene_set) {
            
            setdiff(
                unique(
                    toupper(
                        gene_set
                    )
                ),
                hub_genes
            )
        }
    )
    
    pathways_nohub <- pathways_nohub[
        pathway_order
    ]
    
    figure5c_overlap_records <- list()
    
    for (pathway_name in pathway_order) {
        
        original_set <- unique(
            toupper(
                getElement(
                    pathways,
                    pathway_name
                )
            )
        )
        
        overlap_genes <- intersect(
            original_set,
            hub_genes
        )
        
        nohub_set <- getElement(
            pathways_nohub,
            pathway_name
        )
        
        figure5c_overlap_records <- append(
            figure5c_overlap_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    N_genes_before = length(
                        original_set
                    ),
                    N_overlapping_hub_genes = length(
                        overlap_genes
                    ),
                    Overlapping_hub_genes = paste(
                        overlap_genes,
                        collapse = ";"
                    ),
                    N_genes_after_exclusion = length(
                        nohub_set
                    ),
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    figure5c_overlap <- do.call(
        rbind,
        figure5c_overlap_records
    )
    
    if (
        any(
            figure5c_overlap$N_genes_after_exclusion < 5L
        )
    ) {
        
        stop(
            paste0(
                "At least one pathway contains fewer than five genes ",
                "after removal of all 11 hub genes."
            ),
            call. = FALSE
        )
    }
    
    ssgsea_scores_nohub <- run_ssgsea(
        expression_matrix = expr,
        gene_sets = pathways_nohub
    )
    
    ssgsea_scores_nohub <- ssgsea_scores_nohub[
        pathway_order,
        group_df$SampleID,
        drop = FALSE
    ]
    
    correlation_records <- list()
    
    for (gene_name in hub_genes) {
        
        gene_expression <- as.numeric(
            expr[
                gene_name,
                group_df$SampleID
            ]
        )
        
        for (pathway_name in pathway_order) {
            
            pathway_score <- as.numeric(
                ssgsea_scores_nohub[
                    pathway_name,
                    group_df$SampleID
                ]
            )
            
            complete_rows <- stats::complete.cases(
                gene_expression,
                pathway_score
            )
            
            correlation_test <- suppressWarnings(
                stats::cor.test(
                    gene_expression[
                        complete_rows
                    ],
                    pathway_score[
                        complete_rows
                    ],
                    method = "spearman",
                    exact = FALSE
                )
            )
            
            correlation_records <- append(
                correlation_records,
                list(
                    data.frame(
                        Gene = gene_name,
                        Pathway = pathway_name,
                        rho = unname(
                            correlation_test$estimate
                        ),
                        p_value = correlation_test$p.value,
                        n = sum(
                            complete_rows
                        ),
                        stringsAsFactors = FALSE
                    )
                )
            )
        }
    }
    
    figure5c_correlations <- do.call(
        rbind,
        correlation_records
    )
    
    if (nrow(figure5c_correlations) != 44L) {
        
        stop(
            "Expected 44 correlations, but obtained ",
            nrow(figure5c_correlations),
            ".",
            call. = FALSE
        )
    }
    
    figure5c_correlations$FDR <- stats::p.adjust(
        figure5c_correlations$p_value,
        method = "BH"
    )
    
    figure5c_correlations$significance <- ifelse(
        figure5c_correlations$FDR < 0.001,
        "***",
        ifelse(
            figure5c_correlations$FDR < 0.01,
            "**",
            ifelse(
                figure5c_correlations$FDR < 0.05,
                "*",
                ""
            )
        )
    )
    
    figure5c_correlations$label <- paste0(
        sprintf(
            "%.2f",
            figure5c_correlations$rho
        ),
        figure5c_correlations$significance
    )
    
    figure5c_correlations$Hub_class <- ifelse(
        figure5c_correlations$Gene %in%
            strict_hub_genes,
        "Strict hub",
        "Robust hub"
    )
    
    figure5c_correlations$text_color <- ifelse(
        abs(
            figure5c_correlations$rho
        ) >= 0.55,
        "white",
        "black"
    )
    
    pathway_short_map <- stats::setNames(
        target_pathways$Short,
        target_pathways$Pathway
    )
    
    figure5c_correlations$Pathway_short <- unname(
        pathway_short_map[
            figure5c_correlations$Pathway
        ]
    )
    
    figure5c_correlations$Gene <- factor(
        figure5c_correlations$Gene,
        levels = rev(
            hub_genes
        )
    )
    
    figure5c_correlations$Pathway_short <- factor(
        figure5c_correlations$Pathway_short,
        levels = pathway_short_order
    )
    
    p_C <- ggplot2::ggplot(
        figure5c_correlations,
        ggplot2::aes(
            x = Pathway_short,
            y = Gene,
            fill = rho
        )
    ) +
        ggplot2::geom_tile(
            color = "white",
            linewidth = 0.8
        ) +
        ggplot2::geom_text(
            ggplot2::aes(
                label = label,
                color = text_color
            ),
            size = 3.35
        ) +
        ggplot2::scale_color_identity() +
        ggplot2::scale_fill_gradient2(
            low = "#3B6FB6",
            mid = "white",
            high = "#B73A3A",
            midpoint = 0,
            limits = c(
                -1,
                1
            ),
            breaks = c(
                -1,
                -0.5,
                0,
                0.5,
                1
            ),
            name = "Spearman\nrho"
        ) +
        ggplot2::scale_x_discrete(
            labels = function(x) {
                
                stringr::str_wrap(
                    x,
                    width = 18
                )
            }
        ) +
        ggplot2::labs(
            title = "Hub gene–pathway correlation",
            subtitle = paste0(
                "Eleven-gene analysis after removing all hub genes ",
                "from pathway sets; BH correction across 44 tests"
            ),
            x = NULL,
            y = NULL
        ) +
        ggplot2::theme_minimal(
            base_size = 12.2
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                face = "bold",
                size = 15,
                hjust = 0
            ),
            plot.subtitle = ggplot2::element_text(
                size = 10.2,
                color = "grey35",
                hjust = 0,
                margin = ggplot2::margin(
                    b = 8
                )
            ),
            axis.text.x = ggplot2::element_text(
                angle = 0,
                hjust = 0.5,
                vjust = 0.5,
                color = "black",
                size = 9.4,
                lineheight = 0.92
            ),
            axis.text.y = ggplot2::element_text(
                color = "black",
                size = 10.2
            ),
            panel.grid = ggplot2::element_blank(),
            legend.position = "right",
            legend.title = ggplot2::element_text(
                size = 10
            ),
            legend.text = ggplot2::element_text(
                size = 9.2
            ),
            plot.margin = ggplot2::margin(
                8,
                10,
                8,
                8
            )
        )
    
    write.csv(
        figure5c_overlap,
        file.path(
            figure5_analysis_dir,
            "Figure5C_pathway_hub_overlap_and_sizes.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        ssgsea_scores_nohub,
        file.path(
            figure5_analysis_dir,
            "Figure5C_ssGSEA_scores_after_11hub_exclusion.csv"
        )
    )
    
    write.csv(
        figure5c_correlations,
        file.path(
            figure5_analysis_dir,
            "Figure5C_11gene_44correlation_results.csv"
        ),
        row.names = FALSE
    )
    
    ############################################################
    ## 11. Figure 5D — running enrichment curves
    ############################################################
    
    cat("\n============================================================\n")
    cat("BUILDING FIGURE 5D RUNNING-ENRICHMENT CURVES\n")
    cat("============================================================\n")
    
    get_target_gsea_row <- function(go_id) {
        
        result_rows <- gsea_result[
            gsea_result$ID ==
                go_id,
            ,
            drop = FALSE
        ]
        
        if (nrow(result_rows) == 0L) {
            
            stop(
                "GSEA result not found for ",
                go_id,
                call. = FALSE
            )
        }
        
        result_rows[
            1,
            ,
            drop = FALSE
        ]
    }
    
    target_gsea_records <- list()
    
    for (
        pathway_index in seq_len(
            nrow(target_pathways)
        )
    ) {
        
        current_row <- get_target_gsea_row(
            target_pathways$GO_ID[
                pathway_index
            ]
        )
        
        stored_es <- NA_real_
        
        if (
            "enrichmentScore" %in%
            colnames(current_row)
        ) {
            
            stored_es <- as.numeric(
                current_row$enrichmentScore[1]
            )
        }
        
        stored_set_size <- NA_integer_
        
        if (
            "setSize" %in%
            colnames(current_row)
        ) {
            
            stored_set_size <- as.integer(
                current_row$setSize[1]
            )
        }
        
        target_gsea_records <- append(
            target_gsea_records,
            list(
                data.frame(
                    Pathway = target_pathways$Pathway[
                        pathway_index
                    ],
                    Short = target_pathways$Short[
                        pathway_index
                    ],
                    GO_ID = target_pathways$GO_ID[
                        pathway_index
                    ],
                    NES = as.numeric(
                        current_row$NES[1]
                    ),
                    FDR = as.numeric(
                        current_row$p.adjust[1]
                    ),
                    Stored_ES = stored_es,
                    Stored_setSize = stored_set_size,
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    target_gsea_statistics <- do.call(
        rbind,
        target_gsea_records
    )
    
    calculate_running_es <- function(
        ranked_vector,
        gene_set,
        exponent = 1
    ) {
        
        gene_set <- intersect(
            as.character(gene_set),
            names(ranked_vector)
        )
        
        number_ranked <- length(
            ranked_vector
        )
        
        number_hits <- length(
            gene_set
        )
        
        if (
            number_hits == 0L ||
            number_hits >= number_ranked
        ) {
            
            stop(
                "Invalid gene-set size in running ES calculation.",
                call. = FALSE
            )
        }
        
        hit_indicator <- names(
            ranked_vector
        ) %in%
            gene_set
        
        rank_weights <- abs(
            ranked_vector
        ) ^ exponent
        
        total_hit_weight <- sum(
            rank_weights[
                hit_indicator
            ]
        )
        
        if (
            !is.finite(total_hit_weight) ||
            total_hit_weight <= 0
        ) {
            
            stop(
                "Invalid total hit weight in running ES calculation.",
                call. = FALSE
            )
        }
        
        running_increment <- ifelse(
            hit_indicator,
            rank_weights /
                total_hit_weight,
            -1 /
                (
                    number_ranked -
                        number_hits
                )
        )
        
        running_score <- cumsum(
            running_increment
        )
        
        positive_peak <- max(
            running_score
        )
        
        negative_peak <- min(
            running_score
        )
        
        enrichment_score <- if (
            abs(positive_peak) >=
            abs(negative_peak)
        ) {
            positive_peak
        } else {
            negative_peak
        }
        
        list(
            running_score = running_score,
            hit_positions = which(
                hit_indicator
            ),
            enrichment_score = enrichment_score,
            number_hits = number_hits
        )
    }
    
    rank_table <- data.frame(
        Rank = seq_along(
            rank_vec
        ),
        ENTREZID = names(
            rank_vec
        ),
        rank_metric = as.numeric(
            rank_vec
        ),
        stringsAsFactors = FALSE
    )
    
    figure5d_es_records <- list()
    figure5d_hit_records <- list()
    figure5d_qc_records <- list()
    
    for (
        pathway_index in seq_len(
            nrow(target_pathways)
        )
    ) {
        
        pathway_name <- target_pathways$Pathway[
            pathway_index
        ]
        
        short_name <- target_pathways$Short[
            pathway_index
        ]
        
        current_gene_set <- getElement(
            pathway_entrez,
            pathway_name
        )
        
        running_result <- calculate_running_es(
            ranked_vector = rank_vec,
            gene_set = current_gene_set,
            exponent = 1
        )
        
        figure5d_es_records <- append(
            figure5d_es_records,
            list(
                data.frame(
                    Rank = seq_along(
                        rank_vec
                    ),
                    ES = getElement(
                        running_result,
                        "running_score"
                    ),
                    Pathway = pathway_name,
                    Short = short_name,
                    stringsAsFactors = FALSE
                )
            )
        )
        
        hit_positions <- getElement(
            running_result,
            "hit_positions"
        )
        
        figure5d_hit_records <- append(
            figure5d_hit_records,
            list(
                data.frame(
                    Rank = hit_positions,
                    ENTREZID = names(
                        rank_vec
                    )[
                        hit_positions
                    ],
                    rank_metric = as.numeric(
                        rank_vec[
                            hit_positions
                        ]
                    ),
                    Pathway = pathway_name,
                    Short = short_name,
                    stringsAsFactors = FALSE
                )
            )
        )
        
        stored_statistics <- target_gsea_statistics[
            target_gsea_statistics$Pathway ==
                pathway_name,
            ,
            drop = FALSE
        ]
        
        figure5d_qc_records <- append(
            figure5d_qc_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    Manual_ES = getElement(
                        running_result,
                        "enrichment_score"
                    ),
                    Stored_ES = stored_statistics$Stored_ES[1],
                    Difference = getElement(
                        running_result,
                        "enrichment_score"
                    ) -
                        stored_statistics$Stored_ES[1],
                    Manual_detected_size = getElement(
                        running_result,
                        "number_hits"
                    ),
                    Stored_setSize = stored_statistics$Stored_setSize[1],
                    stringsAsFactors = FALSE
                )
            )
        )
    }
    
    figure5d_es <- do.call(
        rbind,
        figure5d_es_records
    )
    
    figure5d_hits <- do.call(
        rbind,
        figure5d_hit_records
    )
    
    figure5d_qc <- do.call(
        rbind,
        figure5d_qc_records
    )
    
    figure5d_es$Pathway <- factor(
        figure5d_es$Pathway,
        levels = pathway_order
    )
    
    figure5d_hits$Pathway <- factor(
        figure5d_hits$Pathway,
        levels = pathway_order
    )
    
    pathway_colors <- c(
        "Oxidative phosphorylation" = "#4EC21A",
        "Respiratory electron transport chain" = "#2C7FB8",
        "Mitochondrial ATP synthesis" = "#F28E2B",
        "Mitochondrial translation" = "#8E63CE"
    )
    
    legend_labels <- stats::setNames(
        paste0(
            target_gsea_statistics$Short,
            " (NES=",
            sprintf(
                "%.2f",
                target_gsea_statistics$NES
            ),
            ", q",
            ifelse(
                target_gsea_statistics$FDR < 0.001,
                "<0.001",
                paste0(
                    "=",
                    signif(
                        target_gsea_statistics$FDR,
                        2
                    )
                )
            ),
            ")"
        ),
        target_gsea_statistics$Pathway
    )
    
    track_table <- data.frame(
        Pathway = pathway_order,
        Short = pathway_short_order,
        y = rev(
            seq_along(
                pathway_order
            )
        ),
        stringsAsFactors = FALSE
    )
    
    figure5d_hits <- merge(
        figure5d_hits,
        track_table,
        by = c(
            "Pathway",
            "Short"
        ),
        all.x = TRUE,
        sort = FALSE
    )
    
    rank_metric_table <- data.frame(
        Rank = rank_table$Rank,
        rank_metric = rank_table$rank_metric,
        metric_scaled = rank_table$rank_metric /
            max(
                abs(
                    rank_table$rank_metric
                ),
                na.rm = TRUE
            ),
        stringsAsFactors = FALSE
    )
    
    zero_cross <- which.min(
        abs(
            rank_table$rank_metric
        )
    )
    
    number_bins <- 220L
    
    strip_x <- seq(
        from = 1,
        to = length(rank_vec),
        length.out = number_bins
    )
    
    phenotype_strip <- data.frame(
        xmin = strip_x,
        xmax = c(
            strip_x[-1],
            length(rank_vec)
        ),
        value = seq(
            from = 1,
            to = -1,
            length.out = number_bins
        ),
        stringsAsFactors = FALSE
    )
    
    sarcopenia_color <- "#D73027"
    control_color <- "#225EA8"
    
    es_upper <- max(
        0.08,
        max(
            figure5d_es$ES,
            na.rm = TRUE
        ) + 0.06
    )
    
    es_lower <- min(
        -0.95,
        min(
            figure5d_es$ES,
            na.rm = TRUE
        ) - 0.08
    )
    
    p_D_es <- ggplot2::ggplot(
        figure5d_es,
        ggplot2::aes(
            x = Rank,
            y = ES,
            color = Pathway
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0,
            linetype = "dashed",
            linewidth = 0.32,
            color = "grey55"
        ) +
        ggplot2::geom_line(
            linewidth = 0.82,
            alpha = 0.96
        ) +
        ggplot2::scale_color_manual(
            values = pathway_colors,
            breaks = pathway_order,
            labels = legend_labels[
                pathway_order
            ],
            name = NULL
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Enrichment score"
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::comma,
            expand = ggplot2::expansion(
                mult = c(
                    0,
                    0.012
                )
            )
        ) +
        ggplot2::coord_cartesian(
            ylim = c(
                es_lower,
                es_upper
            ),
            clip = "off"
        ) +
        ggplot2::theme_classic(
            base_size = 11.2
        ) +
        ggplot2::theme(
            legend.position = c(
                0.02,
                0.045
            ),
            legend.justification = c(
                0,
                0
            ),
            legend.text = ggplot2::element_text(
                size = 7,
                color = "black",
                lineheight = 0.82
            ),
            legend.background = ggplot2::element_rect(
                fill = "white",
                color = "grey65",
                linewidth = 0.22
            ),
            legend.key.height = grid::unit(
                0.33,
                "lines"
            ),
            legend.key.width = grid::unit(
                1.15,
                "lines"
            ),
            legend.spacing.y = grid::unit(
                0,
                "cm"
            ),
            legend.margin = ggplot2::margin(
                3,
                5,
                3,
                5
            ),
            axis.title.y = ggplot2::element_text(
                face = "bold",
                size = 11.5,
                margin = ggplot2::margin(
                    r = 8
                )
            ),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(
                color = "black",
                size = 9.4
            ),
            plot.margin = ggplot2::margin(
                2,
                10,
                0,
                8
            )
        )
    
    p_D_hits <- ggplot2::ggplot() +
        ggplot2::geom_segment(
            data = figure5d_hits,
            mapping = ggplot2::aes(
                x = Rank,
                xend = Rank,
                y = y - 0.30,
                yend = y + 0.30,
                color = Pathway
            ),
            linewidth = 0.20,
            alpha = 0.78
        ) +
        ggplot2::geom_rect(
            data = phenotype_strip,
            mapping = ggplot2::aes(
                xmin = xmin,
                xmax = xmax,
                ymin = 0.18,
                ymax = 0.42,
                fill = value
            ),
            color = NA
        ) +
        ggplot2::annotate(
            "text",
            x = length(rank_vec) * 0.06,
            y = -0.18,
            label = "Up-ranked in sarcopenia",
            color = sarcopenia_color,
            hjust = 0,
            size = 3,
            fontface = "bold"
        ) +
        ggplot2::annotate(
            "text",
            x = length(rank_vec) * 0.94,
            y = -0.18,
            label = "Up-ranked in control",
            color = control_color,
            hjust = 1,
            size = 3,
            fontface = "bold"
        ) +
        ggplot2::scale_color_manual(
            values = pathway_colors,
            guide = "none"
        ) +
        ggplot2::scale_fill_gradient2(
            low = control_color,
            mid = "white",
            high = sarcopenia_color,
            midpoint = 0,
            guide = "none"
        ) +
        ggplot2::scale_y_continuous(
            breaks = track_table$y,
            labels = track_table$Short,
            limits = c(
                -0.35,
                max(
                    track_table$y
                ) + 0.55
            ),
            expand = c(
                0,
                0
            )
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::comma,
            expand = ggplot2::expansion(
                mult = c(
                    0,
                    0.012
                )
            )
        ) +
        ggplot2::coord_cartesian(
            clip = "off"
        ) +
        ggplot2::theme_classic(
            base_size = 10.5
        ) +
        ggplot2::theme(
            axis.title = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(
                color = "black",
                size = 8,
                face = "bold",
                margin = ggplot2::margin(
                    r = 5
                )
            ),
            plot.margin = ggplot2::margin(
                0,
                10,
                3,
                8
            )
        )
    
    p_D_metric <- ggplot2::ggplot(
        rank_metric_table,
        ggplot2::aes(
            x = Rank,
            y = metric_scaled
        )
    ) +
        ggplot2::geom_area(
            fill = "grey75",
            alpha = 0.95
        ) +
        ggplot2::geom_hline(
            yintercept = 0,
            color = "grey55",
            linewidth = 0.32
        ) +
        ggplot2::geom_vline(
            xintercept = zero_cross,
            linetype = "dashed",
            color = "grey65",
            linewidth = 0.32
        ) +
        ggplot2::annotate(
            "text",
            x = zero_cross,
            y = 0.20,
            label = paste0(
                "Zero cross at ",
                scales::comma(
                    zero_cross
                )
            ),
            hjust = 0.5,
            size = 2.8,
            color = "black"
        ) +
        ggplot2::labs(
            x = "Rank in ordered dataset",
            y = "Ranked list\nmetric"
        ) +
        ggplot2::scale_x_continuous(
            labels = scales::comma,
            expand = ggplot2::expansion(
                mult = c(
                    0,
                    0.012
                )
            )
        ) +
        ggplot2::coord_cartesian(
            ylim = c(
                -1.05,
                1.05
            ),
            clip = "off"
        ) +
        ggplot2::theme_classic(
            base_size = 11.2
        ) +
        ggplot2::theme(
            axis.title.x = ggplot2::element_text(
                face = "bold",
                size = 11.5,
                margin = ggplot2::margin(
                    t = 8
                )
            ),
            axis.title.y = ggplot2::element_text(
                face = "bold",
                size = 10.8,
                margin = ggplot2::margin(
                    r = 8
                )
            ),
            axis.text = ggplot2::element_text(
                color = "black",
                size = 9.3
            ),
            plot.margin = ggplot2::margin(
                8,
                10,
                8,
                8
            )
        )
    
    p_D_title <- ggplot2::ggplot() +
        ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.72,
            label = "Mitochondrial bioenergetic pathways",
            fontface = "bold",
            size = 5.5,
            hjust = 0.5
        ) +
        ggplot2::annotate(
            "text",
            x = 0.5,
            y = 0.24,
            label = paste0(
                "Representative running-enrichment curves for ",
                "four mitochondrial pathways"
            ),
            size = 3.4,
            color = "grey35",
            hjust = 0.5
        ) +
        ggplot2::xlim(
            0,
            1
        ) +
        ggplot2::ylim(
            0,
            1
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(
                0,
                10,
                0,
                8
            )
        )
    
    p_D <- p_D_title /
        p_D_es /
        p_D_hits /
        p_D_metric +
        patchwork::plot_layout(
            heights = c(
                0.11,
                0.45,
                0.25,
                0.32
            )
        )
    
    write.csv(
        target_gsea_statistics,
        file.path(
            figure5_analysis_dir,
            "Figure5D_target_pathway_GSEA_statistics.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        figure5d_es,
        file.path(
            figure5_analysis_dir,
            "Figure5D_running_enrichment_scores.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        figure5d_hits,
        file.path(
            figure5_analysis_dir,
            "Figure5D_gene_set_hit_positions.csv"
        ),
        row.names = FALSE
    )
    
    write.csv(
        figure5d_qc,
        file.path(
            figure5_analysis_dir,
            "Figure5D_manual_ES_vs_stored_ES_QC.csv"
        ),
        row.names = FALSE
    )
    
    ############################################################
    ## 12. Combine Figure 5A–D
    ############################################################
    
    combined_plot <- (
        p_A |
            p_B
    ) /
        (
            p_C |
                patchwork::wrap_elements(
                    full = p_D
                )
        ) +
        patchwork::plot_layout(
            widths = c(
                1.03,
                0.97
            ),
            heights = c(
                1,
                1
            )
        ) +
        patchwork::plot_annotation(
            tag_levels = "A",
            theme = ggplot2::theme(
                plot.tag = ggplot2::element_text(
                    face = "bold",
                    size = 20,
                    color = "black"
                ),
                plot.tag.position = c(
                    0,
                    1
                )
            )
        )
    
    ############################################################
    ## 13. Save helper
    ############################################################
    
    save_png_pdf <- function(
        plot_object,
        file_stem,
        width,
        height
    ) {
        
        png_path <- file.path(
            figure_dir,
            paste0(
                file_stem,
                "_600dpi.png"
            )
        )
        
        pdf_path <- file.path(
            figure_dir,
            paste0(
                file_stem,
                ".pdf"
            )
        )
        
        ggplot2::ggsave(
            filename = png_path,
            plot = plot_object,
            width = width,
            height = height,
            units = "in",
            dpi = 600,
            bg = "white",
            limitsize = FALSE
        )
        
        if (capabilities("cairo")) {
            
            ggplot2::ggsave(
                filename = pdf_path,
                plot = plot_object,
                width = width,
                height = height,
                units = "in",
                bg = "white",
                limitsize = FALSE,
                device = grDevices::cairo_pdf
            )
            
        } else {
            
            ggplot2::ggsave(
                filename = pdf_path,
                plot = plot_object,
                width = width,
                height = height,
                units = "in",
                bg = "white",
                limitsize = FALSE,
                device = grDevices::pdf
            )
        }
        
        required_paths <- c(
            png_path,
            pdf_path
        )
        
        required_sizes <- file.info(
            required_paths
        )$size
        
        if (
            any(!file.exists(required_paths)) ||
            any(is.na(required_sizes)) ||
            any(required_sizes <= 0)
        ) {
            
            stop(
                "Figure output failed for ",
                file_stem,
                ".",
                call. = FALSE
            )
        }
        
        list(
            png = png_path,
            pdf = pdf_path
        )
    }
    
    ############################################################
    ## 14. Save figures
    ############################################################
    
    cat("\n============================================================\n")
    cat("SAVING FIGURE OUTPUTS\n")
    cat("============================================================\n")
    
    figure5a_files <- save_png_pdf(
        plot_object = p_A,
        file_stem = "Figure5A_GSEA_mechanism_dotplot_FINAL",
        width = 10.5,
        height = 6.8
    )
    
    figure5b_files <- save_png_pdf(
        plot_object = p_B,
        file_stem = "Figure5B_mitochondrial_ssGSEA_FINAL",
        width = 8.2,
        height = 6.4
    )
    
    figure5c_files <- save_png_pdf(
        plot_object = p_C,
        file_stem = paste0(
            "Figure5C_11gene_hub_pathway_",
            "correlation_FINAL"
        ),
        width = 8.6,
        height = 6.2
    )
    
    figure5d_files <- save_png_pdf(
        plot_object = p_D,
        file_stem = "Figure5D_mitochondrial_GSEA_curves_FINAL",
        width = 7.8,
        height = 5.45
    )
    
    combined_files <- save_png_pdf(
        plot_object = combined_plot,
        file_stem = "Figure5_ABCD_combined_FINAL",
        width = 17,
        height = 12.4
    )
    
    ############################################################
    ## 15. Freeze canonical Figure 5 core as RData
    ############################################################
    
    canonical_core_path <- file.path(
        canonical_core_dir,
        "M05_Figure5_CANONICAL_CORE_FINAL.RData"
    )
    
    source_bc_input_path <- bc_input_path
    source_bc_input_md5 <- unname(
        tools::md5sum(
            bc_input_path
        )
    )
    
    source_gsea_input_path <- gsea_path
    source_gsea_input_md5 <- unname(
        tools::md5sum(
            gsea_path
        )
    )
    
    save(
        expr,
        group_df,
        rank_vec,
        gsea_result,
        pathway_entrez,
        pathway_symbols_raw,
        pathways,
        pathway_size_table,
        hub_genes,
        strict_hub_genes,
        robust_hub_genes,
        figure5a_display,
        ssgsea_scores,
        figure5b_statistics,
        pathways_nohub,
        ssgsea_scores_nohub,
        figure5c_correlations,
        target_gsea_statistics,
        figure5d_es,
        figure5d_hits,
        figure5d_qc,
        source_bc_input_path,
        source_bc_input_md5,
        source_gsea_input_path,
        source_gsea_input_md5,
        file = canonical_core_path,
        compress = "xz"
    )
    
    if (!file.exists(canonical_core_path)) {
        
        stop(
            "Canonical Figure 5 core was not created.",
            call. = FALSE
        )
    }
    
    ############################################################
    ## 16. Export this complete function as code TXT
    ############################################################
    
    function_code <- deparse(
        run_figure5_abcd_final,
        control = "all",
        width.cutoff = 500L
    )
    
    function_code[1] <- paste0(
        "run_figure5_abcd_final <- ",
        function_code[1]
    )
    
    complete_code <- c(
        "############################################################",
        "## Figure 5A–D complete reconstruction — FINAL FIXED",
        "## Self-contained code exported automatically",
        "############################################################",
        "",
        function_code,
        "",
        "figure5_result <- run_figure5_abcd_final()",
        ""
    )
    
    code_txt_paths <- c(
        file.path(
            figure_dir,
            "Figure5A_GSEA_mechanism_dotplot_FINAL_code.txt"
        ),
        
        file.path(
            figure_dir,
            "Figure5B_mitochondrial_ssGSEA_FINAL_code.txt"
        ),
        
        file.path(
            figure_dir,
            paste0(
                "Figure5C_11gene_hub_pathway_",
                "correlation_FINAL_code.txt"
            )
        ),
        
        file.path(
            figure_dir,
            "Figure5D_mitochondrial_GSEA_curves_FINAL_code.txt"
        ),
        
        file.path(
            figure_dir,
            "Figure5_ABCD_combined_FINAL_code.txt"
        )
    )
    
    for (code_txt_path in code_txt_paths) {
        
        writeLines(
            complete_code,
            con = code_txt_path,
            useBytes = TRUE
        )
    }
    
    code_txt_exists <- file.exists(
        code_txt_paths
    )
    
    if (!all(code_txt_exists)) {
        
        stop(
            "One or more code TXT files were not created.",
            call. = FALSE
        )
    }
    
    ############################################################
    ## 17. Save session information
    ############################################################
    
    session_info_path <- file.path(
        figure5_analysis_dir,
        "Figure5_rebuild_sessionInfo.txt"
    )
    
    capture.output(
        sessionInfo(),
        file = session_info_path
    )
    
    ############################################################
    ## 18. Output manifest
    ############################################################
    
    figure_output_paths <- c(
        getElement(
            figure5a_files,
            "png"
        ),
        getElement(
            figure5a_files,
            "pdf"
        ),
        getElement(
            figure5b_files,
            "png"
        ),
        getElement(
            figure5b_files,
            "pdf"
        ),
        getElement(
            figure5c_files,
            "png"
        ),
        getElement(
            figure5c_files,
            "pdf"
        ),
        getElement(
            figure5d_files,
            "png"
        ),
        getElement(
            figure5d_files,
            "pdf"
        ),
        getElement(
            combined_files,
            "png"
        ),
        getElement(
            combined_files,
            "pdf"
        ),
        code_txt_paths,
        canonical_core_path,
        session_info_path
    )
    
    figure_output_roles <- c(
        "Figure 5A PNG",
        "Figure 5A PDF",
        "Figure 5B PNG",
        "Figure 5B PDF",
        "Figure 5C PNG",
        "Figure 5C PDF",
        "Figure 5D PNG",
        "Figure 5D PDF",
        "Figure 5 combined PNG",
        "Figure 5 combined PDF",
        "Figure 5A code TXT",
        "Figure 5B code TXT",
        "Figure 5C code TXT",
        "Figure 5D code TXT",
        "Figure 5 combined code TXT",
        "Canonical Figure 5 core",
        "Session information"
    )
    
    output_manifest <- data.frame(
        Role = figure_output_roles,
        Path = figure_output_paths,
        Exists = file.exists(
            figure_output_paths
        ),
        Size_bytes = file.info(
            figure_output_paths
        )$size,
        MD5 = unname(
            tools::md5sum(
                figure_output_paths
            )
        ),
        stringsAsFactors = FALSE
    )
    
    manifest_path <- file.path(
        figure5_analysis_dir,
        "Figure5_final_output_manifest.csv"
    )
    
    write.csv(
        output_manifest,
        manifest_path,
        row.names = FALSE
    )
    
    ############################################################
    ## 19. Decision report
    ############################################################
    
    decision_path <- file.path(
        figure5_analysis_dir,
        "Figure5_rebuild_decision_report.txt"
    )
    
    writeLines(
        c(
            "Figure 5A–D complete reconstruction — FINAL FIXED",
            
            "============================================================",
            
            paste0(
                "Run time: ",
                Sys.time()
            ),
            
            paste0(
                "Canonical B/C input: ",
                bc_input_path
            ),
            
            paste0(
                "Canonical B/C MD5: ",
                unname(
                    tools::md5sum(
                        bc_input_path
                    )
                )
            ),
            
            paste0(
                "Canonical GSEA input: ",
                gsea_path
            ),
            
            paste0(
                "Canonical GSEA MD5: ",
                unname(
                    tools::md5sum(
                        gsea_path
                    )
                )
            ),
            
            paste0(
                "Expression dimensions: ",
                nrow(expr),
                " genes x ",
                ncol(expr),
                " samples"
            ),
            
            "Figure 5A: stored GSEA result; no GSEA rerun",
            
            paste0(
                "Figure 5B: ssGSEA for four pathways; ",
                "Wilcoxon tests; BH correction over four pathways"
            ),
            
            paste0(
                "Figure 5C: 11 hub genes; all 11 removed ",
                "simultaneously from each pathway set"
            ),
            
            "Figure 5C: 44 Spearman correlations",
            
            "Figure 5C: BH correction across all 44 tests",
            
            paste0(
                "Figure 5D: running ES calculated from stored ",
                "geneList and stored geneSets; GSEA not rerun"
            ),
            
            paste0(
                "Canonical Figure 5 core: ",
                canonical_core_path
            ),
            
            paste0(
                "Final figure directory: ",
                figure_dir
            ),
            
            "Figure 5E generated: FALSE",
            
            "Active RStudio source tab required: FALSE",
            
            "Backup modified: FALSE"
        ),
        
        con = decision_path,
        useBytes = TRUE
    )
    
    ############################################################
    ## 20. Final validation
    ############################################################
    
    required_final_files <- c(
        figure_output_paths,
        manifest_path,
        decision_path
    )
    
    required_exists <- file.exists(
        required_final_files
    )
    
    required_sizes <- file.info(
        required_final_files
    )$size
    
    failed_files <- required_final_files[
        !required_exists |
            is.na(required_sizes) |
            required_sizes <= 0
    ]
    
    if (length(failed_files) > 0L) {
        
        stop(
            paste0(
                "One or more final files are missing or empty:\n",
                paste(
                    failed_files,
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }
    
    ############################################################
    ## 21. Final report
    ############################################################
    
    cat("\n============================================================\n")
    cat("FIGURE 5A–D GENERATED SUCCESSFULLY\n")
    cat("============================================================\n")
    
    cat("\nExpression dimensions:\n")
    print(dim(expr))
    
    cat("\nGroup table:\n")
    print(table(group_df$Group))
    
    cat("\nFigure 5B statistics:\n")
    
    print(
        figure5b_statistics,
        row.names = FALSE
    )
    
    cat("\nFigure 5C correlation count:\n")
    
    cat(
        nrow(
            figure5c_correlations
        ),
        "\n"
    )
    
    cat("\nFigure 5C FDR-significant correlations:\n")
    
    significant_correlations <- figure5c_correlations[
        figure5c_correlations$FDR < 0.05,
        c(
            "Gene",
            "Pathway",
            "rho",
            "p_value",
            "FDR"
        ),
        drop = FALSE
    ]
    
    print(
        significant_correlations,
        row.names = FALSE
    )
    
    cat("\nFigure 5D target statistics:\n")
    
    print(
        target_gsea_statistics,
        row.names = FALSE
    )
    
    cat("\nOutput files:\n")
    
    print(
        output_manifest[
            ,
            c(
                "Role",
                "Path",
                "Exists"
            ),
            drop = FALSE
        ],
        row.names = FALSE
    )
    
    cat("\nCanonical Figure 5 core:\n")
    cat(canonical_core_path, "\n")
    
    cat("\nCanonical core exists: ")
    cat(file.exists(canonical_core_path), "\n")
    
    cat("\nAll five code TXT files exist:\n")
    
    print(
        stats::setNames(
            file.exists(
                code_txt_paths
            ),
            basename(
                code_txt_paths
            )
        )
    )
    
    cat("\nFinal combined PNG:\n")
    
    cat(
        getElement(
            combined_files,
            "png"
        ),
        "\n"
    )
    
    cat("\nFinal combined PDF:\n")
    
    cat(
        getElement(
            combined_files,
            "pdf"
        ),
        "\n"
    )
    
    cat("\nActive RStudio tab was not used.\n")
    cat("No Figure 5E was generated.\n")
    cat("No files were written to backup.\n")
    cat("============================================================\n")
    
    invisible(
        list(
            plot_A = p_A,
            plot_B = p_B,
            plot_C = p_C,
            plot_D = p_D,
            combined_plot = combined_plot,
            figure5B_statistics = figure5b_statistics,
            figure5C_correlations = figure5c_correlations,
            figure5D_statistics = target_gsea_statistics,
            canonical_core = canonical_core_path,
            combined_png = getElement(
                combined_files,
                "png"
            ),
            combined_pdf = getElement(
                combined_files,
                "pdf"
            )
        )
    )
}
