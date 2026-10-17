rebuild_figure5_bcd <- function() {
    
    options(stringsAsFactors = FALSE)
    
    ############################################################
    ## 0. Paths
    ############################################################
    
    project_root <- "E:/1Data/Cd_Sarcopenia_NHANES"
    
    m05_root <- file.path(
        project_root,
        "new_workspace",
        "03_results",
        "M05_PathwayValidation_GSEA_ssGSEA_Correlation"
    )
    
    core_path <- file.path(
        m05_root,
        "02_analysis_files",
        "canonical_core",
        "M05_Figure5_CANONICAL_CORE_FINAL.RData"
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
        figure_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    if (!file.exists(core_path)) {
        stop(
            "Canonical Figure 5 core was not found:\n",
            core_path,
            call. = FALSE
        )
    }
    
    cat("\n============================================================\n")
    cat("FIGURE 5B–D PUBLICATION-STYLE REVISION\n")
    cat("============================================================\n")
    
    cat("\nCanonical core:\n")
    cat(core_path, "\n")
    
    cat("\nFigure output directory:\n")
    cat(figure_dir, "\n")
    
    cat("\nBackup will not be modified:\n")
    cat(backup_dir, "\n")
    
    ############################################################
    ## 1. Packages
    ############################################################
    
    required_packages <- c(
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
    ## 2. Load frozen canonical Figure 5 objects
    ############################################################
    
    core_environment <- new.env(
        parent = emptyenv()
    )
    
    loaded_objects <- load(
        core_path,
        envir = core_environment
    )
    
    required_objects <- c(
        "group_df",
        "hub_genes",
        "strict_hub_genes",
        "ssgsea_scores",
        "figure5b_statistics",
        "figure5c_correlations",
        "target_gsea_statistics",
        "figure5d_es",
        "figure5d_hits",
        "rank_vec"
    )
    
    missing_objects <- required_objects[
        !vapply(
            required_objects,
            exists,
            logical(1),
            envir = core_environment,
            inherits = FALSE
        )
    ]
    
    if (length(missing_objects) > 0L) {
        stop(
            "Canonical core lacks required object(s): ",
            paste(
                missing_objects,
                collapse = ", "
            ),
            call. = FALSE
        )
    }
    
    group_df <- as.data.frame(
        get(
            "group_df",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    hub_genes <- as.character(
        get(
            "hub_genes",
            envir = core_environment,
            inherits = FALSE
        )
    )
    
    strict_hub_genes <- as.character(
        get(
            "strict_hub_genes",
            envir = core_environment,
            inherits = FALSE
        )
    )
    
    ssgsea_scores <- as.matrix(
        get(
            "ssgsea_scores",
            envir = core_environment,
            inherits = FALSE
        )
    )
    
    figure5b_statistics <- as.data.frame(
        get(
            "figure5b_statistics",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    figure5c_correlations <- as.data.frame(
        get(
            "figure5c_correlations",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    target_gsea_statistics <- as.data.frame(
        get(
            "target_gsea_statistics",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    figure5d_es <- as.data.frame(
        get(
            "figure5d_es",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    figure5d_hits <- as.data.frame(
        get(
            "figure5d_hits",
            envir = core_environment,
            inherits = FALSE
        ),
        stringsAsFactors = FALSE
    )
    
    rank_vec <- get(
        "rank_vec",
        envir = core_environment,
        inherits = FALSE
    )
    
    rm(core_environment)
    invisible(gc())
    
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
    
    stopifnot(
        nrow(group_df) == 40L,
        ncol(ssgsea_scores) == 40L,
        identical(
            colnames(ssgsea_scores),
            group_df$SampleID
        ),
        sum(group_df$Group == "Control") == 20L,
        sum(group_df$Group == "Sarcopenia") == 20L,
        length(hub_genes) == 11L,
        length(strict_hub_genes) == 2L,
        nrow(figure5c_correlations) == 44L
    )
    
    ############################################################
    ## 3. Shared definitions and helpers
    ############################################################
    
    pathway_order <- c(
        "Oxidative phosphorylation",
        "Respiratory electron transport chain",
        "Mitochondrial ATP synthesis",
        "Mitochondrial translation"
    )
    
    pathway_short_order <- c(
        "Oxidative phosphorylation",
        "Respiratory ETC",
        "Mitochondrial ATP synthesis",
        "Mitochondrial translation"
    )
    
    pathway_short_map <- stats::setNames(
        pathway_short_order,
        pathway_order
    )
    
    group_colors <- c(
        "Control" = "#4E9BB8",
        "Sarcopenia" = "#D95F4B"
    )
    
    pathway_colors <- c(
        "Oxidative phosphorylation" = "#4EC21A",
        "Respiratory electron transport chain" = "#2C7FB8",
        "Mitochondrial ATP synthesis" = "#F28E2B",
        "Mitochondrial translation" = "#8E63CE"
    )
    
    sarcopenia_color <- "#D73027"
    control_color <- "#225EA8"
    
    format_fdr <- function(x) {
        
        x <- as.numeric(x)
        
        if (
            length(x) == 0L ||
            is.na(x) ||
            !is.finite(x)
        ) {
            return("FDR = NA")
        }
        
        if (x < 0.001) {
            return("FDR < 0.001")
        }
        
        if (x < 0.10) {
            return(
                sprintf(
                    "FDR = %.3f",
                    x
                )
            )
        }
        
        sprintf(
            "FDR = %.2f",
            x
        )
    }
    
    format_nes <- function(x) {
        
        formatted <- sprintf(
            "%.2f",
            as.numeric(x)
        )
        
        gsub(
            "-",
            "\u2212",
            formatted,
            fixed = TRUE
        )
    }
    
    save_png <- function(
        plot_object,
        filename,
        width,
        height
    ) {
        
        output_path <- file.path(
            figure_dir,
            filename
        )
        
        ggplot2::ggsave(
            filename = output_path,
            plot = plot_object,
            width = width,
            height = height,
            units = "in",
            dpi = 600,
            bg = "white",
            limitsize = FALSE
        )
        
        output_size <- file.info(
            output_path
        )$size
        
        if (
            !file.exists(output_path) ||
            is.na(output_size) ||
            output_size <= 0
        ) {
            stop(
                "PNG output failed:\n",
                output_path,
                call. = FALSE
            )
        }
        
        output_path
    }
    
    ############################################################
    ## 4. Figure 5B — revised sample-level ssGSEA plot
    ############################################################
    
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
    
    figure5b_statistics$Pathway <- as.character(
        figure5b_statistics$Pathway
    )
    
    figure5b_statistics$FDR_label <- vapply(
        figure5b_statistics$FDR,
        format_fdr,
        character(1)
    )
    
    figure5b_annotation_records <- list()
    
    for (pathway_name in pathway_order) {
        
        pathway_plot_data <- figure5b_long[
            as.character(
                figure5b_long$Pathway
            ) == pathway_name,
            ,
            drop = FALSE
        ]
        
        statistics_row <- figure5b_statistics[
            figure5b_statistics$Pathway == pathway_name,
            ,
            drop = FALSE
        ]
        
        if (nrow(statistics_row) != 1L) {
            stop(
                "Figure 5B statistics are not unique for pathway: ",
                pathway_name,
                call. = FALSE
            )
        }
        
        figure5b_annotation_records <- append(
            figure5b_annotation_records,
            list(
                data.frame(
                    Pathway = pathway_name,
                    y_position = max(
                        pathway_plot_data$Score_z,
                        na.rm = TRUE
                    ) + 0.52,
                    label = statistics_row$FDR_label[1],
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
    
    p_B_base <- ggplot2::ggplot(
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
            alpha = 0.24,
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
            size = 3.45
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
        ggplot2::scale_y_continuous(
            expand = ggplot2::expansion(
                mult = c(
                    0.08,
                    0.20
                )
            )
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Pathway activity score, z-scaled"
        ) +
        ggplot2::theme_classic(
            base_size = 12.2
        ) +
        ggplot2::theme(
            legend.position = "top",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(
                size = 10.4
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
                angle = 0,
                hjust = 0.5
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
    
    p_B_with_title <- p_B_base +
        ggplot2::labs(
            title = "Sample-level mitochondrial pathway activity",
            subtitle = "GO-BP-based ssGSEA analysis in GSE111016"
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
            )
        )
    
    p_B_no_title <- p_B_base +
        ggplot2::labs(
            title = NULL,
            subtitle = NULL
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(
                4,
                10,
                8,
                8
            )
        )
    
    figure5b_with_title_path <- save_png(
        plot_object = p_B_with_title,
        filename = "Figure5B_mitochondrial_ssGSEA_with_title_600dpi.png",
        width = 8.2,
        height = 6.4
    )
    
    figure5b_no_title_path <- save_png(
        plot_object = p_B_no_title,
        filename = "Figure5B_mitochondrial_ssGSEA_no_title_600dpi.png",
        width = 8.2,
        height = 5.75
    )
    
    ############################################################
    ## 5. Figure 5C — revised correlation heatmap
    ############################################################
    
    figure5c_correlations$Gene <- as.character(
        figure5c_correlations$Gene
    )
    
    figure5c_correlations$Pathway <- as.character(
        figure5c_correlations$Pathway
    )
    
    figure5c_correlations$rho <- as.numeric(
        figure5c_correlations$rho
    )
    
    figure5c_correlations$FDR <- as.numeric(
        figure5c_correlations$FDR
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
    
    figure5c_correlations$text_color <- ifelse(
        abs(
            figure5c_correlations$rho
        ) >= 0.55,
        "white",
        "black"
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
    
    gene_label_map <- stats::setNames(
        hub_genes,
        hub_genes
    )
    
    gene_label_map[
        strict_hub_genes
    ] <- paste0(
        strict_hub_genes,
        "\u2020"
    )
    
    heatmap_caption <- paste0(
        "\u2020 Strict hub. ",
        "* FDR < 0.05; ** FDR < 0.01; *** FDR < 0.001 ",
        "(Benjamini\u2013Hochberg correction across 44 tests)."
    )
    
    p_C_base <- ggplot2::ggplot(
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
            name = expression(
                "Spearman " * rho
            )
        ) +
        ggplot2::scale_x_discrete(
            labels = function(x) {
                stringr::str_wrap(
                    x,
                    width = 18
                )
            }
        ) +
        ggplot2::scale_y_discrete(
            labels = gene_label_map
        ) +
        ggplot2::labs(
            x = NULL,
            y = NULL,
            caption = heatmap_caption
        ) +
        ggplot2::theme_minimal(
            base_size = 12.2
        ) +
        ggplot2::theme(
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
            plot.caption = ggplot2::element_text(
                size = 8.5,
                color = "grey30",
                hjust = 0,
                margin = ggplot2::margin(
                    t = 7
                )
            ),
            plot.margin = ggplot2::margin(
                8,
                10,
                8,
                8
            )
        )
    
    p_C_with_title <- p_C_base +
        ggplot2::labs(
            title = "Hub gene\u2013mitochondrial pathway correlations",
            subtitle = paste0(
                "Spearman correlations after removing all 11 hub genes ",
                "from pathway sets"
            )
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
            )
        )
    
    p_C_no_title <- p_C_base +
        ggplot2::labs(
            title = NULL,
            subtitle = NULL
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_blank(),
            plot.margin = ggplot2::margin(
                4,
                10,
                8,
                8
            )
        )
    
    figure5c_with_title_path <- save_png(
        plot_object = p_C_with_title,
        filename = paste0(
            "Figure5C_11gene_hub_pathway_correlation_",
            "with_title_600dpi.png"
        ),
        width = 8.8,
        height = 6.25
    )
    
    figure5c_no_title_path <- save_png(
        plot_object = p_C_no_title,
        filename = paste0(
            "Figure5C_11gene_hub_pathway_correlation_",
            "no_title_600dpi.png"
        ),
        width = 8.8,
        height = 5.65
    )
    
    ############################################################
    ## 6. Figure 5D — revised running-enrichment curves
    ############################################################
    
    target_gsea_statistics$Pathway <- as.character(
        target_gsea_statistics$Pathway
    )
    
    target_gsea_statistics$Short <- as.character(
        target_gsea_statistics$Short
    )
    
    target_gsea_statistics$NES <- as.numeric(
        target_gsea_statistics$NES
    )
    
    target_gsea_statistics$FDR <- as.numeric(
        target_gsea_statistics$FDR
    )
    
    target_gsea_statistics <- target_gsea_statistics[
        match(
            pathway_order,
            target_gsea_statistics$Pathway
        ),
        ,
        drop = FALSE
    ]
    
    if (
        anyNA(
            target_gsea_statistics$Pathway
        )
    ) {
        stop(
            "Figure 5D target GSEA statistics do not contain all four pathways.",
            call. = FALSE
        )
    }
    
    legend_labels <- stats::setNames(
        paste0(
            target_gsea_statistics$Short,
            " (NES ",
            vapply(
                target_gsea_statistics$NES,
                format_nes,
                character(1)
            ),
            ", ",
            vapply(
                target_gsea_statistics$FDR,
                format_fdr,
                character(1)
            ),
            ")"
        ),
        target_gsea_statistics$Pathway
    )
    
    figure5d_es$Pathway <- factor(
        as.character(
            figure5d_es$Pathway
        ),
        levels = pathway_order
    )
    
    figure5d_hits$Pathway <- as.character(
        figure5d_hits$Pathway
    )
    
    figure5d_hits$Short <- as.character(
        figure5d_hits$Short
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
    
    if (!"y" %in% colnames(figure5d_hits)) {
        
        figure5d_hits$y <- track_table$y[
            match(
                figure5d_hits$Pathway,
                track_table$Pathway
            )
        ]
    }
    
    if (anyNA(figure5d_hits$y)) {
        stop(
            "Figure 5D hit-track y positions could not be reconstructed.",
            call. = FALSE
        )
    }
    
    rank_metric_table <- data.frame(
        Rank = seq_along(
            rank_vec
        ),
        rank_metric = as.numeric(
            rank_vec
        ),
        metric_scaled = as.numeric(
            rank_vec
        ) /
            max(
                abs(
                    as.numeric(
                        rank_vec
                    )
                ),
                na.rm = TRUE
            ),
        stringsAsFactors = FALSE
    )
    
    zero_cross <- which.min(
        abs(
            rank_metric_table$rank_metric
        )
    )
    
    number_bins <- 220L
    
    strip_x <- seq(
        from = 1,
        to = length(
            rank_vec
        ),
        length.out = number_bins
    )
    
    phenotype_strip <- data.frame(
        xmin = strip_x,
        xmax = c(
            strip_x[-1],
            length(
                rank_vec
            )
        ),
        value = seq(
            from = 1,
            to = -1,
            length.out = number_bins
        ),
        stringsAsFactors = FALSE
    )
    
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
                size = 6.9,
                color = "black",
                lineheight = 0.82
            ),
            legend.background = ggplot2::element_rect(
                fill = "white",
                color = "grey65",
                linewidth = 0.22
            ),
            legend.key.height = grid::unit(
                0.32,
                "lines"
            ),
            legend.key.width = grid::unit(
                1.12,
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
                size = 7.25,
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
    
    p_D_title_panel <- ggplot2::ggplot() +
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
                "Running-enrichment profiles for four prespecified ",
                "mitochondrial pathways"
            ),
            size = 3.35,
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
    
    p_D_with_title <- patchwork::wrap_plots(
        p_D_title_panel,
        p_D_es,
        p_D_hits,
        p_D_metric,
        ncol = 1,
        heights = c(
            0.11,
            0.45,
            0.25,
            0.32
        )
    )
    
    p_D_no_title <- patchwork::wrap_plots(
        p_D_es,
        p_D_hits,
        p_D_metric,
        ncol = 1,
        heights = c(
            0.45,
            0.25,
            0.32
        )
    )
    
    figure5d_with_title_path <- save_png(
        plot_object = p_D_with_title,
        filename = "Figure5D_mitochondrial_GSEA_curves_with_title_600dpi.png",
        width = 7.8,
        height = 5.45
    )
    
    figure5d_no_title_path <- save_png(
        plot_object = p_D_no_title,
        filename = "Figure5D_mitochondrial_GSEA_curves_no_title_600dpi.png",
        width = 7.8,
        height = 4.95
    )
    
    ############################################################
    ## 7. Export complete reproducible code to three TXT files
    ############################################################
    
    function_code <- deparse(
        rebuild_figure5_bcd,
        control = "all",
        width.cutoff = 500L
    )
    
    function_code[1] <- paste0(
        "rebuild_figure5_bcd <- ",
        function_code[1]
    )
    
    common_code <- c(
        "############################################################",
        "## Figure 5B–D publication-style revision — FINAL",
        "## Complete self-contained code exported automatically",
        "############################################################",
        "",
        function_code,
        "",
        "figure5_bcd_result <- rebuild_figure5_bcd()",
        ""
    )
    
    figure5b_code_path <- file.path(
        figure_dir,
        "Figure5B_mitochondrial_ssGSEA_REVISED_code.txt"
    )
    
    figure5c_code_path <- file.path(
        figure_dir,
        "Figure5C_11gene_hub_pathway_correlation_REVISED_code.txt"
    )
    
    figure5d_code_path <- file.path(
        figure_dir,
        "Figure5D_mitochondrial_GSEA_curves_REVISED_code.txt"
    )
    
    writeLines(
        c(
            "## This code generates the titled and title-free Figure 5B PNG files.",
            "",
            common_code
        ),
        con = figure5b_code_path,
        useBytes = TRUE
    )
    
    writeLines(
        c(
            "## This code generates the titled and title-free Figure 5C PNG files.",
            "",
            common_code
        ),
        con = figure5c_code_path,
        useBytes = TRUE
    )
    
    writeLines(
        c(
            "## This code generates the titled and title-free Figure 5D PNG files.",
            "",
            common_code
        ),
        con = figure5d_code_path,
        useBytes = TRUE
    )
    
    ############################################################
    ## 8. Final validation
    ############################################################
    
    expected_outputs <- c(
        figure5b_with_title_path,
        figure5b_no_title_path,
        figure5b_code_path,
        figure5c_with_title_path,
        figure5c_no_title_path,
        figure5c_code_path,
        figure5d_with_title_path,
        figure5d_no_title_path,
        figure5d_code_path
    )
    
    expected_roles <- c(
        "Figure 5B titled PNG",
        "Figure 5B title-free PNG",
        "Figure 5B code TXT",
        "Figure 5C titled PNG",
        "Figure 5C title-free PNG",
        "Figure 5C code TXT",
        "Figure 5D titled PNG",
        "Figure 5D title-free PNG",
        "Figure 5D code TXT"
    )
    
    output_check <- data.frame(
        Role = expected_roles,
        Path = expected_outputs,
        Exists = file.exists(
            expected_outputs
        ),
        Size_bytes = file.info(
            expected_outputs
        )$size,
        stringsAsFactors = FALSE
    )
    
    failed_outputs <- (
        !output_check$Exists |
            is.na(
                output_check$Size_bytes
            ) |
            output_check$Size_bytes <= 0
    )
    
    if (any(failed_outputs)) {
        stop(
            paste0(
                "One or more required outputs are missing or empty:\n",
                paste(
                    output_check$Path[
                        failed_outputs
                    ],
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }
    
    cat("\n============================================================\n")
    cat("FIGURE 5B–D REVISION COMPLETED SUCCESSFULLY\n")
    cat("============================================================\n")
    
    print(
        output_check,
        row.names = FALSE
    )
    
    cat("\nEach panel produced exactly three files:\n")
    cat("1 titled PNG + 1 title-free PNG + 1 code TXT.\n")
    
    cat("\nNo PDF or additional analysis file was generated.\n")
    cat("No files were written to backup.\n")
    cat("============================================================\n")
    
    invisible(
        list(
            plot_B_with_title = p_B_with_title,
            plot_B_no_title = p_B_no_title,
            plot_C_with_title = p_C_with_title,
            plot_C_no_title = p_C_no_title,
            plot_D_with_title = p_D_with_title,
            plot_D_no_title = p_D_no_title,
            outputs = output_check
        )
    )
}
