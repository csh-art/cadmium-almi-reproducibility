combine_figure5_abcd_compact_biglabels <- function() {
    
    options(stringsAsFactors = FALSE)
    
    ############################################################
    ## 0. Paths
    ############################################################
    
    project_root <- "E:/1Data/Cd_Sarcopenia_NHANES"
    
    figure_dir <- file.path(
        project_root,
        "new_workspace",
        "03_results",
        "M05_PathwayValidation_GSEA_ssGSEA_Correlation",
        "03_outputs_figures_tables",
        "figure_re"
    )
    
    figure_a_path <- file.path(
        figure_dir,
        "Figure5A_GSEA_mechanism_dotplot_no_title_600dpi.png"
    )
    
    figure_b_path <- file.path(
        figure_dir,
        "Figure5B_mitochondrial_ssGSEA_no_title_600dpi.png"
    )
    
    figure_c_path <- file.path(
        figure_dir,
        "Figure5C_11gene_hub_pathway_correlation_no_title_600dpi.png"
    )
    
    figure_d_path <- file.path(
        figure_dir,
        "Figure5D_mitochondrial_GSEA_curves_no_title_600dpi.png"
    )
    
    output_png_path <- file.path(
        figure_dir,
        "Figure5_ABCD_combined_no_title_COMPACT_biglabels_600dpi.png"
    )
    
    output_code_path <- file.path(
        figure_dir,
        "Figure5_ABCD_combined_no_title_COMPACT_biglabels_code.txt"
    )
    
    dir.create(
        figure_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    ############################################################
    ## 1. Validate package
    ############################################################
    
    if (!requireNamespace("magick", quietly = TRUE)) {
        stop(
            paste0(
                "Package 'magick' is required.\n",
                "Install it with:\n",
                "install.packages('magick')"
            ),
            call. = FALSE
        )
    }
    
    ############################################################
    ## 2. Validate input files
    ############################################################
    
    input_paths <- c(
        figure_a_path,
        figure_b_path,
        figure_c_path,
        figure_d_path
    )
    
    input_labels <- c(
        "Figure 5A",
        "Figure 5B",
        "Figure 5C",
        "Figure 5D"
    )
    
    input_check <- data.frame(
        Panel = input_labels,
        Path = input_paths,
        Exists = file.exists(input_paths),
        Size_bytes = file.info(input_paths)$size,
        stringsAsFactors = FALSE
    )
    
    missing_inputs <- (
        !input_check$Exists |
            is.na(input_check$Size_bytes) |
            input_check$Size_bytes <= 0
    )
    
    if (any(missing_inputs)) {
        stop(
            paste0(
                "One or more input PNG files are missing or empty:\n",
                paste(
                    input_check$Path[missing_inputs],
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }
    
    cat("\n============================================================\n")
    cat("FIGURE 5A–D COMPACT COMBINATION — BIG LABELS VERSION\n")
    cat("============================================================\n")
    
    cat("\nInput files:\n")
    print(input_check, row.names = FALSE)
    
    ############################################################
    ## 3. Final canvas dimensions
    ############################################################
    
    final_dpi <- 600L
    
    ## same final physical size
    canvas_width <- 10200L
    canvas_height <- 7440L
    
    ## compact layout
    outer_margin_x <- 90L
    outer_margin_y <- 90L
    
    horizontal_gap <- 90L
    vertical_gap <- 90L
    
    ## increased to accommodate larger A/B/C/D letters
    panel_label_band <- 135L
    
    cell_width <- as.integer(
        (
            canvas_width -
                2L * outer_margin_x -
                horizontal_gap
        ) / 2L
    )
    
    cell_height <- as.integer(
        (
            canvas_height -
                2L * outer_margin_y -
                vertical_gap
        ) / 2L
    )
    
    content_height <- cell_height - panel_label_band
    
    inner_width <- 2L * cell_width + horizontal_gap
    inner_height <- 2L * cell_height + vertical_gap
    
    stopifnot(
        cell_width > 0L,
        cell_height > 0L,
        content_height > 0L,
        inner_width <= canvas_width,
        inner_height <= canvas_height
    )
    
    cat("\nLayout settings:\n")
    cat("Canvas:", canvas_width, "x", canvas_height, "px\n")
    cat("Cell:", cell_width, "x", cell_height, "px\n")
    cat("Content height:", content_height, "px\n")
    cat("Panel label band:", panel_label_band, "px\n")
    
    ############################################################
    ## 4. Read source images
    ############################################################
    
    image_a <- magick::image_read(figure_a_path)
    image_b <- magick::image_read(figure_b_path)
    image_c <- magick::image_read(figure_c_path)
    image_d <- magick::image_read(figure_d_path)
    
    ############################################################
    ## 5. Panel preparation helper
    ############################################################
    
    prepare_panel <- function(
        image_object,
        panel_label
    ) {
        
        image_object <- magick::image_background(
            image_object,
            color = "white",
            flatten = TRUE
        )
        
        image_object <- magick::image_trim(
            image_object,
            fuzz = 1
        )
        
        image_object <- magick::image_resize(
            image_object,
            geometry = paste0(
                cell_width,
                "x",
                content_height
            ),
            filter = "Lanczos"
        )
        
        image_object <- magick::image_extent(
            image_object,
            geometry = paste0(
                cell_width,
                "x",
                content_height
            ),
            gravity = "center",
            color = "white"
        )
        
        panel_canvas <- magick::image_extent(
            image_object,
            geometry = paste0(
                cell_width,
                "x",
                cell_height
            ),
            gravity = "south",
            color = "white"
        )
        
        ## enlarged panel label
        panel_canvas <- magick::image_annotate(
            panel_canvas,
            text = panel_label,
            gravity = "northwest",
            location = "+22+6",
            size = 140,
            weight = 700,
            color = "black",
            strokecolor = "transparent"
        )
        
        panel_canvas
    }
    
    ############################################################
    ## 6. Prepare A–D panels
    ############################################################
    
    panel_a <- prepare_panel(image_a, "A")
    panel_b <- prepare_panel(image_b, "B")
    panel_c <- prepare_panel(image_c, "C")
    panel_d <- prepare_panel(image_d, "D")
    
    ############################################################
    ## 7. Create gaps
    ############################################################
    
    horizontal_spacer <- magick::image_blank(
        width = horizontal_gap,
        height = cell_height,
        color = "white"
    )
    
    vertical_spacer <- magick::image_blank(
        width = inner_width,
        height = vertical_gap,
        color = "white"
    )
    
    ############################################################
    ## 8. Assemble rows
    ############################################################
    
    top_row <- magick::image_append(
        c(panel_a, horizontal_spacer, panel_b),
        stack = FALSE
    )
    
    bottom_row <- magick::image_append(
        c(panel_c, horizontal_spacer, panel_d),
        stack = FALSE
    )
    
    combined_inner <- magick::image_append(
        c(top_row, vertical_spacer, bottom_row),
        stack = TRUE
    )
    
    ############################################################
    ## 9. Add outer margins
    ############################################################
    
    combined_figure <- magick::image_extent(
        combined_inner,
        geometry = paste0(canvas_width, "x", canvas_height),
        gravity = "center",
        color = "white"
    )
    
    ############################################################
    ## 10. Write final PNG
    ############################################################
    
    magick::image_write(
        image = combined_figure,
        path = output_png_path,
        format = "png",
        quality = 100,
        depth = 8,
        density = paste0(final_dpi, "x", final_dpi)
    )
    
    output_png_size <- file.info(output_png_path)$size
    
    if (
        !file.exists(output_png_path) ||
        is.na(output_png_size) ||
        output_png_size <= 0
    ) {
        stop(
            "Combined PNG was not created successfully:\n",
            output_png_path,
            call. = FALSE
        )
    }
    
    ############################################################
    ## 11. Verify final dimensions
    ############################################################
    
    final_image_check <- magick::image_read(output_png_path)
    final_image_info <- magick::image_info(final_image_check)
    
    if (
        final_image_info$width[1] != canvas_width ||
        final_image_info$height[1] != canvas_height
    ) {
        stop(
            paste0(
                "Unexpected combined-image dimensions.\n",
                "Expected: ",
                canvas_width,
                " × ",
                canvas_height,
                " pixels\n",
                "Observed: ",
                final_image_info$width[1],
                " × ",
                final_image_info$height[1],
                " pixels"
            ),
            call. = FALSE
        )
    }
    
    ############################################################
    ## 12. Save complete code TXT
    ############################################################
    
    function_code <- deparse(
        combine_figure5_abcd_compact_biglabels,
        control = "all",
        width.cutoff = 500L
    )
    
    function_code[1] <- paste0(
        "combine_figure5_abcd_compact_biglabels <- ",
        function_code[1]
    )
    
    complete_code <- c(
        "############################################################",
        "## Figure 5A–D combined publication figure — COMPACT FINAL",
        "## Larger panel labels version",
        "############################################################",
        "",
        function_code,
        "",
        "figure5_combination_result <- combine_figure5_abcd_compact_biglabels()",
        ""
    )
    
    writeLines(
        complete_code,
        con = output_code_path,
        useBytes = TRUE
    )
    
    output_code_size <- file.info(output_code_path)$size
    
    if (
        !file.exists(output_code_path) ||
        is.na(output_code_size) ||
        output_code_size <= 0
    ) {
        stop(
            "Combined-figure code TXT was not created:\n",
            output_code_path,
            call. = FALSE
        )
    }
    
    ############################################################
    ## 13. Final report
    ############################################################
    
    cat("\n============================================================\n")
    cat("FIGURE 5A–D COMPACT BIG-LABELS COMBINATION SUCCESSFUL\n")
    cat("============================================================\n")
    
    cat("\nPanel arrangement:\n")
    cat("A | B\n")
    cat("C | D\n")
    
    cat("\nFinal dimensions:\n")
    cat(
        final_image_info$width[1],
        " × ",
        final_image_info$height[1],
        " pixels\n",
        sep = ""
    )
    
    cat("\nFinal physical size at 600 dpi:\n")
    cat("17 × 12.4 inches\n")
    
    cat("\nCombined PNG:\n")
    cat(output_png_path, "\n")
    
    cat("\nCombined PNG MD5:\n")
    cat(unname(tools::md5sum(output_png_path)), "\n")
    
    cat("\nCode TXT:\n")
    cat(output_code_path, "\n")
    
    cat("\nA/B/C/D labels were enlarged in this version.\n")
    cat("No overall title was added.\n")
    cat("No files were written to backup.\n")
    cat("============================================================\n")
    
    invisible(
        list(
            combined_png = output_png_path,
            code_txt = output_code_path,
            width_pixels = final_image_info$width[1],
            height_pixels = final_image_info$height[1]
        )
    )
}
