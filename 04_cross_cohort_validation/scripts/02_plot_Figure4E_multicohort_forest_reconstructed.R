# Figure 4E: four-cohort hub-gene forest plot
# Portable reconstruction from the archived 44-row source-data table.
# Historical project files called this analysis "Figure4D_multicohort_forest";
# the final manuscript labels it Figure 4E.
#
# For sensitivity cohorts where lfcSE was not carried into the archived table,
# lfcSE is reconstructed from the DESeq2 Wald identity:
#     Wald_stat = log2FC_unshrunk / lfcSE
# This identity is verified against all rows where both lfcSE and Wald_stat are present.
# Run from the repository root, or set CADMIUM_ALMI_REPO_ROOT.

options(stringsAsFactors = FALSE)

find_repo_root <- function() {
  env_root <- Sys.getenv("CADMIUM_ALMI_REPO_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) return(normalizePath(env_root, winslash = "/"))
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "04_cross_cohort_validation", "source_data", "Figure4D_multicohort_forest_source_data.csv")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Repository root not found. Run from repository root or set CADMIUM_ALMI_REPO_ROOT.", call. = FALSE)
}

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing R package: ggplot2", call. = FALSE)
repo_root <- find_repo_root()
module_dir <- file.path(repo_root, "04_cross_cohort_validation")
source_file <- file.path(module_dir, "source_data", "Figure4D_multicohort_forest_source_data.csv")
output_dir <- file.path(module_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("Dataset", "Gene_symbol", "Hub_order", "Dataset_order", "log2FC_unshrunk", "lfcSE", "Wald_stat", "P_value", "P_adj_target11", "Same_direction", "Available_for_plot")
miss <- setdiff(need, names(d))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
stopifnot(nrow(d) == 44L, length(unique(d$Gene_symbol)) == 11L, length(unique(d$Dataset)) == 4L)

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(z)); out[z %in% c("true","t","1","yes","y")] <- TRUE; out[z %in% c("false","f","0","no","n")] <- FALSE; out
}
d$Available_for_plot <- to_logical(d$Available_for_plot)
d$Same_direction <- to_logical(d$Same_direction)
d$log2FC_unshrunk <- as.numeric(d$log2FC_unshrunk)
d$lfcSE <- as.numeric(d$lfcSE)
d$Wald_stat <- as.numeric(d$Wald_stat)
d$P_value <- as.numeric(d$P_value)
d$P_adj_target11 <- as.numeric(d$P_adj_target11)

d$lfcSE_from_Wald <- ifelse(d$Available_for_plot & is.finite(d$log2FC_unshrunk) & is.finite(d$Wald_stat) & d$Wald_stat != 0, abs(d$log2FC_unshrunk / d$Wald_stat), NA_real_)
both <- is.finite(d$lfcSE) & is.finite(d$lfcSE_from_Wald)
if (any(both)) stopifnot(max(abs(d$lfcSE[both] - d$lfcSE_from_Wald[both])) < 1e-10)
d$lfcSE_final <- ifelse(is.finite(d$lfcSE), d$lfcSE, d$lfcSE_from_Wald)
d$CI_low <- d$log2FC_unshrunk - 1.96 * d$lfcSE_final
d$CI_high <- d$log2FC_unshrunk + 1.96 * d$lfcSE_final
stopifnot(sum(d$Available_for_plot, na.rm = TRUE) == 41L, sum(!d$Available_for_plot, na.rm = TRUE) == 3L, sum(is.finite(d$lfcSE_final[d$Available_for_plot])) == 41L)

d$Evidence_marker <- ""
val <- d$Dataset != "GSE111016" & d$Available_for_plot & d$Same_direction %in% TRUE
d$Evidence_marker[val & is.finite(d$P_adj_target11) & d$P_adj_target11 < 0.05] <- "**"
d$Evidence_marker[val & d$Evidence_marker == "" & is.finite(d$P_value) & d$P_value < 0.05] <- "*"

d$Point_group <- ifelse(d$Dataset == "GSE111016", "Discovery reference", ifelse(!d$Available_for_plot, "Unavailable", ifelse(d$Same_direction %in% TRUE, "Same direction", "Opposite direction")))
d$Gene_symbol <- factor(d$Gene_symbol, levels = rev(unique(d$Gene_symbol[order(d$Hub_order)])))
d$Dataset <- factor(d$Dataset, levels = c("GSE111016","GSE111010","GSE111006","GSE167186"))
facet_labels <- c(GSE111016 = "GSE111016\nDiscovery\n20 vs 20", GSE111010 = "GSE111010\nValidation\n9 vs 14", GSE111006 = "GSE111006\nSensitivity\u00b9\n4 vs 28", GSE167186 = "GSE167186\nSensitivity\u00b2\n24 vs 29")
point_colors <- c("Discovery reference" = "#262626", "Same direction" = "#2C7FB8", "Opposite direction" = "#D64B3C")

available <- d[d$Available_for_plot & is.finite(d$log2FC_unshrunk) & is.finite(d$CI_low) & is.finite(d$CI_high), , drop = FALSE]
unavailable <- d[!d$Available_for_plot, , drop = FALSE]
star_df <- available[nzchar(available$Evidence_marker), , drop = FALSE]
star_df$star_x <- star_df$CI_high + 0.06

p <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.55, color = "#707070") +
  ggplot2::geom_errorbarh(data = available, ggplot2::aes(y = Gene_symbol, xmin = CI_low, xmax = CI_high, color = Point_group), height = 0.15, linewidth = 0.75) +
  ggplot2::geom_point(data = available, ggplot2::aes(x = log2FC_unshrunk, y = Gene_symbol, color = Point_group), size = 2.8) +
  ggplot2::geom_text(data = unavailable, ggplot2::aes(x = -1.20, y = Gene_symbol, label = "NA"), hjust = 0, color = "#8A8A8A", fontface = "italic", size = 3.1) +
  ggplot2::geom_text(data = star_df, ggplot2::aes(x = star_x, y = Gene_symbol, label = Evidence_marker, color = Point_group), hjust = 0, vjust = 0.2, fontface = "bold", size = 3.4) +
  ggplot2::facet_grid(. ~ Dataset, labeller = ggplot2::labeller(Dataset = facet_labels)) +
  ggplot2::scale_color_manual(values = point_colors, breaks = c("Discovery reference","Same direction","Opposite direction"), name = "Cohort result") +
  ggplot2::scale_x_continuous(limits = c(-1.30, 1.55), breaks = seq(-1, 1.5, by = 0.5)) +
  ggplot2::labs(x = "Unshrunk log2 fold-change", y = NULL) +
  ggplot2::theme_minimal(base_size = 10.5) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_line(color = "#ECEFF1", linewidth = 0.4), panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, color = "#535B61", linewidth = 0.65), strip.background = ggplot2::element_rect(fill = "#E9EEF1", color = "#535B61", linewidth = 0.65), strip.text = ggplot2::element_text(face = "bold", size = 10.2, lineheight = 0.95), axis.text.y = ggplot2::element_text(face = "bold.italic", color = "#202020"), axis.title.x = ggplot2::element_text(face = "bold"), legend.position = "bottom", legend.title = ggplot2::element_text(face = "bold"), panel.spacing.x = grid::unit(0.08, "in"))

ggplot2::ggsave(file.path(output_dir, "Figure_04E_multicohort_forest_plain.png"), p, width = 12.3, height = 6.2, units = "in", dpi = 600, bg = "white")
ggplot2::ggsave(file.path(output_dir, "Figure_04E_multicohort_forest_plain.pdf"), p, width = 12.3, height = 6.2, units = "in", bg = "white")
write.csv(d[, c("Dataset","Gene_symbol","Hub_order","log2FC_unshrunk","lfcSE","Wald_stat","lfcSE_from_Wald","lfcSE_final","CI_low","CI_high","P_value","P_adj_target11","Same_direction","Point_group","Evidence_marker","Available_for_plot")], file.path(module_dir, "source_data", "Figure_04E_multicohort_forest_reconstructed_plot_data_from_R.csv"), row.names = FALSE)
cat("Figure 4E completed from: ", source_file, "\n", sep = "")
