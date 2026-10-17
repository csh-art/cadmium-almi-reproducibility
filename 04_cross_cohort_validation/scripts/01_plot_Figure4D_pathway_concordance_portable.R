# Figure 4D: cross-cohort concordance of seven prespecified pathways
# Portable repository reconstruction from the archived final source-data table.
# Run from the repository root, or set CADMIUM_ALMI_REPO_ROOT.

options(stringsAsFactors = FALSE)

find_repo_root <- function() {
  env_root <- Sys.getenv("CADMIUM_ALMI_REPO_ROOT", unset = "")
  if (nzchar(env_root) && dir.exists(env_root)) return(normalizePath(env_root, winslash = "/"))
  current <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  repeat {
    marker <- file.path(current, "04_cross_cohort_validation", "source_data", "Figure4E_pathway_replication_source_data.csv")
    if (file.exists(marker)) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Repository root not found. Run from the repository root or set CADMIUM_ALMI_REPO_ROOT.", call. = FALSE)
}

required <- c("ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)

repo_root <- find_repo_root()
module_dir <- file.path(repo_root, "04_cross_cohort_validation")
source_file <- file.path(module_dir, "source_data", "Figure4E_pathway_replication_source_data.csv")
output_dir <- file.path(module_dir, "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(source_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("Pathway_label", "Dataset", "NES", "P_value", "P_adj_target7")
miss <- setdiff(required_columns, names(d))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

pathway_order <- c("Oxidative phosphorylation", "Respiratory electron transport", "Mitochondrial ATP synthesis", "Mitochondrial translation", "Apoptosis", "Autophagy", "Integrated stress response")
dataset_order <- c("GSE111016", "GSE111010", "GSE111006", "GSE167186")
d <- d[d$Pathway_label %in% pathway_order & d$Dataset %in% dataset_order, , drop = FALSE]
key <- paste(d$Pathway_label, d$Dataset, sep = "::")
stopifnot(nrow(d) == 28L, anyDuplicated(key) == 0L)
if (!"Available" %in% names(d)) d$Available <- is.finite(d$NES)

disc <- d[d$Dataset == "GSE111016", c("Pathway_label", "NES")]
names(disc)[2] <- "Discovery_NES_recalc"
d <- merge(d, disc, by = "Pathway_label", all.x = TRUE, sort = FALSE)
d$Same_direction_recalc <- NA
i <- d$Dataset != "GSE111016" & d$Available & is.finite(d$NES) & is.finite(d$Discovery_NES_recalc) & d$NES != 0 & d$Discovery_NES_recalc != 0
d$Same_direction_recalc[i] <- sign(d$NES[i]) == sign(d$Discovery_NES_recalc[i])
d$Significance_symbol <- ""
d$Significance_symbol[d$Dataset == "GSE111016" & d$Available & d$P_adj_target7 < 0.05] <- "\u2020"
d$Significance_symbol[d$Dataset != "GSE111016" & d$Available & d$Same_direction_recalc %in% TRUE & d$P_adj_target7 < 0.05] <- "**"
d$Significance_symbol[d$Dataset != "GSE111016" & d$Available & d$Same_direction_recalc %in% TRUE & !(d$P_adj_target7 < 0.05) & d$P_value < 0.05] <- "*"
d$Tile_label <- ifelse(!d$Available | !is.finite(d$NES), "NA", paste0(sprintf("%+.2f", d$NES), d$Significance_symbol))

d$Dataset_factor <- factor(d$Dataset, levels = dataset_order)
d$Pathway_factor <- factor(d$Pathway_label, levels = rev(pathway_order))
labels_x <- c(GSE111016 = "GSE111016\nDiscovery", GSE111010 = "GSE111010\nValidation", GSE111006 = "GSE111006\nSensitivity\u00b9", GSE167186 = "GSE167186\nSensitivity\u00b2")
max_abs <- ceiling(max(abs(d$NES[d$Available]), na.rm = TRUE) * 2) / 2

p <- ggplot2::ggplot(d, ggplot2::aes(Dataset_factor, Pathway_factor)) +
  ggplot2::geom_tile(ggplot2::aes(fill = NES), width = 0.965, height = 0.925, color = "white", linewidth = 1.25) +
  ggplot2::geom_vline(xintercept = 1.5, linetype = "dashed", linewidth = 0.55, color = "#666666") +
  ggplot2::geom_text(ggplot2::aes(label = Tile_label), size = 4.2, fontface = "bold") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), oob = scales::squish, name = "Normalized\nenrichment score") +
  ggplot2::scale_x_discrete(labels = labels_x, expand = ggplot2::expansion(add = 0.025)) +
  ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.025)) +
  ggplot2::coord_fixed(ratio = 0.72, clip = "off") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 11.5) +
  ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(fill = NA, color = "#4D4D4D", linewidth = 0.7), axis.text.y = ggplot2::element_text(size = 11.4, face = "bold", color = "#202020"), axis.text.x = ggplot2::element_text(size = 10, color = "#202020"), axis.ticks = ggplot2::element_blank(), legend.position = "right")

ggplot2::ggsave(file.path(output_dir, "Figure_04D_pathway_concordance_plain.png"), p, width = 9.25, height = 7.15, units = "in", dpi = 600, bg = "white")
ggplot2::ggsave(file.path(output_dir, "Figure_04D_pathway_concordance_plain.pdf"), p, width = 9.25, height = 7.15, units = "in", bg = "white")
cat("Figure 4D completed from: ", source_file, "\n", sep = "")
