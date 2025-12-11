#' Compare Two Features with Violin Plots (Side-by-Side)
#'
#' @description
#' Create side-by-side violin plots to compare two features with statistical testing.
#' This matches the web interface "Feature vs Feature" comparison with violin plots.
#'
#' @param seurat_obj Seurat spatial object
#' @param feature1 Character; first feature name
#' @param feature2 Character; second feature name
#' @param feature1_type Character; "gene", "metadata", or "geneset"
#' @param feature2_type Character; "gene", "metadata", or "geneset"
#' @param gene_set1 Character vector; genes for first gene set (if feature1_type = "geneset")
#' @param gene_set2 Character vector; genes for second gene set (if feature2_type = "geneset")
#' @param spots_to_use Character vector; spot IDs to include (default: all spots)
#' @param stat_test Character; "wilcox" (default) or "ttest"
#' @param colors Character vector; colors for features (default: c("#56B4E9", "#0072B2"))
#' @param plot_title Character; custom plot title (optional)
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import Seurat
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two genes in Group 1
#' group1 <- draw_ROI(seurat_obj)
#' plot_feature_violin_comparison(seurat_obj, "Ighg3", "Slc26a3", 
#'                                 "gene", "gene", spots_to_use = group1)
#'
#' # Compare across all spots
#' plot_feature_violin_comparison(seurat_obj, "Cd3d", "Cd8a", "gene", "gene")
#'
#' # With t-test
#' plot_feature_violin_comparison(seurat_obj, "Cd3d", "Cd8a", 
#'                                "gene", "gene", stat_test = "ttest")
#' }

plot_feature_violin_comparison <- function(seurat_obj,
                                           feature1,
                                           feature2,
                                           feature1_type = "gene",
                                           feature2_type = "gene",
                                           gene_set1 = NULL,
                                           gene_set2 = NULL,
                                           spots_to_use = NULL,
                                           stat_test = "wilcox",
                                           colors = c("#56B4E9", "#0072B2"),
                                           plot_title = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
  }
  
  # Helper function to extract values
  extract_values <- function(seurat_obj, feature, feature_type, gene_set = NULL) {
    if (feature_type == "gene") {
      # Case-insensitive gene matching
      if (!feature %in% rownames(seurat_obj)) {
        gene_idx <- grep(paste0("^", feature, "$"), rownames(seurat_obj), ignore.case = TRUE)
        if (length(gene_idx) > 0) {
          gene_match <- rownames(seurat_obj)[gene_idx[1]]
          message("Gene '", feature, "' not found. Using '", gene_match, "' instead.")
          feature <- gene_match
        } else {
          stop("Gene '", feature, "' not found in Seurat object")
        }
      }
      values <- GetAssayData(seurat_obj, slot = "data")[feature, ]
      
    } else if (feature_type == "metadata") {
      if (!feature %in% colnames(seurat_obj@meta.data)) {
        stop("Metadata column '", feature, "' not found")
      }
      values <- seurat_obj@meta.data[[feature]]
      names(values) <- colnames(seurat_obj)
      
    } else if (feature_type == "geneset") {
      if (is.null(gene_set) || length(gene_set) == 0) {
        stop("gene_set is required when feature_type = 'geneset'")
      }
      genes_found <- gene_set[gene_set %in% rownames(seurat_obj)]
      if (length(genes_found) == 0) {
        stop("None of the genes in gene_set found")
      }
      expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_found, , drop = FALSE]
      values <- colMeans(as.matrix(expr_data))
    }
    return(values)
  }
  
  # Fix spot IDs if provided
  if (!is.null(spots_to_use) && !all(spots_to_use %in% colnames(seurat_obj))) {
    spots_to_use <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", spots_to_use)
  }
  
  # If no spots specified, use all
  if (is.null(spots_to_use)) {
    spots_to_use <- colnames(seurat_obj)
  }
  
  # Extract feature 1
  values1 <- extract_values(seurat_obj, feature1, feature1_type, gene_set1)
  label1 <- if (feature1_type == "geneset") {
    paste0(feature1, " (", length(gene_set1[gene_set1 %in% rownames(seurat_obj)]), " genes)")
  } else {
    feature1
  }
  
  # Extract feature 2
  values2 <- extract_values(seurat_obj, feature2, feature2_type, gene_set2)
  label2 <- if (feature2_type == "geneset") {
    paste0(feature2, " (", length(gene_set2[gene_set2 %in% rownames(seurat_obj)]), " genes)")
  } else {
    feature2
  }
  
  # Filter to selected spots
  values1 <- values1[spots_to_use]
  values2 <- values2[spots_to_use]
  
  # Create long-format data frame
  plot_data <- data.frame(
    value = c(values1, values2),
    feature = rep(c(label1, label2), each = length(spots_to_use)),
    spot_id = rep(spots_to_use, 2)
  )
  
  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$value), ]
  
  # Make feature a factor
  plot_data$feature <- factor(plot_data$feature, levels = c(label1, label2))
  
  # Statistical test
  if (stat_test == "wilcox") {
    test_result <- wilcox.test(values1, values2, na.rm = TRUE)
    test_name <- "Wilcoxon"
  } else {
    test_result <- t.test(values1, values2, na.rm = TRUE)
    test_name <- "t-test"
  }
  
  p_value <- test_result$p.value
  
  # Format p-value
  if (p_value < 0.001) {
    p_label <- "p < 0.001"
  } else {
    p_label <- format(p_value, scientific = TRUE, digits = 2)
  }
  
  # Add significance stars
  if (p_value < 0.001) {
    sig_label <- "***"
  } else if (p_value < 0.01) {
    sig_label <- "**"
  } else if (p_value < 0.05) {
    sig_label <- "*"
  } else {
    sig_label <- "ns"
  }
  
  stat_label <- sprintf("%s test p = %s (n = %d spots)",
                        test_name, p_label, length(spots_to_use))
  
  # Create plot title
  if (is.null(plot_title)) {
    plot_title <- "Two Features Comparison"
  }
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = feature, y = value, fill = feature)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = colors) +
    labs(
      title = plot_title,
      subtitle = stat_label,
      x = "",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkgray"),
      legend.position = "none",
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  # Add significance bracket if p < 0.05
  if (p_value < 0.05) {
    y_max <- max(plot_data$value, na.rm = TRUE)
    y_min <- min(plot_data$value, na.rm = TRUE)
    y_pos <- y_max + 0.1 * (y_max - y_min)
    
    p <- p +
      annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos, size = 0.5) +
      annotate("text", x = 1.5, y = y_pos * 1.02, label = sig_label, size = 5)
  }
  
  print(p)
  return(invisible(p))
}


#' Quick Two-Feature Violin Comparison
#'
#' @description
#' Simplified wrapper to compare two genes with violin plots
#'
#' @param seurat_obj Seurat object
#' @param gene1 First gene name
#' @param gene2 Second gene name
#' @param spots Optional; spot IDs to include
#' @param stat_test "wilcox" or "ttest"
#'
#' @export

quick_feature_violin <- function(seurat_obj, gene1, gene2, spots = NULL, stat_test = "wilcox") {
  plot_feature_violin_comparison(
    seurat_obj,
    feature1 = gene1,
    feature2 = gene2,
    feature1_type = "gene",
    feature2_type = "gene",
    spots_to_use = spots,
    stat_test = stat_test
  )
}