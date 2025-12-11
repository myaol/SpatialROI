#' Compare Features Between Groups (Violin Plot)
#'
#' @description
#' Create violin plots to compare gene expression or metadata between groups
#' with statistical testing
#'
#' @param seurat_obj Seurat spatial object
#' @param feature Character; gene name, metadata column, or "custom" for gene set
#' @param feature_type Character; "gene", "metadata", or "geneset"
#' @param group1_spots Character vector of spot IDs for group 1 (optional)
#' @param group2_spots Character vector of spot IDs for group 2 (optional)
#' @param comparison Character; "group1_vs_rest", "group2_vs_rest", or "group1_vs_group2"
#' @param gene_set Character vector; gene names for gene set scoring (if feature_type = "geneset")
#' @param stat_test Character; "wilcox" (default) or "ttest"
#' @param colors Character vector; colors for groups (default: c("#E41A1C", "#377EB8", "#999999"))
#' @param plot_title Character; custom plot title (optional)
#' @param return_data Logical; return data instead of plot (default: FALSE)
#'
#' @return ggplot2 object or data frame (if return_data = TRUE)
#' @import ggplot2
#' @import Seurat
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare gene expression between groups
#' group1 <- draw_ROI(seurat_obj)
#' group2 <- draw_ROI(seurat_obj)
#' plot_group_comparison(seurat_obj, "Cd3d", "gene",
#'                       group1_spots = group1, group2_spots = group2,
#'                       comparison = "group1_vs_group2")
#'
#' # Compare gene vs rest
#' plot_group_comparison(seurat_obj, "Cd8a", "gene",
#'                       group1_spots = group1,
#'                       comparison = "group1_vs_rest")
#'
#' # Compare metadata
#' plot_group_comparison(seurat_obj, "nCount_Spatial", "metadata",
#'                       group1_spots = group1, group2_spots = group2,
#'                       comparison = "group1_vs_group2")
#' }

plot_group_comparison <- function(seurat_obj,
                                  feature,
                                  feature_type = "gene",
                                  group1_spots = NULL,
                                  group2_spots = NULL,
                                  comparison = "group1_vs_group2",
                                  gene_set = NULL,
                                  stat_test = "wilcox",
                                  colors = c("#E41A1C", "#377EB8", "#999999"),
                                  plot_title = NULL,
                                  return_data = FALSE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
  }
  
  # Validate inputs
  if (!feature_type %in% c("gene", "metadata", "geneset")) {
    stop("feature_type must be 'gene', 'metadata', or 'geneset'")
  }
  
  if (!comparison %in% c("group1_vs_rest", "group2_vs_rest", "group1_vs_group2")) {
    stop("comparison must be 'group1_vs_rest', 'group2_vs_rest', or 'group1_vs_group2'")
  }
  
  # Fix spot IDs if needed (remove prefixes)
  if (!is.null(group1_spots) && !all(group1_spots %in% colnames(seurat_obj))) {
    group1_spots <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", group1_spots)
  }
  if (!is.null(group2_spots) && !all(group2_spots %in% colnames(seurat_obj))) {
    group2_spots <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", group2_spots)
  }
  
  # Validate spots
  if (comparison == "group1_vs_rest" && (is.null(group1_spots) || length(group1_spots) == 0)) {
    stop("group1_spots is required for group1_vs_rest comparison")
  }
  if (comparison == "group2_vs_rest" && (is.null(group2_spots) || length(group2_spots) == 0)) {
    stop("group2_spots is required for group2_vs_rest comparison")
  }
  if (comparison == "group1_vs_group2" &&
      (is.null(group1_spots) || length(group1_spots) == 0 ||
       is.null(group2_spots) || length(group2_spots) == 0)) {
    stop("Both group1_spots and group2_spots are required for group1_vs_group2 comparison")
  }
  
  # Extract feature values
  if (feature_type == "gene") {
    if (!feature %in% rownames(seurat_obj)) {
      stop("Gene '", feature, "' not found in Seurat object")
    }
    values <- GetAssayData(seurat_obj, slot = "data")[feature, ]
    feature_label <- feature
    
  } else if (feature_type == "metadata") {
    if (!feature %in% colnames(seurat_obj@meta.data)) {
      stop("Metadata column '", feature, "' not found in Seurat object")
    }
    values <- seurat_obj@meta.data[[feature]]
    names(values) <- colnames(seurat_obj)
    feature_label <- feature
    
  } else if (feature_type == "geneset") {
    if (is.null(gene_set) || length(gene_set) == 0) {
      stop("gene_set is required when feature_type = 'geneset'")
    }
    genes_found <- gene_set[gene_set %in% rownames(seurat_obj)]
    if (length(genes_found) == 0) {
      stop("None of the genes in gene_set found in Seurat object")
    }
    expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_found, , drop = FALSE]
    values <- colMeans(as.matrix(expr_data))
    feature_label <- paste0(feature, " (", length(genes_found), " genes)")
  }
  
  # Create group labels
  group_labels <- rep("Other", ncol(seurat_obj))
  names(group_labels) <- colnames(seurat_obj)
  
  if (!is.null(group1_spots)) {
    group_labels[group1_spots] <- "Group 1"
  }
  if (!is.null(group2_spots)) {
    group_labels[group2_spots] <- "Group 2"
  }
  
  # Filter based on comparison
  if (comparison == "group1_vs_rest") {
    keep_spots <- names(group_labels)
    group_labels <- factor(group_labels, levels = c("Group 1", "Other"))
  } else if (comparison == "group2_vs_rest") {
    keep_spots <- names(group_labels)
    group_labels <- factor(group_labels, levels = c("Group 2", "Other"))
  } else if (comparison == "group1_vs_group2") {
    keep_spots <- c(group1_spots, group2_spots)
    group_labels <- group_labels[keep_spots]
    values <- values[keep_spots]
    group_labels <- factor(group_labels, levels = c("Group 1", "Group 2"))
  }
  
  # Create data frame
  plot_data <- data.frame(
    value = values[keep_spots],
    group = group_labels
  )
  
  # Remove NAs
  plot_data <- plot_data[!is.na(plot_data$value), ]
  
  # Statistical test
  groups <- levels(plot_data$group)
  if (length(groups) == 2) {
    group1_vals <- plot_data$value[plot_data$group == groups[1]]
    group2_vals <- plot_data$value[plot_data$group == groups[2]]
    
    if (stat_test == "wilcox") {
      test_result <- wilcox.test(group1_vals, group2_vals)
      test_name <- "Wilcoxon"
    } else {
      test_result <- t.test(group1_vals, group2_vals)
      test_name <- "t-test"
    }
    
    p_value <- test_result$p.value
    if (p_value < 0.001) {
      p_label <- "p < 0.001"
    } else if (p_value < 0.01) {
      p_label <- sprintf("p = %.3f", p_value)
    } else {
      p_label <- sprintf("p = %.2f", p_value)
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
    
    stat_label <- paste0(test_name, ": ", p_label, " ", sig_label)
  } else {
    stat_label <- ""
  }
  
  # Return data if requested
  if (return_data) {
    result <- list(
      data = plot_data,
      p_value = if (exists("p_value")) p_value else NA,
      test = stat_test
    )
    return(result)
  }
  
  # Create plot title
  if (is.null(plot_title)) {
    if (comparison == "group1_vs_rest") {
      plot_title <- paste(feature_label, "- Group 1 vs Rest")
    } else if (comparison == "group2_vs_rest") {
      plot_title <- paste(feature_label, "- Group 2 vs Rest")
    } else {
      plot_title <- paste(feature_label, "- Group 1 vs Group 2")
    }
  }
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    scale_fill_manual(values = colors[1:length(groups)]) +
    labs(
      title = plot_title,
      subtitle = stat_label,
      x = "Group",
      y = feature_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkgray"),
      legend.position = "none",
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Add significance bracket if p < 0.05
  if (exists("p_value") && p_value < 0.05) {
    y_max <- max(plot_data$value, na.rm = TRUE)
    y_pos <- y_max + 0.1 * (y_max - min(plot_data$value, na.rm = TRUE))
    
    p <- p +
      annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos, size = 0.5) +
      annotate("text", x = 1.5, y = y_pos * 1.02, label = sig_label, size = 5)
  }
  
  print(p)
  return(invisible(p))
}


#' Compare Two Features (Scatter Plot)
#'
#' @description
#' Create scatter plot to compare two features with correlation analysis
#'
#' @param seurat_obj Seurat spatial object
#' @param feature1 Character; first feature name
#' @param feature2 Character; second feature name
#' @param feature1_type Character; "gene", "metadata", or "geneset"
#' @param feature2_type Character; "gene", "metadata", or "geneset"
#' @param gene_set1 Character vector; genes for first gene set (if feature1_type = "geneset")
#' @param gene_set2 Character vector; genes for second gene set (if feature2_type = "geneset")
#' @param spots_to_use Character vector; spot IDs to include (default: all spots)
#' @param cor_method Character; "spearman" (default) or "pearson"
#' @param add_smooth Logical; add smoothing line (default: TRUE)
#' @param point_size Numeric; size of points (default: 1.5)
#' @param point_alpha Numeric; transparency of points (default: 0.6)
#' @param plot_title Character; custom plot title (optional)
#'
#' @return ggplot2 object
#' @import ggplot2
#' @import Seurat
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two genes
#' plot_feature_comparison(seurat_obj, "Cd3d", "Cd8a", "gene", "gene")
#'
#' # Gene vs metadata
#' plot_feature_comparison(seurat_obj, "Cd3d", "nCount_Spatial",
#'                        "gene", "metadata")
#'
#' # Compare within specific spots
#' roi_spots <- draw_ROI(seurat_obj)
#' plot_feature_comparison(seurat_obj, "Cd4", "Cd8a", "gene", "gene",
#'                        spots_to_use = roi_spots)
#' }

plot_feature_comparison <- function(seurat_obj,
                                    feature1,
                                    feature2,
                                    feature1_type = "gene",
                                    feature2_type = "gene",
                                    gene_set1 = NULL,
                                    gene_set2 = NULL,
                                    spots_to_use = NULL,
                                    cor_method = "spearman",
                                    add_smooth = TRUE,
                                    point_size = 1.5,
                                    point_alpha = 0.6,
                                    plot_title = NULL) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required")
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
  values1 <- extract_feature_values(seurat_obj, feature1, feature1_type, gene_set1)
  label1 <- if (feature1_type == "geneset") {
    paste0(feature1, " (", length(gene_set1[gene_set1 %in% rownames(seurat_obj)]), " genes)")
  } else {
    feature1
  }
  
  # Extract feature 2
  values2 <- extract_feature_values(seurat_obj, feature2, feature2_type, gene_set2)
  label2 <- if (feature2_type == "geneset") {
    paste0(feature2, " (", length(gene_set2[gene_set2 %in% rownames(seurat_obj)]), " genes)")
  } else {
    feature2
  }
  
  # Filter to selected spots
  values1 <- values1[spots_to_use]
  values2 <- values2[spots_to_use]
  
  # Create data frame
  plot_data <- data.frame(
    feature1 = values1,
    feature2 = values2
  )
  
  # Remove NAs
  plot_data <- plot_data[complete.cases(plot_data), ]
  
  # Calculate correlation
  cor_result <- cor.test(plot_data$feature1, plot_data$feature2, method = cor_method)
  cor_value <- cor_result$estimate
  p_value <- cor_result$p.value
  
  # Format correlation label
  if (p_value < 0.001) {
    cor_label <- sprintf("%s r = %.3f, p < 0.001",
                         tools::toTitleCase(cor_method), cor_value)
  } else {
    cor_label <- sprintf("%s r = %.3f, p = %.3f",
                         tools::toTitleCase(cor_method), cor_value, p_value)
  }
  
  # Create plot title
  if (is.null(plot_title)) {
    plot_title <- paste(label1, "vs", label2)
  }
  
  # Create scatter plot
  p <- ggplot(plot_data, aes(x = feature1, y = feature2)) +
    geom_point(size = point_size, alpha = point_alpha, color = "#377EB8") +
    labs(
      title = plot_title,
      subtitle = cor_label,
      x = label1,
      y = label2
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "darkgray"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold")
    )
  
  # Add smoothing line
  if (add_smooth) {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "#E41A1C", alpha = 0.2)
  }
  
  print(p)
  return(invisible(p))
}


#' Helper function to extract feature values
#' @keywords internal

extract_feature_values <- function(seurat_obj, feature, feature_type, gene_set = NULL) {
  
  if (feature_type == "gene") {
    if (!feature %in% rownames(seurat_obj)) {
      stop("Gene '", feature, "' not found in Seurat object")
    }
    values <- GetAssayData(seurat_obj, slot = "data")[feature, ]
    
  } else if (feature_type == "metadata") {
    if (!feature %in% colnames(seurat_obj@meta.data)) {
      stop("Metadata column '", feature, "' not found in Seurat object")
    }
    values <- seurat_obj@meta.data[[feature]]
    names(values) <- colnames(seurat_obj)
    
  } else if (feature_type == "geneset") {
    if (is.null(gene_set) || length(gene_set) == 0) {
      stop("gene_set is required when feature_type = 'geneset'")
    }
    genes_found <- gene_set[gene_set %in% rownames(seurat_obj)]
    if (length(genes_found) == 0) {
      stop("None of the genes in gene_set found in Seurat object")
    }
    expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_found, , drop = FALSE]
    values <- colMeans(as.matrix(expr_data))
  }
  
  return(values)
}


#' Quick Violin Plot Wrapper
#'
#' @description
#' Simplified function to quickly compare a gene between two groups
#'
#' @param seurat_obj Seurat object
#' @param gene Gene name
#' @param group1_spots Spot IDs for group 1
#' @param group2_spots Spot IDs for group 2 (optional)
#'
#' @export

quick_violin <- function(seurat_obj, gene, group1_spots, group2_spots = NULL) {
  
  if (is.null(group2_spots)) {
    # Group 1 vs Rest
    plot_group_comparison(
      seurat_obj,
      feature = gene,
      feature_type = "gene",
      group1_spots = group1_spots,
      comparison = "group1_vs_rest"
    )
  } else {
    # Group 1 vs Group 2
    plot_group_comparison(
      seurat_obj,
      feature = gene,
      feature_type = "gene",
      group1_spots = group1_spots,
      group2_spots = group2_spots,
      comparison = "group1_vs_group2"
    )
  }
}


#' Quick Scatter Plot Wrapper
#'
#' @description
#' Simplified function to quickly compare two genes
#'
#' @param seurat_obj Seurat object
#' @param gene1 First gene name
#' @param gene2 Second gene name
#' @param spots Optional; spot IDs to include
#'
#' @export

quick_scatter <- function(seurat_obj, gene1, gene2, spots = NULL) {
  plot_feature_comparison(
    seurat_obj,
    feature1 = gene1,
    feature2 = gene2,
    feature1_type = "gene",
    feature2_type = "gene",
    spots_to_use = spots
  )
}