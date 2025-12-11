#' Perform Differential Expression Analysis on Spatial Transcriptomics Data
#'
#' @description
#' Run differential expression analysis comparing different groups of spots.
#' Supports group vs rest, group vs group, and cluster-based comparisons.
#'
#' @param seurat_obj A Seurat spatial object (full dataset)
#' @param group1_spots Character vector of spot IDs for group 1
#' @param group2_spots Character vector of spot IDs for group 2 (optional)
#' @param comparison Character; type of comparison:
#'   \itemize{
#'     \item "group1_vs_rest": Group 1 vs all other spots
#'     \item "group2_vs_rest": Group 2 vs all other spots
#'     \item "group1_vs_group2": Group 1 vs Group 2 only
#'     \item "clusters_in_group1": Find markers for each cluster in Group 1
#'     \item "clusters_in_group2": Find markers for each cluster in Group 2
#'   }
#' @param clustering_result Optional; result from run_spatial_clustering() (required for cluster comparisons)
#' @param min_pct Numeric; minimum fraction of cells expressing gene (default: 0.1)
#' @param logfc_threshold Numeric; minimum log2 fold change (default: 0.25)
#' @param top_n Integer; number of top genes to display (default: 50)
#' @param export_csv Logical; whether to export results to CSV (default: FALSE)
#' @param output_file Character; filename for CSV export (default: "DEG_results.csv")
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A data frame of DEG results with columns:
#'   \itemize{
#'     \item gene: Gene name
#'     \item avg_log2FC: Average log2 fold change
#'     \item pct.1: Percentage of cells expressing in group 1
#'     \item pct.2: Percentage of cells expressing in group 2
#'     \item p_val: P-value
#'     \item p_val_adj: Adjusted p-value
#'   }
#'
#' @import Seurat
#' @importFrom Seurat FindMarkers FindAllMarkers Idents
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Group 1 vs Rest
#' group1_spots <- draw_ROI(seurat_obj, sample_name = "Group1")
#' degs <- run_deg_analysis(seurat_obj, group1_spots, comparison = "group1_vs_rest")
#'
#' # Example 2: Group 1 vs Group 2
#' group1 <- draw_ROI(seurat_obj, sample_name = "Group1")
#' group2 <- draw_ROI(seurat_obj, sample_name = "Group2")
#' degs <- run_deg_analysis(seurat_obj, group1, group2, comparison = "group1_vs_group2")
#'
#' # Example 3: Clusters within Group 1
#' result <- run_spatial_clustering(subset_obj)
#' degs <- run_deg_analysis(seurat_obj, group1_spots,
#'                          comparison = "clusters_in_group1",
#'                          clustering_result = result)
#' }

run_deg_analysis <- function(seurat_obj,
                             group1_spots = NULL,
                             group2_spots = NULL,
                             comparison = "group1_vs_rest",
                             clustering_result = NULL,
                             min_pct = 0.1,
                             logfc_threshold = 0.25,
                             top_n = 50,
                             export_csv = FALSE,
                             output_file = "DEG_results.csv",
                             verbose = TRUE) {

  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  valid_comparisons <- c("group1_vs_rest", "group2_vs_rest", "group1_vs_group2",
                         "clusters_in_group1", "clusters_in_group2")
  if (!comparison %in% valid_comparisons) {
    stop("comparison must be one of: ", paste(valid_comparisons, collapse = ", "))
  }

  # Check if cluster-based analysis
  is_cluster_analysis <- grepl("clusters_in_group", comparison)

  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("🔬 Running Differential Expression Analysis")
    message(paste(rep("=", 70), collapse = ""))
    message("Comparison type: ", comparison)
  }


  # GROUP-BASED ANALYSIS

  if (!is_cluster_analysis) {

    # Validate group inputs
    if (comparison == "group1_vs_rest" && (is.null(group1_spots) || length(group1_spots) == 0)) {
      stop("group1_spots cannot be empty for group1_vs_rest comparison")
    }

    if (comparison == "group2_vs_rest" && (is.null(group2_spots) || length(group2_spots) == 0)) {
      stop("group2_spots cannot be empty for group2_vs_rest comparison")
    }

    if (comparison == "group1_vs_group2") {
      if (is.null(group1_spots) || length(group1_spots) == 0) {
        stop("group1_spots cannot be empty for group1_vs_group2 comparison")
      }
      if (is.null(group2_spots) || length(group2_spots) == 0) {
        stop("group2_spots cannot be empty for group1_vs_group2 comparison")
      }
    }

    if (verbose) {
      if (!is.null(group1_spots)) message("Group 1: ", length(group1_spots), " spots")
      if (!is.null(group2_spots)) message("Group 2: ", length(group2_spots), " spots")
      message("")
    }

    # Create temporary identities
    temp_idents <- rep("Other", ncol(seurat_obj))
    names(temp_idents) <- colnames(seurat_obj)

    if (!is.null(group1_spots)) {
      # Ensure group1_spots are valid cell names in the object
      valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]
      if (length(valid_group1) == 0) {
        # Provide helpful diagnostic information
        cat("\n=== DIAGNOSTIC INFORMATION ===\n")
        cat("group1_spots format (first 5):\n")
        print(head(group1_spots, 5))
        cat("\nSeurat cell names format (first 5):\n")
        print(head(colnames(seurat_obj), 5))
        cat("\n")

        # Try automatic fix
        cat("Attempting automatic fix...\n")

        # Fix 1: Convert to character
        group1_spots <- as.character(group1_spots)
        valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]

        if (length(valid_group1) == 0) {
          # Fix 2: Strip sample name prefix (e.g., "P12_N_AAACAAGTATCTCCCA-1" -> "AAACAAGTATCTCCCA-1")
          if (any(grepl("_", group1_spots))) {
            cat("Detected prefixes, attempting to strip sample names...\n")
            # Try to find the pattern: PREFIX_BARCODE
            # Extract everything after the last underscore before the barcode
            group1_spots_stripped <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", group1_spots)
            valid_group1 <- group1_spots_stripped[group1_spots_stripped %in% colnames(seurat_obj)]

            if (length(valid_group1) > 0) {
              cat("✓ Fixed by stripping sample name prefix\n")
              cat("  Example: ", head(group1_spots, 1), " -> ", head(group1_spots_stripped, 1), "\n\n")
              group1_spots <- group1_spots_stripped
            }
          }
        } else {
          cat("✓ Fixed by converting to character\n\n")
        }

        if (length(valid_group1) == 0) {
          # Fix 3: Check if numeric indices
          if (all(grepl("^[0-9]+$", group1_spots))) {
            cat("Detected numeric indices, converting to cell names...\n")
            indices <- as.numeric(group1_spots)
            if (all(indices > 0 & indices <= ncol(seurat_obj))) {
              group1_spots <- colnames(seurat_obj)[indices]
              valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]
              if (length(valid_group1) > 0) {
                cat("✓ Fixed by converting indices to cell names\n\n")
              }
            }
          }
        }

        if (length(valid_group1) == 0) {
          cat("Could not automatically fix spot ID mismatch.\n")
          cat("\nSuggested manual fix:\n")
          cat("# Strip the sample name prefix:\n")
          cat("group1_spots_fixed <- gsub('^.*_([ACGT]{16}-[0-9]+)$', '\\\\1', group1_spots)\n")
          cat("# Then re-run with the fixed IDs\n\n")
          stop("None of the group1_spots are found in the Seurat object")
        }
      }
      if (length(valid_group1) < length(group1_spots)) {
        warning(length(group1_spots) - length(valid_group1), " spots from group1_spots not found in object")
      }
      temp_idents[valid_group1] <- "Group1"
      group1_spots <- valid_group1  # Update to valid spots only
    }

    if (!is.null(group2_spots)) {
      # Ensure group2_spots are valid cell names in the object
      valid_group2 <- group2_spots[group2_spots %in% colnames(seurat_obj)]
      if (length(valid_group2) == 0) {
        # Try automatic fix
        cat("Attempting to fix group2_spots...\n")

        # Convert to character
        group2_spots <- as.character(group2_spots)
        valid_group2 <- group2_spots[group2_spots %in% colnames(seurat_obj)]

        if (length(valid_group2) == 0) {
          # Strip sample name prefix
          if (any(grepl("_", group2_spots))) {
            group2_spots_stripped <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", group2_spots)
            valid_group2 <- group2_spots_stripped[group2_spots_stripped %in% colnames(seurat_obj)]
            if (length(valid_group2) > 0) {
              cat("✓ Fixed group2_spots by stripping prefix\n\n")
              group2_spots <- group2_spots_stripped
            }
          }
        }

        if (length(valid_group2) == 0) {
          stop("None of the group2_spots are found in the Seurat object")
        }
      }
      if (length(valid_group2) < length(group2_spots)) {
        warning(length(group2_spots) - length(valid_group2), " spots from group2_spots not found in object")
      }
      temp_idents[valid_group2] <- "Group2"
      group2_spots <- valid_group2  # Update to valid spots only
    }

    # Determine which cells to use
    if (comparison == "group1_vs_group2") {
      cells_to_use <- c(group1_spots, group2_spots)
    } else {
      cells_to_use <- colnames(seurat_obj)
    }

    if (verbose) message("Running FindMarkers...")

    # Subset and set identities
    temp_seurat <- subset(seurat_obj, cells = cells_to_use)
    temp_seurat$temp_ident <- factor(temp_idents[cells_to_use])
    Idents(temp_seurat) <- temp_seurat$temp_ident

    # Verify identities were set correctly
    if (verbose) {
      message("  Identity levels: ", paste(levels(Idents(temp_seurat)), collapse = ", "))
      message("  Cells per identity: ", paste(table(Idents(temp_seurat)), collapse = ", "))
    }

    # Run DEG analysis
    if (comparison == "group1_vs_group2") {
      markers <- FindMarkers(temp_seurat,
                             ident.1 = "Group1",
                             ident.2 = "Group2",
                             verbose = FALSE,
                             min.pct = min_pct,
                             logfc.threshold = logfc_threshold)
    } else if (comparison == "group1_vs_rest") {
      markers <- FindMarkers(temp_seurat,
                             ident.1 = "Group1",
                             ident.2 = "Other",
                             verbose = FALSE,
                             min.pct = min_pct,
                             logfc.threshold = logfc_threshold)
    } else if (comparison == "group2_vs_rest") {
      markers <- FindMarkers(temp_seurat,
                             ident.1 = "Group2",
                             ident.2 = "Other",
                             verbose = FALSE,
                             min.pct = min_pct,
                             logfc.threshold = logfc_threshold)
    }

    # Format results
    markers$gene <- rownames(markers)
    markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
    markers <- markers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]

    if (verbose) {
      message("  ✓ Found ", nrow(markers), " differentially expressed genes")
      message("")
      message("Top 10 DEGs by adjusted p-value:")
      print_table <- head(markers, 10)
      for (i in 1:nrow(print_table)) {
        message(sprintf("  %2d. %s (log2FC: %.3f, p_adj: %.2e)",
                        i, print_table$gene[i], print_table$avg_log2FC[i], print_table$p_val_adj[i]))
      }
    }


    # CLUSTER-BASED ANALYSIS

  } else {

    # Validate clustering result
    if (is.null(clustering_result)) {
      stop("clustering_result is required for cluster-based comparison. ",
           "Run run_spatial_clustering() first.")
    }

    if (!all(c("clusters", "seurat_obj") %in% names(clustering_result))) {
      stop("clustering_result must be output from run_spatial_clustering()")
    }

    # Determine which group to analyze
    target_group <- if (comparison == "clusters_in_group1") {
      "group1"
    } else {
      "group2"
    }

    target_spots <- if (target_group == "group1") {
      group1_spots
    } else {
      group2_spots
    }

    if (is.null(target_spots) || length(target_spots) == 0) {
      group_name <- if (target_group == "group1") "Group 1" else "Group 2"
      stop(group_name, " spots cannot be empty for cluster-based analysis")
    }

    if (verbose) {
      message("Target group: ", if (target_group == "group1") "Group 1" else "Group 2")
      message("Number of spots: ", length(target_spots))
    }

    # Get clusters
    clusters <- clustering_result$clusters

    # Check if clustering was done on the target spots
    clustered_spots <- names(clusters)
    if (!all(clustered_spots %in% target_spots)) {
      warning("Clustering appears to be done on different spots than specified group. ",
              "Proceeding with clustered spots...")
    }

    if (verbose) message("Number of clusters: ", clustering_result$n_clusters, "\n")

    # Subset to clustered spots
    temp_seurat <- subset(seurat_obj, cells = clustered_spots)
    temp_seurat$cluster_ident <- as.factor(clusters)
    Idents(temp_seurat) <- "cluster_ident"

    if (verbose) message("Running FindAllMarkers (one vs rest for each cluster)...")

    # Find markers for all clusters
    markers <- FindAllMarkers(temp_seurat,
                              only.pos = FALSE,
                              min.pct = min_pct,
                              logfc.threshold = logfc_threshold,
                              verbose = FALSE)

    # Format results
    markers$gene <- rownames(markers)
    markers <- markers[order(markers$cluster, markers$p_val_adj, -abs(markers$avg_log2FC)), ]
    markers <- markers[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]

    if (verbose) {
      message("  ✓ Found markers for ", clustering_result$n_clusters, " clusters")
      message("")
      message("Summary by cluster:")
      cluster_counts <- table(markers$cluster)
      for (cl in names(cluster_counts)) {
        top_gene <- markers[markers$cluster == cl, ][1, ]
        message(sprintf("  Cluster %s: %d DEGs (top: %s, log2FC: %.3f)",
                        cl, cluster_counts[cl], top_gene$gene, top_gene$avg_log2FC))
      }
    }
  }

  # Export to CSV if requested
  if (export_csv) {
    write.csv(markers, file = output_file, row.names = FALSE)
    if (verbose) message("\n✓ Results exported to: ", output_file)
  }

  if (verbose) {
    message("")
    message(paste(rep("=", 70), collapse = ""))
    message("✅ DEG analysis complete!")
    message("Total DEGs: ", nrow(markers))
    message("Access results: View top results with head(result, ", top_n, ")")
    message(paste(rep("=", 70), collapse = ""))
    message("")
  }

  return(markers)
}


#' Quick DEG Analysis Wrapper
#'
#' @description
#' Simplified wrapper for common DEG comparisons between two Seurat subsets
#'
#' @param subset1 Seurat object for group 1
#' @param subset2 Seurat object for group 2
#' @param parent_seurat Optional; full Seurat object (if not provided, uses subset1)
#' @param comparison Character; "group1_vs_group2" or "group1_vs_rest" (default: "group1_vs_group2")
#'
#' @return Data frame of DEG results
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare two ROI subsets
#' roi1 <- draw_ROI(seurat_obj, sample_name = "ROI1")
#' roi2 <- draw_ROI(seurat_obj, sample_name = "ROI2")
#' subset1 <- subset(seurat_obj, cells = roi1)
#' subset2 <- subset(seurat_obj, cells = roi2)
#'
#' degs <- compare_subsets(subset1, subset2, parent_seurat = seurat_obj)
#' }

compare_subsets <- function(subset1,
                            subset2 = NULL,
                            parent_seurat = NULL,
                            comparison = "group1_vs_group2") {

  # Get spot IDs
  group1_spots <- colnames(subset1)

  if (is.null(parent_seurat)) {
    parent_seurat <- subset1
  }

  if (is.null(subset2)) {
    # Group 1 vs Rest
    return(run_deg_analysis(
      parent_seurat,
      group1_spots = group1_spots,
      comparison = "group1_vs_rest"
    ))
  } else {
    # Group 1 vs Group 2
    group2_spots <- colnames(subset2)
    return(run_deg_analysis(
      parent_seurat,
      group1_spots = group1_spots,
      group2_spots = group2_spots,
      comparison = "group1_vs_group2"
    ))
  }
}


#' Compare Specific Clusters
#'
#' @description
#' Compare specific clusters from clustering results
#'
#' @param seurat_obj Full Seurat object
#' @param clustering_result Result from run_spatial_clustering()
#' @param cluster1 Character or numeric; first cluster ID
#' @param cluster2 Character or numeric; second cluster ID
#'
#' @return Data frame of DEG results
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_spatial_clustering(seurat_obj)
#' degs <- compare_clusters(seurat_obj, result, cluster1 = "0", cluster2 = "1")
#' }

compare_clusters <- function(seurat_obj,
                             clustering_result,
                             cluster1,
                             cluster2) {

  if (!all(c("cluster_ids", "clusters") %in% names(clustering_result))) {
    stop("clustering_result must be output from run_spatial_clustering()")
  }

  cluster1 <- as.character(cluster1)
  cluster2 <- as.character(cluster2)

  # Get spot IDs for each cluster
  cluster1_spots <- clustering_result$cluster_ids[[cluster1]]
  cluster2_spots <- clustering_result$cluster_ids[[cluster2]]

  if (is.null(cluster1_spots) || length(cluster1_spots) == 0) {
    stop("Cluster ", cluster1, " not found in clustering results")
  }

  if (is.null(cluster2_spots) || length(cluster2_spots) == 0) {
    stop("Cluster ", cluster2, " not found in clustering results")
  }

  message("\nComparing Cluster ", cluster1, " (", length(cluster1_spots), " spots) vs ",
          "Cluster ", cluster2, " (", length(cluster2_spots), " spots)\n")

  # Run DEG analysis
  return(run_deg_analysis(
    seurat_obj,
    group1_spots = cluster1_spots,
    group2_spots = cluster2_spots,
    comparison = "group1_vs_group2"
  ))
}


#' Get Top DEGs
#'
#' @description
#' Extract top N differentially expressed genes from DEG results
#'
#' @param deg_results Data frame from run_deg_analysis()
#' @param n Integer; number of top genes to return (default: 20)
#' @param by Character; column to sort by - "p_val_adj", "log2FC", or "both" (default: "p_val_adj")
#' @param direction Character; "up", "down", or "both" (default: "both")
#'
#' @return Data frame of top DEGs
#' @export

get_top_degs <- function(deg_results,
                         n = 20,
                         by = "p_val_adj",
                         direction = "both") {

  if (!by %in% c("p_val_adj", "log2FC", "both")) {
    stop("by must be 'p_val_adj', 'log2FC', or 'both'")
  }

  if (!direction %in% c("up", "down", "both")) {
    stop("direction must be 'up', 'down', or 'both'")
  }

  # Filter by direction if specified
  if (direction == "up") {
    deg_results <- deg_results[deg_results$avg_log2FC > 0, ]
  } else if (direction == "down") {
    deg_results <- deg_results[deg_results$avg_log2FC < 0, ]
  }

  # Sort by specified column
  if (by == "p_val_adj") {
    deg_results <- deg_results[order(deg_results$p_val_adj), ]
  } else if (by == "log2FC") {
    deg_results <- deg_results[order(-abs(deg_results$avg_log2FC)), ]
  } else if (by == "both") {
    deg_results <- deg_results[order(deg_results$p_val_adj, -abs(deg_results$avg_log2FC)), ]
  }

  return(head(deg_results, n))
}


#' Visualize DEG Results
#'
#' @description
#' Create volcano plot or bar plot of DEG results
#'
#' @param deg_results Data frame from run_deg_analysis()
#' @param plot_type Character; "volcano" or "bar" (default: "volcano")
#' @param top_n Integer; for bar plot, number of top genes to show (default: 20)
#' @param fc_threshold Numeric; log2 fold change threshold for coloring (default: 0.5)
#' @param pval_threshold Numeric; -log10(p_val_adj) threshold for coloring (default: 1.3, i.e., p=0.05)
#'
#' @import ggplot2
#' @export

visualize_degs <- function(deg_results,
                           plot_type = "volcano",
                           top_n = 20,
                           fc_threshold = 0.5,
                           pval_threshold = 1.3) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }

  if (plot_type == "volcano") {
    # Volcano plot
    deg_results$neg_log10_p <- -log10(deg_results$p_val_adj)
    deg_results$significant <- ifelse(
      abs(deg_results$avg_log2FC) > fc_threshold & deg_results$neg_log10_p > pval_threshold,
      "Significant", "Not Significant"
    )

    p <- ggplot(deg_results, aes(x = avg_log2FC, y = neg_log10_p, color = significant)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
      geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "blue") +
      labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)",
           title = "Volcano Plot of DEGs") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "top"
      )

  } else if (plot_type == "bar") {
    # Bar plot of top genes
    top_genes <- get_top_degs(deg_results, n = top_n, by = "p_val_adj")
    top_genes$gene <- factor(top_genes$gene, levels = rev(top_genes$gene))

    p <- ggplot(top_genes, aes(x = gene, y = avg_log2FC, fill = avg_log2FC > 0)) +
      geom_col() +
      scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                        labels = c("Down", "Up")) +
      coord_flip() +
      labs(x = "Gene", y = "Log2 Fold Change",
           title = paste("Top", top_n, "DEGs"),
           fill = "Direction") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "top"
      )

  } else {
    stop("plot_type must be 'volcano' or 'bar'")
  }

  print(p)
  return(invisible(p))
}
