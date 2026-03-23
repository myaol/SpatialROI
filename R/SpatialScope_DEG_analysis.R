#' Perform Differential Expression Analysis on Spatial Transcriptomics Data
#' (with optional Moran's I spatial autocorrelation)
#'
#' @description
#' Run differential expression analysis comparing different groups of spots.
#' Supports group vs rest, group vs group, and cluster-based comparisons.
#' Optionally computes Moran's I spatial autocorrelation (permutation-based,
#' BH-corrected) for the top DEGs to identify spatially structured genes.
#'
#' @param seurat_obj        A Seurat spatial object (full dataset)
#' @param group1_spots      Character vector of spot IDs for group 1
#' @param group2_spots      Character vector of spot IDs for group 2 (optional)
#' @param comparison        Character; type of comparison:
#'   \itemize{
#'     \item "group1_vs_rest":    Group 1 vs all other spots
#'     \item "group2_vs_rest":    Group 2 vs all other spots
#'     \item "group1_vs_group2":  Group 1 vs Group 2 only
#'     \item "clusters_in_group1": Find markers for each cluster in Group 1
#'     \item "clusters_in_group2": Find markers for each cluster in Group 2
#'   }
#' @param clustering_result Optional; result from run_spatial_clustering()
#'                          (required for cluster comparisons)
#' @param min_pct           Numeric; minimum fraction of cells expressing gene (default: 0.1)
#' @param logfc_threshold   Numeric; minimum log2 fold change (default: 0.25)
#' @param top_n             Integer; number of top genes to display (default: 50)
#' @param export_csv        Logical; whether to export results to CSV (default: FALSE)
#' @param output_file       Character; filename for CSV export (default: "DEG_results.csv")
#' @param verbose           Logical; whether to print progress messages (default: TRUE)
#'
#' @param add_morans        Logical; compute Moran's I for top DEGs (default: FALSE)
#'                          High Moran_I = spatially clustered gene ("blob" pattern)
#'                          Low  Moran_I = randomly scattered gene ("salt-and-pepper")
#' @param morans_top_n      Integer; number of top DEGs to score with Moran's I (default: 100)
#'                          Capped to avoid excessive compute time on large datasets
#' @param morans_k          Integer; number of nearest neighbours for spatial weight matrix (default: 6)
#'                          k=6 approximates the hexagonal Visium grid neighbourhood; use 4-8
#' @param morans_subset_obj Logical; compute Moran's I only on the relevant group spots
#'                          rather than the full object (default: TRUE — recommended for ROI analyses)
#'                          Rationale: global coordinates dilute local ROI structure
#'
#' @return A data frame of DEG results with columns:
#'   \itemize{
#'     \item gene:          Gene name
#'     \item avg_log2FC:    Average log2 fold change
#'     \item pct.1:         Percentage of cells expressing in group 1
#'     \item pct.2:         Percentage of cells expressing in group 2
#'     \item p_val:         P-value
#'     \item p_val_adj:     Adjusted p-value
#'   }
#'   If add_morans = TRUE, additional columns:
#'   \itemize{
#'     \item Moran_I:       Moran's I statistic (-1 to 1; >0 = spatially clustered)
#'     \item Moran_pval:    Permutation-based p-value (999 simulations)
#'     \item Moran_padj:    BH-corrected p-value across all tested genes
#'     \item spatial_sig:   Logical; TRUE if Moran_padj < 0.05
#'     \item spatial_class: Factor; "Spatially structured" / "Not structured" / "Not tested"
#'   }
#'   Sort order when Moran's I is included:
#'     Primary p_val_adj ascending, secondary Moran_I descending.
#'     Top rows = strong DEGs that are ALSO spatially structured (highest priority for pathology).
#'
#' @import Seurat
#' @importFrom Seurat FindMarkers FindAllMarkers Idents
#' @importFrom spdep knearneigh knn2nb nb2listw moran.mc
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Group 1 vs Rest (no Moran's I)
#' group1_spots <- draw_ROI(seurat_obj, sample_name = "Group1")
#' degs <- run_deg_analysis(seurat_obj, group1_spots, comparison = "group1_vs_rest")
#'
#' # Example 2: Group 1 vs Rest WITH Moran's I
#' degs <- run_deg_analysis(seurat_obj, group1_spots,
#'                          comparison   = "group1_vs_rest",
#'                          add_morans   = TRUE,
#'                          morans_top_n = 100,
#'                          morans_k     = 6)
#'
#' # Example 3: Group 1 vs Group 2
#' group1 <- draw_ROI(seurat_obj, sample_name = "Group1")
#' group2 <- draw_ROI(seurat_obj, sample_name = "Group2")
#' degs <- run_deg_analysis(seurat_obj, group1, group2,
#'                          comparison = "group1_vs_group2",
#'                          add_morans = TRUE)
#'
#' # Example 4: Clusters within Group 1
#' result <- run_spatial_clustering(subset_obj)
#' degs <- run_deg_analysis(seurat_obj, group1_spots,
#'                          comparison        = "clusters_in_group1",
#'                          clustering_result = result)
#'
#' # Filter high-priority genes (significant DEG + spatially structured)
#' high_priority <- degs[degs$p_val_adj < 0.05 &
#'                       degs$spatial_class == "Spatially structured", ]
#' }

run_deg_analysis <- function(seurat_obj,
                             group1_spots       = NULL,
                             group2_spots       = NULL,
                             comparison         = "group1_vs_rest",
                             clustering_result  = NULL,
                             min_pct            = 0.1,
                             logfc_threshold    = 0.25,
                             top_n              = 50,
                             export_csv         = FALSE,
                             output_file        = "DEG_results.csv",
                             verbose            = TRUE,
                             # ── Moran's I options ──────────────────────────────────────────
                             add_morans         = FALSE,
                             morans_top_n       = 100,
                             morans_k           = 6,
                             morans_subset_obj  = TRUE) {

  # ══════════════════════════════════════════════════════════════════════════
  # INPUT VALIDATION
  # ══════════════════════════════════════════════════════════════════════════

  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  valid_comparisons <- c("group1_vs_rest", "group2_vs_rest", "group1_vs_group2",
                         "clusters_in_group1", "clusters_in_group2")
  if (!comparison %in% valid_comparisons) {
    stop("comparison must be one of: ", paste(valid_comparisons, collapse = ", "))
  }

  is_cluster_analysis <- grepl("clusters_in_group", comparison)

  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("🔬 Running Differential Expression Analysis")
    message(paste(rep("=", 70), collapse = ""))
    message("Comparison type: ", comparison)
  }

  # ══════════════════════════════════════════════════════════════════════════
  # GROUP-BASED ANALYSIS
  # ══════════════════════════════════════════════════════════════════════════

  if (!is_cluster_analysis) {

    # Validate group inputs
    if (comparison == "group1_vs_rest" && (is.null(group1_spots) || length(group1_spots) == 0)) {
      stop("group1_spots cannot be empty for group1_vs_rest comparison")
    }
    if (comparison == "group2_vs_rest" && (is.null(group2_spots) || length(group2_spots) == 0)) {
      stop("group2_spots cannot be empty for group2_vs_rest comparison")
    }
    if (comparison == "group1_vs_group2") {
      if (is.null(group1_spots) || length(group1_spots) == 0)
        stop("group1_spots cannot be empty for group1_vs_group2 comparison")
      if (is.null(group2_spots) || length(group2_spots) == 0)
        stop("group2_spots cannot be empty for group1_vs_group2 comparison")
    }

    if (verbose) {
      if (!is.null(group1_spots)) message("Group 1: ", length(group1_spots), " spots")
      if (!is.null(group2_spots)) message("Group 2: ", length(group2_spots), " spots")
      message("")
    }

    # ── Create temporary identities ────────────────────────────────────────
    temp_idents <- rep("Other", ncol(seurat_obj))
    names(temp_idents) <- colnames(seurat_obj)

    # ── Validate & fix group1_spots ────────────────────────────────────────
    if (!is.null(group1_spots)) {
      valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]

      if (length(valid_group1) == 0) {
        cat("\n=== DIAGNOSTIC INFORMATION ===\n")
        cat("group1_spots format (first 5):\n"); print(head(group1_spots, 5))
        cat("\nSeurat cell names format (first 5):\n"); print(head(colnames(seurat_obj), 5))
        cat("\nAttempting automatic fix...\n")

        # Fix 1: Convert to character
        group1_spots <- as.character(group1_spots)
        valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]

        if (length(valid_group1) == 0) {
          # Fix 2: Strip sample name prefix
          if (any(grepl("_", group1_spots))) {
            cat("Detected prefixes, attempting to strip sample names...\n")
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
          # Fix 3: Numeric indices
          if (all(grepl("^[0-9]+$", group1_spots))) {
            cat("Detected numeric indices, converting to cell names...\n")
            indices <- as.numeric(group1_spots)
            if (all(indices > 0 & indices <= ncol(seurat_obj))) {
              group1_spots <- colnames(seurat_obj)[indices]
              valid_group1 <- group1_spots[group1_spots %in% colnames(seurat_obj)]
              if (length(valid_group1) > 0) cat("✓ Fixed by converting indices to cell names\n\n")
            }
          }
        }

        if (length(valid_group1) == 0) {
          cat("Could not automatically fix spot ID mismatch.\n")
          cat("\nSuggested manual fix:\n")
          cat("group1_spots_fixed <- gsub('^.*_([ACGT]{16}-[0-9]+)$', '\\\\1', group1_spots)\n\n")
          stop("None of the group1_spots are found in the Seurat object")
        }
      }

      if (length(valid_group1) < length(group1_spots)) {
        warning(length(group1_spots) - length(valid_group1),
                " spots from group1_spots not found in object")
      }
      temp_idents[valid_group1] <- "Group1"
      group1_spots <- valid_group1
    }

    # ── Validate & fix group2_spots ────────────────────────────────────────
    if (!is.null(group2_spots)) {
      valid_group2 <- group2_spots[group2_spots %in% colnames(seurat_obj)]

      if (length(valid_group2) == 0) {
        cat("Attempting to fix group2_spots...\n")
        group2_spots <- as.character(group2_spots)
        valid_group2 <- group2_spots[group2_spots %in% colnames(seurat_obj)]

        if (length(valid_group2) == 0 && any(grepl("_", group2_spots))) {
          group2_spots_stripped <- gsub("^.*_([ACGT]{16}-[0-9]+)$", "\\1", group2_spots)
          valid_group2 <- group2_spots_stripped[group2_spots_stripped %in% colnames(seurat_obj)]
          if (length(valid_group2) > 0) {
            cat("✓ Fixed group2_spots by stripping prefix\n\n")
            group2_spots <- group2_spots_stripped
          }
        }

        if (length(valid_group2) == 0)
          stop("None of the group2_spots are found in the Seurat object")
      }

      if (length(valid_group2) < length(group2_spots)) {
        warning(length(group2_spots) - length(valid_group2),
                " spots from group2_spots not found in object")
      }
      temp_idents[valid_group2] <- "Group2"
      group2_spots <- valid_group2
    }

    # ── Subset & run FindMarkers ───────────────────────────────────────────
    cells_to_use <- if (comparison == "group1_vs_group2") {
      c(group1_spots, group2_spots)
    } else {
      colnames(seurat_obj)
    }

    if (verbose) message("Running FindMarkers...")

    temp_seurat <- subset(seurat_obj, cells = cells_to_use)
    temp_seurat$temp_ident <- factor(temp_idents[cells_to_use])
    Idents(temp_seurat) <- temp_seurat$temp_ident

    if (verbose) {
      message("  Identity levels: ", paste(levels(Idents(temp_seurat)), collapse = ", "))
      message("  Cells per identity: ", paste(table(Idents(temp_seurat)), collapse = ", "))
    }

    if (comparison == "group1_vs_group2") {
      markers <- FindMarkers(temp_seurat, ident.1 = "Group1", ident.2 = "Group2",
                             verbose = FALSE, min.pct = min_pct, logfc.threshold = logfc_threshold)
    } else if (comparison == "group1_vs_rest") {
      markers <- FindMarkers(temp_seurat, ident.1 = "Group1", ident.2 = "Other",
                             verbose = FALSE, min.pct = min_pct, logfc.threshold = logfc_threshold)
    } else if (comparison == "group2_vs_rest") {
      markers <- FindMarkers(temp_seurat, ident.1 = "Group2", ident.2 = "Other",
                             verbose = FALSE, min.pct = min_pct, logfc.threshold = logfc_threshold)
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
      for (i in seq_len(nrow(print_table))) {
        message(sprintf("  %2d. %s (log2FC: %.3f, p_adj: %.2e)",
                        i, print_table$gene[i],
                        print_table$avg_log2FC[i],
                        print_table$p_val_adj[i]))
      }
    }

  # ══════════════════════════════════════════════════════════════════════════
  # CLUSTER-BASED ANALYSIS
  # ══════════════════════════════════════════════════════════════════════════

  } else {

    if (is.null(clustering_result)) {
      stop("clustering_result is required for cluster-based comparison. ",
           "Run run_spatial_clustering() first.")
    }
    if (!all(c("clusters", "seurat_obj") %in% names(clustering_result))) {
      stop("clustering_result must be output from run_spatial_clustering()")
    }

    target_group <- if (comparison == "clusters_in_group1") "group1" else "group2"
    target_spots <- if (target_group == "group1") group1_spots else group2_spots

    if (is.null(target_spots) || length(target_spots) == 0) {
      stop(if (target_group == "group1") "Group 1" else "Group 2",
           " spots cannot be empty for cluster-based analysis")
    }

    if (verbose) {
      message("Target group: ", if (target_group == "group1") "Group 1" else "Group 2")
      message("Number of spots: ", length(target_spots))
    }

    clusters        <- clustering_result$clusters
    clustered_spots <- names(clusters)

    if (!all(clustered_spots %in% target_spots)) {
      warning("Clustering appears to be done on different spots than specified group. ",
              "Proceeding with clustered spots...")
    }

    if (verbose) message("Number of clusters: ", clustering_result$n_clusters, "\n")

    temp_seurat <- subset(seurat_obj, cells = clustered_spots)
    temp_seurat$cluster_ident <- as.factor(clusters)
    Idents(temp_seurat) <- "cluster_ident"

    if (verbose) message("Running FindAllMarkers (one vs rest for each cluster)...")

    markers <- FindAllMarkers(temp_seurat,
                              only.pos = FALSE,
                              min.pct = min_pct,
                              logfc.threshold = logfc_threshold,
                              verbose = FALSE)

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

  # ══════════════════════════════════════════════════════════════════════════
  # MORAN'S I SPATIAL AUTOCORRELATION (optional)
  # Permutation-based (999 simulations) + BH correction — matches Squidpy standard
  # Interpretation: Moran_I > 0 = "blob" (spatially clustered, tissue-structure-like)
  #                 Moran_I ~ 0 = "salt-and-pepper" (spatially random)
  # ══════════════════════════════════════════════════════════════════════════

  if (add_morans) {

    # ── 1. Dependency check ────────────────────────────────────────────────
    if (!requireNamespace("spdep", quietly = TRUE)) {
      warning(
        "Package 'spdep' is required for Moran's I. ",
        "Install with: install.packages('spdep')\n",
        "Skipping spatial autocorrelation step."
      )
    } else {

      if (verbose) {
        message("")
        message(paste(rep("-", 70), collapse = ""))
        message("📍 Computing Moran's I spatial autocorrelation")
        message(paste(rep("-", 70), collapse = ""))
      }

      # ── 2. Decide scope of spots for Moran's I ────────────────────────
      #       Default: relevant group spots only — keeps analysis within the ROI
      if (morans_subset_obj) {
        if (comparison %in% c("group1_vs_rest", "group1_vs_group2") && !is.null(group1_spots)) {
          moran_spots <- group1_spots
          if (verbose) message("  Moran scope: Group 1 spots (n=", length(moran_spots), ")")
        } else if (comparison == "group2_vs_rest" && !is.null(group2_spots)) {
          moran_spots <- group2_spots
          if (verbose) message("  Moran scope: Group 2 spots (n=", length(moran_spots), ")")
        } else {
          moran_spots <- colnames(seurat_obj)
          if (verbose) message("  Moran scope: All spots (cluster/rest comparison)")
        }
      } else {
        moran_spots <- colnames(seurat_obj)
        if (verbose) message("  Moran scope: Full object (n=", length(moran_spots), ")")
      }

      # ── 3. Minimum spot guard ──────────────────────────────────────────
      MIN_SPOTS_MORAN <- 30
      if (length(moran_spots) < MIN_SPOTS_MORAN) {
        warning(
          "Only ", length(moran_spots), " spots available for Moran's I ",
          "(minimum recommended: ", MIN_SPOTS_MORAN, "). ",
          "Results may be unreliable. Skipping."
        )
      } else {

        # ── 4. Get spatial coordinates ──────────────────────────────────
        coords_full <- Seurat::GetTissueCoordinates(seurat_obj)
        coords      <- coords_full[rownames(coords_full) %in% moran_spots, ]
        coords      <- coords[, c("imagerow", "imagecol")]  # row = y, col = x
        coords      <- as.matrix(coords)

        if (nrow(coords) == 0) {
          warning("Could not match spot IDs to tissue coordinates. Skipping Moran's I.")
          coords <- NULL
        }

        if (!is.null(coords)) {

          # ── 5. Build spatial weight matrix (kNN) ─────────────────────
          if (verbose) message("  Building spatial weight matrix (k=", morans_k, " neighbours)...")

          effective_k <- min(morans_k, nrow(coords) - 1)
          if (effective_k < morans_k && verbose)
            message("  (k reduced to ", effective_k, " due to small spot count)")

          knn_obj   <- spdep::knearneigh(coords, k = effective_k)
          nb_obj    <- spdep::knn2nb(knn_obj)
          listw_obj <- spdep::nb2listw(nb_obj, style = "W", zero.policy = TRUE)

          # ── 6. Select candidate genes ─────────────────────────────────
          all_genes       <- markers$gene
          candidate_genes <- head(
            intersect(all_genes, rownames(seurat_obj@assays$Spatial$data)),
            morans_top_n
          )

          if (length(candidate_genes) == 0) {
            warning("No DEG genes found in Spatial assay data. Skipping Moran's I.")
          } else {

            if (verbose)
              message("  Testing ", length(candidate_genes),
                      " genes (top ", morans_top_n, " DEGs by p_val_adj)...")

            # ── 7. Extract expression (target spots only) ──────────────
            expr_matrix <- as.matrix(
              seurat_obj@assays$Spatial$data[candidate_genes, rownames(coords), drop = FALSE]
            )

            # ── 8. Permutation-based Moran's I per gene ───────────────
            #       nsim=999: standard for permutation tests
            #       alternative="greater": test for positive autocorrelation (blob)
            moran_results <- lapply(candidate_genes, function(gene) {
              x <- expr_matrix[gene, ]

              # Zero-variance genes → Moran's I undefined
              if (var(x) == 0) {
                return(data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_))
              }

              tryCatch({
                mt <- spdep::moran.mc(
                  x,
                  listw       = listw_obj,
                  nsim        = 999,
                  zero.policy = TRUE,
                  alternative = "greater"
                )
                data.frame(
                  gene       = gene,
                  Moran_I    = mt$statistic,
                  Moran_pval = mt$p.value
                )
              }, error = function(e) {
                data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_)
              })
            })

            moran_df <- do.call(rbind, moran_results)

            # ── 9. BH correction across all tested genes ───────────────
            #       Matches Squidpy's FDR-correction standard
            moran_df$Moran_padj <- p.adjust(moran_df$Moran_pval, method = "BH")

            # ── 10. Merge into markers table ───────────────────────────
            markers <- merge(markers, moran_df, by = "gene", all.x = TRUE)

            # ── 11. Annotate spatial significance ──────────────────────
            markers$spatial_sig <- ifelse(
              !is.na(markers$Moran_padj) & markers$Moran_padj < 0.05,
              TRUE, FALSE
            )

            markers$spatial_class <- dplyr::case_when(
              is.na(markers$Moran_I)       ~ "Not tested",
              markers$spatial_sig == TRUE  ~ "Spatially structured",
              TRUE                         ~ "Not structured"
            )
            markers$spatial_class <- factor(
              markers$spatial_class,
              levels = c("Spatially structured", "Not structured", "Not tested")
            )

            # ── 12. Re-sort: best DEGs AND most spatially structured first
            #        Primary:   p_val_adj ↑ (most significant DEG first)
            #        Secondary: Moran_I   ↓ (most spatially clustered first)
            markers <- markers[
              order(markers$p_val_adj,
                    -replace(markers$Moran_I, is.na(markers$Moran_I), -Inf)), ]

            # ── 13. Verbose summary ────────────────────────────────────
            if (verbose) {
              n_spatial <- sum(markers$spatial_sig, na.rm = TRUE)
              n_both    <- sum(markers$p_val_adj < 0.05 & markers$spatial_sig == TRUE,
                               na.rm = TRUE)
              message("")
              message("  Moran's I summary:")
              message("  • Genes tested:                   ", length(candidate_genes))
              message("  • Spatially structured (padj<0.05): ", n_spatial)
              message("  • Sig DEG + spatially structured:  ", n_both,
                      "  ← HIGH PRIORITY for pathology")
              message("")
              message("  Top 5 by Moran's I (among tested genes):")
              top_moran <- moran_df[!is.na(moran_df$Moran_I), ]
              top_moran <- head(top_moran[order(-top_moran$Moran_I), ], 5)
              for (i in seq_len(nrow(top_moran))) {
                message(sprintf("    %d. %-12s  I = %+.3f  p_adj = %.2e",
                                i, top_moran$gene[i],
                                top_moran$Moran_I[i],
                                top_moran$Moran_padj[i]))
              }
            }
          }
        }
      }
    }
  }

  # ══════════════════════════════════════════════════════════════════════════
  # EXPORT + RETURN
  # ══════════════════════════════════════════════════════════════════════════

  if (export_csv) {
    write.csv(markers, file = output_file, row.names = FALSE)
    if (verbose) message("\n✓ Results exported to: ", output_file)
  }

  if (verbose) {
    message("")
    message(paste(rep("=", 70), collapse = ""))
    message("✅ DEG analysis complete!")
    message("Total DEGs: ", nrow(markers))
    if (add_morans && "Moran_I" %in% colnames(markers)) {
      message("Columns: gene | avg_log2FC | p_val_adj | Moran_I | Moran_pval | Moran_padj | spatial_class")
      message("Sort order: p_val_adj ↑  then  Moran_I ↓")
    }
    message("Access results: View top results with head(result, ", top_n, ")")
    message(paste(rep("=", 70), collapse = ""))
    message("")
  }

  return(markers)
}


# ══════════════════════════════════════════════════════════════════════════════
#' Quick DEG Analysis Wrapper
#'
#' @description
#' Simplified wrapper for common DEG comparisons between two Seurat subsets
#'
#' @param subset1       Seurat object for group 1
#' @param subset2       Seurat object for group 2 (optional)
#' @param parent_seurat Full Seurat object (if not provided, uses subset1)
#' @param comparison    Character; "group1_vs_group2" or "group1_vs_rest" (default: "group1_vs_group2")
#'
#' @return Data frame of DEG results
#' @export
# ══════════════════════════════════════════════════════════════════════════════

compare_subsets <- function(subset1,
                            subset2       = NULL,
                            parent_seurat = NULL,
                            comparison    = "group1_vs_group2") {

  group1_spots <- colnames(subset1)
  if (is.null(parent_seurat)) parent_seurat <- subset1

  if (is.null(subset2)) {
    return(run_deg_analysis(parent_seurat,
                            group1_spots = group1_spots,
                            comparison   = "group1_vs_rest"))
  } else {
    group2_spots <- colnames(subset2)
    return(run_deg_analysis(parent_seurat,
                            group1_spots = group1_spots,
                            group2_spots = group2_spots,
                            comparison   = "group1_vs_group2"))
  }
}


# ══════════════════════════════════════════════════════════════════════════════
#' Compare Specific Clusters
#'
#' @description
#' Compare specific clusters from clustering results
#'
#' @param seurat_obj        Full Seurat object
#' @param clustering_result Result from run_spatial_clustering()
#' @param cluster1          Character or numeric; first cluster ID
#' @param cluster2          Character or numeric; second cluster ID
#'
#' @return Data frame of DEG results
#' @export
# ══════════════════════════════════════════════════════════════════════════════

compare_clusters <- function(seurat_obj,
                             clustering_result,
                             cluster1,
                             cluster2) {

  if (!all(c("cluster_ids", "clusters") %in% names(clustering_result))) {
    stop("clustering_result must be output from run_spatial_clustering()")
  }

  cluster1 <- as.character(cluster1)
  cluster2 <- as.character(cluster2)

  cluster1_spots <- clustering_result$cluster_ids[[cluster1]]
  cluster2_spots <- clustering_result$cluster_ids[[cluster2]]

  if (is.null(cluster1_spots) || length(cluster1_spots) == 0)
    stop("Cluster ", cluster1, " not found in clustering results")
  if (is.null(cluster2_spots) || length(cluster2_spots) == 0)
    stop("Cluster ", cluster2, " not found in clustering results")

  message("\nComparing Cluster ", cluster1, " (", length(cluster1_spots), " spots) vs ",
          "Cluster ", cluster2, " (", length(cluster2_spots), " spots)\n")

  return(run_deg_analysis(seurat_obj,
                          group1_spots = cluster1_spots,
                          group2_spots = cluster2_spots,
                          comparison   = "group1_vs_group2"))
}


# ══════════════════════════════════════════════════════════════════════════════
#' Get Top DEGs
#'
#' @description
#' Extract top N differentially expressed genes from DEG results
#'
#' @param deg_results Data frame from run_deg_analysis()
#' @param n           Integer; number of top genes to return (default: 20)
#' @param by          Character; "p_val_adj", "log2FC", or "both" (default: "p_val_adj")
#' @param direction   Character; "up", "down", or "both" (default: "both")
#'
#' @return Data frame of top DEGs
#' @export
# ══════════════════════════════════════════════════════════════════════════════

get_top_degs <- function(deg_results,
                         n         = 20,
                         by        = "p_val_adj",
                         direction = "both") {

  if (!by %in% c("p_val_adj", "log2FC", "both"))
    stop("by must be 'p_val_adj', 'log2FC', or 'both'")
  if (!direction %in% c("up", "down", "both"))
    stop("direction must be 'up', 'down', or 'both'")

  if (direction == "up")   deg_results <- deg_results[deg_results$avg_log2FC > 0, ]
  if (direction == "down") deg_results <- deg_results[deg_results$avg_log2FC < 0, ]

  if (by == "p_val_adj") {
    deg_results <- deg_results[order(deg_results$p_val_adj), ]
  } else if (by == "log2FC") {
    deg_results <- deg_results[order(-abs(deg_results$avg_log2FC)), ]
  } else {
    deg_results <- deg_results[order(deg_results$p_val_adj, -abs(deg_results$avg_log2FC)), ]
  }

  return(head(deg_results, n))
}


# ══════════════════════════════════════════════════════════════════════════════
#' Visualize DEG Results
#'
#' @description
#' Create volcano plot or bar plot of DEG results
#'
#' @param deg_results    Data frame from run_deg_analysis()
#' @param plot_type      Character; "volcano" or "bar" (default: "volcano")
#' @param top_n          Integer; number of top genes for bar plot (default: 20)
#' @param fc_threshold   Numeric; log2FC threshold for colouring (default: 0.5)
#' @param pval_threshold Numeric; -log10(p_val_adj) threshold (default: 1.3, i.e. p=0.05)
#'
#' @import ggplot2
#' @export
# ══════════════════════════════════════════════════════════════════════════════

visualize_degs <- function(deg_results,
                           plot_type      = "volcano",
                           top_n          = 20,
                           fc_threshold   = 0.5,
                           pval_threshold = 1.3) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 package is required for visualization")

  if (plot_type == "volcano") {
    deg_results$neg_log10_p <- -log10(deg_results$p_val_adj)
    deg_results$significant <- ifelse(
      abs(deg_results$avg_log2FC) > fc_threshold & deg_results$neg_log10_p > pval_threshold,
      "Significant", "Not Significant"
    )
    p <- ggplot2::ggplot(deg_results,
                         ggplot2::aes(x = avg_log2FC, y = neg_log10_p, color = significant)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::scale_color_manual(
        values = c("Not Significant" = "grey", "Significant" = "red")) +
      ggplot2::geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "blue") +
      ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                          linetype = "dashed", color = "blue") +
      ggplot2::labs(x = "Log2 Fold Change", y = "-Log10(Adjusted P-value)",
                    title = "Volcano Plot of DEGs") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
                     legend.position = "top")

  } else if (plot_type == "bar") {
    top_genes       <- get_top_degs(deg_results, n = top_n, by = "p_val_adj")
    top_genes$gene  <- factor(top_genes$gene, levels = rev(top_genes$gene))
    p <- ggplot2::ggplot(top_genes,
                         ggplot2::aes(x = gene, y = avg_log2FC, fill = avg_log2FC > 0)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                                 labels = c("Down", "Up")) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Gene", y = "Log2 Fold Change",
                    title = paste("Top", top_n, "DEGs"), fill = "Direction") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
                     legend.position = "top")
  } else {
    stop("plot_type must be 'volcano' or 'bar'")
  }

  print(p)
  return(invisible(p))
}