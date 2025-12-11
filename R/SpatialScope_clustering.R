#' Perform Clustering on Spatial Transcriptomics Data
#'
#' @description
#' Run clustering analysis on Seurat spatial objects with customizable parameters.
#' Returns cluster assignments and generates visualization plots.
#'
#' @param seurat_obj A Seurat spatial object (can be full dataset or subset)
#' @param resolution Numeric; clustering resolution (default: 0.8). Higher values = more clusters
#' @param n_pcs Integer; number of principal components to use (default: 30)
#' @param run_umap Logical; whether to compute UMAP for visualization (default: TRUE)
#' @param plot_umap Logical; whether to display UMAP plot (default: TRUE)
#' @param plot_spatial Logical; whether to display spatial plot if coordinates available (default: TRUE)
#' @param verbose Logical; whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \item{seurat_obj}{The Seurat object with clustering results added}
#'   \item{clusters}{Named factor vector of cluster assignments}
#'   \item{cluster_ids}{List of spot IDs for each cluster}
#'   \item{n_clusters}{Number of clusters found}
#'   \item{resolution}{Resolution used}
#'   \item{n_pcs}{Number of PCs used}
#'
#' @import Seurat
#' @importFrom Seurat DefaultAssay NormalizeData FindVariableFeatures ScaleData
#' @importFrom Seurat RunPCA FindNeighbors FindClusters RunUMAP Idents DimPlot
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage - cluster entire dataset
#' result <- run_spatial_clustering(seurat_obj)
#'
#' # Cluster with custom parameters
#' result <- run_spatial_clustering(seurat_obj, resolution = 0.5, n_pcs = 20)
#'
#' # Use selected spots from draw_ROI
#' selected_spots <- draw_ROI(seurat_obj, sample_name = "Sample1")
#' subset_obj <- subset(seurat_obj, cells = selected_spots)
#' result <- run_spatial_clustering(subset_obj, resolution = 1.0)
#'
#' # Access results
#' print(result$clusters)              # View all cluster assignments
#' table(result$clusters)              # Count spots per cluster
#' cluster_0_spots <- result$cluster_ids[["0"]]  # Get spot IDs for cluster 0
#'
#' # Use clustered object for downstream analysis
#' clustered_seurat <- result$seurat_obj
#' }

run_spatial_clustering <- function(seurat_obj,
                                   resolution = 0.8,
                                   n_pcs = 30,
                                   run_umap = TRUE,
                                   plot_umap = TRUE,
                                   plot_spatial = TRUE,
                                   verbose = TRUE) {

  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (resolution <= 0 || resolution > 5) {
    stop("Resolution must be between 0 and 5")
  }

  if (n_pcs < 5 || n_pcs > 100) {
    stop("n_pcs must be between 5 and 100")
  }

  n_spots <- ncol(seurat_obj)
  if (verbose) {
    message("\n", paste(rep("=", 70), collapse = ""))
    message("🔬 Running Spatial Clustering")
    message(paste(rep("=", 70), collapse = ""))
    message("Spots to cluster: ", n_spots)
    message("Resolution: ", resolution)
    message("PCs: ", n_pcs)
    message("")
  }

  # Get default assay
  default_assay <- DefaultAssay(seurat_obj)
  if (verbose) message("Using assay: ", default_assay)

  # Check if data needs preprocessing
  needs_preprocessing <- FALSE

  if (!"SCT" %in% names(seurat_obj@assays)) {
    current_assay <- seurat_obj@assays[[default_assay]]
    scale_data <- slot(current_assay, "scale.data")
    var_features <- slot(current_assay, "var.features")

    if (is.null(scale_data) || length(scale_data) == 0 ||
        is.null(var_features) || length(var_features) == 0) {
      needs_preprocessing <- TRUE
    }
  }

  # Preprocessing if needed
  if (needs_preprocessing) {
    if (verbose) message("Preprocessing data...")

    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    if (verbose) message("  ✓ Normalized")

    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    if (verbose) message("  ✓ Found variable features: ", length(VariableFeatures(seurat_obj)))

    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    if (verbose) message("  ✓ Scaled data")
  }

  # Run PCA if not present
  if (is.null(seurat_obj@reductions$pca)) {
    if (verbose) message("Running PCA...")
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    if (verbose) message("  ✓ PCA complete")
  } else {
    if (verbose) message("Using existing PCA")
  }

  # Check if we have enough PCs
  available_pcs <- ncol(seurat_obj@reductions$pca)
  if (n_pcs > available_pcs) {
    warning("Requested ", n_pcs, " PCs but only ", available_pcs, " available. Using ", available_pcs, " PCs.")
    n_pcs <- available_pcs
  }

  # Find neighbors and clusters
  if (verbose) message("Finding neighbors...")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
  if (verbose) message("  ✓ Neighbors found")

  if (verbose) message("Finding clusters...")
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  clusters <- Idents(seurat_obj)
  n_clusters <- length(unique(clusters))
  if (verbose) message("  ✓ Found ", n_clusters, " clusters")

  # Run UMAP if requested
  if (run_umap) {
    if (is.null(seurat_obj@reductions$umap)) {
      if (verbose) message("Running UMAP...")
      seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
      if (verbose) message("  ✓ UMAP complete")
    } else {
      if (verbose) message("Using existing UMAP")
    }
  }

  # Create cluster_ids list (spot IDs for each cluster)
  cluster_ids <- list()
  for (cluster_name in levels(clusters)) {
    cluster_ids[[cluster_name]] <- names(clusters)[clusters == cluster_name]
  }

  # Print summary
  if (verbose) {
    message("")
    message("Cluster Summary:")
    cluster_table <- table(clusters)
    for (i in seq_along(cluster_table)) {
      message(sprintf("  Cluster %s: %d spots (%.1f%%)",
                      names(cluster_table)[i],
                      cluster_table[i],
                      100 * cluster_table[i] / n_spots))
    }
    message("")
  }

  # Generate plots
  if (plot_umap && !is.null(seurat_obj@reductions$umap)) {
    if (verbose) message("Generating UMAP plot...")

    p_umap <- DimPlot(seurat_obj,
                      reduction = "umap",
                      label = TRUE,
                      label.size = 5,
                      pt.size = 1.2) +
      ggtitle(paste0("UMAP Clustering (Resolution: ", resolution, ", PCs: ", n_pcs, ")")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right"
      )

    print(p_umap)
  }

  if (plot_spatial && length(seurat_obj@images) > 0) {
    if (verbose) message("Generating spatial plot...")

    # Try to create spatial plot
    tryCatch({
      p_spatial <- DimPlot(seurat_obj,
                           reduction = NULL,
                           group.by = "seurat_clusters",
                           label = TRUE,
                           label.size = 4,
                           pt.size.factor = 1.5) +
        ggtitle(paste0("Spatial Clustering (Resolution: ", resolution, ")")) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right"
        )

      print(p_spatial)
    }, error = function(e) {
      if (verbose) message("  Note: Could not generate spatial plot (", e$message, ")")
    })
  }

  # Prepare return object
  result <- list(
    seurat_obj = seurat_obj,
    clusters = clusters,
    cluster_ids = cluster_ids,
    n_clusters = n_clusters,
    resolution = resolution,
    n_pcs = n_pcs
  )

  if (verbose) {
    message(paste(rep("=", 70), collapse = ""))
    message("✅ Clustering complete!")
    message("Access results with: result$clusters, result$cluster_ids, result$seurat_obj")
    message(paste(rep("=", 70), collapse = ""))
    message("")
  }

  return(invisible(result))
}


#' Quick Clustering Wrapper
#'
#' @description
#' Simplified wrapper for quick clustering with fewer parameters.
#' Useful for rapid exploration.
#'
#' @param seurat_obj A Seurat spatial object
#' @param resolution Numeric; clustering resolution (default: 0.8)
#'
#' @return Cluster assignments as a named factor vector
#' @export
#'
#' @examples
#' \dontrun{
#' # Quick clustering
#' clusters <- cluster_spots(seurat_obj)
#' table(clusters)
#'
#' # With custom resolution
#' clusters <- cluster_spots(seurat_obj, resolution = 1.2)
#' }

cluster_spots <- function(seurat_obj, resolution = 0.8) {
  result <- run_spatial_clustering(seurat_obj,
                                   resolution = resolution,
                                   plot_umap = TRUE,
                                   plot_spatial = FALSE,
                                   verbose = TRUE)
  return(result$clusters)
}


#' Extract Spots by Cluster
#'
#' @description
#' Extract spot IDs for specific clusters from clustering results
#'
#' @param clustering_result Result object from run_spatial_clustering()
#' @param clusters Character or numeric vector of cluster IDs to extract
#'
#' @return Character vector of spot IDs
#' @export
#'
#' @examples
#' \dontrun{
#' result <- run_spatial_clustering(seurat_obj)
#'
#' # Get spots from cluster 0
#' cluster_0_spots <- get_cluster_spots(result, clusters = "0")
#'
#' # Get spots from multiple clusters
#' cluster_0_1_spots <- get_cluster_spots(result, clusters = c("0", "1"))
#'
#' # Create subset
#' subset_obj <- subset(seurat_obj, cells = cluster_0_spots)
#' }

get_cluster_spots <- function(clustering_result, clusters) {
  if (!all(c("cluster_ids", "clusters") %in% names(clustering_result))) {
    stop("Input must be a result object from run_spatial_clustering()")
  }

  clusters <- as.character(clusters)
  available_clusters <- names(clustering_result$cluster_ids)

  invalid_clusters <- setdiff(clusters, available_clusters)
  if (length(invalid_clusters) > 0) {
    stop("Invalid cluster IDs: ", paste(invalid_clusters, collapse = ", "),
         "\nAvailable clusters: ", paste(available_clusters, collapse = ", "))
  }

  # Extract spot IDs for requested clusters
  spot_ids <- unlist(clustering_result$cluster_ids[clusters], use.names = FALSE)

  return(spot_ids)
}


#' Compare Clustering Resolutions
#'
#' @description
#' Run clustering at multiple resolutions and visualize results side-by-side
#'
#' @param seurat_obj A Seurat spatial object
#' @param resolutions Numeric vector of resolutions to test (default: c(0.3, 0.5, 0.8, 1.0, 1.5))
#' @param n_pcs Integer; number of PCs to use (default: 30)
#'
#' @return A list of clustering results for each resolution
#' @export
#'
#' @examples
#' \dontrun{
#' # Test multiple resolutions
#' results <- compare_resolutions(seurat_obj, resolutions = c(0.5, 0.8, 1.2))
#'
#' # Access specific resolution result
#' res_08 <- results[["0.8"]]
#' }

compare_resolutions <- function(seurat_obj,
                                resolutions = c(0.3, 0.5, 0.8, 1.0, 1.5),
                                n_pcs = 30) {

  message("\n", paste(rep("=", 70), collapse = ""))
  message("🔍 Comparing Clustering Resolutions")
  message(paste(rep("=", 70), collapse = ""))
  message("Resolutions to test: ", paste(resolutions, collapse = ", "))
  message("")

  # Preprocess once
  default_assay <- DefaultAssay(seurat_obj)
  needs_preprocessing <- FALSE

  if (!"SCT" %in% names(seurat_obj@assays)) {
    current_assay <- seurat_obj@assays[[default_assay]]
    scale_data <- slot(current_assay, "scale.data")
    var_features <- slot(current_assay, "var.features")

    if (is.null(scale_data) || length(scale_data) == 0 ||
        is.null(var_features) || length(var_features) == 0) {
      needs_preprocessing <- TRUE
    }
  }

  if (needs_preprocessing) {
    message("Preprocessing data...")
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    message("  ✓ Preprocessing complete\n")
  }

  if (is.null(seurat_obj@reductions$pca)) {
    message("Running PCA...")
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    message("  ✓ PCA complete\n")
  }

  if (is.null(seurat_obj@reductions$umap)) {
    message("Running UMAP...")
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
    message("  ✓ UMAP complete\n")
  }

  # Run clustering at each resolution
  results <- list()
  plots <- list()

  for (res in resolutions) {
    message("Testing resolution: ", res)

    temp_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs, verbose = FALSE)
    temp_obj <- FindClusters(temp_obj, resolution = res, verbose = FALSE)
    clusters <- Idents(temp_obj)
    n_clusters <- length(unique(clusters))

    message("  → Found ", n_clusters, " clusters\n")

    # Store result
    results[[as.character(res)]] <- list(
      seurat_obj = temp_obj,
      clusters = clusters,
      n_clusters = n_clusters,
      resolution = res
    )

    # Create UMAP plot
    p <- DimPlot(temp_obj,
                 reduction = "umap",
                 label = TRUE,
                 label.size = 4,
                 pt.size = 0.8) +
      ggtitle(paste0("Resolution: ", res, " (", n_clusters, " clusters)")) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "none"
      )

    plots[[as.character(res)]] <- p
  }

  # Display all plots together
  message("Generating comparison plot...")
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- patchwork::wrap_plots(plots, ncol = min(3, length(plots)))
    print(combined_plot)
  } else {
    # If patchwork not available, print plots sequentially
    for (p in plots) {
      print(p)
    }
  }

  message("\n", paste(rep("=", 70), collapse = ""))
  message("✅ Comparison complete!")
  message("Summary:")
  for (res in names(results)) {
    message(sprintf("  Resolution %s: %d clusters", res, results[[res]]$n_clusters))
  }
  message(paste(rep("=", 70), collapse = ""), "\n")

  return(invisible(results))
}
