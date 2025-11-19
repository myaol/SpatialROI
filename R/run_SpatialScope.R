#' Launch SpatialScope Shiny App
#'
#' Wrapper for run_spatial_selector() so users can call
#' run_SpatialScope("demo") or use their own Seurat object.
#'
#' @param seurat_input Seurat object, path to .rds file, or "demo" for example data
#' @param sample_name Optional dataset name
#' @param show_image Logical; whether to display H&E image
#' @export
run_SpatialScope <- function(seurat_input = "demo", sample_name = "sample", show_image = TRUE) {
  # Use the function defined in the current environment (development)
  run_spatial_selector(seurat_input,
                       sample_name = sample_name,
                       show_image = show_image)
}
