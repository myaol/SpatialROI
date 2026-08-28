library(shiny)
library(shinyjs)
library(leaflet)
library(dplyr)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(SpatialROI)

options(shiny.maxRequestSize = 500 * 1024^2)
data_path <- getOption("SpatialROI.data_path", default = NULL)

if (!is.null(data_path) && file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  SpatialROI::run_spatial_selector(seurat_obj, basename(data_path), show_image = TRUE)
} else {
  # Launch the complete application with the packaged example. Users can replace
  # it from the session-local upload panel; no nested/second Shiny app is started.
  SpatialROI::run_spatial_selector("demo")
}
