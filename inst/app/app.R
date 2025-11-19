# ---- Load dependencies ----
library(shiny)
library(shinyjs)
library(leaflet)
library(leaflet.extras)
library(dplyr)
library(sf)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(SpatialScope)

data_path <- getOption("SpatialScope.data_path", default = NULL)

# If data path provided, load and run directly
if (!is.null(data_path) && file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  SpatialScope::run_spatial_selector(seurat_obj, basename(data_path), show_image = TRUE)
} else {
  # Otherwise just run with default (you'd need example data in the package)
  # Or show a simple loader UI
  ui <- fluidPage(
    shinyjs::useShinyjs(),
    titlePanel("Load your data"),
    fileInput("file", "Upload RDS:")
  )
  server <- function(input, output, session) {
    observeEvent(input$file, {
      seurat_obj <- readRDS(input$file$datapath)
      stopApp()
      SpatialScope::run_spatial_selector(seurat_obj, input$file$name)
    })
  }
  shinyApp(ui, server)
}
