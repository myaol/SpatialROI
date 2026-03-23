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
library(SpatialScopeDev)

data_path <- getOption("SpatialScope.data_path", default = NULL)  # ← this must match run_SpatialScope

if (!is.null(data_path) && file.exists(data_path)) {
  seurat_obj <- readRDS(data_path)
  SpatialScopeDev::run_spatial_selector(seurat_obj, basename(data_path), show_image = TRUE)
} else {
  ui <- fluidPage(
    shinyjs::useShinyjs(),
    titlePanel("Load your data"),
    fileInput("file", "Upload RDS:")
  )
  server <- function(input, output, session) {
    observeEvent(input$file, {
      seurat_obj <- readRDS(input$file$datapath)
      stopApp()
      SpatialScopeDev::run_spatial_selector(seurat_obj, input$file$name)
    })
  }
  shinyApp(ui, server)
}