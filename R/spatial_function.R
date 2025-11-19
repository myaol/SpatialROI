#' Run the SpatialScope Spatial Selector
#'
#' @description
#' Launch an interactive ROI selection and analysis tool for spatial transcriptomics data.
#' Accepts a Seurat object, a data list, or `"demo"` to load the included example dataset.
#'
#' @param seurat_input Seurat object, data list, or `"demo"` to load example data.
#' @param sample_name Character label for the dataset.
#' @param show_image Logical; whether to show the H&E image (if available).
#'
#' @details
#' This function powers the interactive module for selecting regions of interest (ROIs),
#' visualizing gene expression, running clustering, and performing differential expression
#' directly within a spatial transcriptomics context.
#'
#' @import ggplot2
#'
#' @importFrom Seurat
#'   AddModuleScore
#'   DefaultAssay
#'   DimPlot
#'   FetchData
#'   FindAllMarkers
#'   FindMarkers
#'   FindVariableFeatures
#'   GetAssayData
#'   GetTissueCoordinates
#'   Idents
#'   "Idents<-"
#'   RunPCA
#'   RunUMAP
#'   FindNeighbors
#'   FindClusters
#'   NormalizeData
#'   ScaleData
#'
#' @importFrom dplyr %>% filter
#' @importFrom grDevices as.raster colorRampPalette dev.off png rainbow
#' @importFrom graphics par
#' @importFrom methods slot
#' @importFrom stats cor.test setNames t.test wilcox.test
#' @importFrom utils capture.output head write.csv
#'
#' @return A Shiny application object; launching the GUI as a side effect.
#' @export

run_spatial_selector <- function(seurat_input, sample_name = "sample", show_image = TRUE) {

  # Handle demo data
  if (is.character(seurat_input) && seurat_input == "demo") {
    demo_file <- system.file("extdata", "example_visium.rds", package = "SpatialScope")  # ← Changed here

    if (demo_file == "" || !file.exists(demo_file)) {
      stop("Demo data not found. Please ensure the package is installed correctly.\n",
           "The example data should be located in inst/extdata/example_visium.rds")
    }

    message("Loading example Visium dataset from inst/extdata...")
    seurat_input <- list(seurat = readRDS(demo_file))
    sample_name <- "Example_Visium"
  }

  # Prepare data (this replaces your first ~60 lines)
  prepared_data <- prepare_seurat_data(seurat_input, sample_name, show_image)


  # Extract prepared components
  seurat_obj <- prepared_data$seurat_obj
  coords <- prepared_data$coords
  image_name <- prepared_data$image_name
  spots_sf <- prepared_data$spots_sf
  he_image_base64 <- prepared_data$he_image_base64
  he_image_bounds <- prepared_data$he_image_bounds
  signature_library_human <- prepared_data$signature_library_human
  signature_library_mouse <- prepared_data$signature_library_mouse


  all_genes <- rownames(seurat_obj)

  all_metadata <- colnames(seurat_obj@meta.data)

  x_range <- range(spots_sf$x)
  y_range <- range(spots_sf$y)
  x_buffer <- diff(x_range) * 0.1
  y_buffer <- diff(y_range) * 0.1

  # H&E image processing
  he_image_base64 <- NULL
  he_image_bounds <- NULL
  image_obj <- seurat_obj@images[[image_name]]

  if (show_image) {
    tryCatch({
      if (class(image_obj)[1] == "VisiumV1") {
        he_image_data <- image_obj@image
        coords_full <- GetTissueCoordinates(seurat_obj, image = image_name)
        if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
          pixel_x <- coords_full$pxl_col_in_fullres
          pixel_y <- coords_full$pxl_row_in_fullres
        } else {
          pixel_x <- spots_sf$x
          pixel_y <- spots_sf$y
        }
        pixel_x_min <- min(pixel_x, na.rm = TRUE)
        pixel_x_max <- max(pixel_x, na.rm = TRUE)
        pixel_y_min <- min(pixel_y, na.rm = TRUE)
        pixel_y_max <- max(pixel_y, na.rm = TRUE)
        x_buffer_px <- (pixel_x_max - pixel_x_min) * 0.1
        y_buffer_px <- (pixel_y_max - pixel_y_min) * 0.1
        crop_x_min <- max(1, floor(pixel_x_min - x_buffer_px))
        crop_x_max <- min(dim(he_image_data)[2], ceiling(pixel_x_max + x_buffer_px))
        crop_y_min <- max(1, floor(pixel_y_min - y_buffer_px))
        crop_y_max <- min(dim(he_image_data)[1], ceiling(pixel_y_max + y_buffer_px))
        he_image_cropped <- he_image_data[crop_y_min:crop_y_max, crop_x_min:crop_x_max, ]
        temp_file <- tempfile(fileext = ".png")
        png(temp_file, width = dim(he_image_cropped)[2], height = dim(he_image_cropped)[1])
        par(mar = c(0,0,0,0))
        plot(as.raster(he_image_cropped), axes = FALSE)
        dev.off()
        he_image_base64 <- paste0("data:image/png;base64,", base64enc::base64encode(temp_file))
        unlink(temp_file)
        he_image_bounds <- list(
          south = crop_y_max,
          west = crop_x_min,
          north = crop_y_min,
          east = crop_x_max
        )
      }
    }, error = function(e) {
      print(paste("Could not extract H&E image:", e$message))
    })
  }



  # UI
  ui <- fluidPage(
    useShinyjs(),  # Enable shinyjs
    tags$head(
      tags$style(HTML("
        body { margin: 0; padding: 0; overflow: hidden; }
        .main-container { display: flex; height: 100vh; }

        /* Initial loading screen - shown before everything */
        .initial-loading {
          position: fixed;
          top: 0;
          left: 0;
          right: 0;
          bottom: 0;
          background: linear-gradient(135deg, #0072B5 0%, #E18727 100%);
          z-index: 9999;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          color: white;
        }
        .initial-loading.loaded {
          opacity: 0;
          visibility: hidden;
          transition: opacity 0.5s, visibility 0.5s;
        }
        .loading-spinner-large {
          border: 12px solid rgba(255,255,255,0.3);
          border-radius: 50%;
          border-top: 12px solid white;
          width: 100px;
          height: 100px;
          animation: spin 1s linear infinite;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        .loading-title {
          font-size: 48px;
          font-weight: bold;
          margin-top: 30px;
          margin-bottom: 10px;
        }
        .loading-message {
          font-size: 24px;
          margin-top: 10px;
        }

        /* Top header bar */
        .top-header {
          position: fixed;
          top: 0;
          left: 0;
          right: 0;
          height: 80px;
          background: linear-gradient(135deg, #0072B5 0%, #E18727 100%);
          color: white;
          display: flex;
          align-items: center;
          padding: 0 20px;
          z-index: 2000;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          visibility: hidden;
          opacity: 0;
          transition: opacity 0.3s, visibility 0.3s;
        }
        .top-header.active {
          visibility: visible;
          opacity: 1;
        }
        .top-header h1 {
          margin: 0;
          font-size: 28px;
          font-weight: bold;
          text-shadow: 1px 1px 2px rgba(0,0,0,0.2);
        }

        /* Sidebar styling */
        .sidebar-left {
          width: 120px;
          background: white;
          display: flex;
          flex-direction: column;
          padding: 15px 0;
          margin-top: 80px;
          height: calc(100vh - 80px);
          box-shadow: 2px 0 10px rgba(0,0,0,0.1);
        }
        .sidebar-button {
          width: 90px;
          height: 70px;
          margin: 12px auto;
          background: transparent;
          border: none;
          color: #2c3e50;
          cursor: pointer;
          border-radius: 10px;
          transition: all 0.3s;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          font-size: 28px;
        }
        .sidebar-button:hover {
          background: #f0f0f0;
        }
        .sidebar-button.active {
          background: #0072B5;
          color: white;
        }
        .sidebar-button-label {
          font-size: 12px;
          margin-top: 6px;
          font-weight: 500;
        }

        /* Map container */
        .map-container {
          flex: 1;
          position: relative;
          background: #f5f5f5;
          margin-top: 80px;
          height: calc(100vh - 80px);
        }

        /* Control panel */
        .control-panel {
          width: 0;
          background: white;
          box-shadow: -2px 0 10px rgba(0,0,0,0.1);
          overflow-y: auto;
          overflow-x: hidden;
          transition: width 0.3s;
          position: relative;
          margin-top: 80px;
          height: calc(100vh - 80px);
        }
        .control-panel.open {
          width: 400px;
        }
        .control-content {
          padding: 20px;
          display: none;
        }
        .control-content.active {
          display: block;
        }

        /* Responsive Design */
        @media (max-width: 1200px) {
          .control-panel.open {
            width: 350px;
          }
        }

        @media (max-width: 992px) {
          .sidebar-left {
            width: 80px;
          }
          .sidebar-button {
            width: 60px;
            height: 60px;
            font-size: 24px;
          }
          .sidebar-button-label {
            font-size: 10px;
          }
          .control-panel.open {
            width: 300px;
          }
          .top-header h1 {
            font-size: 20px;
          }
        }

        @media (max-width: 768px) {
          .main-container {
            flex-direction: column;
          }
          .sidebar-left {
            width: 100%;
            height: 60px;
            flex-direction: row;
            padding: 0;
            margin-top: 80px;
            overflow-x: auto;
            overflow-y: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
          }
          .sidebar-button {
            width: 70px;
            height: 50px;
            margin: 5px;
            font-size: 20px;
          }
          .sidebar-button-label {
            font-size: 9px;
          }
          .map-container {
            margin-top: 140px;
            height: calc(100vh - 140px);
          }
          .control-panel {
            position: fixed;
            bottom: 0;
            left: 0;
            right: 0;
            width: 100% !important;
            height: 0;
            margin-top: 0;
            z-index: 3000;
            transition: height 0.3s;
          }
          .control-panel.open {
            height: 60vh;
            width: 100% !important;
          }
          .top-header h1 {
            font-size: 18px;
          }
        }

        @media (max-width: 576px) {
          .top-header {
            height: 50px;
            padding: 0 10px;
          }
          .top-header h1 {
            font-size: 16px;
          }
          .sidebar-left {
            margin-top: 50px;
            height: 55px;
          }
          .sidebar-button {
            width: 60px;
            height: 45px;
            font-size: 18px;
          }
          .sidebar-button-label {
            font-size: 8px;
          }
          .map-container {
            margin-top: 105px;
            height: calc(100vh - 105px);
          }
          .control-panel.open {
            height: 70vh;
          }
        }

        /* Full screen content for Home */
        .control-content.full-screen-content.active {
          position: fixed;
          top: 80px;
          left: 135px;
          right: 0;
          bottom: 0;
          z-index: 3000;
          background: #f5f7fa;
          padding: 40px 80px;
          display: flex;
          justify-content: center;
          overflow-y: auto;
        }
        .panel-header {
          font-size: 20px;
          font-weight: bold;
          margin-bottom: 20px;
          color: #2c3e50;
          border-bottom: 2px solid #0072B5;
          padding-bottom: 10px;
        }
        .control-section {
          margin-bottom: 20px;
          padding: 15px;
          background: #f8f9fa;
          border-radius: 8px;
        }
        .control-section h4 {
          margin-top: 0;
          color: #34495e;
        }
        #map { height: 100% !important; }
        .leaflet-container { background: #f5f5f5; }

        /* Map controls */
        .map-controls {
          position: absolute;
          top: 20px;
          right: 20px;
          z-index: 1000;
          background: white;
          padding: 15px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }

        /* Clear selection button */
        .clear-selection-btn {
          position: absolute;
          bottom: 20px;
          right: 20px;
          z-index: 1000;
          background: #e74c3c;
          color: white;
          border: none;
          padding: 12px 24px;
          border-radius: 8px;
          cursor: pointer;
          font-weight: bold;
          font-size: 14px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          transition: all 0.3s;
        }
        .clear-selection-btn:hover {
          background: #c0392b;
          transform: translateY(-2px);
          box-shadow: 0 4px 15px rgba(0,0,0,0.3);
        }

        /* Spot count */
        .spot-count {
          position: absolute;
          top: 20px;
          left: 20px;
          z-index: 1000;
          background: white;
          padding: 10px 15px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          font-weight: bold;
        }
        .close-panel-btn {
          position: absolute;
          top: 15px;
          right: 15px;
          background: transparent;
          border: none;
          font-size: 24px;
          cursor: pointer;
          color: #7f8c8d;
        }
        .close-panel-btn:hover {
          color: #e74c3c;
        }

        /* Custom styling for show_groups checkbox */
        .checkbox#show_groups {
          background: white;
          padding: 10px 15px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          margin: 0 !important;
          display: inline-block;
        }
        .checkbox#show_groups label {
          margin: 0 !important;
          display: flex;
          align-items: center;
          gap: 8px;
        }
        .checkbox#show_groups input {
          margin: 0 !important;
          position: relative;
          top: 0;
        }

        /* Landing page styles */
        .landing-container {
          position: fixed;
          top: 0;
          left: 0;
          right: 0;
          bottom: 0;
          background: linear-gradient(135deg, #0072B5 0%, #E18727 100%);
          z-index: 5000;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          color: white;
          transition: opacity 0.5s, visibility 0.5s;
        }
        .landing-container.hidden {
          opacity: 0;
          visibility: hidden;
          pointer-events: none;
        }

        /* Hide entire app content until started */
        .app-content {
          visibility: hidden;
          opacity: 0;
          transition: opacity 0.3s, visibility 0.3s;
        }
        .app-content.active {
          visibility: visible;
          opacity: 1;
        }

        /* Loading overlay */
        .loading-overlay {
          position: absolute;
          top: 0;
          left: 0;
          right: 0;
          bottom: 0;
          background: rgba(0, 114, 181, 0.95);
          z-index: 3500;
          display: none;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          color: white;
        }
        .loading-overlay.active {
          display: flex;
        }
        .loading-spinner {
          border: 8px solid rgba(255,255,255,0.3);
          border-radius: 50%;
          border-top: 8px solid white;
          width: 80px;
          height: 80px;
          animation: spin 1s linear infinite;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        .loading-text {
          margin-top: 30px;
          font-size: 24px;
          font-weight: bold;
        }
        .landing-title {
          font-size: 64px;
          font-weight: bold;
          margin-bottom: 20px;
          text-shadow: 2px 2px 4px rgba(0,0,0,0.3);
        }
        .landing-subtitle {
          font-size: 24px;
          font-weight: 300;
          margin-bottom: 40px;
        }
        .landing-description {
          font-size: 18px;
          max-width: 600px;
          text-align: center;
          line-height: 1.6;
          margin-bottom: 40px;
        }
        .feature-grid {
          display: grid;
          grid-template-columns: repeat(4, 1fr);
          gap: 30px;
          max-width: 1200px;
          margin: 40px 0;
        }
        .feature-card {
          background: rgba(255,255,255,0.1);
          backdrop-filter: blur(10px);
          padding: 20px;
          border-radius: 12px;
          text-align: center;
          transition: transform 0.3s;
        }
        .feature-card:hover {
          transform: translateY(-5px);
          background: rgba(255,255,255,0.15);
        }
        .feature-icon {
          font-size: 48px;
          margin-bottom: 10px;
        }
        .feature-title {
          font-size: 16px;
          font-weight: bold;
          margin-bottom: 5px;
        }
        .feature-desc {
          font-size: 12px;
          opacity: 0.9;
        }
        .start-button {
          background: white;
          color: #0072B5;
          border: none;
          padding: 15px 40px;
          font-size: 20px;
          font-weight: bold;
          border-radius: 50px;
          cursor: pointer;
          box-shadow: 0 4px 15px rgba(0,0,0,0.3);
          transition: all 0.3s;
        }
        .start-button:hover {
          transform: scale(1.05);
          box-shadow: 0 6px 20px rgba(0,0,0,0.4);
        }
      "))
    ),

    # Initial Loading Screen (before landing page)
    tags$div(class = "initial-loading", id = "initial_loading",
             div(class = "loading-spinner-large"),
             div(class = "loading-title", "🔬 SpatialScope"),
             div(class = "loading-message", "Loading spatial data..."),
             div(style = "margin-top: 20px; font-size: 16px; opacity: 0.8;",
                 "Please wait while we prepare your analysis")
    ),

    # Landing Page
    tags$div(class = "landing-container", id = "landing_page",
             div(class = "landing-title", "🔬 SpatialScope"),
             div(class = "landing-subtitle", "Interactive Spatial Transcriptomics Analysis Platform"),
             div(class = "landing-description",
                 "Explore, analyze, and visualize spatial gene expression data with powerful interactive tools"
             ),
             div(class = "feature-grid",
                 div(class = "feature-card",
                     div(class = "feature-icon", "🗺️"),
                     div(class = "feature-title", "Interactive Map"),
                     div(class = "feature-desc", "Select regions with drawing tools")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "🎨"),
                     div(class = "feature-title", "Visualization"),
                     div(class = "feature-desc", "Display gene expression patterns")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "🧬"),
                     div(class = "feature-title", "Gene Sets"),
                     div(class = "feature-desc", "Calculate multi-gene signatures")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "📊"),
                     div(class = "feature-title", "Clustering"),
                     div(class = "feature-desc", "Identify spatial domains")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "📈"),
                     div(class = "feature-title", "DEG Analysis"),
                     div(class = "feature-desc", "Find marker genes")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "⚖️"),
                     div(class = "feature-title", "Comparisons"),
                     div(class = "feature-desc", "Compare features & groups")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "📚"),
                     div(class = "feature-title", "Statistics"),
                     div(class = "feature-desc", "Perform statistical tests")
                 ),
                 div(class = "feature-card",
                     div(class = "feature-icon", "💾"),
                     div(class = "feature-title", "Export"),
                     div(class = "feature-desc", "Download results & subsets")
                 )
             ),
             tags$button(class = "start-button",
                         onclick = "$('#landing_page').addClass('hidden'); $('#loading_overlay').addClass('active'); Shiny.setInputValue('start_analysis_from_landing', Math.random());",
                         "Start Analysis →")
    ),

    # Loading Overlay
    tags$div(class = "loading-overlay", id = "loading_overlay",
             div(class = "loading-spinner"),
             div(class = "loading-text", "Please wait a moment..."),
             div(id = "loading_message", style = "margin-top: 10px; font-size: 16px;", "Loading spatial data")
    ),

    # Top Header Bar
    tags$div(class = "top-header",
             h1(textOutput("header_title", inline = TRUE))
    ),

    tags$div(class = "main-container app-content",
             # Left sidebar with icon buttons
             tags$div(class = "sidebar-left",
                      actionButton("btn_home", HTML("<div style='font-size:24px;'>🏠</div><div class='sidebar-button-label'>Home</div>"),
                                   class = "sidebar-button active"),
                      actionButton("btn_viz", HTML("<div style='font-size:24px;'>🎨</div><div class='sidebar-button-label'>Visualization</div>"),
                                   class = "sidebar-button"),
                      actionButton("btn_geneset", HTML("<div style='font-size:24px;'>🧬</div><div class='sidebar-button-label'>Gene Sets</div>"),
                                   class = "sidebar-button"),
                      actionButton("btn_cluster", HTML("<div style='font-size:24px;'>📊</div><div class='sidebar-button-label'>Clusters</div>"),
                                   class = "sidebar-button"),
                      actionButton("btn_deg", HTML("<div style='font-size:24px;'>📈</div><div class='sidebar-button-label'>DEGs</div>"),
                                   class = "sidebar-button"),
                      actionButton("btn_compare", HTML("<div style='font-size:24px;'>⚖️</div><div class='sidebar-button-label'>Compare</div>"),
                                   class = "sidebar-button")
             ),

             # Map container
             tags$div(class = "map-container",
                      leafletOutput("map", height = "100%"),

                      # Clear selection button (bottom right)
                      tags$button(class = "clear-selection-btn",
                                  onclick = "Shiny.setInputValue('clear_selection_click', Math.random());",
                                  "🗑️ Clear Selection"),

                      # Map controls (top right)
                      tags$div(class = "map-controls",
                               tags$details(
                                 open = "open",  # Starts expanded
                                 tags$summary(
                                   tags$span(style = "font-weight: bold; font-size: 14px;", "▼ Spot Size")
                                 ),
                                 sliderInput("spot_size", NULL, min = 0, max = 15, value = 4, step = 0.5, width = "200px")
                               ),
                               tags$details(
                                 open = "open",  # Starts expanded
                                 tags$summary(
                                   tags$span(style = "font-weight: bold; font-size: 14px;", "▼ H&E Opacity")
                                 ),
                                 sliderInput("image_opacity", NULL, min = 0, max = 1, value = 0.6, step = 0.05, width = "200px")
                               )
                      ),

                      # Spot count (top middle)
                      tags$div(
                        style = "position: absolute; top: 20px; left: 50%; transform: translateX(-50%); z-index: 1000; background: white; padding: 12px 18px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.2); font-weight: bold; font-size: 15px;",
                        textOutput("spot_count_display", inline = TRUE)
                      ),

                      # Group buttons (bottom center) - MODIFIED: Added download buttons
                      tags$div(
                        style = "position: absolute; bottom: 20px; left: 50%; transform: translateX(-50%); z-index: 1000; display: flex; gap: 15px; align-items: center;",
                        tags$div(style = "display: flex; flex-direction: column; gap: 5px;",
                                 tags$button(
                                   style = "background: white; color: #E41A1C; border: 3px solid #E41A1C; padding: 12px 24px; border-radius: 8px; cursor: pointer; font-weight: bold; font-size: 14px; box-shadow: 0 2px 10px rgba(0,0,0,0.2);",
                                   onclick = "Shiny.setInputValue('save_group1_map', Math.random());",
                                   onmouseover = "this.style.background='#E41A1C'; this.style.color='white';",
                                   onmouseout = "this.style.background='white'; this.style.color='#E41A1C';",
                                   "Group 1"
                                 ),
                                 downloadButton("dl_group1", "⬇ Download Group 1", class = "btn btn-success btn-sm",
                                                style = "font-size: 11px; padding: 4px 12px;")
                        ),
                        tags$div(style = "display: flex; flex-direction: column; gap: 5px;",
                                 tags$button(
                                   style = "background: white; color: #377EB8; border: 3px solid #377EB8; padding: 12px 24px; border-radius: 8px; cursor: pointer; font-weight: bold; font-size: 14px; box-shadow: 0 2px 10px rgba(0,0,0,0.2);",
                                   onclick = "Shiny.setInputValue('save_group2_map', Math.random());",
                                   onmouseover = "this.style.background='#377EB8'; this.style.color='white';",
                                   onmouseout = "this.style.background='white'; this.style.color='#377EB8';",
                                   "Group 2"
                                 ),
                                 downloadButton("dl_group2", "⬇ Download Group 2", class = "btn btn-success btn-sm",
                                                style = "font-size: 11px; padding: 4px 12px;")
                        ),
                        checkboxInput("show_groups", "Show Groups on Map", value = FALSE)
                      )
             ),

             # Right control panel
             tags$div(class = "control-panel", id = "control_panel",
                      tags$button(class = "close-panel-btn", onclick = "Shiny.setInputValue('close_panel', Math.random());", "×"),

                      # Home content - FULL SCREEN WITHIN APP
                      tags$div(class = "control-content full-screen-content", id = "content_home",
                               tags$div(style = "width: 100%; max-width: 1300px; height: 100vh; overflow-y: auto; background: #f5f7fa; padding: 0px 20px 40px 20px;",
                                        tags$div(style = "max-width: 1100px; margin: 0 auto;",
                                                 tags$div(style = "background: linear-gradient(135deg, #0072B5 0%, #E18727 100%); color: white; padding: 30px 40px; border-radius: 12px; margin-bottom: 30px;",
                                                          tags$h1(style = "margin: 0; font-size: 32px;", "🔬 SpatialScope - User Guide & Tutorial")
                                                 ),

                                                 tags$div(style = "background: white; border-radius: 12px; padding: 35px; box-shadow: 0 2px 10px rgba(0,0,0,0.1);",
                                                          tags$div(style = "margin-bottom: 35px;",
                                                                   tags$h2(style = "color: #0072B5; font-size: 26px; margin-bottom: 15px; border-bottom: 3px solid #E18727; padding-bottom: 10px;", "Welcome to SpatialScope"),
                                                                   tags$p(style = "font-size: 15px; line-height: 1.7; color: #333;",
                                                                          "SpatialScope is an interactive platform for analyzing and visualizing spatial transcriptomics data. This tool allows you to explore gene expression patterns, identify spatial domains, and perform comprehensive statistical analyses on your Seurat spatial objects.")
                                                          ),

                                                          tags$div(style = "margin-bottom: 35px;",
                                                                   tags$h2(style = "color: #0072B5; font-size: 26px; margin-bottom: 20px; border-bottom: 3px solid #E18727; padding-bottom: 10px;", "📖 Quick Start Tutorial"),
                                                                   tags$div(style = "display: flex; margin-bottom: 25px; align-items: flex-start;",
                                                                            tags$div(style = "background: #0072B5; color: white; width: 45px; height: 45px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 22px; font-weight: bold; flex-shrink: 0; margin-right: 18px;", "1"),
                                                                            tags$div(
                                                                              tags$h3(style = "margin: 0 0 8px 0; color: #2c3e50; font-size: 18px;", "Upload Your Data"),
                                                                              tags$p(style = "margin: 0; font-size: 14px; color: #555; line-height: 1.6;", "Navigate to the 📤 Upload panel to load your own Seurat spatial object (.rds file), or use the example dataset to get started immediately.")
                                                                            )
                                                                   ),
                                                                   tags$div(style = "display: flex; margin-bottom: 25px; align-items: flex-start;",
                                                                            tags$div(style = "background: #0072B5; color: white; width: 45px; height: 45px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 22px; font-weight: bold; flex-shrink: 0; margin-right: 18px;", "2"),
                                                                            tags$div(
                                                                              tags$h3(style = "margin: 0 0 8px 0; color: #2c3e50; font-size: 18px;", "Select Regions of Interest"),
                                                                              tags$p(style = "margin: 0; font-size: 14px; color: #555; line-height: 1.6;", "Use the freehand drawing tool (✏️ pencil icon) on the map to draw around regions you want to analyze. Click and drag to create custom selection shapes.")
                                                                            )
                                                                   ),
                                                                   tags$div(style = "display: flex; margin-bottom: 25px; align-items: flex-start;",
                                                                            tags$div(style = "background: #0072B5; color: white; width: 45px; height: 45px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 22px; font-weight: bold; flex-shrink: 0; margin-right: 18px;", "3"),
                                                                            tags$div(
                                                                              tags$h3(style = "margin: 0 0 8px 0; color: #2c3e50; font-size: 18px;", "Save Groups for Comparison"),
                                                                              tags$p(style = "margin: 0; font-size: 14px; color: #555; line-height: 1.6;", "After selecting spots, click 'Group 1' or 'Group 2' buttons at the bottom of the map to save your selections for comparative analysis.")
                                                                            )
                                                                   ),
                                                                   tags$div(style = "display: flex; margin-bottom: 25px; align-items: flex-start;",
                                                                            tags$div(style = "background: #0072B5; color: white; width: 45px; height: 45px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 22px; font-weight: bold; flex-shrink: 0; margin-right: 18px;", "4"),
                                                                            tags$div(
                                                                              tags$h3(style = "margin: 0 0 8px 0; color: #2c3e50; font-size: 18px;", "Explore & Analyze"),
                                                                              tags$p(style = "margin: 0; font-size: 14px; color: #555; line-height: 1.6;", "Use the sidebar tools to visualize gene expression, calculate gene set scores, perform clustering, find differentially expressed genes, and export your results.")
                                                                            )
                                                                   )
                                                          ),

                                                          tags$div(style = "margin-bottom: 35px;",
                                                                   tags$h2(style = "color: #0072B5; font-size: 26px; margin-bottom: 20px; border-bottom: 3px solid #E18727; padding-bottom: 10px;", "🛠️ Available Tools"),
                                                                   tags$div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 15px;",
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "🏠 Home"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Access this documentation and tutorial anytime")
                                                                            ),
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "📤 Upload & Visualize"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Upload Seurat objects and visualize gene expression, metadata, or gene set scores")
                                                                            ),
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "🧬 Gene Sets"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Calculate multi-gene signatures using pre-defined libraries (human/mouse) or custom lists")
                                                                            ),
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "📊 Clustering"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Identify spatial domains using graph-based clustering")
                                                                            ),
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "📈 DEG Analysis"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Find differentially expressed genes between saved groups")
                                                                            ),
                                                                            tags$div(style = "background: #f8f9fa; padding: 15px; border-radius: 6px; border-left: 4px solid #0072B5;",
                                                                                     tags$h3(style = "margin: 0 0 8px 0; color: #0072B5; font-size: 16px;", "⚖️ Compare"),
                                                                                     tags$p(style = "margin: 0; font-size: 13px; color: #555; line-height: 1.5;", "Generate plots to compare features with statistics")
                                                                            )
                                                                   )
                                                          ),

                                                          tags$div(style = "margin-bottom: 35px;",
                                                                   tags$h2(style = "color: #0072B5; font-size: 26px; margin-bottom: 20px; border-bottom: 3px solid #E18727; padding-bottom: 10px;", "💡 Tips & Best Practices"),
                                                                   tags$ul(style = "font-size: 14px; line-height: 1.8; color: #333; padding-left: 25px; margin: 0;",
                                                                           tags$li(tags$strong("Selection Strategy:"), " Use the freehand tool (✏️) for precise region selection. You can draw multiple regions - they will all be combined into your selection."),
                                                                           tags$li(tags$strong("Group Management:"), " The 'Show Groups on Map' checkbox is enabled by default. Save selections to groups before clearing to preserve them."),
                                                                           tags$li(tags$strong("Gene Sets:"), " Choose the correct species (Human/Mouse) before loading pre-defined signatures to ensure proper gene symbol matching."),
                                                                           tags$li(tags$strong("Clustering:"), " Start with default resolution (0.8) and adjust based on your tissue complexity. Higher resolution = more clusters."),
                                                                           tags$li(tags$strong("Export:"), " Download your spot IDs and Seurat subsets for further analysis in R or other tools."),
                                                                           tags$li(tags$strong("Clear Everything:"), " Use the 🗑️ Clear Selection button to reset all selections and saved groups.")
                                                                   )
                                                          ),

                                                          tags$div(
                                                            style = "background: linear-gradient(135deg, #0072B5 0%, #E18727 100%);
           color: white; padding: 25px; border-radius: 8px; text-align: center;",

                                                            tags$h2(style = "margin: 0 0 12px 0; font-size: 22px;", "🚀 Ready to Begin?"),

                                                            tags$p(style = "margin: 0; font-size: 15px; line-height: 1.6;",
                                                                   "Click on the other tools in the sidebar (📤 Upload, 🧬 Gene Sets, etc.) to start your spatial analysis! The map will be visible when you switch to other tools."),

                                                            tags$hr(style = "border: 0; height: 1px; background-color: rgba(255,255,255,0.5); margin: 15px 0;"),

                                                            tags$p(style = "margin: 0; font-size: 14px; line-height: 1.6;",
                                                                   HTML("Contact: Aodong Qiu (<a href='mailto:qiuaodon@pitt.edu' style='color: #fff; text-decoration: underline;'>qiuaodon@pitt.edu</a>),
               Mengyao Lu (<a href='mailto:my.lu@pitt.edu' style='color: #fff; text-decoration: underline;'>mel373@pitt.edu</a>),
               Lujia Chen (<a href='mailto:luc17@pitt.edu' style='color: #fff; text-decoration: underline;'>luc17@pitt.edu</a>)"))
                                                          )
                                                 )
                                        )
                               )
                      ),

                      # Visualization content
                      tags$div(class = "control-content", id = "content_viz",
                               div(class = "panel-header", "🎨 Upload & Visualize"),
                               div(class = "control-section",
                                   h4("Data Source"),
                                   div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                                       actionButton("use_example_data", "📊 Use Example Data",
                                                    class = "btn btn-primary",
                                                    style = "flex: 1;"),
                                       actionButton("show_upload_panel", "📤 Upload Data",
                                                    class = "btn btn-info",
                                                    style = "flex: 1;")
                                   ),
                                   conditionalPanel(
                                     condition = "input.show_upload_panel % 2 == 1",
                                     fileInput("upload_seurat", "Select Seurat Object (.rds)",
                                               accept = c(".rds")),
                                     actionButton("load_uploaded_seurat", "Load Uploaded Data",
                                                  class = "btn btn-success btn-block"),
                                     tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                            "⚠️ Loading new data will replace current analysis")
                                   ),
                                   verbatimTextOutput("upload_status")
                               ),
                               div(class = "control-section",
                                   h4("Feature Selection"),
                                   selectInput("feature_type", "Feature Type:",
                                               choices = c("None", "Gene Expression", "Metadata")),
                                   conditionalPanel(
                                     condition = "input.feature_type == 'Gene Expression'",
                                     selectInput("gene_select", "Gene:", choices = c("", all_genes))
                                   ),
                                   conditionalPanel(
                                     condition = "input.feature_type == 'Metadata'",
                                     selectInput("meta_select", "Metadata:", choices = c("", all_metadata))
                                   )
                               ),
                               div(class = "control-section",
                                   h4("Color Scheme"),
                                   conditionalPanel(
                                     condition = "input.feature_type != 'None'",
                                     selectInput("color_scheme", "Palette:",
                                                 choices = c("Grey to Red" = "greyred",
                                                             "Blue-White-Red" = "bwr",
                                                             "Rainbow" = "rainbow")),
                                     plotOutput("color_legend", height = "150px")
                                   )
                               )
                      ),

                      # Gene Set content - MODIFIED: Added species selection
                      tags$div(class = "control-content", id = "content_geneset",
                               div(class = "panel-header", "🧬 Gene Set Analysis"),

                               # CellMarker 2.0 Citation Panel
                               div(class = "control-section",
                                   style = "background-color: #e8f4f8; border-left: 4px solid #0072B5; padding: 10px; margin-bottom: 15px;",
                                   tags$div(
                                     tags$p(style = "margin: 0; font-size: 13px; font-weight: bold; color: #0072B5;",
                                            "📚 Cell Marker Database"),
                                     tags$p(style = "margin: 5px 0; font-size: 12px; line-height: 1.5;",
                                            "Pre-defined signatures are curated from CellMarker 2.0, a manually curated database of ",
                                            tags$b("26,915 cell markers"), " across ", tags$b("2,578 cell types"), " and ", tags$b("656 tissues.")),
                                     tags$p(style = "margin: 5px 0 0 0; font-size: 11px; line-height: 1.4;",
                                            tags$b("Citation:"), " Hu C, Li T, Xu Y, et al.",
                                            tags$i("Nucleic Acids Res."), " 2023;51(D1):D870-D876."),
                                     tags$div(style = "margin-top: 8px;",
                                              tags$a(href = "https://academic.oup.com/nar/article/51/D1/D870/6775381",
                                                     target = "_blank",
                                                     style = "font-size: 11px; color: #0072B5; text-decoration: none; margin-right: 10px;",
                                                     "📄 Read Paper"),
                                              tags$a(href = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/",
                                                     target = "_blank",
                                                     style = "font-size: 11px; color: #0072B5; text-decoration: none; margin-right: 10px;",
                                                     "🌐 Visit Database"),
                                              tags$a(href = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html",
                                                     target = "_blank",
                                                     style = "font-size: 11px; color: #0072B5; text-decoration: none;",
                                                     "⬇️ Download Data")
                                     )
                                   )
                               ),

                               div(class = "control-section",
                                   h4("Species Selection"),
                                   selectInput("species_select", "Select Species:",
                                               choices = c("Human" = "human", "Mouse" = "mouse"),
                                               selected = "human"),
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                          "💡 Gene symbols will be updated based on species")
                               ),
                               div(class = "control-section",
                                   h4("Select Signature"),
                                   selectInput("signature_library", "Pre-defined Signatures:",
                                               choices = names(signature_library_human),
                                               selected = "Custom"),
                                   conditionalPanel(
                                     condition = "input.signature_library != 'Custom'",
                                     actionButton("load_signature", "Load Selected Signature",
                                                  class = "btn btn-info btn-block")
                                   ),
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 10px;",
                                          "💡 Select a pre-defined signature or enter custom genes below")
                               ),
                               div(class = "control-section",
                                   h4("Gene Input"),
                                   textAreaInput("gene_set_input", NULL,
                                                 placeholder = "Enter genes (one per line or comma-separated):\nCD3D\nCD3E\nCD8A",
                                                 height = "150px")
                               ),
                               div(class = "control-section",
                                   h4("Parameters"),
                                   selectInput("gene_set_method", "Method:",
                                               choices = c("Mean" = "mean",
                                                           "AddModuleScore" = "addmodulescore",
                                                           "GSVA (rank-based)" = "gsva")),
                                   tags$p(style = "font-size: 11px; color: #7f8c8d; margin-top: -8px; margin-bottom: 10px;",
                                          HTML("💡 <b>Mean</b>: Fastest, simple average<br>
                                                <b>AddModuleScore</b>: Fast, with control features<br>
                                                <b>GSVA</b>: Rank-based enrichment (may take 10-30s)")),
                                   textInput("gene_set_name", "Name:", value = "GeneSet1"),
                                   selectInput("geneset_color_scheme", "Color Scheme:",
                                               choices = c("Grey to Red" = "greyred",
                                                           "Blue-White-Red" = "bwr",
                                                           "Rainbow" = "rainbow"
                                               ),
                                               selected = "greyred"),
                                   actionButton("calculate_gene_set", "Calculate Score", class = "btn btn-primary btn-block"),
                                   br(),
                                   actionButton("save_gene_set", "Save to gene set", class = "btn btn-success btn-block")
                               ),
                               conditionalPanel(
                                 condition = "output.geneset_calculated",
                                 div(class = "control-section",
                                     h4("Color Legend"),
                                     plotOutput("geneset_color_legend", height = "150px")
                                 )
                               )
                      ),


                      # Clustering content
                      tags$div(class = "control-content", id = "content_cluster",
                               div(class = "panel-header", "📊 Clustering Analysis"),
                               div(class = "control-section",
                                   h4("Spot Selection"),
                                   selectInput("cluster_spot_selection", "Cluster which spots:",
                                               choices = c("All Spots" = "all",
                                                           "Group 1 Only" = "group1",
                                                           "Group 2 Only" = "group2"),
                                               selected = "all"),
                                   textOutput("cluster_spot_count"),
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                          "💡 Tip: Save spots to Group 1 or Group 2 first using the map buttons")
                               ),
                               div(class = "control-section",
                                   h4("Parameters"),
                                   numericInput("cluster_resolution", "Resolution:", value = 0.8, min = 0.1, max = 2, step = 0.1),
                                   numericInput("cluster_dims", "PCs:", value = 30, min = 5, max = 50, step = 5),
                                   actionButton("run_clustering", "Run Clustering", class = "btn btn-primary btn-block"),
                                   br(),
                                   checkboxInput("show_clusters", "Show on Map", value = TRUE)
                               ),
                               div(class = "control-section",
                                   h4("Results"),
                                   verbatimTextOutput("cluster_info"),
                                   conditionalPanel(
                                     condition = "output.clustering_done",
                                     plotOutput("cluster_umap", height = "300px")
                                   )
                               )
                      ),

                      # DEG content
                      tags$div(class = "control-content", id = "content_deg",
                               div(class = "panel-header", "📈 Differential Expression"),
                               div(class = "control-section",
                                   h4("Group Information"),
                                   div(style = "display: flex; align-items: center; margin: 10px 0;",
                                       div(style = "width: 20px; height: 20px; background: #E41A1C; border-radius: 50%; margin-right: 10px;"),
                                       div(textOutput("group1_info", inline = TRUE))
                                   ),
                                   div(style = "display: flex; align-items: center; margin: 10px 0;",
                                       div(style = "width: 20px; height: 20px; background: #377EB8; border-radius: 50%; margin-right: 10px;"),
                                       div(textOutput("group2_info", inline = TRUE))
                                   ),
                                   checkboxInput("show_groups", "Show on Map",  value = FALSE),
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 10px;",
                                          "💡 Tip: Use the group buttons at the bottom of the map to save selections and download spot IDs.")
                               ),
                               div(class = "control-section",
                                   h4("Analysis"),
                                   selectInput("deg_comparison", "Compare:",
                                               choices = c("Group 1 vs Rest",
                                                           "Group 2 vs Rest",
                                                           "Group 1 vs Group 2",
                                                           "Clusters in Group 1 (one vs rest)",
                                                           "Clusters in Group 2 (one vs rest)")),
                                   conditionalPanel(
                                     condition = "input.deg_comparison == 'Clusters in Group 1 (one vs rest)' || input.deg_comparison == 'Clusters in Group 2 (one vs rest)'",
                                     tags$p(style = "font-size: 12px; color: #e67e22; margin: 10px 0;",
                                            "⚠️ Clustering must be run first on the selected group")
                                   ),
                                   actionButton("run_deg", "Find DEGs", class = "btn btn-danger btn-block")
                               ),
                               conditionalPanel(
                                 condition = "output.deg_results_available",
                                 div(class = "control-section",
                                     h4("Top DEGs"),
                                     div(style = "max-height: 400px; overflow-y: auto;",
                                         tableOutput("deg_table")),
                                     br(),
                                     downloadButton("dl_deg", "Download DEG Results", class = "btn btn-warning btn-block")
                                 )
                               )
                      ),

                      # Compare content
                      tags$div(class = "control-content", id = "content_compare",
                               div(class = "panel-header", "⚖️ Feature or Group Comparison"),
                               div(class = "control-section",
                                   h4("Group vs Group"),
                                   selectInput("violin_feature_type", "Feature:",
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Gene Set" = "geneset")),
                                   conditionalPanel(
                                     condition = "input.violin_feature_type == 'gene'",
                                     selectInput("violin_gene", "Gene:", choices = c("", all_genes))
                                   ),
                                   conditionalPanel(
                                     condition = "input.violin_feature_type == 'metadata'",
                                     selectInput("violin_metadata", "Metadata:", choices = c("", all_metadata))
                                   ),
                                   conditionalPanel(
                                     condition = "input.violin_feature_type == 'geneset'",
                                     selectInput("violin_geneset", "Gene Set:", choices = character(0))
                                   ),
                                   selectInput("violin_comparison", "Compare:",
                                               choices = c("Group 1 vs Rest", "Group 2 vs Rest", "Group 1 vs Group 2")),
                                   selectInput("violin_stat_test", "Test:",
                                               choices = c("Wilcoxon" = "wilcox", "t-test" = "ttest"),
                                               selected = "wilcox"),
                                   actionButton("plot_violin", "Generate Plot", class = "btn btn-info btn-block"),
                                   conditionalPanel(
                                     condition = "output.violin_available",
                                     plotOutput("violin_plot", height = "300px")
                                   )
                               ),
                               div(class = "control-section",
                                   h4("Feature vs Feature"),
                                   selectInput("compare_type1", "Feature 1:",
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Gene Set" = "geneset")),
                                   conditionalPanel(
                                     condition = "input.compare_type1 == 'gene'",
                                     selectInput("compare_gene1", "Gene:", choices = c("", all_genes))
                                   ),
                                   conditionalPanel(
                                     condition = "input.compare_type1 == 'metadata'",
                                     selectInput("compare_meta1", "Metadata:", choices = c("", all_metadata))
                                   ),
                                   conditionalPanel(
                                     condition = "input.compare_type1 == 'geneset'",
                                     selectInput("compare_geneset1", "Gene Set:", choices = character(0))
                                   ),
                                   selectInput("compare_type2", "Feature 2:",
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Gene Set" = "geneset")),
                                   conditionalPanel(
                                     condition = "input.compare_type2 == 'gene'",
                                     selectInput("compare_gene2", "Gene:", choices = c("", all_genes))
                                   ),
                                   conditionalPanel(
                                     condition = "input.compare_type2 == 'metadata'",
                                     selectInput("compare_meta2", "Metadata:", choices = c("", all_metadata))
                                   ),
                                   conditionalPanel(
                                     condition = "input.compare_type2 == 'geneset'",
                                     selectInput("compare_geneset2", "Gene Set:", choices = character(0))
                                   ),
                                   selectInput("compare_spots_selection", "Use Spots:",
                                               choices = c("All Spots" = "all", "Group 1" = "group1", "Group 2" = "group2")),
                                   selectInput("stat_test", "Test:", choices = c("Wilcoxon" = "wilcox", "t-test" = "ttest")),
                                   # selectInput("cor_method", "Correlation:", choices = c("Spearman" = "spearman", "Pearson" = "pearson")),
                                   actionButton("plot_compare", "Compare", class = "btn btn-info btn-block"),
                                   conditionalPanel(
                                     condition = "output.compare_available",
                                     plotOutput("compare_plot", height = "400px")
                                   )
                               )
                      )
             )
    ),

    # JavaScript for panel management
    tags$script(HTML("
      $(document).ready(function() {
        // Button click handlers
        $('.sidebar-button').click(function() {
          var btnId = $(this).attr('id');
          var contentId = 'content_' + btnId.replace('btn_', '');

          // Toggle active button
          $('.sidebar-button').removeClass('active');
          $(this).addClass('active');

          // Show corresponding content
          $('.control-content').removeClass('active');
          $('#' + contentId).addClass('active');

          // Open panel
          $('#control_panel').addClass('open');
        });
      });

      Shiny.addCustomMessageHandler('closePanel', function(message) {
        $('#control_panel').removeClass('open');
        $('.sidebar-button').removeClass('active');
        $('#btn_home').addClass('active');
      });
    "))
  )

  server <- function(input, output, session) {

    # Print CellMarker 2.0 citation to console
    message("\n", paste(rep("=", 80), collapse = ""))
    message("Cell Marker Gene Signatures")
    message(paste(rep("=", 80), collapse = ""))
    message("Pre-defined signatures are curated from CellMarker 2.0 database")
    message("Database: 26,915 cell markers | 2,578 cell types | 656 tissues")
    message("")
    message("Citation:")
    message("  Hu C, Li T, Xu Y, et al. (2023)")
    message("  CellMarker 2.0: an updated database of manually curated cell markers")
    message("  in human/mouse and web tools based on scRNA-seq data.")
    message("  Nucleic Acids Res. 51(D1):D870-D876. doi: 10.1093/nar/gkac947")
    message("")
    message("Resources:")
    message("  Paper:    https://academic.oup.com/nar/article/51/D1/D870/6775381")
    message("  Database: http://bio-bigdata.hrbmu.edu.cn/CellMarker/")
    message("  Download: http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html")
    message(paste(rep("=", 80), collapse = ""), "\n")

    # Increase file upload size limit (default is 5MB, set to 500MB)
    options(shiny.maxRequestSize = 500*1024^2)  # 500MB in bytes

    # Hide initial loading screen after app is ready
    observe({
      # Wait a moment for everything to load
      Sys.sleep(1)

      # Hide the initial loading screen
      shinyjs::runjs("$('#initial_loading').addClass('loaded');")

      print("App fully loaded - showing landing page")
    }) %>%
      bindEvent(once = TRUE, {TRUE})

    # Global error handler - suppress termination errors when user stops the app
    options(shiny.error = function() {
      # Silently catch errors during app termination
      # This prevents the "stopifnot(length(class2) == 1L)" error from showing
      return(invisible(NULL))
    })

    # Suppress warnings during session to prevent termination messages
    options(warn = -1)

    # Reactive values
    current_sample_name <- reactiveVal(sample_name)  # Make sample_name reactive
    drawn_feats <- reactiveVal(list())
    selected_spots <- reactiveVal(character(0))
    current_values <- reactiveVal(NULL)
    is_categorical <- reactiveVal(FALSE)  # Track if current data is categorical
    category_labels <- reactiveVal(NULL)
    group1_spots <- reactiveVal(character(0))
    group2_spots <- reactiveVal(character(0))
    output_group1 <- reactiveVal(character(0))  # Stored copy for export
    output_group2 <- reactiveVal(character(0))  # Stored copy for export
    deg_results <- reactiveVal(NULL)
    violin_data <- reactiveVal(NULL)
    compare_data <- reactiveVal(NULL)
    gene_set_scores <- reactiveVal(list())
    current_gene_set_score <- reactiveVal(NULL)
    cluster_results <- reactiveVal(NULL)
    showing_gene_set <- reactiveVal(FALSE)  # Track if we're displaying gene set scores

    # Dynamic header title
    output$header_title <- renderText({
      paste("🔬 SpatialScope -", current_sample_name())
    })

    # FIXED: Observer to enable/disable download buttons based on group content
    observe({
      has_group1 <- length(output_group1()) > 0
      has_group2 <- length(output_group2()) > 0

      # Enable/disable Group 1 buttons
      if (has_group1) {
        shinyjs::enable("dl_group1")
        shinyjs::enable("dl_seurat_group1")
      } else {
        shinyjs::disable("dl_group1")
        shinyjs::disable("dl_seurat_group1")
      }

      # Enable/disable Group 2 buttons
      if (has_group2) {
        shinyjs::enable("dl_group2")
        shinyjs::enable("dl_seurat_group2")
      } else {
        shinyjs::disable("dl_group2")
        shinyjs::disable("dl_seurat_group2")
      }
    })

    # MODIFIED: Reactive for current signature library based on species
    current_signature_library <- reactive({
      if (input$species_select == "human") {
        signature_library_human
      } else {
        signature_library_mouse
      }
    })

    # MODIFIED: Update signature choices when species changes
    observe({
      updateSelectInput(session, "signature_library",
                        choices = names(current_signature_library()),
                        selected = "Custom")
    })

    # Start analysis from landing page - go directly to showing Home
    observeEvent(input$start_analysis_from_landing, {
      Sys.sleep(0.5)
      shinyjs::runjs("$('#loading_overlay').removeClass('active'); $('.app-content').addClass('active'); $('.top-header').addClass('active');")

      # Make sure Home is the active panel by default
      shinyjs::runjs("$('.sidebar-button').removeClass('active'); $('#btn_home').addClass('active'); $('.control-content').removeClass('active'); $('#content_home').addClass('active');")

      # Clear all reactive values
      drawn_feats(list())
      selected_spots(character(0))
      current_values(NULL)

      # Clear all map layers
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        clearGroup("selected") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")

      # Send clear message to JavaScript
      session$sendCustomMessage("clearFreehandDrawings", list())

      print("Analysis started - Home documentation shown first")
      showNotification("Welcome! Review the guide below, then click other sidebar tools to begin analysis.",
                       type = "message", duration = 5)
    })

    # Initialize when "Start Analysis" is clicked
    observeEvent(input$start_analysis_clicked, {
      print("Start Analysis clicked - initializing...")

      # Clear all reactive values immediately
      drawn_feats(list())
      selected_spots(character(0))
      current_values(NULL)

      # Clear all map layers immediately
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        clearGroup("selected") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")

      # Send clear message to JavaScript
      session$sendCustomMessage("clearFreehandDrawings", list())

      # Wait a moment for map to fully initialize
      Sys.sleep(0.5)

      # Hide loading overlay
      session$sendCustomMessage("hideLoading", list())

      print("Analysis started - all selections cleared and initialized")
      showNotification("Welcome to SpatialScope! Draw a region on the map to begin.",
                       type = "message", duration = 5)
    })

    # Initialize when "Start Analysis" is clicked
    observeEvent(input$start_analysis_clicked, {
      # Clear all reactive values
      drawn_feats(list())
      selected_spots(character(0))
      current_values(NULL)

      # Small delay to ensure map is loaded
      Sys.sleep(0.3)

      # Clear all map layers
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        clearGroup("selected") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")

      # Send clear message to JavaScript
      session$sendCustomMessage("clearFreehandDrawings", list())

      print("Analysis started - all selections cleared and initialized")
      showNotification("Welcome to SpatialScope! Draw a region on the map to begin.",
                       type = "message", duration = 5)
    })

    # Close panel handler
    observeEvent(input$close_panel, {
      session$sendCustomMessage("closePanel", list())
    })

    # Spot count display
    output$spot_count_display <- renderText({
      paste("Selected:", length(selected_spots()), "spots")
    })

    # Upload status display
    output$upload_status <- renderText({
      if (is.null(input$upload_seurat)) {
        "No file uploaded. Using current data."
      } else {
        paste("File ready:", input$upload_seurat$name)
      }
    })

    # Load uploaded Seurat object
    observeEvent(input$load_uploaded_seurat, {
      req(input$upload_seurat)

      # Show loading overlay with custom message
      shinyjs::runjs("
        $('#loading_message').text('Uploading and processing new dataset...');
        $('#loading_overlay').addClass('active');
      ")

      showNotification("Loading new Seurat object...", type = "message", duration = NULL, id = "load_seurat")

      # Add a small delay to ensure loading screen shows
      Sys.sleep(0.3)

      tryCatch({
        # Load the new Seurat object
        new_seurat <- readRDS(input$upload_seurat$datapath)

        # Update sample name from uploaded file
        uploaded_filename <- tools::file_path_sans_ext(input$upload_seurat$name)
        current_sample_name(uploaded_filename)

        # Validate it's a Seurat object
        if (!inherits(new_seurat, "Seurat")) {
          stop("Uploaded file is not a valid Seurat object")
        }

        # Check for spatial images
        if (length(new_seurat@images) == 0) {
          stop("No spatial images found in the uploaded Seurat object")
        }

        # Replace the global seurat object (dangerous but necessary for this use case)
        seurat_obj <<- new_seurat

        # Update loading message
        shinyjs::runjs("$('#loading_message').text('Updating gene and metadata lists...');")

        # Update all gene and metadata lists
        all_genes <- rownames(seurat_obj)
        all_metadata <- colnames(seurat_obj@meta.data)

        # Update the selectInputs with new gene/metadata lists
        updateSelectInput(session, "gene_select", choices = c("", all_genes))
        updateSelectInput(session, "violin_gene", choices = c("", all_genes))
        updateSelectInput(session, "compare_gene1", choices = c("", all_genes))
        updateSelectInput(session, "compare_gene2", choices = c("", all_genes))
        updateSelectInput(session, "meta_select", choices = c("", all_metadata))
        updateSelectInput(session, "violin_metadata", choices = c("", all_metadata))
        updateSelectInput(session, "compare_meta1", choices = c("", all_metadata))
        updateSelectInput(session, "compare_meta2", choices = c("", all_metadata))

        # Update loading message
        shinyjs::runjs("$('#loading_message').text('Extracting spatial coordinates...');")

        # Extract new coordinates
        image_name <- names(new_seurat@images)[1]  # Use the first image
        tryCatch({
          coords <- GetTissueCoordinates(new_seurat, image = image_name)
        }, error = function(e) {
          coords <- new_seurat@images[[image_name]]@coordinates
        })

        spots_df <- data.frame(spot_id = rownames(coords), stringsAsFactors = FALSE)
        if ("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
          spots_df$x <- coords$imagecol
          spots_df$y <- coords$imagerow
        } else if ("row" %in% colnames(coords) && "col" %in% colnames(coords)) {
          spots_df$x <- coords$col
          spots_df$y <- coords$row
        } else {
          spots_df$x <- coords[,1]
          spots_df$y <- coords[,2]
        }

        spots_sf <<- st_as_sf(spots_df, coords = c("x","y"), crs = NA)
        coords_matrix <- do.call(rbind, st_geometry(spots_sf)) %>% as.matrix()
        spots_sf$x <<- coords_matrix[, 1]
        spots_sf$y <<- coords_matrix[, 2]
        spots_sf$y <<- max(spots_sf$y) - spots_sf$y + min(spots_sf$y)

        # Update loading message
        shinyjs::runjs("$('#loading_message').text('Extracting H&E image...');")

        # Extract and update H&E image from new data
        tryCatch({
          if (show_image) {
            image_obj <- new_seurat@images[[image_name]]
            he_image_data <- image_obj@image
            coords_full <- GetTissueCoordinates(new_seurat, image = image_name)

            if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
              pixel_x <- coords_full$pxl_col_in_fullres
              pixel_y <- coords_full$pxl_row_in_fullres
            } else {
              pixel_x <- spots_sf$x
              pixel_y <- spots_sf$y
            }

            pixel_x_min <- min(pixel_x, na.rm = TRUE)
            pixel_x_max <- max(pixel_x, na.rm = TRUE)
            pixel_y_min <- min(pixel_y, na.rm = TRUE)
            pixel_y_max <- max(pixel_y, na.rm = TRUE)
            x_buffer_px <- (pixel_x_max - pixel_x_min) * 0.1
            y_buffer_px <- (pixel_y_max - pixel_y_min) * 0.1
            crop_x_min <- max(1, floor(pixel_x_min - x_buffer_px))
            crop_x_max <- min(dim(he_image_data)[2], ceiling(pixel_x_max + x_buffer_px))
            crop_y_min <- max(1, floor(pixel_y_min - y_buffer_px))
            crop_y_max <- min(dim(he_image_data)[1], ceiling(pixel_y_max + y_buffer_px))
            he_image_cropped <- he_image_data[crop_y_min:crop_y_max, crop_x_min:crop_x_max, ]

            temp_file <- tempfile(fileext = ".png")
            png(temp_file, width = dim(he_image_cropped)[2], height = dim(he_image_cropped)[1])
            par(mar = c(0,0,0,0))
            plot(as.raster(he_image_cropped), axes = FALSE)
            dev.off()

            he_image_base64 <<- paste0("data:image/png;base64,", base64enc::base64encode(temp_file))
            unlink(temp_file)

            he_image_bounds <<- list(
              south = crop_y_max,
              west = crop_x_min,
              north = crop_y_min,
              east = crop_x_max
            )

            # Update H&E image on map
            session$sendCustomMessage("updateHEImage", list(
              imageUrl = he_image_base64,
              bounds = he_image_bounds
            ))

            # Wait for H&E image to load
            Sys.sleep(0.5)
          }
        }, error = function(e) {
          print(paste("Could not extract H&E image from new data:", e$message))
          he_image_base64 <<- NULL
        })

        # Clear all selections and results
        drawn_feats(list())
        selected_spots(character(0))
        group1_spots(character(0))
        group2_spots(character(0))
        current_values(NULL)
        deg_results(NULL)
        violin_data(NULL)
        compare_data(NULL)
        gene_set_scores(list())
        current_gene_set_score(NULL)
        cluster_results(NULL)

        # Update loading message
        shinyjs::runjs("$('#loading_message').text('Rendering map with new spots...');")

        # Clear map and redraw with new data
        x_range <- range(spots_sf$x)
        y_range <- range(spots_sf$y)
        x_buffer <- diff(x_range) * 0.1
        y_buffer <- diff(y_range) * 0.1

        leafletProxy("map") %>%
          clearGroup("spots") %>%
          clearGroup("drawn") %>%
          clearGroup("selected") %>%
          clearGroup("group1_display") %>%
          clearGroup("group2_display") %>%
          fitBounds(
            lng1 = x_range[1] - x_buffer, lat1 = y_range[1] - y_buffer,
            lng2 = x_range[2] + x_buffer, lat2 = y_range[2] + y_buffer
          ) %>%
          addCircleMarkers(
            lng = spots_sf$x, lat = spots_sf$y,
            radius = 5, stroke = TRUE, color = "black", weight = 0.5,
            fillColor = "lightblue", fillOpacity = 0.8, group = "spots"
          )

        # Update loading message for final step
        shinyjs::runjs("$('#loading_message').text('Finalizing visualization...');")

        removeNotification(id = "load_seurat")

        # Use JavaScript setTimeout with longer delay to ensure map fully renders
        # This gives the browser enough time to render all spots and H&E image
        shinyjs::runjs("
          setTimeout(function() {
            $('#loading_overlay').removeClass('active');
          }, 5000);
        ")

        showNotification(paste("Successfully loaded", current_sample_name(), "with",
                               nrow(spots_sf), "spots. All analyses will now use this new dataset."),
                         type = "message", duration = 7)

      }, error = function(e) {
        removeNotification(id = "load_seurat")

        # Hide loading overlay on error
        shinyjs::runjs("$('#loading_overlay').removeClass('active');")

        showNotification(paste("Error loading Seurat object:", e$message),
                         type = "error", duration = 10)
      })
    })

    # Use Example Data button handler
    observeEvent(input$use_example_data, {
      showNotification("You are currently using the example dataset. The example data is loaded by default.",
                       type = "message", duration = 5)
    })

    # Selection summary
    output$selection_summary <- renderText({
      paste("Total selected spots:", length(selected_spots()))
    })

    output$selected_spots_table <- renderTable({
      sel <- selected_spots()
      if (length(sel) > 0) {
        data.frame(`Spot ID` = head(sel, 50), check.names = FALSE)
      } else {
        data.frame(`Spot ID` = character(0), check.names = FALSE)
      }
    }, rownames = FALSE)

    # Clear selection handler (from the map button)
    observeEvent(input$clear_selection_click, {
      drawn_feats(list())
      selected_spots(character(0))
      group1_spots(character(0))  # Clear Group 1
      group2_spots(character(0))  # Clear Group 2
      session$sendCustomMessage("clearFreehandDrawings", list())
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        clearGroup("selected") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")
      showNotification("Selection and all groups cleared", type = "message", duration = 2)
    })

    # Clustering
    observeEvent(input$run_clustering, {
      showNotification("Running clustering...", type = "message", duration = NULL, id = "clustering_run")

      tryCatch({
        # Determine which spots to cluster based on selection
        spots_to_cluster <- switch(input$cluster_spot_selection,
                                   "all" = spots_sf$spot_id,
                                   "group1" = output_group1(),
                                   "group2" = output_group2()
        )

        # Validate that we have spots to cluster
        if (length(spots_to_cluster) == 0) {
          stop(paste("No spots in",
                     switch(input$cluster_spot_selection,
                            "all" = "dataset",
                            "group1" = "Group 1. Please save spots to Group 1 first.",
                            "group2" = "Group 2. Please save spots to Group 2 first.")))
        }

        spatial_obj <- subset(seurat_obj, cells = spots_to_cluster)

        # Get the default assay name dynamically
        default_assay <- DefaultAssay(spatial_obj)

        # Check if data needs preprocessing (works with any assay name: Spatial, RNA, etc.)
        if (!"SCT" %in% names(spatial_obj@assays)) {
          current_assay <- spatial_obj@assays[[default_assay]]
          scale_data <- slot(current_assay, "scale.data")
          var_features <- slot(current_assay, "var.features")

          # Check if scale.data is empty or variable features are missing
          if (is.null(scale_data) || length(scale_data) == 0 ||
              is.null(var_features) || length(var_features) == 0) {
            spatial_obj <- NormalizeData(spatial_obj, verbose = FALSE)
            spatial_obj <- Seurat::FindVariableFeatures(spatial_obj, verbose = FALSE)
            spatial_obj <- ScaleData(spatial_obj, verbose = FALSE)
          }
        }

        if (is.null(spatial_obj@reductions$pca)) {
          spatial_obj <- Seurat::RunPCA(spatial_obj, verbose = FALSE)
        }

        spatial_obj <- FindNeighbors(spatial_obj, dims = 1:input$cluster_dims, verbose = FALSE)
        spatial_obj <- FindClusters(spatial_obj, resolution = input$cluster_resolution, verbose = FALSE)

        if (is.null(spatial_obj@reductions$umap)) {
          spatial_obj <- RunUMAP(spatial_obj, dims = 1:input$cluster_dims, verbose = FALSE)
        }


        clusters <- Idents(spatial_obj)

        # Store clusters in metadata with info about which spots were clustered
        cluster_col_name <- paste0("seurat_clusters_res", input$cluster_resolution)
        seurat_obj@meta.data[[cluster_col_name]] <<- NA  # Initialize with NA
        seurat_obj@meta.data[names(clusters), cluster_col_name] <<- as.character(clusters)

        cluster_results(list(
          seurat = spatial_obj,
          clusters = clusters,
          resolution = input$cluster_resolution,
          n_clusters = length(unique(clusters)),
          spot_selection = input$cluster_spot_selection,
          n_spots = length(spots_to_cluster)
        ))

        removeNotification(id = "clustering_run")

        selection_text <- switch(input$cluster_spot_selection,
                                 "all" = "all spots",
                                 "group1" = "Group 1 spots",
                                 "group2" = "Group 2 spots"
        )

        showNotification(paste("Found", length(unique(clusters)), "clusters in",
                               length(spots_to_cluster), selection_text),
                         type = "message", duration = 5)

      }, error = function(e) {
        removeNotification(id = "clustering_run")
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

    output$cluster_info <- renderText({
      cluster_res <- cluster_results()
      if (is.null(cluster_res)) {
        "No clustering results yet"
      } else {
        selection_text <- switch(cluster_res$spot_selection,
                                 "all" = "All spots",
                                 "group1" = "Group 1 spots only",
                                 "group2" = "Group 2 spots only"
        )

        paste0("Spot Selection: ", selection_text, "\n",
               "Number of spots: ", cluster_res$n_spots, "\n",
               "Resolution: ", cluster_res$resolution, "\n",
               "Clusters: ", cluster_res$n_clusters, "\n\n",
               paste(capture.output(table(cluster_res$clusters)), collapse = "\n"))
      }
    })

    output$clustering_done <- reactive({
      !is.null(cluster_results())
    })
    outputOptions(output, "clustering_done", suspendWhenHidden = FALSE)

    # Show spot count for selected clustering option
    output$cluster_spot_count <- renderText({
      spot_count <- switch(input$cluster_spot_selection,
                           "all" = nrow(spots_sf),
                           "group1" = length(output_group1()),
                           "group2" = length(output_group2())
      )

      if (spot_count == 0 && input$cluster_spot_selection != "all") {
        paste("⚠️ No spots saved in",
              ifelse(input$cluster_spot_selection == "group1", "Group 1", "Group 2"))
      } else {
        paste("✓", spot_count, "spots will be clustered")
      }
    })

    output$cluster_umap <- renderPlot({
      cluster_res <- cluster_results()
      if (!is.null(cluster_res)) {
        DimPlot(cluster_res$seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
          ggtitle(paste("UMAP - Resolution:", cluster_res$resolution))
      }
    })

    observe({
      if (input$show_clusters) {
        cluster_res <- cluster_results()

        if (is.null(cluster_res)) {
          # showNotification("Run clustering first!", type = "warning")
          updateCheckboxInput(session, "show_clusters",  value = FALSE)
          return()
        }

        # Clear groups when showing clusters
        leafletProxy("map") %>%
          clearGroup("group1_display") %>%
          clearGroup("group2_display")
        updateCheckboxInput(session, "show_groups", value = FALSE)

        clusters <- cluster_res$clusters
        n_clusters <- cluster_res$n_clusters

        cluster_levels <- levels(clusters)
        if (n_clusters <= 12) {
          cluster_colors <- RColorBrewer::brewer.pal(min(12, max(3, n_clusters)), "Set3")
        } else {
          cluster_colors <- rainbow(n_clusters, s = 1, v = 0.9)
        }
        names(cluster_colors) <- cluster_levels

        spot_colors <- rep("lightgrey", nrow(spots_sf))
        for (i in 1:nrow(spots_sf)) {
          spot_id <- spots_sf$spot_id[i]
          if (spot_id %in% names(clusters)) {
            cluster_id <- as.character(clusters[spot_id])
            if (cluster_id %in% names(cluster_colors)) {
              spot_colors[i] <- cluster_colors[cluster_id]
            }
          }
        }

        leafletProxy("map") %>%
          clearGroup("spots") %>%
          addCircleMarkers(
            lng = spots_sf$x, lat = spots_sf$y,
            radius = input$spot_size, stroke = TRUE, color = "#333333", weight = 0.5,
            fillColor = spot_colors, fillOpacity = 1, group = "spots"
          )
      } else {
        update_map_colors()
      }
    })

    # Load signature from library - MODIFIED: Use current_signature_library()
    observeEvent(input$load_signature, {
      selected_sig <- input$signature_library
      sig_library <- current_signature_library()

      if (selected_sig != "Custom" && selected_sig %in% names(sig_library)) {
        genes <- sig_library[[selected_sig]]$genes

        # Update gene input area
        updateTextAreaInput(session, "gene_set_input",
                            value = paste(genes, collapse = "\n"))

        # Update gene set name
        clean_name <- gsub(":", "", selected_sig)
        clean_name <- gsub(" ", "_", clean_name)
        updateTextInput(session, "gene_set_name", value = clean_name)

        showNotification(paste("Loaded signature:", selected_sig, "with", length(genes), "genes",
                               "(", input$species_select, ")"),
                         type = "message", duration = 3)
      }
    })

    # Gene set functions
    calculate_gene_set_score <- function(genes, method = "mean") {
      genes_in_data <- intersect(genes, rownames(seurat_obj))
      genes_missing <- setdiff(genes, rownames(seurat_obj))

      if (length(genes_in_data) == 0) {
        stop("None of the provided genes are found in the dataset")
      }

      if (length(genes_in_data) < length(genes)) {
        showNotification(paste("Only", length(genes_in_data), "out of", length(genes),
                               "genes found"), type = "warning", duration = 5)
      }

      if (method == "mean") {
        expr_data <- FetchData(seurat_obj, vars = genes_in_data, slot = "data")
        scores <- rowMeans(expr_data, na.rm = TRUE)

      } else if (method == "addmodulescore") {
        temp_seurat <- AddModuleScore(seurat_obj, features = list(genes_in_data),
                                      name = "GeneSet", assay = DefaultAssay(seurat_obj))
        scores <- temp_seurat$GeneSet1

      } else if (method == "gsva") {
        # VECTORIZED GSVA ALTERNATIVE: Fast rank-based enrichment score
        # This uses matrix operations to rank all cells at once (10-50x faster than loop)

        # Get expression matrix
        expr_mat <- GetAssayData(seurat_obj, slot = "data", assay = DefaultAssay(seurat_obj))

        # VECTORIZED RANKING: Rank all genes for all cells simultaneously
        # Wrap in tryCatch to handle interruptions gracefully
        rank_mat <- tryCatch({
          # apply(expr_mat, 2, rank) ranks each column (cell) independently
          apply(expr_mat, 2, rank, ties.method = "average")
        }, error = function(e) {
          # If interrupted, return NULL to trigger proper error handling upstream
          stop("Calculation interrupted")
        })

        # Get ranks for gene set genes across all cells
        gene_set_ranks <- rank_mat[genes_in_data, , drop = FALSE]

        # Calculate mean rank of gene set for each cell (vectorized)
        mean_rank_geneset <- colMeans(gene_set_ranks)

        # Calculate overall mean rank for each cell (should be constant: (n+1)/2)
        mean_rank_all <- colMeans(rank_mat)

        # Calculate enrichment scores (vectorized across all cells)
        # Normalize to [-1, 1] scale like GSVA
        scores <- (mean_rank_geneset - mean_rank_all) / (nrow(expr_mat) / 2)

        names(scores) <- colnames(expr_mat)
      }

      return(list(
        scores = scores,
        genes_found = genes_in_data,
        genes_missing = genes_missing,
        n_found = length(genes_in_data),
        n_missing = length(genes_missing),
        n_total = length(genes)
      ))
    }

    observeEvent(input$calculate_gene_set, {
      req(input$gene_set_input)

      genes <- unlist(strsplit(input$gene_set_input, "[,\n\r\t ]+"))
      genes <- genes[genes != ""]

      if (length(genes) == 0) {
        showNotification("Please enter at least one gene", type = "warning")
        return()
      }

      showNotification("Calculating gene set scores...", type = "message", id = "calc_gene_set")

      tryCatch({
        # Suppress all warnings and messages during calculation to prevent termination errors
        result <- suppressWarnings(suppressMessages(
          calculate_gene_set_score(genes, input$gene_set_method)
        ))
        current_gene_set_score(result$scores)

        removeNotification(id = "calc_gene_set")
        # Create detailed notification message
        msg <- paste0("✓ Scores calculated using ", input$gene_set_method, "\n",
                      "Found: ", result$n_found, "/", result$n_total, " genes")

        if (result$n_missing > 0) {
          missing_genes_str <- paste(result$genes_missing, collapse = ", ")
          if (nchar(missing_genes_str) > 100) {
            missing_genes_str <- paste0(paste(head(result$genes_missing, 10), collapse = ", "),
                                        ", ... (", result$n_missing - 10, " more)")
          }
          msg <- paste0(msg, "\n\nMissing genes (", result$n_missing, "):\n", missing_genes_str)
        }

        showNotification(msg, type = "message", duration = 10)

        # Also print to console for reference
        if (result$n_missing > 0) {
          cat("\n=== Gene Set Calculation ===\n")
          cat("Total genes provided:", result$n_total, "\n")
          cat("Genes found:", result$n_found, "\n")
          cat("Genes missing:", result$n_missing, "\n")
          cat("\nMissing genes:\n")
          cat(paste(result$genes_missing, collapse = ", "), "\n")
          cat("===========================\n\n")
        }

        # CLEAR visualization settings first to fully separate from visualization function
        updateSelectInput(session, "feature_type", selected = "None")
        updateSelectInput(session, "gene_select", selected = "")
        updateSelectInput(session, "meta_select", selected = "")

        # Set gene set data as current display
        current_values(result$scores)
        is_categorical(FALSE)  # Gene set scores are always continuous
        showing_gene_set(TRUE)  # Set flag to indicate we're showing gene set
        update_map_colors()


      }, error = function(e) {
        removeNotification(id = "calc_gene_set")
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

    observeEvent(input$save_gene_set, {
      req(current_gene_set_score())
      req(input$gene_set_name)

      gene_set_name <- make.names(input$gene_set_name)
      scores <- current_gene_set_score()

      current_list <- gene_set_scores()
      current_list[[gene_set_name]] <- scores
      gene_set_scores(current_list)

      seurat_obj@meta.data[[gene_set_name]] <<- scores

      # Update all gene set selectors
      updateSelectInput(session, "geneset_select",
                        choices = names(gene_set_scores()),
                        selected = gene_set_name)
      updateSelectInput(session, "violin_geneset",
                        choices = names(gene_set_scores()),
                        selected = gene_set_name)
      updateSelectInput(session, "compare_geneset1",
                        choices = names(gene_set_scores()),
                        selected = gene_set_name)
      updateSelectInput(session, "compare_geneset2",
                        choices = names(gene_set_scores()),
                        selected = gene_set_name)

      showNotification(paste("Gene set saved as:", gene_set_name), type = "message")
    })

    # Color functions
    get_color_palette <- function(scheme, n = 100) {
      if (scheme == "greyred") {
        colorRampPalette(c("grey90", "red"))(n)
      } else if (scheme == "bwr") {
        colorRampPalette(c("blue", "white", "red"))(n)
      } else if (scheme == "rainbow") {
        colorRampPalette(c("darkblue", "cyan", "green", "yellow", "orange", "red"))(n)
      }
    }

    map_to_colors <- function(scheme, values) {
      if (all(is.na(values))) return(rep("lightgrey", length(values)))

      # Check if values are categorical (stored as character/factor labels)
      if (is_categorical()) {
        # Use discrete colors for categorical data
        unique_cats <- unique(values[!is.na(values)])
        n_cats <- length(unique_cats)

        # Handle case with no categories
        if (n_cats == 0) {
          return(rep("lightgrey", length(values)))
        }

        # Generate distinct colors for each category
        if (n_cats <= 8) {
          # Use a standard discrete palette for small number of categories
          cat_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                          "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[1:n_cats]
        } else {
          # For more categories, generate colors from rainbow palette
          cat_colors <- rainbow(n_cats)
        }

        # Create a named vector mapping categories to colors
        names(cat_colors) <- as.character(unique_cats)

        # Map colors to values - convert both to character for matching
        colors <- cat_colors[as.character(values)]

        # Replace NA colors (from failed matches) with lightgrey
        colors[is.na(colors) | is.na(values)] <- "lightgrey"

        # Ensure all colors are valid (convert to character and remove names)
        colors <- as.character(colors)

        return(colors)
      } else {
        # Continuous data - use gradient
        pal <- get_color_palette(scheme, 100)
        val_range <- range(values, na.rm = TRUE)

        if (val_range[1] == val_range[2]) {
          norm_vals <- rep(0.5, length(values))
        } else {
          norm_vals <- (values - val_range[1]) / (val_range[2] - val_range[1])
        }

        color_indices <- pmax(1, pmin(100, ceiling(norm_vals * 100)))
        colors <- pal[color_indices]
        colors[is.na(values)] <- "lightgrey"

        return(colors)
      }
    }
    observeEvent({
      input$feature_type
      input$gene_select
      input$meta_select
      input$geneset_select
      input$color_scheme
      input$geneset_color_scheme
      input$spot_size
    }, {
      if (input$feature_type == "Gene Expression" && !is.null(input$gene_select) && input$gene_select != "") {
        gene_expr <- FetchData(seurat_obj, vars = input$gene_select, slot = "data")
        values <- gene_expr[spots_sf$spot_id, 1]
        is_categorical(FALSE)  # Gene expression is continuous
        current_values(values)
        showing_gene_set(FALSE)  # Not showing gene set


      } else if (input$feature_type == "Metadata" && !is.null(input$meta_select) && input$meta_select != "") {
        meta_vals <- seurat_obj@meta.data[spots_sf$spot_id, input$meta_select]
        if (is.factor(meta_vals) || is.character(meta_vals)) {
          # Keep as categorical - don't convert to numeric
          is_categorical(TRUE)
          current_values(as.character(meta_vals))
        } else {
          # Numeric metadata
          is_categorical(FALSE)
          current_values(meta_vals)
        }
        showing_gene_set(FALSE)  # Not showing gene set

      } else if (input$feature_type == "Gene Set" && !is.null(input$geneset_select) && input$geneset_select != "") {
        scores <- gene_set_scores()[[input$geneset_select]]
        is_categorical(FALSE)  # Gene set scores are continuous
        current_values(scores)
        showing_gene_set(FALSE)  # This is old gene set from removed feature

      } else if (input$feature_type == "Clustering") {
        cluster_res <- cluster_results()
        if (!is.null(cluster_res)) {
          clusters <- cluster_res$clusters
          # Keep clusters as categorical labels
          cluster_vals <- as.character(clusters)[match(spots_sf$spot_id, names(clusters))]
          is_categorical(TRUE)  # Clusters are categorical
          current_values(cluster_vals)
          showing_gene_set(FALSE)  # Not showing gene set
        } else {
          showNotification("Run clustering first!", type = "warning")
          current_values(NULL)
          showing_gene_set(FALSE)
        }

      } else {
        # If geneset_color_scheme changed while showing gene set, just update colors
        if (showing_gene_set() && !is.null(current_values())) {
          # Don't change current_values, just update colors by calling update_map_colors
        } else {
          current_values(NULL)
          showing_gene_set(FALSE)
        }
      }

      # Always update map colors when any color scheme or feature changes
      update_map_colors()
    }, ignoreInit = TRUE)


    update_map_colors <- function() {
      values <- current_values()
      spot_radius <- input$spot_size

      # Clear groups when updating visualization
      leafletProxy("map") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")

      updateCheckboxInput(session, "show_groups", value = FALSE)

      if (is.null(values)) {
        colors <- rep("lightblue", nrow(spots_sf))
      } else {
        # Use geneset_color_scheme if showing gene set, otherwise use color_scheme
        scheme <- if (showing_gene_set()) {
          input$geneset_color_scheme
        } else {
          input$color_scheme
        }
        colors <- map_to_colors(scheme, values)
      }

      # Create popup text based on whether data is categorical or continuous
      popup_text <- if (is.null(values)) {
        paste0("Spot: ", spots_sf$spot_id)
      } else if (is_categorical()) {
        # For categorical data, show the category label directly
        paste0("Spot: ", spots_sf$spot_id, "<br>Category: ", values)
      } else {
        # For continuous data, round to 3 decimal places
        paste0("Spot: ", spots_sf$spot_id, "<br>Value: ", round(values, 3))
      }

      leafletProxy("map") %>%
        clearGroup("spots") %>%
        addCircleMarkers(
          lng = spots_sf$x, lat = spots_sf$y,
          radius = spot_radius, stroke = TRUE, color = "black", weight = 0.5,
          fillColor = colors, fillOpacity = 0.8, group = "spots",
          popup = popup_text
        )
    }

    output$color_legend <- renderPlot({
      values <- current_values()
      if (is.null(values)) return(NULL)

      # Check if data is categorical
      if (is_categorical()) {
        # Create categorical legend
        unique_cats <- unique(values[!is.na(values)])
        n_cats <- length(unique_cats)

        # Generate same colors as map_to_colors function
        if (n_cats <= 8) {
          cat_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                          "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[1:n_cats]
        } else {
          cat_colors <- rainbow(n_cats)
        }

        # Create legend data
        legend_data <- data.frame(
          category = factor(unique_cats, levels = unique_cats),
          color = cat_colors,
          x = 1,
          y = seq_along(unique_cats)
        )

        ggplot(legend_data, aes(x = x, y = y, fill = category)) +
          geom_tile(width = 0.8, height = 0.8) +
          scale_fill_manual(values = setNames(cat_colors, unique_cats),
                            name = "",
                            labels = unique_cats) +
          theme_void() +
          theme(
            legend.position = "bottom",
            legend.text = element_text(size = 9),
            legend.title = element_blank(),
            legend.key.size = unit(0.5, "cm")
          ) +
          guides(fill = guide_legend(nrow = ceiling(n_cats / 4)))
      } else {
        # Continuous data - use gradient
        val_range <- range(values, na.rm = TRUE)
        n_colors <- 100

        # Use geneset_color_scheme if showing gene set, otherwise use color_scheme
        scheme <- if (showing_gene_set()) {
          input$geneset_color_scheme
        } else {
          input$color_scheme
        }
        color_pal <- get_color_palette(scheme, n_colors)

        # Create a horizontal gradient bar
        legend_data <- data.frame(
          x = 1:n_colors,
          y = 1,
          value = seq(val_range[1], val_range[2], length.out = n_colors)
        )

        ggplot(legend_data, aes(x = x, y = y, fill = value)) +
          geom_tile() +
          scale_fill_gradientn(
            colors = color_pal,
            name = "",
            breaks = seq(val_range[1], val_range[2], length.out = 5),
            labels = round(seq(val_range[1], val_range[2], length.out = 5), 2)
          ) +
          theme_void() +
          theme(
            legend.position = "bottom",
            legend.key.width = unit(3, "cm"),
            legend.key.height = unit(0.5, "cm"),
            legend.text = element_text(size = 9),
            legend.title = element_blank()
          ) +
          guides(fill = guide_colorbar(
            barwidth = 15,
            barheight = 0.8,
            title.position = "top",
            label.position = "bottom"
          ))
      }
    }, bg = "transparent")

    # Gene set color legend (separate from main color_legend)
    output$geneset_color_legend <- renderPlot({
      # Add explicit dependency on geneset_color_scheme
      input$geneset_color_scheme

      scores <- current_gene_set_score()
      if (is.null(scores)) return(NULL)

      val_range <- range(scores, na.rm = TRUE)
      n_colors <- 100
      color_pal <- get_color_palette(input$geneset_color_scheme, n_colors)

      # Create a horizontal gradient bar
      legend_data <- data.frame(
        x = 1:n_colors,
        y = 1,
        value = seq(val_range[1], val_range[2], length.out = n_colors)
      )

      ggplot(legend_data, aes(x = x, y = y, fill = value)) +
        geom_tile() +
        scale_fill_gradientn(
          colors = color_pal,
          name = "",
          breaks = seq(val_range[1], val_range[2], length.out = 5),
          labels = round(seq(val_range[1], val_range[2], length.out = 5), 2)
        ) +
        theme_void() +
        theme(
          legend.position = "bottom",
          legend.key.width = unit(3, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          plot.margin = margin(5, 5, 5, 5)
        ) +
        guides(fill = guide_colorbar(
          barwidth = 15,
          barheight = 0.8,
          title.position = "top",
          label.position = "bottom"
        ))
    }, bg = "transparent")

    # Reactive to control visibility of geneset color legend
    output$geneset_calculated <- reactive({
      !is.null(current_gene_set_score())
    })
    outputOptions(output, "geneset_calculated", suspendWhenHidden = FALSE)

    # Map rendering
    output$map <- renderLeaflet({
      m <- leaflet(options = leafletOptions(
        crs = leafletCRS(crsClass = "L.CRS.Simple"),
        zoomControl = TRUE
      )) %>%
        fitBounds(
          lng1 = x_range[1] - x_buffer, lat1 = y_range[1] - y_buffer,
          lng2 = x_range[2] + x_buffer, lat2 = y_range[2] + y_buffer
        )

      # Add freehand drawing
      m <- m %>%
        htmlwidgets::onRender("
          function(el, x) {
            var map = this;

            if (!map.drawnItems) {
              map.drawnItems = new L.FeatureGroup();
              map.addLayer(map.drawnItems);
            }

            L.Control.FreehandDraw = L.Control.extend({
              onAdd: function(map) {
                var container = L.DomUtil.create('div', 'leaflet-bar leaflet-control');
                var button = L.DomUtil.create('a', 'leaflet-draw-draw-polygon', container);
                button.href = '#';
                button.title = 'Draw freehand polygon';
                button.innerHTML = '✏️';
                button.style.fontSize = '18px';

                L.DomEvent.on(button, 'click', function(e) {
                  L.DomEvent.stopPropagation(e);
                  L.DomEvent.preventDefault(e);
                  startFreehandDraw();
                });

                return container;
              }
            });

            if (!map.freehandControl) {
              map.freehandControl = new L.Control.FreehandDraw({ position: 'topleft' });
              map.addControl(map.freehandControl);
            }

            var isDrawing = false;
            var freehandPoints = [];
            var tempPolyline = null;

            function startFreehandDraw() {
              map.dragging.disable();
              map.getContainer().style.cursor = 'crosshair';
              isDrawing = true;
              freehandPoints = [];
            }

            map.on('mousedown', function(e) {
              if (!isDrawing) return;
              freehandPoints = [e.latlng];
              tempPolyline = L.polyline(freehandPoints, {
                color: '#ff0000',
                weight: 2
              }).addTo(map);
            });

            map.on('mousemove', function(e) {
              if (!isDrawing || freehandPoints.length === 0) return;
              freehandPoints.push(e.latlng);
              tempPolyline.setLatLngs(freehandPoints);
            });

            map.on('mouseup', function(e) {
              if (!isDrawing || freehandPoints.length < 3) {
                if (tempPolyline) map.removeLayer(tempPolyline);
                isDrawing = false;
                map.dragging.enable();
                map.getContainer().style.cursor = '';
                return;
              }

              freehandPoints.push(freehandPoints[0]);
              var polygon = L.polygon(freehandPoints, {
                color: '#ff0000',
                weight: 2,
                fillOpacity: 0.3
              });

              map.drawnItems.addLayer(polygon);
              if (tempPolyline) map.removeLayer(tempPolyline);

              var feature = polygon.toGeoJSON();
              feature.properties = feature.properties || {};
              feature.properties._leaflet_id = polygon._leaflet_id;

              Shiny.setInputValue('map_draw_new_feature', feature, {priority: 'event'});

              isDrawing = false;
              freehandPoints = [];
              map.dragging.enable();
              map.getContainer().style.cursor = '';
            });

            Shiny.addCustomMessageHandler('clearFreehandDrawings', function(message) {
              if (map.drawnItems) {
                map.drawnItems.clearLayers();
              }
            });

            Shiny.addCustomMessageHandler('toggleHEImage', function(message) {
              if (window.heImageOverlay && window.leafletMap) {
                if (message.show) {
                  if (!window.leafletMap.hasLayer(window.heImageOverlay)) {
                    window.heImageOverlay.addTo(window.leafletMap);
                    window.heImageOverlay.bringToBack();
                  }
                } else {
                  if (window.leafletMap.hasLayer(window.heImageOverlay)) {
                    window.leafletMap.removeLayer(window.heImageOverlay);
                  }
                }
              }
            });

            Shiny.addCustomMessageHandler('updateHEOpacity', function(message) {
              if (window.heImageOverlay) {
                window.heImageOverlay.setOpacity(message.opacity);
              }
            });

            Shiny.addCustomMessageHandler('updateHEImage', function(message) {
              if (window.leafletMap) {
                // Remove old overlay if exists
                if (window.heImageOverlay && window.leafletMap.hasLayer(window.heImageOverlay)) {
                  window.leafletMap.removeLayer(window.heImageOverlay);
                }

                // Create new overlay with updated image and bounds
                window.heImageData = message.imageUrl;
                window.heImageBounds = [[message.bounds.south, message.bounds.west],
                                       [message.bounds.north, message.bounds.east]];

                window.heImageOverlay = L.imageOverlay(window.heImageData, window.heImageBounds, {
                  opacity: 0.6,
                  interactive: false
                });

                // Add to map if show_he is checked
                // Always add to map (no checkbox needed)
                window.heImageOverlay.addTo(window.leafletMap);
                window.heImageOverlay.bringToBack();
              }
            });

            Shiny.addCustomMessageHandler('hideLoading', function(message) {
              $('#loading_overlay').removeClass('active');
            });
          }
        ")

      # Add H&E image
      if (!is.null(he_image_base64)) {
        m <- m %>%
          htmlwidgets::onRender(paste0("
            function(el, x) {
              var map = this;
              window.heImageData = '", he_image_base64, "';
              window.heImageBounds = [[", he_image_bounds$south, ", ", he_image_bounds$west, "],
                                      [", he_image_bounds$north, ", ", he_image_bounds$east, "]];

              if (!window.heImageOverlay) {
                window.heImageOverlay = L.imageOverlay(window.heImageData, window.heImageBounds, {
                  opacity: 0.6,
                  interactive: false
                });
                window.heImageOverlay.addTo(map);
                window.heImageOverlay.bringToBack();
              }

              window.leafletMap = map;
            }
          "))
      }

      m %>% addCircleMarkers(
        lng = spots_sf$x, lat = spots_sf$y,
        radius = 5, stroke = TRUE, color = "black", weight = 0.5,
        fillColor = "lightblue", fillOpacity = 0.8, group = "spots"
      )
    })

    # Drawing handlers
    recompute_selection <- function() {
      feats <- drawn_feats()
      if (length(feats) == 0) {
        selected_spots(character(0))
        return()
      }

      all_selected <- character(0)

      for (feat in feats) {
        if (feat$type == "Feature") {
          geom_type <- feat$geometry$type
          coords <- feat$geometry$coordinates

          if (geom_type == "Polygon") {
            poly_coords <- coords[[1]]
            poly_x <- sapply(poly_coords, function(p) p[1])
            poly_y <- sapply(poly_coords, function(p) p[2])

            inside <- sp::point.in.polygon(spots_sf$x, spots_sf$y, poly_x, poly_y) > 0
            inside_spots <- spots_sf$spot_id[inside]
            all_selected <- c(all_selected, inside_spots)

          } else if (geom_type == "Rectangle" || (geom_type == "Polygon" && length(coords[[1]]) == 5)) {
            poly_coords <- coords[[1]]
            poly_x <- sapply(poly_coords, function(p) p[1])
            poly_y <- sapply(poly_coords, function(p) p[2])

            x_min <- min(poly_x)
            x_max <- max(poly_x)
            y_min <- min(poly_y)
            y_max <- max(poly_y)

            inside <- spots_sf$x >= x_min & spots_sf$x <= x_max &
              spots_sf$y >= y_min & spots_sf$y <= y_max
            inside_spots <- spots_sf$spot_id[inside]
            all_selected <- c(all_selected, inside_spots)
          }
        }
      }

      unique_selected <- unique(all_selected)
      selected_spots(unique_selected)

      # Print for debugging
      print(paste("Selected", length(unique_selected), "spots"))
    }

    observeEvent(input$map_draw_new_feature, {
      feature <- input$map_draw_new_feature
      current_feats <- drawn_feats()
      drawn_feats(c(current_feats, list(feature)))

      # Force immediate recomputation
      isolate({
        recompute_selection()
      })
    })

    observeEvent(input$map_draw_edited_features, {
      edited <- input$map_draw_edited_features
      if (!is.null(edited) && length(edited$features) > 0) {
        drawn_feats(edited$features)
        isolate({
          recompute_selection()
        })
      }
    })

    observeEvent(input$map_draw_deleted_features, {
      deleted <- input$map_draw_deleted_features
      if (!is.null(deleted) && length(deleted$features) > 0) {
        drawn_feats(list())
        isolate({
          recompute_selection()
        })
      }
    })

    observe({
      sel_spots <- selected_spots()
      proxy <- leafletProxy("map") %>% clearGroup("selected")

      if (length(sel_spots) > 0) {
        sel_indices <- which(spots_sf$spot_id %in% sel_spots)
        proxy <- proxy %>%
          addCircleMarkers(
            lng = spots_sf$x[sel_indices],
            lat = spots_sf$y[sel_indices],
            radius = 4, stroke = TRUE, color = "red", weight = 2,
            fillColor = "yellow", fillOpacity = 1, group = "selected"
          )
      }
    })

    # Group management (from map buttons)
    observeEvent(input$save_group1_map, {
      sel <- selected_spots()
      if (length(sel) > 0) {
        group1_spots(sel)
        output_group1(sel)  # Save to export variable
        showNotification(paste("Saved", length(sel), "spots to Group 1"), type = "message")

        # Clear selection after saving
        drawn_feats(list())
        selected_spots(character(0))
        session$sendCustomMessage("clearFreehandDrawings", list())
        leafletProxy("map") %>%
          clearGroup("drawn") %>%
          clearGroup("selected")

        # Immediately show the group on map if checkbox is checked
        if (input$show_groups) {
          Sys.sleep(0.1)  # Small delay to ensure state is updated
          g1_indices <- which(spots_sf$spot_id %in% sel)
          leafletProxy("map") %>%
            clearGroup("group1_display") %>%
            addCircleMarkers(
              lng = spots_sf$x[g1_indices],
              lat = spots_sf$y[g1_indices],
              radius = 4, stroke = TRUE, color = "darkred", weight = 2,
              fillColor = "red", fillOpacity = 0.8, group = "group1_display"
            )
        }
      } else {
        showNotification("No spots selected!", type = "warning")
      }
    })

    observeEvent(input$save_group2_map, {
      sel <- selected_spots()
      if (length(sel) > 0) {
        group2_spots(sel)
        output_group2(sel)  # Save to export variable
        showNotification(paste("Saved", length(sel), "spots to Group 2"), type = "message")

        # Clear selection after saving
        drawn_feats(list())
        selected_spots(character(0))
        session$sendCustomMessage("clearFreehandDrawings", list())
        leafletProxy("map") %>%
          clearGroup("drawn") %>%
          clearGroup("selected")

        # Immediately show the group on map if checkbox is checked
        if (input$show_groups) {
          Sys.sleep(0.1)  # Small delay to ensure state is updated
          g2_indices <- which(spots_sf$spot_id %in% sel)
          leafletProxy("map") %>%
            clearGroup("group2_display") %>%
            addCircleMarkers(
              lng = spots_sf$x[g2_indices],
              lat = spots_sf$y[g2_indices],
              radius = 4, stroke = TRUE, color = "darkblue", weight = 2,
              fillColor = "blue", fillOpacity = 0.8, group = "group2_display"
            )
        }
      } else {
        showNotification("No spots selected!", type = "warning")
      }
    })

    output$group1_info <- renderText({
      g1 <- group1_spots()
      if (length(g1) == 0) "Group 1: Empty" else paste("Group 1:", length(g1), "spots")
    })

    output$group2_info <- renderText({
      g2 <- group2_spots()
      if (length(g2) == 0) "Group 2: Empty" else paste("Group 2:", length(g2), "spots")
    })

    observe({
      if (input$show_groups) {
        g1 <- group1_spots()
        g2 <- group2_spots()

        proxy <- leafletProxy("map") %>%
          clearGroup("group1_display") %>%
          clearGroup("group2_display")

        if (length(g1) > 0) {
          g1_indices <- which(spots_sf$spot_id %in% g1)
          proxy <- proxy %>%
            addCircleMarkers(
              lng = spots_sf$x[g1_indices],
              lat = spots_sf$y[g1_indices],
              radius = 4, stroke = TRUE, color = "darkred", weight = 2,
              fillColor = "red", fillOpacity = 0.8, group = "group1_display"
            )
        }

        if (length(g2) > 0) {
          g2_indices <- which(spots_sf$spot_id %in% g2)
          proxy <- proxy %>%
            addCircleMarkers(
              lng = spots_sf$x[g2_indices],
              lat = spots_sf$y[g2_indices],
              radius = 4, stroke = TRUE, color = "darkblue", weight = 2,
              fillColor = "blue", fillOpacity = 0.8, group = "group2_display"
            )
        }
      } else {
        leafletProxy("map") %>%
          clearGroup("group1_display") %>%
          clearGroup("group2_display")
      }
    }, priority = 1000)


    ###1
    # DEG Analysis
    observeEvent(input$run_deg, {
      g1 <- group1_spots()
      g2 <- group2_spots()

      # Check if we're doing cluster-based analysis
      is_cluster_analysis <- grepl("Clusters in Group", input$deg_comparison)

      if (!is_cluster_analysis) {
        # Original group-based analysis
        if (input$deg_comparison == "Group 1 vs Rest" && length(g1) == 0) {
          showNotification("Group 1 is empty!", type = "error")
          return()
        }

        if (input$deg_comparison == "Group 2 vs Rest" && length(g2) == 0) {
          showNotification("Group 2 is empty!", type = "error")
          return()
        }

        if (input$deg_comparison == "Group 1 vs Group 2" && (length(g1) == 0 || length(g2) == 0)) {
          showNotification("Both groups must have spots!", type = "error")
          return()
        }

        showNotification("Running DEG analysis...", type = "message", duration = NULL, id = "deg_running")

        tryCatch({
          temp_idents <- rep("Other", ncol(seurat_obj))
          names(temp_idents) <- colnames(seurat_obj)
          temp_idents[g1] <- "Group1"
          temp_idents[g2] <- "Group2"

          if (input$deg_comparison == "Group 1 vs Group 2") {
            cells_to_use <- c(g1, g2)
          } else {
            cells_to_use <- colnames(seurat_obj)
          }

          temp_seurat <- subset(seurat_obj, cells = cells_to_use)
          temp_seurat$temp_ident <- temp_idents[cells_to_use]
          Idents(temp_seurat) <- "temp_ident"

          if (input$deg_comparison == "Group 1 vs Group 2") {
            markers <- Seurat::FindMarkers(temp_seurat, ident.1 = "Group1", ident.2 = "Group2",
                                   verbose = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
          } else if (input$deg_comparison == "Group 1 vs Rest") {
            markers <- Seurat::FindMarkers(temp_seurat, ident.1 = "Group1", ident.2 = "Other",
                                   verbose = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
          } else if (input$deg_comparison == "Group 2 vs Rest") {
            markers <- Seurat::FindMarkers(temp_seurat, ident.1 = "Group2", ident.2 = "Other",
                                   verbose = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
          }

          markers$gene <- rownames(markers)
          markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]

          deg_results(markers)

          removeNotification(id = "deg_running")
          showNotification(paste("Found", nrow(markers), "DEGs"), type = "message")

        }, error = function(e) {
          removeNotification(id = "deg_running")
          showNotification(paste("Error:", e$message), type = "error", duration = 10)
        })

      } else {
        # Cluster-based analysis
        cluster_res <- cluster_results()

        if (is.null(cluster_res)) {
          showNotification("Please run clustering first!", type = "error")
          return()
        }

        # Determine which group to analyze
        target_group <- if (input$deg_comparison == "Clusters in Group 1 (one vs rest)") {
          "group1"
        } else {
          "group2"
        }

        # Check if clustering was done on the target group
        if (cluster_res$spot_selection != target_group) {
          group_name <- if (target_group == "group1") "Group 1" else "Group 2"
          showNotification(paste("Clustering was not performed on", group_name, "spots!"),
                           type = "error", duration = 8)
          return()
        }

        target_spots <- if (target_group == "group1") g1 else g2

        if (length(target_spots) == 0) {
          group_name <- if (target_group == "group1") "Group 1" else "Group 2"
          showNotification(paste(group_name, "is empty!"), type = "error")
          return()
        }

        showNotification("Running cluster marker analysis...", type = "message", duration = NULL, id = "deg_running")

        tryCatch({
          # Get clusters for the target spots
          clusters <- cluster_res$clusters

          # Subset seurat object to only clustered spots
          temp_seurat <- subset(seurat_obj, cells = names(clusters))
          temp_seurat$cluster_ident <- as.factor(clusters)
          Idents(temp_seurat) <- "cluster_ident"

          # Find markers for all clusters (one vs rest)
          all_markers <- FindAllMarkers(temp_seurat,
                                        only.pos = FALSE,
                                        min.pct = 0.1,
                                        logfc.threshold = 0.25,
                                        verbose = FALSE)

          # Add cluster column and reorganize
          all_markers$gene <- rownames(all_markers)
          all_markers <- all_markers[order(all_markers$cluster, all_markers$p_val_adj, -abs(all_markers$avg_log2FC)), ]

          # Reorder columns for better readability
          all_markers <- all_markers[, c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]

          deg_results(all_markers)

          removeNotification(id = "deg_running")
          showNotification(paste("Found markers for", cluster_res$n_clusters, "clusters"),
                           type = "message", duration = 5)

        }, error = function(e) {
          removeNotification(id = "deg_running")
          showNotification(paste("Error:", e$message), type = "error", duration = 10)
        })
      }
    })

    output$deg_results_available <- reactive({
      !is.null(deg_results())
    })
    outputOptions(output, "deg_results_available", suspendWhenHidden = FALSE)

    output$deg_table <- renderTable({
      deg <- deg_results()
      if (!is.null(deg)) {
        top_deg <- head(deg, 50)  # Show more results for cluster analysis

        # Check if this is cluster-based analysis (has 'cluster' column)
        if ("cluster" %in% colnames(top_deg)) {
          data.frame(
            Cluster = top_deg$cluster,
            Gene = top_deg$gene,
            log2FC = round(top_deg$avg_log2FC, 3),
            pct.1 = round(top_deg$pct.1, 3),
            pct.2 = round(top_deg$pct.2, 3),
            p_adj = format(top_deg$p_val_adj, scientific = TRUE, digits = 3)
          )
        } else {
          data.frame(
            Gene = top_deg$gene,
            log2FC = round(top_deg$avg_log2FC, 3),
            p_adj = format(top_deg$p_val_adj, scientific = TRUE, digits = 3)
          )
        }
      }
    }, rownames = FALSE)


    # Violin plot
    observeEvent(input$plot_violin, {
      if (input$violin_feature_type == "gene") {
        feature <- input$violin_gene
      } else if (input$violin_feature_type == "metadata") {
        feature <- input$violin_metadata
      } else {
        feature <- input$violin_geneset
      }

      if (is.null(feature) || feature == "") {
        showNotification("Select a feature!", type = "warning")
        return()
      }

      g1 <- group1_spots()
      g2 <- group2_spots()

      # Validate based on comparison type
      if (input$violin_comparison == "Group 1 vs Rest" && length(g1) == 0) {
        showNotification("Group 1 is empty!", type = "error")
        return()
      }

      if (input$violin_comparison == "Group 2 vs Rest" && length(g2) == 0) {
        showNotification("Group 2 is empty!", type = "error")
        return()
      }

      if (input$violin_comparison == "Group 1 vs Group 2" && (length(g1) == 0 || length(g2) == 0)) {
        showNotification("Both groups must have spots!", type = "error")
        return()
      }

      tryCatch({
        if (input$violin_feature_type == "gene") {
          all_values <- FetchData(seurat_obj, vars = feature, slot = "data")
          is_numeric_data <- TRUE
        } else if (input$violin_feature_type == "metadata") {
          meta_col <- seurat_obj@meta.data[[feature]]
          all_values <- data.frame(value = meta_col)
          colnames(all_values) <- feature
          rownames(all_values) <- rownames(seurat_obj@meta.data)
          is_numeric_data <- is.numeric(meta_col)
        } else {
          scores <- gene_set_scores()[[feature]]
          all_values <- data.frame(value = scores)
          colnames(all_values) <- feature
          rownames(all_values) <- names(scores)
          is_numeric_data <- TRUE
        }

        if (input$violin_comparison == "Group 1 vs Group 2") {
          cells_to_plot <- c(g1, g2)
          group_label <- c(rep("Group 1", length(g1)), rep("Group 2", length(g2)))
          names(group_label) <- cells_to_plot
        } else if (input$violin_comparison == "Group 1 vs Rest") {
          cells_to_plot <- spots_sf$spot_id
          group_label <- rep("Rest", length(cells_to_plot))
          names(group_label) <- cells_to_plot
          group_label[g1] <- "Group 1"
        } else if (input$violin_comparison == "Group 2 vs Rest") {
          cells_to_plot <- spots_sf$spot_id
          group_label <- rep("Rest", length(cells_to_plot))
          names(group_label) <- cells_to_plot
          group_label[g2] <- "Group 2"
        }

        valid_cells <- intersect(cells_to_plot, rownames(all_values))

        plot_data <- data.frame(
          value = all_values[valid_cells, 1],
          group = factor(group_label[valid_cells])
        )

        plot_data$is_numeric <- is_numeric_data
        # Calculate statistics for group comparison
        pvalue_test <- NA

        if (is_numeric_data) {
          group_levels <- unique(plot_data$group)

          if (length(group_levels) == 2) {
            group1_values <- plot_data$value[plot_data$group == group_levels[1]]
            group2_values <- plot_data$value[plot_data$group == group_levels[2]]

            group1_values <- group1_values[!is.na(group1_values)]
            group2_values <- group2_values[!is.na(group2_values)]

            if (length(group1_values) >= 3 && length(group2_values) >= 3) {
              if (input$violin_stat_test == "ttest") {
                test_result <- t.test(group1_values, group2_values)
              } else {
                test_result <- wilcox.test(group1_values, group2_values)
              }
              pvalue_test <- test_result$p.value
            }
          }
        }

        attr(plot_data, "pvalue_test") <- pvalue_test
        attr(plot_data, "stat_test") <- input$violin_stat_test

        violin_data(plot_data)

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })

    output$violin_available <- reactive({
      !is.null(violin_data())
    })
    outputOptions(output, "violin_available", suspendWhenHidden = FALSE)

    output$violin_plot <- renderPlot({
      plot_data <- violin_data()
      if (!is.null(plot_data) && plot_data$is_numeric[1]) {
        # Get comparison type for title
        comparison_text <- if(input$violin_comparison == "Group 1 vs Group 2") {
          "Group 1 vs Group 2"
        } else if(input$violin_comparison == "Group 2 vs Rest") {
          "Group 2 vs Rest"
        } else {
          "Group 1 vs Rest"
        }

        # Calculate n for each group
        group_counts <- table(plot_data$group)
        n_text <- paste(names(group_counts), "n =", group_counts, collapse = ", ")

        # Extract statistical test information
        pvalue_test <- attr(plot_data, "pvalue_test")
        stat_test <- attr(plot_data, "stat_test")

        # Create subtitle with test results
        subtitle_text <- if (!is.na(pvalue_test)) {
          test_name <- if(stat_test == "ttest") "t-test" else "Wilcoxon test"
          paste0(n_text, "\n", test_name, " p = ", format(pvalue_test, scientific = TRUE, digits = 3))
        } else {
          n_text
        }

        ggplot(plot_data, aes(x = group, y = value, fill = group)) +
          geom_violin(trim = FALSE, alpha = 0.5) +
          geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
          labs(title = comparison_text,
               subtitle = subtitle_text,
               x = "", y = "Value") +
          theme_classic(base_size = 12) +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5)) +
          scale_fill_manual(values = c("Group 1" = "#E41A1C", "Group 2" = "#377EB8", "Rest" = "#999999"))
      }
    })

    # Feature comparison
    observeEvent(input$plot_compare, {
      # Get feature 1
      if (input$compare_type1 == "gene") {
        feature1 <- input$compare_gene1
      } else if (input$compare_type1 == "metadata") {
        feature1 <- input$compare_meta1
      } else {
        feature1 <- input$compare_geneset1
      }

      # Get feature 2
      if (input$compare_type2 == "gene") {
        feature2 <- input$compare_gene2
      } else if (input$compare_type2 == "metadata") {
        feature2 <- input$compare_meta2
      } else {
        feature2 <- input$compare_geneset2
      }

      if (is.null(feature1) || feature1 == "" || is.null(feature2) || feature2 == "") {
        showNotification("Please select both features to compare!", type = "warning")
        return()
      }

      if (input$compare_type1 == "geneset") {
        if (is.null(gene_set_scores()[[feature1]])) {
          showNotification("Gene set for Feature 1 not found. Please calculate it first!", type = "error")
          return()
        }
      }

      if (input$compare_type2 == "geneset") {
        if (is.null(gene_set_scores()[[feature2]])) {
          showNotification("Gene set for Feature 2 not found. Please calculate it first!", type = "error")
          return()
        }
      }

      # Use spots based on selection
      if (input$compare_spots_selection == "group1") {
        sel_spots <- group1_spots()
        if (length(sel_spots) == 0) {
          showNotification("No spots in Group 1! Please select spots for Group 1 first.",
                           type = "warning", duration = 5)
          return()
        }
      } else if (input$compare_spots_selection == "group2") {
        sel_spots <- group2_spots()
        if (length(sel_spots) == 0) {
          showNotification("No spots in Group 2! Please select spots for Group 2 first.",
                           type = "warning", duration = 5)
          return()
        }
      } else {
        # Use all spots
        sel_spots <- spots_sf$spot_id
      }

      showNotification(paste("Using", length(sel_spots), "spots from",
                             switch(input$compare_spots_selection,
                                    "all" = "all spots",
                                    "group1" = "Group 1",
                                    "group2" = "Group 2")),
                       type = "message", duration = 3)

      tryCatch({
        # Get ALL values for feature 1
        if (input$compare_type1 == "gene") {
          all_values1 <- FetchData(seurat_obj, vars = feature1, slot = "data")[, 1]
          if (is.null(names(all_values1))) {
            names(all_values1) <- colnames(seurat_obj)
          }
          is_numeric1 <- TRUE
        } else if (input$compare_type1 == "metadata") {
          all_values1 <- seurat_obj@meta.data[[feature1]]
          names(all_values1) <- rownames(seurat_obj@meta.data)
          is_numeric1 <- is.numeric(all_values1)
        } else {
          all_values1 <- gene_set_scores()[[feature1]]
          if (is.null(names(all_values1))) {
            names(all_values1) <- colnames(seurat_obj)
          }
          is_numeric1 <- TRUE
        }

        # Get ALL values for feature 2
        if (input$compare_type2 == "gene") {
          all_values2 <- FetchData(seurat_obj, vars = feature2, slot = "data")[, 1]
          if (is.null(names(all_values2))) {
            names(all_values2) <- colnames(seurat_obj)
          }
          is_numeric2 <- TRUE
        } else if (input$compare_type2 == "metadata") {
          all_values2 <- seurat_obj@meta.data[[feature2]]
          names(all_values2) <- rownames(seurat_obj@meta.data)
          is_numeric2 <- is.numeric(all_values2)
        } else {
          all_values2 <- gene_set_scores()[[feature2]]
          if (is.null(names(all_values2))) {
            names(all_values2) <- colnames(seurat_obj)
          }
          is_numeric2 <- TRUE
        }

        # Filter for selected spots
        valid_spots <- intersect(sel_spots, names(all_values1))
        valid_spots <- intersect(valid_spots, names(all_values2))

        if (length(valid_spots) == 0) {
          valid_spots_sf <- intersect(sel_spots, spots_sf$spot_id)

          if (length(valid_spots_sf) > 0) {
            valid_spots <- intersect(valid_spots_sf, names(all_values1))
            valid_spots <- intersect(valid_spots, names(all_values2))
          }
        }

        if (length(valid_spots) == 0) {
          showNotification(paste("No valid spots found in selection!",
                                 "\nSelected:", length(sel_spots), "spots"),
                           type = "error", duration = 10)
          return()
        }

        values1_sel <- all_values1[valid_spots]
        values2_sel <- all_values2[valid_spots]

        # CRITICAL FIX: Force numeric conversion for proper scatter plot display
        # This handles cases where metadata might be stored as factors or characters
        if (!is.numeric(values1_sel)) {
          numeric_attempt1 <- suppressWarnings(as.numeric(as.character(values1_sel)))
          if (!all(is.na(numeric_attempt1))) {
            values1_sel <- numeric_attempt1
            is_numeric1 <- TRUE
          }
        }

        if (!is.numeric(values2_sel)) {
          numeric_attempt2 <- suppressWarnings(as.numeric(as.character(values2_sel)))
          if (!all(is.na(numeric_attempt2))) {
            values2_sel <- numeric_attempt2
            is_numeric2 <- TRUE
          }
        }

        if (length(values1_sel) != length(values2_sel)) {
          showNotification("Feature values have different lengths!", type = "error")
          return()
        }

        # Calculate statistics if both are numeric
        cor_result <- NULL
        pvalue_cor <- NA
        pvalue_comp <- NA

        if (is_numeric1 && is_numeric2) {
          complete_idx <- which(!is.na(values1_sel) & !is.na(values2_sel))

          if (length(complete_idx) >= 3) {
            cor_result <- cor.test(values1_sel[complete_idx], values2_sel[complete_idx],
                                   method = input$cor_method)
            pvalue_cor <- cor_result$p.value

            if (input$stat_test == "ttest") {
              comp_test <- t.test(values1_sel[complete_idx], values2_sel[complete_idx],
                                  paired = TRUE)
            } else {
              comp_test <- wilcox.test(values1_sel[complete_idx], values2_sel[complete_idx],
                                       paired = TRUE)
            }
            pvalue_comp <- comp_test$p.value
          }
        }

        # Create combined dataframe for plotting
        plot_data <- data.frame(
          value = c(values1_sel, values2_sel),
          feature = factor(rep(c(feature1, feature2), each = length(valid_spots)),
                           levels = c(feature1, feature2)),
          stringsAsFactors = FALSE
        )

        # Store additional data as attributes
        attr(plot_data, "pvalue_cor") <- pvalue_cor
        attr(plot_data, "pvalue_comp") <- pvalue_comp
        attr(plot_data, "is_numeric1") <- is_numeric1
        attr(plot_data, "is_numeric2") <- is_numeric2
        #attr(plot_data, "correlation") <- if(!is.null(cor_result)) cor_result$estimate else NA
        attr(plot_data, "feature1") <- feature1
        attr(plot_data, "feature2") <- feature2
        attr(plot_data, "values1") <- values1_sel
        attr(plot_data, "values2") <- values2_sel
        attr(plot_data, "cor_method") <- input$cor_method
        attr(plot_data, "stat_test") <- input$stat_test
        attr(plot_data, "n_spots") <- length(valid_spots)

        compare_data(plot_data)

        showNotification(paste("Comparison complete for", length(valid_spots), "spots"),
                         type = "message")

      }, error = function(e) {
        showNotification(paste("Error comparing features:", e$message),
                         type = "error", duration = 10)
      })
    })

    output$compare_available <- reactive({
      !is.null(compare_data())
    })
    outputOptions(output, "compare_available", suspendWhenHidden = FALSE)

    output$compare_plot <- renderPlot({
      plot_data <- compare_data()
      if (!is.null(plot_data)) {

        # Extract attributes
        feature1 <- attr(plot_data, "feature1")
        feature2 <- attr(plot_data, "feature2")
        is_numeric1 <- attr(plot_data, "is_numeric1")
        is_numeric2 <- attr(plot_data, "is_numeric2")
        pvalue_comp <- attr(plot_data, "pvalue_comp")
        pvalue_cor <- attr(plot_data, "pvalue_cor")
        #correlation <- attr(plot_data, "correlation")
        values1 <- attr(plot_data, "values1")
        values2 <- attr(plot_data, "values2")
        cor_method <- attr(plot_data, "cor_method")
        stat_test <- attr(plot_data, "stat_test")
        n_spots <- attr(plot_data, "n_spots")

        if (is_numeric1 && is_numeric2) {
          # Violin plot - Two features comparison (use teal/cyan color scheme)
          p1 <- ggplot(plot_data, aes(x = feature, y = value, fill = feature)) +
            geom_violin(trim = FALSE, alpha = 0.5) +
            geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
            labs(title = "Two Features Comparison",
                 subtitle = if(!is.na(pvalue_comp)) {
                   test_name <- if(stat_test == "ttest") "Paired t-test" else "Wilcoxon test"
                   paste0(test_name, " p = ", format(pvalue_comp, scientific = TRUE, digits = 3),
                          " (n = ", n_spots, " spots)")
                 } else {
                   paste0("n = ", n_spots, " spots")
                 },
                 x = "", y = "Value") +
            theme_classic(base_size = 12) +
            theme(legend.position = "none",
                  plot.title = element_text(hjust = 0.5, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5)) +
            scale_fill_manual(values = c("#17BECF", "#1F77B4"))  # Teal and blue

          print(p1)

        } else {
          # Handle non-numeric comparisons
          p <- ggplot(plot_data, aes(x = feature, y = value, fill = feature)) +
            geom_boxplot() +
            labs(title = "Feature Comparison", x = "", y = "Value") +
            theme_classic() +
            theme(legend.position = "none")

          print(p)
        }
      }
    })

    # # H&E image controls
    # observeEvent(input$show_he, {
    #   if (!is.null(he_image_base64)) {
    #     session$sendCustomMessage("toggleHEImage", list(show = input$show_he))
    #   }
    # })

    observeEvent(input$image_opacity, {
      if (!is.null(he_image_base64)) {
        session$sendCustomMessage("updateHEOpacity", list(opacity = input$image_opacity))
      }
    })

    # Download handlers for groups
    output$dl_group1 <- downloadHandler(
      filename = function() {
        paste0(current_sample_name(), "_group1_spots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        g1 <- isolate(output_group1())
        print(paste("Downloading", length(g1), "Group 1 spots"))

        if (length(g1) > 0) {
          writeLines(g1, file)
        } else {
          writeLines("No spots in Group 1", file)
        }
      },
      contentType = "text/plain"
    )

    output$dl_group2 <- downloadHandler(
      filename = function() {
        paste0(current_sample_name(), "_group2_spots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
      },
      content = function(file) {
        g2 <- isolate(output_group2())
        print(paste("Downloading", length(g2), "Group 2 spots"))

        if (length(g2) > 0) {
          writeLines(g2, file)
        } else {
          writeLines("No spots in Group 2", file)
        }
      },
      contentType = "text/plain"
    )

    output$dl_seurat_group1 <- downloadHandler(
      filename = function() {
        paste0(current_sample_name(), "_group1_subset_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        g1 <- isolate(output_group1())
        print(paste("Downloading Group 1 Seurat subset with", length(g1), "spots"))

        if (length(g1) > 0) {
          subset_obj <- subset(seurat_obj, cells = g1)
          saveRDS(subset_obj, file)
        } else {
          showNotification("No spots in Group 1!", type = "warning", duration = 5)
        }
      },
      contentType = "application/octet-stream"
    )

    output$dl_seurat_group2 <- downloadHandler(
      filename = function() {
        paste0(current_sample_name(), "_group2_subset_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      },
      content = function(file) {
        g2 <- isolate(output_group2())
        print(paste("Downloading Group 2 Seurat subset with", length(g2), "spots"))

        if (length(g2) > 0) {
          subset_obj <- subset(seurat_obj, cells = g2)
          saveRDS(subset_obj, file)
        } else {
          showNotification("No spots in Group 2!", type = "warning", duration = 5)
        }
      },
      contentType = "application/octet-stream"
    )

    output$dl_deg <- downloadHandler(
      filename = function() paste0(current_sample_name(), "_DEGs.csv"),
      content = function(file) {
        if (!is.null(deg_results())) {
          write.csv(deg_results(), file, row.names = FALSE)
        }
      },
      contentType = "text/csv"
    )
  }

  shinyApp(ui, server)
}
