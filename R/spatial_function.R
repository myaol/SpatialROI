future::plan("sequential")
options(future.globals.maxSize = 2000 * 1024^2)

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
#' @importFrom DT datatable
#' @importFrom dplyr %>% filter
#' @importFrom grDevices as.raster colorRampPalette dev.off png rainbow
#' @importFrom graphics par
#' @importFrom methods slot
#' @importFrom stats cor.test setNames t.test wilcox.test
#' @importFrom utils capture.output head write.csv
#'
#' @return A Shiny application object; launching the GUI as a side effect.
#' @export



#'
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
  hallmark_library_human <- prepared_data$hallmark_library_human
  hallmark_library_mouse <- prepared_data$hallmark_library_mouse

  cellmarker_db <- lapply(signature_library_human, function(x) x$genes)
  cellmarker_db_mouse <- lapply(signature_library_mouse, function(x) x$genes)

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

# ── ADD THESE DEBUG LINES ──
      cat("Image class:", class(image_obj)[1], "\n")
      cat("Slot names:", paste(slotNames(image_obj), collapse=", "), "\n")
      tryCatch(cat("@image dim:", paste(dim(image_obj@image), collapse="x"), "\n"),
               error = function(e) cat("@image failed:", e$message, "\n"))
      tryCatch({
        img <- GetImage(seurat_obj, image = image_name, mode = "raster")
        cat("GetImage worked, dim:", paste(dim(img), collapse="x"), "\n")
      }, error = function(e) cat("GetImage failed:", e$message, "\n"))
      # ── END DEBUG ──


      he_image_data <- image_obj@image
      scale_factor  <- image_obj@scale.factors$lowres
      H <- dim(he_image_data)[1]
      W <- dim(he_image_data)[2]

      coords_full <- GetTissueCoordinates(seurat_obj, image = image_name)

      # Always produce lowres-pixel coordinates
      if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
        pixel_x <- coords_full$pxl_col_in_fullres * scale_factor
        pixel_y <- coords_full$pxl_row_in_fullres  * scale_factor
      } else if ("imagecol" %in% colnames(coords_full)) {
        pixel_x <- coords_full$imagecol
        pixel_y <- coords_full$imagerow
      } else {
        # spots_sf$y is already flipped — un-flip before using as pixel row
        pixel_x <- spots_sf$x
        pixel_y <- max(spots_sf$y) + min(spots_sf$y) - spots_sf$y
      }

      pixel_x_min <- min(pixel_x, na.rm = TRUE); pixel_x_max <- max(pixel_x, na.rm = TRUE)
      pixel_y_min <- min(pixel_y, na.rm = TRUE); pixel_y_max <- max(pixel_y, na.rm = TRUE)
      x_buffer_px <- (pixel_x_max - pixel_x_min) * 0.1
      y_buffer_px <- (pixel_y_max - pixel_y_min) * 0.1

      crop_x_min <- max(1, floor(pixel_x_min - x_buffer_px))
      crop_x_max <- min(W, ceiling(pixel_x_max + x_buffer_px))
      crop_y_min <- max(1, floor(pixel_y_min - y_buffer_px))
      crop_y_max <- min(H, ceiling(pixel_y_max + y_buffer_px))

      # Guard against invalid ranges
      if (crop_x_min >= crop_x_max) { crop_x_min <- 1; crop_x_max <- W }
      if (crop_y_min >= crop_y_max) { crop_y_min <- 1; crop_y_max <- H }

      he_image_cropped <- he_image_data[crop_y_min:crop_y_max, crop_x_min:crop_x_max, ]
      temp_file <- tempfile(fileext = ".png")
      png(temp_file, width = dim(he_image_cropped)[2], height = dim(he_image_cropped)[1])
      par(mar = c(0, 0, 0, 0))
      plot(as.raster(he_image_cropped), axes = FALSE)
      dev.off()
      he_image_base64 <- paste0("data:image/png;base64,", base64enc::base64encode(temp_file))
      unlink(temp_file)

      # Bounds in flipped leaflet coordinate system
      y_sum <- max(spots_sf$y) + min(spots_sf$y)
      he_image_bounds <- list(
        north = y_sum - crop_y_min,
        south = y_sum - crop_y_max,
        west  = crop_x_min,
        east  = crop_x_max
      )
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
          z-index: 1100;
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
        .leaflet-bottom.leaflet-right .leaflet-control {
          margin-bottom: 70px !important;
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
                        # ── Two-tier ROI / Group panel (Reviewer 1, item 1) ────────────
                        # ① Draw a region → name it → Save region (named ROI).
                        # ② Optionally group several named ROIs for group-vs-group.
                        tags$div(style = "background:white; padding:10px 14px; border-radius:10px; box-shadow:0 3px 14px rgba(0,0,0,0.25); display:flex; flex-direction:row; gap:16px; align-items:flex-start;",
                                 # Column ① — draw & save a named region
                                 tags$div(style = "display:flex; flex-direction:column; gap:6px; width:270px;",
                                   tags$div(style = "font-weight:700; font-size:13px; color:#2c3e50;", "① Draw a region"),
                                   div(style = "width:100%;",
                                       textInput("roi_name", NULL, value = "ROI 1",
                                                 placeholder = "e.g. Tumor edge", width = "100%")),
                                   actionButton("save_roi_btn", "➕ Save region",
                                                class = "btn btn-primary", style = "width:100%; font-weight:700;"),
                                   uiOutput("roi_chips")
                                 ),
                                 tags$div(style = "width:1px; align-self:stretch; background:#e5e7eb;"),
                                 # Column ② — group ROIs for comparison
                                 tags$div(style = "display:flex; flex-direction:column; gap:6px; width:270px;",
                                   tags$div(style = "font-weight:700; font-size:13px; color:#2c3e50;", "② Group ROIs"),
                                   div(style = "width:100%;",
                                       selectizeInput("group_member_rois", NULL, choices = NULL, multiple = TRUE,
                                                      options = list(placeholder = "pick ROIs", dropdownParent = "body"),
                                                      width = "100%")),
                                   div(style = "display:flex; gap:6px;",
                                       div(style = "flex:1;",
                                           textInput("group_name", NULL, placeholder = "group name", width = "100%")),
                                       actionButton("create_group_btn", "＋ Group",
                                                    class = "btn btn-success", style = "white-space:nowrap; font-weight:700;")),
                                   uiOutput("group_chips")
                                 )
                        ),
                        tags$div(style = "display: flex; flex-direction: column; gap: 6px; min-width:180px;",
                          # Focus selector sits ABOVE the checkbox; includes ROIs AND groups.
                          div(style = "font-size:12px;",
                              selectizeInput("roi_show_filter", "Show on map:",
                                             choices = c("All ROIs" = "__all__"),
                                             selected = "__all__",
                                             options = list(dropdownParent = "body"), width = "100%")),
                          checkboxInput("show_groups", "Show ROIs on Map", value = FALSE),
                          checkboxInput("transparent_groups", "Transparent ROI Display", value = FALSE),
                          # Reviewer 1, item 3: high-contrast outline of each saved ROI so the
                          # boundary stays visible even over bright/high expression. Fixed
                          # white-cased black dashed style (no color/width options, by design).
                          checkboxInput("show_roi_contours", "Show ROI contours", value = FALSE)
                        )
                      )
             ),

             # Right control panel
             tags$div(class = "control-panel", id = "control_panel",
                      tags$button(class = "close-panel-btn", onclick = "Shiny.setInputValue('close_panel', Math.random());", "×"),

                      # Home content - FULL SCREEN WITHIN APP
                      tags$div(class = "control-content full-screen-content", id = "content_home",
                               
                               tags$div(style = "width: 100%; max-width: 100%; height: 100vh; overflow-y: auto; background: #f5f7fa; padding: 0px 8px 40px 8px;",
                                tags$div(style = "max-width: 100%; margin: 0 auto; padding: 20px 16px;",

                                  # Hero

                                  tags$div(style = "background: linear-gradient(135deg, #0072B5 0%, #E18727 100%); color: white; padding: 28px 24px; border-radius: 12px; margin-bottom: 24px;",
                                    tags$h1(style = "margin: 0 0 10px 0; font-size: 32px; font-weight: bold;", "🔬 SpatialScope"),
                                    tags$p(style = "margin: 0 0 8px 0; font-size: 17px; line-height: 1.7; opacity: 0.97;",
                                      "An interactive platform for spatial transcriptomics analysis. ",
                                      tags$strong("Draw freehand regions of interest"), 
                                      " directly on your tissue image, compare spatial domains, and perform statistical analyses — with no coding required."),
                                    tags$p(style = "margin: 0; font-size: 15px; opacity: 0.85;",
                                      "✅ Accepts: Seurat object (.rds); 10x Visium SpaceRanger raw output")
                                  ),

                                  # Quick Start - vertical rows
                                  tags$div(style = "background: white; border-radius: 12px; padding: 24px 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 20px;",
                                    tags$h2(style = "color: #0072B5; font-size: 22px; margin: 0 0 20px 0; padding-bottom: 8px; border-bottom: 2px solid #E18727;", "📖 Quick Start"),
                                    
                                    tags$div(style = "display: flex; flex-direction: column; gap: 14px;",
                                      
                                      tags$div(style = "display: flex; align-items: flex-start; gap: 18px;",
                                        tags$div(style = "background: #0072B5; color: white; min-width: 32px; height: 32px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 15px; font-weight: bold;", "1"),
                                        tags$div(
                                          tags$div(style = "font-weight: bold; color: #2c3e50; font-size: 16px; margin-bottom: 3px;", "Upload Your Data"),
                                          tags$div(style = "font-size: 14px; color: #555; line-height: 1.5;", "Navigate to the 🎨 Visualization panel to load a Seurat .rds file, a 10x SpaceRanger output folder, or use the built-in example dataset to get started immediately.")
                                        )
                                      ),

                                      tags$div(style = "display: flex; align-items: flex-start; gap: 18px;",
                                        tags$div(style = "background: #0072B5; color: white; min-width: 32px; height: 32px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 15px; font-weight: bold;", "2"),
                                        tags$div(
                                          tags$div(style = "font-weight: bold; color: #2c3e50; font-size: 16px; margin-bottom: 3px;", "Select Regions of Interest"),
                                          tags$div(style = "font-size: 14px; color: #555; line-height: 1.5;", "Use the ✏️ freehand drawing tool on the map to draw around any tissue regions you want to analyze, then assign the selected spots to Group 1 or Group 2 for comparison.")
                                        )
                                      ),

                                      tags$div(style = "display: flex; align-items: flex-start; gap: 18px;",
                                        tags$div(style = "background: #0072B5; color: white; min-width: 32px; height: 32px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 15px; font-weight: bold;", "3"),
                                        tags$div(
                                          tags$div(style = "font-weight: bold; color: #2c3e50; font-size: 16px; margin-bottom: 3px;", "Save Groups for Comparison"),
                                          tags$div(style = "font-size: 14px; color: #555; line-height: 1.5;", 
                                            " Selected groups remain available throughout the session. Use ", 
                                            tags$strong("⬇ Download"), " beneath each button to export the corresponding spots as a new Seurat object for reuse or downstream analysis.")
                                        )
                                      ),

                                      tags$div(style = "display: flex; align-items: flex-start; gap: 18px;",
                                        tags$div(style = "background: #0072B5; color: white; min-width: 32px; height: 32px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-size: 15px; font-weight: bold;", "4"),
                                        tags$div(
                                          tags$div(style = "font-weight: bold; color: #2c3e50; font-size: 16px; margin-bottom: 3px;", "Explore & Analyze"),
                                          tags$div(style = "font-size: 14px; color: #555; line-height: 1.5;", "Use the sidebar tools to visualize gene expression, score gene signatures, perform clustering, differential expression analysis, ligand–receptor colocalization, cell-type deconvolution, and export results.")
                                        )
                                      )
                                    )
                                  ),

                                  # Tips - bullet list style
                                  tags$div(style = "background: white; border-radius: 12px; padding: 24px 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin-bottom: 20px;",
                                    tags$h2(style = "color: #0072B5; font-size: 22px; margin: 0 0 15px 0; padding-bottom: 8px; border-bottom: 2px solid #E18727;", "💡 Tips & Best Practices"),
                                    tags$ul(style = "font-size: 14px; line-height: 1.9; color: #333; padding-left: 20px; margin: 0;",
                                      tags$li(tags$strong("Selection:"), " Draw a region of interest using the freehand tool, then assign it to Group 1 or Group 2 for downstream comparison."),
                                      tags$li(tags$strong("Show Groups:"), " Enable “Show Groups on Map” to visualize saved group selections overlaid on the tissue image."),
                                      tags$li(tags$strong("Species:"), " Select the correct species (Human/Mouse) before using built-in gene signatures or pathway gene sets."),
                                      tags$li(tags$strong("Clustering:"), " Start with the default resolution (0.8) and increase it for finer subgroup identification"),
                                      tags$li(tags$strong("Export:"), " Export selected ROIs as Seurat subsets for reuse in SpatialScope or downstream analysis in external tools.")
                                    )
                                  ),

                                  # Ready to Begin
                                  tags$div(
                                    style = "background: linear-gradient(135deg, #0072B5 0%, #E18727 100%); color: white; padding: 25px; border-radius: 8px; text-align: center; margin-bottom: 20px;",
                                    tags$h2(style = "margin: 0 0 12px 0; font-size: 22px;", "🚀 Ready to Begin?"),
                                    tags$p(style = "margin: 0; font-size: 15px; line-height: 1.6;",
                                      "Click on the tools in the sidebar (🎨 Visualization, 🧬 Gene Sets, etc.) to start your spatial analysis. The tissue map will appear when you switch to any analysis tool.")
                                  ),

                                  # Contact
                                  tags$div(style = "text-align: center; padding: 16px; background: #f0f4f8; border-radius: 8px; font-size: 14px; color: #555;",
                                    HTML("Contact: Mengyao Lu (<a href='mailto:mel373@pitt.edu'>mel373@pitt.edu</a>) &nbsp;·&nbsp;
                                          Aodong Qiu (<a href='mailto:qiuaodon@pitt.edu'>qiuaodon@pitt.edu</a>) &nbsp;·&nbsp;
                                          Lujia Chen (<a href='mailto:luc17@pitt.edu'>luc17@pitt.edu</a>)")
                                  )
                                )

                              ) 
                      

                      ),

                      # Visualization content
                      tags$div(class = "control-content", id = "content_viz",
                               div(class = "panel-header", "🎨 Upload & Visualize"),
                                div(class = "control-section",
                                    h4("Data Source"),
                                    radioButtons("data_input_type", NULL,
                                                choices = c("Seurat Object (.rds)" = "rds",
                                                            "10x Visium Raw Output" = "raw"),
                                                selected = "rds"),

                                    # ── Existing RDS upload ──────────────────────────────────────────
                                    conditionalPanel(
                                      condition = "input.data_input_type == 'rds'",
                                      div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                                          actionButton("use_example_data", "📊 Use Example Data",
                                                      class = "btn btn-primary", style = "flex: 1;"),
                                          actionButton("show_upload_panel", "📤 Upload .rds",
                                                      class = "btn btn-info", style = "flex: 1;")
                                      ),
                                      conditionalPanel(
                                        condition = "input.show_upload_panel % 2 == 1",
                                        fileInput("upload_seurat", "Select Seurat Object (.rds)",
                                                  accept = c(".rds")),
                                        actionButton("load_uploaded_seurat", "Load Uploaded Data",
                                                    class = "btn btn-success btn-block"),
                                        tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                              "⚠️ Loading new data will replace current analysis")
                                      )
                                    ),

                                    # ── NEW: 10x raw upload ──────────────────────────────────────────
                                    conditionalPanel(
                                      condition = "input.data_input_type == 'raw'",
                                      tags$div(
                                        style = "background: #e8f4f8; border-left: 4px solid #0072B5; padding: 10px; margin-bottom: 12px;",
                                        tags$p(style = "font-size: 12px; font-weight: bold;", "📋 How to prepare your zip:"),
                                        tags$ol(style = "font-size: 11px; color: #555; padding-left: 18px;",
                                          tags$li("Locate your Space Ranger output folder"),
                                          tags$li("Make sure all files are uncompressed (no .gz files)"),
                                          tags$li("The folder must contain: .h5 file + spatial/ subfolder"),
                                          tags$li("Compress the entire folder as a .zip and upload")
                                        )
                                      ),

                                      fileInput("upload_visium_zip", "Upload Space Ranger Output (.zip)",
                                                accept = c(".zip")),
                                      # ⚙ Advanced settings (collapsed) — QC spot filters (Reviewer 2, item 3).
                                      tags$details(style = "margin:6px 0;",
                                        tags$summary(style = "cursor:pointer; font-weight:600; color:#2c3e50;",
                                                     "⚙ QC settings"),
                                        tags$div(style = "padding:8px 4px;",
                                          numericInput("qc_min_features", "Min. genes per spot:", value = 200, min = 0, step = 50),
                                          numericInput("qc_min_counts", "Min. counts per spot:", value = 500, min = 0, step = 100),
                                          numericInput("qc_max_mt", "Max. % mitochondrial:", value = 30, min = 0, max = 100, step = 5),
                                          tags$p(style = "font-size:11px; color:#7f8c8d;",
                                                 "Applied to raw SpaceRanger output only. Uploaded Seurat objects are used as provided.")
                                        )
                                      ),
                                      actionButton("load_raw_visium", "Load 10x Visium Data",
                                                  class = "btn btn-success btn-block",
                                                  style = "margin-top: 8px;")
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

                               ),

                              # ── NEW: Independent LR Score Section (NO RCTD needed) ────────────────
                              div(class = "control-section",
                                  h4("🔗 L-R Colocalization Score"),
                                  tags$p(style = "font-size: 12px; color: #7f8c8d; margin-bottom: 10px;",
                                        "Compute ligand-receptor geometric mean scores in selected region."),

                                  selectInput("lr_species", "Species:",
                                              choices = c("Human" = "human", "Mouse" = "mouse"),
                                              selected = "human"),

                                  radioButtons("lr_entity_type", "Analyze by:",
                                               choices = c("Groups" = "Groups", "ROIs" = "ROIs"),
                                               selected = "Groups", inline = TRUE),
                                  selectizeInput("lr_solo_target", "Region:", choices = NULL,
                                                 options = list(placeholder = "select an ROI or group")),

                                  checkboxGroupInput("lr_solo_db_filter", "Include databases:",
                                                    choices  = c("CellChat"              = "cellchat",
                                                                  "KEGG"                 = "kegg",
                                                                  "Guide to Pharmacology" = "guide2pharmacology",
                                                                  "Ramilowski"            = "ramilowski"),
                                                    selected = c("cellchat", "kegg", "guide2pharmacology", "ramilowski")),

                                  # ⚙ Advanced settings (collapsed) — Gaussian kernel bandwidth (Reviewer 1, item 4).
                                  tags$details(style = "margin:6px 0;",
                                    tags$summary(style = "cursor:pointer; font-weight:600; color:#2c3e50;",
                                                 "⚙ Advanced settings"),
                                    tags$div(style = "padding:8px 4px;",
                                      numericInput("lr_bandwidth_mult", "Kernel bandwidth × (Gaussian σ multiplier):",
                                                   value = 1.0, min = 0.25, max = 4, step = 0.25),
                                      tags$p(style = "font-size:11px; color:#7f8c8d; margin-top:4px;",
                                             "σ adapts to local spot spacing (median distance to the 12 nearest spots); this multiplies it. Larger = smoother, wider neighborhood.")
                                    )
                                  ),

                                  actionButton("run_lr_solo", " Run LR Scoring",
                                              class = "btn btn-primary btn-block",
                                              style = "margin-top: 10px;"),

                                  verbatimTextOutput("lr_solo_status"),

                                  conditionalPanel(
                                    condition = "output.lr_solo_results_available",
                                    tags$hr(),
                                    h5("Top L-R Pairs by Mean Score"),
                                    tags$p(style = "font-size: 12px; color: #7f8c8d; margin-bottom: 6px;",
                                      "👉 Click a row to show that pair on the map (grey → red by ligand × receptor)."),
                                    DT::dataTableOutput("lr_simple_table"),
                                    hr(),
                                    downloadButton("dl_lr_solo", "Download LR Enrichment Table",
                                                  class = "btn btn-warning btn-block")
                                  )
                              ),

                              # ── existing deconvolution section follows below ───────────────────────



                               # Cell Type Deconvolution Section

                                div(class = "control-section",
                                    h4("🔬 Cell Type Deconvolution"),
                                    tags$p(style = "font-size: 12px; color: #7f8c8d; margin-bottom: 10px;",
                                          "Estimate cell type proportions in ROI using RCTD."),

                                    # ── Reference source selector ──────────────────────────────────────────
                                    h5("Reference Data"),
                                    radioButtons("ref_source", NULL,
                                                choices = c("Use built-in reference" = "builtin",
                                                            "Upload my own reference" = "upload"),
                                                selected = "builtin"),

                                    # Built-in reference picker
                                    conditionalPanel(
                                      condition = "input.ref_source == 'builtin'",
                                      selectInput("builtin_ref_choice",
                                                  "Select built-in reference:",
                                                  choices = c(
                                                    "CRC — Colorectal Cancer (SMC cohort)" = "crc"
                                                    # Add more here as you package new references, e.g.:
                                                    # "BRCA — Breast Cancer" = "brca",
                                                    # "LUAD — Lung Adenocarcinoma" = "luad"
                                                  )),
                                      
                                      tags$p(style = "font-size: 11px; color: #7f8c8d; margin-top: 4px;",
                                            "Additional references available at our ", 
                                            tags$a("GitHub Releases", href = "https://github.com/YOUR_REPO/releases", target = "_blank"), "."),

                                      actionButton("load_builtin_ref", "Load Built-in Reference",
                                                  class = "btn btn-info btn-block",
                                                  style = "margin-top: 6px;")
                                    ),

                                    # User upload
                                    conditionalPanel(
                                      condition = "input.ref_source == 'upload'",
                                      div(style = "display: flex; gap: 10px; margin-bottom: 10px;",
                                          actionButton("show_ref_upload", "📤 Upload scRNA-seq Reference",
                                                      class = "btn btn-info",
                                                      style = "flex: 1;")
                                      ),
                                      conditionalPanel(
                                        condition = "input.show_ref_upload % 2 == 1",
                                        fileInput("upload_reference", "Select scRNA-seq Reference (.rds)",
                                                  accept = c(".rds")),
                                        tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                              "⚠️ Loading a new reference will reset deconvolution results.")
                                      )
                                    ),


                                    # Only show cell type column selector for user-uploaded references


                                    verbatimTextOutput("ref_status"),

                                    # ── Run button ─────────────────────────────────────────────────────────
                                    actionButton("run_deconv", "🔬 Estimate Cell Composition",
                                                class = "btn btn-primary btn-block",
                                                style = "margin-top: 10px;"),
                                    tags$p(style = "font-size: 11px; color: #e67e22; margin-top: 5px;",
                                          "💡 Save spots to Group 1 or Group 2 on the map first."),

                                    # ── Results ────────────────────────────────────────────────────────────
                                    conditionalPanel(
                                      condition = "output.deconv_results_available",
                                      tags$hr(),
                                      h5("Cell Type Proportions"),
                                      verbatimTextOutput("deconv_gene_overlap_msg"),
                                      selectInput("deconv_group", "Show results for:",
                                                  choices = c("Group 1" = "group1", "Group 2" = "group2", "All Spots" = "all")),
                                      plotOutput("deconv_barplot", height = "280px"),
                                      div(style = "max-height: 300px; overflow-y: auto;",
                                          tableOutput("deconv_table")),
                                      br(),
                                      downloadButton("dl_deconv", "Download Cell Type Proportions",
                                                    class = "btn btn-warning btn-block")
                                    )
                                ),


                                # ── L-R Colocalization Section ────────────────────────────────────────────────
                                # Shows only once RCTD has produced results (same conditional gate)
                                # conditionalPanel(
                                #   condition = "output.deconv_results_available",
                                #   checkboxInput("use_ewa", "🔬 Rare Signal Recovery (reference-anchored attribution)", value = TRUE),
                                #   tags$p(style = "font-size: 11px; color: #7f8c8d;",
                                #         "Uses reference gene signatures to correct dominant cell-type bias in L/R attribution."),
                                #   div(class = "control-section",
                                #       h4("🧬 L-R Cell Communication"),
                                #       tags$p(style = "font-size: 12px; color: #7f8c8d; margin-bottom: 10px;",
                                #             "L-R colocalization with inferred sender/receiver cell types using RCTD + reference-anchored attribution."),
                                #   div(style = "background-color: #fff8e1; border-left: 4px solid #f39c12; padding: 10px; margin-bottom: 10px;",
                                #       tags$p(style = "margin: 0; font-size: 12px; color: #7f8c8d;",
                                #             "⚠️ ", tags$b("Reference data quality matters."),
                                #             " Cell type attribution accuracy depends heavily on your reference. ",
                                #             "We provide a built-in CRC reference for the demo dataset, but ",
                                #             tags$b("we strongly encourage uploading your own tissue-matched scRNA-seq reference"),
                                #             " for best results. Mismatched references may lead to incorrect sender/receiver inference.")
                                #   ),                                                                                  

                                #       # ── Group selector ──────────────────────────────────────────────────────
                                #       selectInput("lr_group", "Analyze group:",
                                #                   choices = c("Group 1" = "group1", "Group 2" = "group2")),

                                #       # # ── L-R database source ─────────────────────────────────────────────────
                                #       # h5("L-R Database"),
                                #       # radioButtons("lr_db_source", NULL,
                                #       #             choices = c("Use built-in database" = "builtin",
                                #       #                         "Upload my own (.tsv)"  = "upload"),
                                #       #             selected = "builtin"),

                                #       # conditionalPanel(
                                #       #   condition = "input.lr_db_source == 'upload'",
                                #       #   fileInput("upload_lr_db", "Select lr_network .tsv file",
                                #       #             accept = c(".tsv", ".txt")),
                                #       #   tags$p(style = "font-size: 11px; color: #7f8c8d;",
                                #       #         "Required columns: from, to, database")
                                #       # ),

                                #       # ── Database filter ─────────────────────────────────────────────────────
                                #       checkboxGroupInput("lr_db_filter", "Include databases:",
                                #                         choices  = c("KEGG"                 = "kegg",
                                #                                       "Guide to Pharmacology" = "guide2pharmacology",
                                #                                       "Ramilowski"            = "ramilowski"),
                                #                         selected = c("kegg", "guide2pharmacology", "ramilowski")),

                                #       # ── Spatial lag neighbors ───────────────────────────────────────────────
                                    
                                #       tags$p(style = "font-size: 11px; color: #7f8c8d;",
                                #               "Spatial smoothing neighborhood size: 12 spots"),            

                                #       # ── Run button ──────────────────────────────────────────────────────────
                                #       actionButton("run_lr", "🧬 Run L-R Cell Communication",
                                #                   class = "btn btn-primary btn-block",
                                #                   style = "margin-top: 10px;"),
                                #       tags$p(style = "font-size: 11px; color: #e67e22; margin-top: 5px;",
                                #             "⚠️ Run RCTD deconvolution first — correlations need cell-type proportions."),

                                #       verbatimTextOutput("lr_status"),

                                #       # ── Results ─────────────────────────────────────────────────────────────
                                #       conditionalPanel(
                                #         condition = "output.lr_results_available",
                                #         tags$hr(),
                                #         h5("Top L-R Pairs by Mean Score"),
                                #         tags$p(style = "font-size: 11px; color: #7f8c8d; margin-bottom: 8px;",
                                #               "Ranked by interaction strength in selected ROI. "),
                                #              # "Score = geometric mean of spatially-lagged Ligand × Receptor. ",
                                #              # "Correlation = Spearman r between spot-level score and top RCTD cell type."),

                                #         # numericInput("lr_top_n", "Show top N pairs:",
                                #         #             value = 20, min = 5, max = 200, step = 5),

                              
                                #         DT::dataTableOutput("lr_table"),

                                #         downloadButton("dl_lr_results", "⬇ Download Full Results",
                                #                       class = "btn btn-success btn-sm",
                                #                       style = "margin-top: 8px;"),
                                #       )
                                #   )
                                # ),  # end LR section

                         
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
                                  h4("Pathway Signatures (MSigDB Hallmark)"),
                                  selectInput("pathway_library", "Hallmark Pathways:",
                                              choices = c("None", names(hallmark_library_human)),
                                              selected = "None"),
                                  conditionalPanel(
                                    condition = "input.pathway_library != 'None'",
                                    actionButton("load_pathway", "Load Pathway Genes",
                                                class = "btn btn-warning btn-block")
                                  ),
                                  tags$p(style = "font-size: 11px; color: #7f8c8d; margin-top: 5px;",
                                        "💡 Loads into Gene Input below — same scoring and visualization applies")
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
                                   selectizeInput("cluster_region", "Cluster which spots:",
                                                  choices = c("All spots" = "__all__"),
                                                  selected = "__all__",
                                                  options = list(placeholder = "All spots, an ROI, or a group")),
                                   textOutput("cluster_spot_count"),
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 5px;",
                                          "💡 Tip: draw & save ROIs (or make a group) to cluster within a region")
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
                                     plotOutput("cluster_umap", height = "300px"),
                                     downloadButton("dl_cluster_umap", "Download UMAP Figure (PDF)",
                                                    class = "btn btn-info btn-block"),
                                     br(),
                                      downloadButton("dl_cluster", "Download Cluster Assignments",
                                                    class = "btn btn-warning btn-block")
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
                                   # Removed the redundant "Show on Map" checkbox (Reviewer 1, item 7):
                                   # it shared the same inputId ("show_groups") as the central
                                   # "Show Groups on Map" control, which made the two boxes drift out
                                   # of visual sync. Group highlighting is now driven by that single
                                   # central checkbox; this panel keeps the Group 1 / Group 2 readouts.
                                   tags$p(style = "font-size: 12px; color: #7f8c8d; margin-top: 10px;",
                                          "💡 Tip: Toggle \"Show Groups on Map\" in the visualization panel to highlight the saved groups. Use the group buttons at the bottom of the map to save selections and download spot IDs.")
                               ),
                               div(class = "control-section",
                                   h4("Analysis"),
                                   # Compare ROI-vs-ROI or Group-vs-Group (Reviewer 1, item 1).
                                   # Side A / Side B can each be any ROI or Group (cross-type allowed).
                                   selectizeInput("deg_side_a", "Side A:", choices = NULL,
                                                  options = list(placeholder = "select an ROI or group")),
                                   selectizeInput("deg_side_b", "Side B:", choices = NULL,
                                                  options = list(placeholder = "ROI/group, or Rest of tissue")),
                                   tags$p(style = "font-size:11px; color:#7f8c8d; margin:2px 0;",
                                          "Positive log2FC = higher in Side A."),
                                   # ⚙ Advanced settings (collapsed) — DE test + thresholds (Reviewer 2, item 3).
                                   tags$details(style = "margin:6px 0;",
                                     tags$summary(style = "cursor:pointer; font-weight:600; color:#2c3e50;",
                                                  "⚙ Advanced settings"),
                                     tags$div(style = "padding:8px 4px;",
                                       selectInput("deg_test", "DE test:",
                                                   choices = c("Wilcoxon" = "wilcox", "t-test" = "t",
                                                               "Logistic regression" = "LR", "ROC" = "roc"),
                                                   selected = "wilcox"),
                                       numericInput("deg_logfc", "|log2FC| ≥", value = 0.25, min = 0, max = 5, step = 0.05),
                                       numericInput("deg_minpct", "Min. spot fraction (min.pct):", value = 0.05, min = 0, max = 1, step = 0.01),
                                       sliderInput("volcano_fdr", "Adjusted p (BH) threshold:", min = 0.01, max = 0.2,
                                                   value = 0.05, step = 0.01)
                                     )
                                   ),
                                   actionButton("run_deg", "Find DEGs", class = "btn btn-danger btn-block"),
                                   
                               ),
                               conditionalPanel(
                                 condition = "output.deg_results_available",
                                 div(class = "control-section",
                                     h4("Top DEGs"),
                                     tags$p(style = "font-size:12px; color:#2c3e50; font-weight:600;",
                                            textOutput("deg_direction", inline = TRUE)),
                                      div(style = "max-height: 400px; overflow-y: auto;",
                                          tableOutput("deg_table")),
                                  plotOutput("deg_volcano", height = "450px"),
                                  downloadButton("dl_deg_volcano", "Download Volcano Figure (PDF)",
                                                 class = "btn btn-warning btn-block"),
                                  plotOutput("deg_moran_volcano", height = "450px"),

                                      tags$p(style = "font-size:12px; color:#7f8c8d; margin-top: 6px;",
                                        "💡 Tip: ", tags$strong("p_adj"), " reflects differential expression; ",
                                        tags$strong("Moran's I"), " (with adjusted p-values) quantifies spatial autocorrelation."
                                      ),                                         
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
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Pathway" = "pathway", "Cell Signature" = "cellsig", "Saved Gene Set" = "geneset")),
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
                                    # ← INSERT HERE
                                    conditionalPanel(
                                      condition = "input.violin_feature_type == 'pathway'",
                                      selectInput("violin_pathway", "Pathway:", choices = c("", names(hallmark_library_human)))
                                    ),        
                                    conditionalPanel(
                                      condition = "input.violin_feature_type == 'cellsig'",
                                      selectInput("violin_cellsig", "Cell Type:", choices = character(0))
                                    ),                                                               
                                   # Side A / Side B can each be any ROI or Group (cross-type allowed).
                                   selectizeInput("violin_side_a", "Side A:", choices = NULL,
                                                  options = list(placeholder = "select an ROI or group")),
                                   selectizeInput("violin_side_b", "Side B:", choices = NULL,
                                                  options = list(placeholder = "ROI/group, or Rest of tissue")),
                                   selectInput("violin_stat_test", "Test:",
                                               choices = c("Wilcoxon" = "wilcox", "t-test" = "ttest"),
                                               selected = "wilcox"),
                                   actionButton("plot_violin", "Generate Plot", class = "btn btn-info btn-block"),
                                   conditionalPanel(
                                     condition = "output.violin_available",
                                     plotOutput("violin_plot", height = "300px"),
                                      br(),
                                        downloadButton("dl_compare_plot_group", "Download Comparison Plot",
                                                      class = "btn btn-warning btn-block"))
                               ),
                               div(class = "control-section",
                                   h4("Feature vs Feature"),
                                   selectInput("compare_type1", "Feature 1:",
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Pathway" = "pathway", "Cell Signature" = "cellsig", "Saved Gene Set" = "geneset")),
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

                                    conditionalPanel(
                                      condition = "input.compare_type1 == 'pathway'",
                                      selectInput("compare_pathway1", "Pathway:", choices = c("", names(hallmark_library_human)))
                                    ),    

                                    conditionalPanel(
                                      condition = "input.compare_type1 == 'cellsig'",
                                      selectInput("compare_cellsig1", "Cell Type:", choices = character(0))
                                    ),                                                                   
                                   selectInput("compare_type2", "Feature 2:",
                                               choices = c("Gene" = "gene", "Metadata" = "metadata", "Pathway" = "pathway", "Cell Signature" = "cellsig", "Saved Gene Set" = "geneset")),
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
                                    conditionalPanel(
                                      condition = "input.compare_type2 == 'pathway'",
                                      selectInput("compare_pathway2", "Pathway:", choices = c("", names(hallmark_library_human)))
                                    ),
                                    conditionalPanel(
                                      condition = "input.compare_type2 == 'cellsig'",
                                      selectInput("compare_cellsig2", "Cell Type:", choices = character(0))
                                    ),                                                                                                  
                                   selectizeInput("compare_spots_selection", "Use spots from:",
                                                  choices = c("All spots" = "__all__"), selected = "__all__",
                                                  options = list(placeholder = "All spots, an ROI, or a group")),
                                   selectInput("stat_test", "Test:", choices = c("Wilcoxon" = "wilcox", "t-test" = "ttest")),
                                   # selectInput("cor_method", "Correlation:", choices = c("Spearman" = "spearman", "Pearson" = "pearson")),
                                   actionButton("plot_compare", "Compare", class = "btn btn-info btn-block"),
                                   conditionalPanel(
                                     condition = "output.compare_available",
                                     plotOutput("compare_plot", height = "400px"),
                                     br(),
                                      downloadButton("dl_compare_plot", "Download Comparison Plot",
                                                    class = "btn btn-warning btn-block")
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

    # ── Session isolation ─────────────────────────────────────────────────────
    # Create session-local copies of every mutable data object. These shadow the
    # outer-scope objects, so the `<<-` writes in the upload / clustering /
    # gene-set handlers below resolve to THIS session's copy instead of the
    # shared app-level binding. Without this, one user's upload overwrites the
    # data for every concurrent session (Reviewer 3, item 1: session isolation /
    # data privacy). R is copy-on-modify, so this is a cheap reference copy.
    seurat_obj      <- seurat_obj
    spots_sf        <- spots_sf
    he_image_base64 <- he_image_base64
    he_image_bounds <- he_image_bounds
    all_genes       <- all_genes
    all_metadata    <- all_metadata
    # ──────────────────────────────────────────────────────────────────────────

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
    options(future.globals.maxSize = 8000 * 1024^2) 
    # ── Fix for Seurat FindMarkers future globals error ──


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
    # ── N-group model (Reviewer 1, item 1) ────────────────────────────────────
    # Single source of truth: a named list of groups, each holding its member
    # spot IDs and its drawn ROI boundary rings. Supports any number of groups and
    # multiple ROIs per group. Two default groups exist so the app behaves exactly
    # as before out of the box. The legacy group1_*/group2_* accessors below derive
    # from the first two groups, so existing analysis code keeps working unchanged.
    # Palette for distinguishing ROIs on the map (cycles if > length).
    .roi_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                      "#A65628", "#F781BF", "#1B9E77", "#666666", "#66A61E")
    group_color <- function(i) .roi_palette[((i - 1) %% length(.roi_palette)) + 1]

    # ROIs: first-class named regions, each with its own spots + boundary rings.
    rois        <- reactiveVal(list())   # "ROI 1" = list(spots=chr, rings=list())
    roi_counter <- reactiveVal(0)        # for default names "ROI 1", "ROI 2", ...
    # Groups: named collections of ROI names (each ROI keeps its identity inside).
    groups      <- reactiveVal(list())   # "Tumor" = list(members=chr)
    roi_names   <- reactive(names(rois()))
    group_names <- reactive(names(groups()))

    # Save the current drawing as a named ROI.
    save_roi <- function(name, spots, rings) {
      r <- rois(); r[[name]] <- list(spots = spots, rings = rings); rois(r)
    }
    # Combined spots / boundary rings for a group = union over its member ROIs.
    group_spots_of <- function(gname) {
      gg <- groups()[[gname]]
      if (is.null(gg) || length(gg$members) == 0) return(character(0))
      r <- rois()
      unique(unlist(lapply(gg$members, function(m) if (!is.null(r[[m]])) r[[m]]$spots else NULL),
                    use.names = FALSE))
    }
    group_rings_of <- function(gname) {
      gg <- groups()[[gname]]
      if (is.null(gg) || length(gg$members) == 0) return(list())
      r <- rois(); out <- list()
      for (m in gg$members) if (!is.null(r[[m]])) out <- c(out, r[[m]]$rings)
      out
    }
    reset_groups <- function() { rois(list()); groups(list()); roi_counter(0) }

    # Legacy accessors: map the first two GROUPS to the old group1/group2 slots so
    # existing DEG / plotting / download code keeps working until A3 wires N groups.
    .grp_spots <- function(i) { nm <- names(groups()); if (length(nm) >= i) group_spots_of(nm[i]) else character(0) }
    .grp_rings <- function(i) { nm <- names(groups()); if (length(nm) >= i) group_rings_of(nm[i]) else list() }
    group1_spots <- reactive(.grp_spots(1))
    group2_spots <- reactive(.grp_spots(2))
    group1_roi   <- reactive(.grp_rings(1))
    group2_roi   <- reactive(.grp_rings(2))

    # ── Unified region picker (Reviewer feedback) ─────────────────────────────
    # One choice list covering every ROI and Group, plus a resolver to spot IDs,
    # so any analysis can compare across types (e.g. a Group vs a single ROI).
    # Keys are prefixed "roi:" / "grp:"; "__rest__" means the rest of the tissue.
    region_choices <- reactive({
      rn <- roi_names(); gn <- group_names()
      c(if (length(rn)) setNames(paste0("roi:", rn), paste0("ROI: ", rn)),
        if (length(gn)) setNames(paste0("grp:", gn), paste0("Group: ", gn)))
    })
    region_spots <- function(key) {
      if (is.null(key) || length(key) == 0 || key == "") return(character(0))
      if (startsWith(key, "roi:")) { nm <- sub("^roi:", "", key); rr <- rois()[[nm]]; if (is.null(rr)) character(0) else rr$spots }
      else if (startsWith(key, "grp:")) group_spots_of(sub("^grp:", "", key))
      else character(0)
    }
    region_label <- function(key) {
      if (is.null(key) || length(key) == 0 || key == "") return("")
      if (identical(key, "__rest__")) return("Rest")
      if (startsWith(key, "roi:")) sub("^roi:", "", key)
      else if (startsWith(key, "grp:")) sub("^grp:", "", key)
      else key
    }
    cluster_colors_palette <- reactiveVal(NULL) 
    last_group_plot <- reactiveVal(NULL) 

    lr_db_path <- reactive({
      species <- if (is.null(input$lr_species) || input$lr_species == "") "human" else input$lr_species
      if (species == "mouse") {
        system.file("extdata", "lr_network_combined_mouse.tsv", package = "SpatialScope")
      } else {
        system.file("extdata", "lr_network_combined_human.tsv", package = "SpatialScope")
      }
    })


    lr_network <- reactive({
      read.table(lr_db_path(), sep = "\t", header = TRUE)
    })    
    output_group1 <- reactive(.grp_spots(1))  # bridged: first group's combined spots
    output_group2 <- reactive(.grp_spots(2))  # bridged: second group's combined spots
    deg_results <- reactiveVal(NULL)
    violin_data <- reactiveVal(NULL)
    compare_data <- reactiveVal(NULL)
    gene_set_scores <- reactiveVal(list())
    current_gene_set_score <- reactiveVal(NULL)
    cluster_results <- reactiveVal(NULL)
    showing_gene_set <- reactiveVal(FALSE)  # Track if we're displaying gene set scores

    # ── Figure downloads (Reviewer 3, item 5) ─────────────────────────────────
    # Each plot output captures its ggplot into one of these so it can be saved.
    deg_volcano_rv  <- reactiveVal(NULL)
    violin_rv       <- reactiveVal(NULL)
    cluster_umap_rv <- reactiveVal(NULL)
    # Shared publication-quality PDF download handler for a captured ggplot.
    make_plot_download <- function(rv, stem, w = 8, h = 6) {
      downloadHandler(
        filename = function() paste0(stem, "_", format(Sys.time(), "%Y%m%d"), ".pdf"),
        content  = function(file) {
          p <- rv()
          if (is.null(p)) { showNotification("Generate the plot first.", type = "warning"); return() }
          ggplot2::ggsave(file, plot = p, width = w, height = h, dpi = 300, device = "pdf")
        }
      )
    }
    output$dl_deg_volcano  <- make_plot_download(deg_volcano_rv,  "DEG_volcano")
    output$dl_violin_plot  <- make_plot_download(violin_rv,       "violin", w = 6, h = 6)
    output$dl_cluster_umap <- make_plot_download(cluster_umap_rv, "UMAP_clusters", w = 7, h = 6)

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

    updateSelectInput(session, "compare_cellsig1", choices = c("", names(cellmarker_db)))
    updateSelectInput(session, "compare_cellsig2", choices = c("", names(cellmarker_db)))
    updateSelectInput(session, "violin_cellsig",   choices = c("", names(cellmarker_db)))

    current_hallmark_library <- reactive({
      if (input$species_select == "human") hallmark_library_human else hallmark_library_mouse
    })
    # Populate the Hallmark dropdowns after load / on species change. If the
    # library is unexpectedly empty, show an explicit message instead of a blank
    # list so the failure is visible rather than silent (Reviewer 3, item 4).
    observe({
      lib <- current_hallmark_library()
      pw_names <- names(lib)
      if (length(pw_names) == 0) {
        updateSelectInput(session, "pathway_library",
                          choices = "No Hallmark pathways available on this server",
                          selected = NULL)
      } else {
        updateSelectInput(session, "pathway_library",
                          choices = c("None", pw_names), selected = "None")
        # Keep the DEG/violin/comparison pathway pickers in sync with species too.
        updateSelectInput(session, "violin_pathway",   choices = c("", pw_names))
        updateSelectInput(session, "compare_pathway1", choices = c("", pw_names))
        updateSelectInput(session, "compare_pathway2", choices = c("", pw_names))
      }
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
    # ── Shared dataset loader ─────────────────────────────────────────────
    # Applies an already-read Seurat object to the current session: swaps in
    # the session-local objects, refreshes dropdowns, rebuilds coordinates and
    # the H&E image, clears results, and redraws the map. Called by BOTH the
    # upload handler and the "Use Example Data" button so the example can be
    # reloaded on demand (Reviewer 3, item 1).
    apply_loaded_seurat <- function(new_seurat) {
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
            image_obj     <- new_seurat@images[[image_name]]
            he_image_data <- image_obj@image
            scale_factor  <- image_obj@scale.factors$lowres
            H <- dim(he_image_data)[1]
            W <- dim(he_image_data)[2]

            coords_full <- GetTissueCoordinates(new_seurat, image = image_name)

            if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
              pixel_x <- coords_full$pxl_col_in_fullres * scale_factor
              pixel_y <- coords_full$pxl_row_in_fullres  * scale_factor
            } else if ("imagecol" %in% colnames(coords_full)) {
              pixel_x <- coords_full$imagecol
              pixel_y <- coords_full$imagerow
            } else {
              # spots_sf$y is already flipped — un-flip before using as pixel row
              pixel_x <- spots_sf$x
              pixel_y <- max(spots_sf$y) + min(spots_sf$y) - spots_sf$y
            }

            pixel_x_min <- min(pixel_x, na.rm = TRUE); pixel_x_max <- max(pixel_x, na.rm = TRUE)
            pixel_y_min <- min(pixel_y, na.rm = TRUE); pixel_y_max <- max(pixel_y, na.rm = TRUE)
            x_buffer_px <- (pixel_x_max - pixel_x_min) * 0.1
            y_buffer_px <- (pixel_y_max - pixel_y_min) * 0.1

            crop_x_min <- max(1, floor(pixel_x_min - x_buffer_px))
            crop_x_max <- min(W, ceiling(pixel_x_max + x_buffer_px))
            crop_y_min <- max(1, floor(pixel_y_min - y_buffer_px))
            crop_y_max <- min(H, ceiling(pixel_y_max + y_buffer_px))

            # Guard against invalid ranges
            if (crop_x_min >= crop_x_max) { crop_x_min <- 1; crop_x_max <- W }
            if (crop_y_min >= crop_y_max) { crop_y_min <- 1; crop_y_max <- H }

            he_image_cropped <- he_image_data[crop_y_min:crop_y_max, crop_x_min:crop_x_max, ]
            temp_file <- tempfile(fileext = ".png")
            png(temp_file, width = dim(he_image_cropped)[2], height = dim(he_image_cropped)[1])
            par(mar = c(0, 0, 0, 0))
            plot(as.raster(he_image_cropped), axes = FALSE)
            dev.off()
            he_image_base64 <<- paste0("data:image/png;base64,", base64enc::base64encode(temp_file))
            unlink(temp_file)

            # Bounds in flipped leaflet coordinate system
            y_sum <- max(spots_sf$y) + min(spots_sf$y)
            he_image_bounds <<- list(
              north = y_sum - crop_y_min,
              south = y_sum - crop_y_max,
              west  = crop_x_min,
              east  = crop_x_max
            )

            session$sendCustomMessage("updateHEImage", list(
              imageUrl = he_image_base64,
              bounds   = he_image_bounds
            ))
            Sys.sleep(0.5)
          }
        }, error = function(e) {
          print(paste("Could not extract H&E image from new data:", e$message))
          he_image_base64 <<- NULL
        })

        # Clear all selections and results
        drawn_feats(list())
        selected_spots(character(0))
        reset_groups()
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
    }

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

        new_seurat <- tryCatch(
          UpdateSeuratObject(new_seurat),
          error = function(e) {
            message("UpdateSeuratObject failed, using original: ", e$message)
            readRDS(input$upload_seurat$datapath)  # reload original
          }
        )

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

        apply_loaded_seurat(new_seurat)

      }, error = function(e) {
        removeNotification(id = "load_seurat")

        # Hide loading overlay on error
        shinyjs::runjs("$('#loading_overlay').removeClass('active');")

        showNotification(paste("Error loading Seurat object:", e$message),
                         type = "error", duration = 10)
      })
    })

### allow user to upload the whole raw data with a zip file
    observeEvent(input$load_raw_visium, {
      req(input$upload_visium_zip)

      tmp_dir <- tempfile()
      dir.create(tmp_dir)

      showNotification("Unzipping and loading 10x Visium data...", id = "raw_loading", duration = NULL)
      shinyjs::runjs("
        $('#loading_message').text('Loading 10x Visium data...');
        $('#loading_overlay').addClass('active');
      ")

      tryCatch({
        # ── 1. Unzip ──────────────────────────────────────────────────────────
        unzip(input$upload_visium_zip$datapath, exdir = tmp_dir)

        # ── 2. Find and normalize spatial/ directory ──────────────────────────
        spatial_dirs <- list.dirs(tmp_dir, recursive = TRUE, full.names = TRUE)
        spatial_dirs <- spatial_dirs[basename(spatial_dirs) == "spatial"]

        spatial_dir <- NULL
        for (d in spatial_dirs) {
          if (length(list.files(d)) > 0) { spatial_dir <- d; break }
        }

        if (!is.null(spatial_dir)) {
          cat("Using spatial dir:", spatial_dir, "\n")
          all_files <- list.files(spatial_dir, full.names = TRUE)
          name_patterns <- list(
            "tissue_lowres_image.png" = "tissue_lowres_image\\.png",
            "tissue_hires_image.png"  = "tissue_hires_image\\.png",
            "tissue_positions_list.csv" = "tissue_positions(_list)?\\.csv",
            "scalefactors_json.json"  = "scalefactors_json\\.json"
          )
          for (std_name in names(name_patterns)) {
            matches <- all_files[grepl(name_patterns[[std_name]], basename(all_files))]
            if (length(matches) > 0 && basename(matches[1]) != std_name) {
              file.rename(matches[1], file.path(spatial_dir, std_name))
            }
          }
          cat("Files after renaming:\n"); print(list.files(spatial_dir))
       
          lowres_path <- file.path(spatial_dir, "tissue_lowres_image.png")
          hires_path  <- file.path(spatial_dir, "tissue_hires_image.png")
          if (!file.exists(lowres_path) && file.exists(hires_path)) {
            json_path <- file.path(spatial_dir, "scalefactors_json.json")
            sf <- jsonlite::fromJSON(json_path)
            lf <- sf$tissue_lowres_scalef
            hf <- sf$tissue_hires_scalef
            
            hires_img <- magick::image_read(hires_path)
            info <- magick::image_info(hires_img)
            target_w <- round(info$width  * lf / hf)
            target_h <- round(info$height * lf / hf)
            
            magick::image_scale(hires_img, paste0(target_w, "x", target_h, "!")) %>%
              magick::image_write(lowres_path)
            cat("Resized hires to lowres dimensions:", target_h, "x", target_w, "\n")
          }                  
        } else {
            cat("WARNING: No non-empty spatial directory found\n")
          }  

        # ── 3. Locate .h5 and spatial/ ────────────────────────────────────────
        h5_files <- list.files(tmp_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
        if (length(h5_files) == 0) stop("No .h5 file found in zip.")
        h5_file  <- h5_files[1]
        data_dir <- dirname(h5_file)

        spatial_dirs2 <- list.dirs(tmp_dir, recursive = TRUE, full.names = TRUE)
        spatial_dir2  <- spatial_dirs2[basename(spatial_dirs2) == "spatial"]
        spatial_dir2  <- spatial_dir2[sapply(spatial_dir2, function(d) length(list.files(d)) > 0)][1]
        if (dirname(spatial_dir2) != data_dir) {
          file.copy(spatial_dir2, data_dir, recursive = TRUE)
        }

        # ── 4. Load with Seurat ───────────────────────────────────────────────
        new_obj <- Seurat::Load10X_Spatial(data.dir = data_dir, filename = basename(h5_file))
        new_obj <- tryCatch(UpdateSeuratObject(new_obj), error = function(e) new_obj)
        new_obj <- NormalizeData(new_obj)

        new_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(new_obj, pattern = "^MT-|^mt-")
        new_obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(new_obj, pattern = "^RP[SL]|^Rp[sl]")
        n_before <- ncol(new_obj)
        # QC thresholds from the UI (Reviewer 2, item 3), with the original defaults.
        nf_thr <- if (is.null(input$qc_min_features)) 200 else input$qc_min_features
        nc_thr <- if (is.null(input$qc_min_counts))   500 else input$qc_min_counts
        mt_thr <- if (is.null(input$qc_max_mt))         30 else input$qc_max_mt
        .md <- new_obj@meta.data
        keep_cells <- rownames(.md)[.md$nFeature_Spatial > nf_thr &
                                    .md$nCount_Spatial   > nc_thr &
                                    .md$percent.mt        < mt_thr]
        new_obj  <- subset(new_obj, cells = keep_cells)
        n_after <- ncol(new_obj)
        showNotification(paste0("QC: kept ", n_after, "/", n_before, " spots"),
                        type = "message", duration = 5)

        seurat_obj   <<- new_obj
        all_genes    <<- rownames(seurat_obj)
        all_metadata <<- colnames(seurat_obj@meta.data)

        updateSelectInput(session, "gene_select",      choices = c("", all_genes))
        updateSelectInput(session, "violin_gene",      choices = c("", all_genes))
        updateSelectInput(session, "compare_gene1",    choices = c("", all_genes))
        updateSelectInput(session, "compare_gene2",    choices = c("", all_genes))
        updateSelectInput(session, "meta_select",      choices = c("", all_metadata))
        updateSelectInput(session, "violin_metadata",  choices = c("", all_metadata))
        updateSelectInput(session, "compare_meta1",    choices = c("", all_metadata))
        updateSelectInput(session, "compare_meta2",    choices = c("", all_metadata))

        # ── 5. Extract coordinates → lowres pixel space ───────────────────────
        image_name   <- names(seurat_obj@images)[1]
        image_obj    <- seurat_obj@images[[image_name]]      # ← saved for H&E below
        scale_factor <- image_obj@scale.factors$lowres

        tryCatch({
          coords <- GetTissueCoordinates(seurat_obj, image = image_name)
        }, error = function(e) {
          coords <- seurat_obj@images[[image_name]]@coordinates
        })

        spots_df <- data.frame(spot_id = rownames(coords), stringsAsFactors = FALSE)

        if ("pxl_col_in_fullres" %in% colnames(coords)) {
          spots_df$x <- coords$pxl_col_in_fullres * scale_factor
          spots_df$y <- coords$pxl_row_in_fullres  * scale_factor
        } else if ("imagecol" %in% colnames(coords)) {
          spots_df$x <- coords$imagecol
          spots_df$y <- coords$imagerow
        } else if ("col" %in% colnames(coords)) {
          spots_df$x <- coords$col * scale_factor
          spots_df$y <- coords$row * scale_factor
        } else {
          spots_df$x <- coords[, 1] * scale_factor
          spots_df$y <- coords[, 2] * scale_factor
        }

        spots_sf <<- st_as_sf(spots_df, coords = c("x", "y"), crs = NA)
        coords_matrix <- do.call(rbind, st_geometry(spots_sf)) %>% as.matrix()
        spots_sf$x <<- coords_matrix[, 1]
        spots_sf$y <<- coords_matrix[, 2]
        spots_sf$y <<- max(spots_sf$y) - spots_sf$y + min(spots_sf$y)   # flip Y

        # ── 6. H&E image ──────────────────────────────────────────────────────
        # NOTE: spots_sf$y is already flipped here, so bounds formula is correct
        tryCatch({
          if (show_image) {
            he_image_data <- image_obj@image
            H <- dim(he_image_data)[1]
            W <- dim(he_image_data)[2]
            cat("H&E image dim:", H, "x", W, "\n")

            temp_file <- tempfile(fileext = ".png")
            png(temp_file, width = W, height = H)
            par(mar = c(0, 0, 0, 0))
            plot(as.raster(he_image_data), axes = FALSE)
            dev.off()

            he_image_base64 <<- paste0(
              "data:image/png;base64,",
              base64enc::base64encode(temp_file)
            )
            unlink(temp_file)

            # Align image to flipped coordinate system used by spots
            y_sum <- max(spots_sf$y) + min(spots_sf$y)
            he_image_bounds <<- list(
              north = y_sum,       # pixel row 0   → top of image
              south = y_sum - H,   # pixel row H   → bottom of image
              west  = 0,           # pixel col 0   → left edge
              east  = W            # pixel col W   → right edge
            )

            session$sendCustomMessage("updateHEImage", list(
              imageUrl = he_image_base64,
              bounds   = he_image_bounds
            ))
            Sys.sleep(0.5)
          }
        }, error = function(e) {
          cat("Could not extract H&E image:", e$message, "\n")
          he_image_base64 <<- NULL
        })

        # ── 7. Clear state and redraw map ─────────────────────────────────────
        drawn_feats(list())
        selected_spots(character(0))
        reset_groups()
        current_values(NULL)
        deg_results(NULL)
        gene_set_scores(list())
        current_gene_set_score(NULL)
        cluster_results(NULL)

        x_range <- range(spots_sf$x)
        y_range <- range(spots_sf$y)
        x_buf   <- diff(x_range) * 0.1
        y_buf   <- diff(y_range) * 0.1

        leafletProxy("map") %>%
          clearGroup("spots") %>%
          clearGroup("drawn") %>%
          clearGroup("selected") %>%
          clearGroup("group1_display") %>%
          clearGroup("group2_display") %>%
          fitBounds(x_range[1] - x_buf, y_range[1] - y_buf,
                    x_range[2] + x_buf, y_range[2] + y_buf) %>%
          addCircleMarkers(
            lng = spots_sf$x, lat = spots_sf$y,
            radius = 5, stroke = TRUE, color = "black", weight = 0.5,
            fillColor = "lightblue", fillOpacity = 0.8,
            group = "spots"
          )

        # ── 8. Cleanup ────────────────────────────────────────────────────────
        unlink(tmp_dir, recursive = TRUE)
        removeNotification(id = "raw_loading")
        output$upload_status <- renderText("✓ 10x Visium data loaded successfully")
        shinyjs::runjs("setTimeout(function() { $('#loading_overlay').removeClass('active'); }, 3000);")
        showNotification(paste0("✓ Loaded ", n_after, " spots from zip"),
                        type = "message", duration = 5)

      }, error = function(e) {
        removeNotification(id = "raw_loading")
        shinyjs::runjs("$('#loading_overlay').removeClass('active');")
        output$upload_status <- renderText(paste("Error:", e$message))
        showNotification(paste("Load error:", e$message), type = "error", duration = 15)
        unlink(tmp_dir, recursive = TRUE)
      })
    })




    # Use Example Data button handler
    observeEvent(input$use_example_data, {
      # Reviewer 3, item 1: actually (re)load the bundled server-side example
      # dataset instead of just showing a message. This makes the example load
      # reliably and lets a user who uploaded their own data return to the demo.
      # Reuses the shared apply_loaded_seurat() loader.
      shinyjs::runjs("
        $('#loading_message').text('Loading example dataset...');
        $('#loading_overlay').addClass('active');
      ")
      showNotification("Loading example dataset...", type = "message",
                       duration = NULL, id = "load_seurat")
      Sys.sleep(0.3)

      tryCatch({
        example_path <- system.file("extdata", "example_visium.rds", package = "SpatialScope")
        if (example_path == "" || !file.exists(example_path)) {
          # source / dev fallback when the package isn't installed
          example_path <- file.path("inst", "extdata", "example_visium.rds")
        }
        if (!file.exists(example_path)) {
          stop("Example dataset not found on the server (expected inst/extdata/example_visium.rds).")
        }

        new_seurat <- readRDS(example_path)
        new_seurat <- tryCatch(UpdateSeuratObject(new_seurat), error = function(e) new_seurat)

        if (!inherits(new_seurat, "Seurat")) stop("Example dataset is not a valid Seurat object.")
        if (length(new_seurat@images) == 0)  stop("Example dataset has no spatial images.")

        current_sample_name("Example_Visium")
        apply_loaded_seurat(new_seurat)   # swaps in session data, refreshes UI, redraws map
      }, error = function(e) {
        removeNotification(id = "load_seurat")
        shinyjs::runjs("$('#loading_overlay').removeClass('active');")
        showNotification(paste("Error loading example data:", e$message),
                         type = "error", duration = 10)
      })
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
      reset_groups()  # Clear all ROI groups
      session$sendCustomMessage("clearFreehandDrawings", list())
      leafletProxy("map") %>%
        clearGroup("drawn") %>%
        clearGroup("selected") %>%
        clearGroup("group1_display") %>%
        clearGroup("group2_display")
      showNotification("Selection and all groups cleared", type = "message", duration = 2)
    })


    ## celltype deconvolution


    # ── Built-in reference paths ───────────────────────────────────────────────────
    # Place your .rds files in a data/ subfolder next to app.R
    
    builtin_refs <- list(
      crc = {
        # Works when installed as a package
        pkg_path <- system.file("extdata", "CRC_reference_RCTD.rds", package = "SpatialScope")
        if (nchar(pkg_path) > 0 && file.exists(pkg_path)) {
          pkg_path
        } else {
          # Fallback for running app.R directly during development
          # app.R is at inst/app/app.R, so extdata is one level up
          file.path(dirname(getwd()), "extdata", "CRC_reference_RCTD.rds")
        }
      }
    )

    # ── Reactive state ─────────────────────────────────────────────────────────────
    ref_seurat   <- reactiveVal(NULL)   # raw Seurat reference
    rctd_ref_cache <- reactiveVal(NULL) # cached spacexr::Reference object
                                        # invalidated when ref_seurat changes

    # Helper flag: is any reference currently loaded?

    output$ref_loaded <- reactive({ 
      !is.null(ref_seurat()) || !is.null(rctd_ref_cache()) 
    })
    outputOptions(output, "ref_loaded", suspendWhenHidden = FALSE)

    # ── Shared helper: ingest any loaded reference Seurat object ───────────────────
    ingest_reference <- function(ref, source_label) {
      if (!inherits(ref, "Seurat")) {
        output$ref_status <- renderText("✗ Error: File is not a Seurat object.")
        return()
      }

      ref_seurat(ref)
      rctd_ref_cache(NULL)   # invalidate cache whenever reference changes

      meta_cols <- colnames(ref@meta.data)
      # updateSelectInput(session, "ref_celltype_col",
      #                   choices = c("Idents (default)" = "idents", meta_cols))

      n_cells  <- ncol(ref)
      n_types  <- length(unique(Idents(ref)))
      output$ref_status <- renderText(
        paste0("✓ ", source_label, " loaded: ",
              n_cells, " cells, ", n_types, " cell types (by Idents)")
      )
    }

    # ── Load built-in reference ────────────────────────────────────────────────────
    # ── Load built-in reference (already a spacexr::Reference object) ─────────────
    observeEvent(input$load_builtin_ref, {
      path <- builtin_refs[[input$builtin_ref_choice]]

      if (is.null(path) || !file.exists(path)) {
        output$ref_status <- renderText(
          paste0("✗ Built-in reference file not found: ", path)
        )
        return()
      }

      # Diagnostic guard (Reviewer 1 #6): a valid reference is ~20 MB. A tiny file
      # here means the .rds did not deploy correctly — typically a Git-LFS pointer
      # stub or a truncated copy — which makes readRDS fail with the cryptic
      # "unknown input format". Catch that case and report the real cause.
      if (file.size(path) < 1e6) {
        output$ref_status <- renderText(paste0(
          "✗ Reference file looks incomplete (", file.size(path),
          " bytes; expected ~20 MB). It may not have deployed correctly ",
          "(e.g., a Git-LFS pointer stub). Re-deploy the full CRC_reference_RCTD.rds."
        ))
        return()
      }

      output$ref_status <- renderText("⏳ Loading built-in reference...")
      tryCatch({
        ref <- readRDS(path)

        if (inherits(ref, "Reference")) {
          # Already a spacexr::Reference — cache it directly, skip the build step
          rctd_ref_cache(ref)

          n_cells <- ncol(ref@counts)
          n_types <- length(unique(ref@cell_types))
          output$ref_status <- renderText(
            paste0("✓ CRC reference loaded: ", n_cells, " cells, ",
                  n_types, " cell types")
          )
          # Populate the cell type dropdown for display purposes
          updateSelectInput(session, "ref_celltype_col",
                            choices = c("Pre-built (spacexr::Reference)" = "prebuilt"))

        } else if (inherits(ref, "Seurat")) {
          # It's a Seurat object — go through normal ingest
          ingest_reference(ref, "Built-in CRC reference")

        } else {
          output$ref_status <- renderText(
            paste0("✗ Unrecognized format: ", class(ref),
                  ". Expected a Seurat or spacexr::Reference object.")
          )
        }
      }, error = function(e) {
        # Make the classic readRDS failure ("unknown input format") actionable.
        detail <- if (grepl("unknown input format|magic number|not a|corrupt",
                            e$message, ignore.case = TRUE)) {
          " — the reference file failed to load; it may not have deployed correctly (expected a ~20 MB RDS, not an LFS pointer)."
        } else ""
        output$ref_status <- renderText(paste0("✗ Error: ", e$message, detail))
      })
    })

    # ── Load user-uploaded reference ───────────────────────────────────────────────
    observeEvent(input$upload_reference, {
      req(input$upload_reference)
      output$ref_status <- renderText("⏳ Loading uploaded reference...")
      tryCatch({
        ref <- readRDS(input$upload_reference$datapath)
        ingest_reference(ref, "Uploaded reference")
      }, error = function(e) {
        output$ref_status <- renderText(paste("✗ Error:", e$message))
      })
    })

    # ── Run RCTD ───────────────────────────────────────────────────────────────────
    deconv_results       <- reactiveVal(NULL)
    deconv_overlap_msg   <- reactiveVal(NULL)   # gene overlap info for the UI

    output$deconv_results_available <- reactive({ !is.null(deconv_results()) })
    outputOptions(output, "deconv_results_available", suspendWhenHidden = FALSE)

    observeEvent(input$run_deconv, {
      if (is.null(ref_seurat()) && is.null(rctd_ref_cache())) {
        showNotification("Please load a reference first (built-in or upload).", type = "error")
        return()
      }

      g1 <- group1_spots()
      g2 <- group2_spots()


      if (!requireNamespace("spacexr", quietly = TRUE)) {
        showNotification(
          "spacexr not installed. Run: devtools::install_github('dmcable/spacexr')",
          type = "error", duration = 15
        )
        return()
      }

      # Diagnostic guard (Reviewer 1 #6 / Reviewer 3 #3): the deployed spacexr can
      # be a wrong/old version or a partial install that does not export the
      # functions we call. Detect that up front and report the actual problem
      # instead of surfacing a cryptic "'Reference' is not an exported object"
      # crash mid-run. The real fix is pinning spacexr on the server (renv).
      spacexr_exports <- getNamespaceExports("spacexr")
      missing_api <- setdiff(c("Reference", "SpatialRNA", "create.RCTD", "run.RCTD"),
                             spacexr_exports)
      if (length(missing_api) > 0) {
        removeNotification(id = "rctd_running")
        msg <- paste0(
          "Incompatible spacexr version on this server (", packageVersion("spacexr"),
          "): missing ", paste(missing_api, collapse = ", "),
          ". Reinstall a compatible build, e.g. ",
          "remotes::install_github('dmcable/spacexr')."
        )
        showNotification(msg, type = "error", duration = 20)
        output$ref_status <- renderText(paste("✗", msg))
        return()
      }

      showNotification("Running RCTD deconvolution... this may take a minute.",
                      type = "message", duration = NULL, id = "rctd_running")

      tryCatch({

        # ── Build (or reuse cached) spacexr::Reference ─────────────────────────
        if (is.null(rctd_ref_cache())) {
          # Only runs for user-uploaded Seurat objects
          ref <- ref_seurat()

          if (is.null(ref)) {
            stop("No reference loaded. Please load a built-in or upload a reference first.")
          }
          cell_types <- Idents(ref)

          cell_types <- as.factor(cell_types)

          # Drop cell types with < 25 cells (RCTD hard requirement)
          type_counts <- table(cell_types)
          valid_types <- names(type_counts[type_counts >= 25])
          keep_cells  <- names(cell_types)[cell_types %in% valid_types]

          if (length(valid_types) < 2) {
            stop("Need ≥2 cell types with ≥25 cells each in the reference.")
          }

          dropped <- setdiff(levels(cell_types), valid_types)
          if (length(dropped) > 0) {
            showNotification(
              paste0("Dropped ", length(dropped), " rare cell type(s) with <25 cells: ",
                    paste(dropped, collapse = ", ")),
              type = "warning", duration = 10
            )
          }

          ref_counts <- Seurat::GetAssayData(ref[, keep_cells], layer = "counts")
          ref_counts <- round(ref_counts)
          cell_types <- droplevels(cell_types[keep_cells])

          rctd_ref_cache(spacexr::Reference(ref_counts, cell_types))
        }

        # At this point rctd_ref_cache() is always populated
        # (either just built from Seurat, or pre-loaded from built-in)
        rctd_ref <- rctd_ref_cache()

        # ── Determine spatial assay ────────────────────────────────────────────
        spatial_assay <- if ("Spatial" %in% names(seurat_obj@assays)) {
          "Spatial"
        } else {
          DefaultAssay(seurat_obj)
        }

        # ── Run per group ──────────────────────────────────────────────────────
        results      <- list()
        overlap_msgs <- c()

        for (grp in c("group1", "group2", "all")) {
          spots <- if (grp == "group1") g1 else if (grp == "group2") g2 else spots_sf$spot_id
          if (length(spots) == 0) next

          grp_label <- if (grp == "group1") "Group 1" else if (grp == "group2") "Group 2" else "All Spots"

          if (length(spots) > 800) {
            showNotification(
              paste0(grp_label, ": ", length(spots),
                    " spots selected — RCTD may take 2-4 minutes."),
              type = "warning", duration = 8
            )
          }

          # Spatial counts
          sp_counts <- Seurat::GetAssayData(seurat_obj,
                                            assay = spatial_assay,
                                            layer = "counts")[, spots, drop = FALSE]
          sp_counts <- round(sp_counts)

          # ── FIX 1: Gene intersection ───────────────────────────────────────
          ref_genes    <- rownames(rctd_ref@counts)
          common_genes <- intersect(ref_genes, rownames(sp_counts))

          overlap_pct <- round(length(common_genes) / length(ref_genes) * 100, 1)
          overlap_msgs <- c(overlap_msgs,
                            paste0(grp_label, ": ", length(common_genes),
                                  " / ", length(ref_genes),
                                  " reference genes matched (", overlap_pct, "%)"))

          if (length(common_genes) < 100) {
            stop(paste0(grp_label, ": Only ", length(common_genes),
                        " genes overlap between reference and spatial data. ",
                        "Check that both use the same species and gene symbol convention."))
          }

          # Subset both to common genes
          sp_counts_sub  <- sp_counts[common_genes, , drop = FALSE]
          ref_counts_sub <- rctd_ref@counts[common_genes, , drop = FALSE]
          # DEBUG - check console output


          # Rebuild reference with matched genes only
          rctd_ref_sub <- spacexr::Reference(ref_counts_sub, rctd_ref@cell_types)
          message("Reference built OK")

          # ── FIX 2: Coordinates with nUMI ──────────────────────────────────
          coords_all <- Seurat::GetTissueCoordinates(seurat_obj)
          x_col <- intersect(c("imagecol", "x"), colnames(coords_all))[1]
          y_col <- intersect(c("imagerow", "y"), colnames(coords_all))[1]
          coords <- coords_all[spots, c(x_col, y_col)]
          colnames(coords) <- c("x", "y")

          coords <- as.data.frame(coords)
          coords$x <- as.numeric(coords$x)
          coords$y <- as.numeric(coords$y)
          rownames(coords) <- spots

                 
          # Force alignment
          sp_counts_sub <- sp_counts_sub[, rownames(coords), drop = FALSE]

          # Convert to plain matrix — some spacexr versions don't handle dgCMatrix well
          sp_counts_mat <- as.matrix(sp_counts_sub)
          nUMI_vec <- colSums(sp_counts_mat)

          sp_obj <- spacexr::SpatialRNA(coords,
                                        sp_counts_mat,
                                        nUMI = nUMI_vec)
          message("SpatialRNA built OK")

          # Run RCTD
          rctd_obj <- spacexr::create.RCTD(sp_obj, rctd_ref_sub,
                                            max_cores = 1,
                                            CELL_MIN_INSTANCE = 25)
                              
          rctd_obj <- spacexr::run.RCTD(rctd_obj, doublet_mode = "full")



          weights_mat <- rctd_obj@results$weights

          if (!is.matrix(weights_mat)) {
            weights_mat <- as.matrix(weights_mat)
          }
          props <- spacexr::normalize_weights(weights_mat)
          results[[grp]] <- as.data.frame(as.matrix(props))
        }

        # Warn if gene overlap is low
        for (msg in overlap_msgs) {
          pct <- as.numeric(gsub(".*\\((.*)%\\).*", "\\1", msg))
          if (!is.na(pct) && pct < 20) {
            showNotification(
              paste0("Low gene overlap (", pct, "%) — results may be less accurate."),
              type = "warning", duration = 8
            )
          }
        }

        deconv_results(results)
        deconv_overlap_msg(paste(overlap_msgs, collapse = "\n"))

        removeNotification(id = "rctd_running")
        showNotification("✓ Deconvolution complete!", type = "message")

      }, error = function(e) {
        removeNotification(id = "rctd_running")
        showNotification(paste("RCTD Error:", e$message), type = "error", duration = 15)
      })
    })

    # ── Gene overlap info display ──────────────────────────────────────────────────
    output$deconv_gene_overlap_msg <- renderText({
      req(deconv_overlap_msg())
      deconv_overlap_msg()
    })

    # ── Barplot ────────────────────────────────────────────────────────────────────

    output$deconv_barplot <- renderPlot({
      req(deconv_results())
      grp   <- input$deconv_group
      props <- deconv_results()[[grp]]
      if (is.null(props)) {
        plot.new()
        text(0.5, 0.5, "No results for this group.", cex = 1.2, col = "grey50")
        return()
      }

      # Build per-spot stacked bar (sample up to 50 spots for readability)
      props_mat <- as.matrix(props)
      if (nrow(props_mat) > 50) {
        props_mat <- props_mat[sample(nrow(props_mat), 50), ]
      }

      df <- reshape2::melt(props_mat, varnames = c("Spot", "CellType"),
                          value.name = "Proportion")

      ggplot(df, aes(x = Spot, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(
          title = paste("Cell Type Proportions —",
              switch(grp, "group1" = "Group 1", "group2" = "Group 2", "all" = "All Spots")),              
          x = NULL, y = "Proportion", fill = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom",
              legend.text = element_text(size = 9),
              plot.title = element_text(size = 12, face = "bold"))
    })


    # ── Table ──────────────────────────────────────────────────────────────────────
    output$deconv_table <- renderTable({
      req(deconv_results())
      grp   <- input$deconv_group
      props <- deconv_results()[[grp]]
      if (is.null(props)) return(NULL)

      avg_props <- sort(colMeans(props, na.rm = TRUE), decreasing = TRUE)
      data.frame(
        Cell_Type  = names(avg_props),
        Proportion = paste0(round(avg_props * 100, 1), "%")
      )
    }, rownames = FALSE)

    # ── Download ───────────────────────────────────────────────────────────────────
    output$dl_deconv <- downloadHandler(
      filename = function() paste0("deconvolution_", input$deconv_group, ".csv"),
      content  = function(file) {
        props <- deconv_results()[[input$deconv_group]]
        req(props)
        write.csv(props, file)
      }
    )



        # ══════════════════════════════════════════════════════════════════════
    ## Independent LR Scoring (no RCTD)
    # ══════════════════════════════════════════════════════════════════════
    lr_solo_results      <- reactiveVal(NULL)
    lr_solo_score_matrix <- reactiveVal(NULL)
    lr_solo_status_msg   <- reactiveVal(NULL)

    output$lr_solo_results_available <- reactive({ !is.null(lr_solo_results()) })
    outputOptions(output, "lr_solo_results_available", suspendWhenHidden = FALSE)


    # Populate the LR region picker from ROIs or Groups (Reviewer 1, item 1).
    observe({
      choices <- if (identical(input$lr_entity_type, "ROIs")) roi_names() else group_names()
      updateSelectizeInput(session, "lr_solo_target", choices = choices,
                           selected = isolate(input$lr_solo_target))
    })

    observeEvent(input$run_lr_solo, {
      req(seurat_obj)

      grp   <- input$lr_solo_target
      spots <- if (identical(input$lr_entity_type, "ROIs")) {
        rr <- rois()[[grp]]; if (is.null(rr)) character(0) else rr$spots
      } else {
        group_spots_of(grp)
      }

      if (length(spots) == 0) {
        showNotification(paste0("'", grp, "' has no spots — pick a non-empty ROI/group."), type = "error")
        return()
      }

      # load LR database

      lrpair <- lr_network()
            if (is.null(lrpair) || nrow(lrpair) == 0) {
              showNotification("L-R database not found or empty.", type = "error")
              return()
            }     

      showNotification("Running LR scoring...", id = "lr_solo_running", duration = NULL)

      withProgress(message = "LR Scoring", value = 0, {
        tryCatch({

          # ── Subset Seurat ──────────────────────────────────────────────
          spatial_assay <- if ("Spatial" %in% names(seurat_obj@assays)) "Spatial" else
                          if ("SCT"     %in% names(seurat_obj@assays)) "SCT" else
                          DefaultAssay(seurat_obj)

          # Avoid full subset which triggers VisiumV1 image validation
          # Instead subset assay and coordinates manually
          seurat_sub <- tryCatch({
            subset(seurat_obj, cells = spots)
          }, error = function(e) {
            # Fallback: update object first then subset
            tryCatch({
              updated <- UpdateSeuratObject(seurat_obj)
              subset(updated, cells = spots)
            }, error = function(e2) {
              # Last resort: subset only the assay, skip image
              seurat_obj[, spots]
            })
          })


          DefaultAssay(seurat_sub) <- spatial_assay
          # FIXED — works for both Assay (v4) and Assay5 (v5)
          counts_mat <- tryCatch(
            GetAssayData(seurat_sub, assay = spatial_assay, layer = "counts"),
            error = function(e)
              GetAssayData(seurat_sub, assay = spatial_assay, slot  = "counts")
          )
          counts_mat <- as(counts_mat, "CsparseMatrix")  # ensure sparse matrix format

          avail_genes <- rownames(counts_mat)
          n_spots     <- ncol(counts_mat)
          spot_names  <- colnames(counts_mat)

          # ── Spatial coordinates ────────────────────────────────────────
          coords     <- GetTissueCoordinates(seurat_sub)
          cc         <- colnames(coords)
          row_col    <- if ("imagerow" %in% cc) "imagerow" else cc[1]
          col_col    <- if ("imagecol" %in% cc) "imagecol" else cc[2]
          coords_mat <- as.matrix(coords[spot_names, c(row_col, col_col)])

          # ── Filter LR pairs ────────────────────────────────────────────
          lrpair_avail <- dplyr::filter(lrpair, from %in% avail_genes, to %in% avail_genes)
          incProgress(0.1, detail = paste(nrow(lrpair_avail), "valid L-R pairs"))

          if (nrow(lrpair_avail) == 0) {
            removeNotification(id = "lr_solo_running")
            showNotification("No L-R pairs found in this ROI.", type = "warning")
            return()
          }

          # ── KNN ────────────────────────────────────────────────────────
          k_nbrs <- 12
          knn_nb <- spdep::knn2nb(spdep::knearneigh(coords_mat, k = k_nbrs))
          # Gaussian-kernel bandwidth multiplier on the adaptive sigma (R1 #4).
          bw_mult <- if (is.null(input$lr_bandwidth_mult) || input$lr_bandwidth_mult <= 0) 1 else input$lr_bandwidth_mult

          # ── Spatial lag for ligands ────────────────────────────────────
          unique_L   <- unique(lrpair_avail$from)
          incProgress(0.2, detail = paste("Spatial lag for", length(unique_L), "ligands"))

          L_lag_cache <- list()
          for (gene in unique_L) {
            expr <- as.numeric(counts_mat[gene, ])
            lag  <- numeric(n_spots)
            for (i in seq_len(n_spots)) {
              nb_idx <- knn_nb[[i]]
              if (length(nb_idx) > 0) {
                dists  <- sqrt(rowSums((coords_mat[nb_idx, ] - coords_mat[i, ])^2))
                sigma  <- median(dists) * bw_mult
                w      <- exp(-dists^2 / (2 * sigma^2))
                lag[i] <- (expr[i] + sum(expr[nb_idx] * w)) / (1 + sum(w))
              } else {
                lag[i] <- expr[i]
              }
            }
            L_lag_cache[[gene]] <- lag
          }

          # ── Geometric mean score matrix ────────────────────────────────
          n_pairs    <- nrow(lrpair_avail)
          pair_names <- paste0(lrpair_avail$from, "_x_", lrpair_avail$to)
          incProgress(0.5, detail = paste("Scoring", n_pairs, "pairs"))

          score_mat <- matrix(0, nrow = n_spots, ncol = n_pairs,
                              dimnames = list(spot_names, pair_names))
          for (j in seq_len(n_pairs)) {
            score_mat[, j] <- sqrt(L_lag_cache[[ lrpair_avail$from[j] ]] *
                                  as.numeric(counts_mat[ lrpair_avail$to[j], ]))
          }

          # ── Summary table: LR_Pair + Mean_Score only ──────────────────
          mean_scores <- colMeans(score_mat, na.rm = TRUE)
          summary_df  <- data.frame(
            LR_Pair    = pair_names,
            Ligand     = lrpair_avail$from,
            Receptor   = lrpair_avail$to,
            Mean_Score = round(mean_scores, 4),
            stringsAsFactors = FALSE
          ) |> dplyr::arrange(dplyr::desc(Mean_Score))

          lr_solo_results(summary_df)
          lr_solo_score_matrix(score_mat)

          incProgress(1, detail = "Done")
          removeNotification(id = "lr_solo_running")
          lr_solo_status_msg(paste0("✓ ", nrow(summary_df), " pairs | ",
                                    n_spots, " spots | Top: ",
                                    summary_df$LR_Pair[1]))
          showNotification(paste0("✓ LR scoring complete: ", nrow(summary_df), " pairs"),
                          type = "message")

        }, error = function(e) {
          removeNotification(id = "lr_solo_running")
          showNotification(paste("Error:", e$message), type = "error", duration = 15)
        })
      })
    })

    output$lr_solo_status <- renderText({
      req(lr_solo_status_msg())
      lr_solo_status_msg()
    })

    # ── Simple table ────────────────────────────────────────────────────────
    output$lr_simple_table <- DT::renderDataTable({
      req(lr_solo_results())
      df <- lr_solo_results()[, c("LR_Pair", "Ligand", "Receptor", "Mean_Score")]
      df$Mean_Score <- round(df$Mean_Score, 4)
      df
    }, rownames = FALSE, filter = "none", selection = "single",
      options = list(
        pageLength = 10,
        lengthChange = FALSE,
        searching = TRUE,
        dom = 'ftp',
        scrollX = TRUE,
        columnDefs = list(list(width = '40%', targets = 0),
                          list(width = '20%', targets = 1),
                          list(width = '20%', targets = 2),
                          list(width = '20%', targets = 3))
      ))  



    # ── Show an LR pair on the main map (shared by the button and row-click) ──
    plot_lr_pair <- function(pair) {
      score_mat <- lr_solo_score_matrix()
      if (is.null(score_mat) || is.null(pair) || !pair %in% colnames(score_mat)) {
        showNotification("LR pair not found.", type = "error"); return(invisible(NULL))
      }
      scores    <- score_mat[, pair]
      all_spots <- rownames(seurat_obj@meta.data)
      full_scores <- setNames(rep(NA_real_, length(all_spots)), all_spots)
      full_scores[names(scores)] <- scores

      updateSelectInput(session, "feature_type", selected = "None")
      showing_gene_set(FALSE)
      current_values(full_scores)
      is_categorical(FALSE)
      updateSelectInput(session, "color_scheme", selected = "greyred")
      update_map_colors()
      showNotification(paste0("Showing: ", pair), type = "message", duration = 3)
      invisible(NULL)
    }

    # Reviewer 1, item 5: clicking a row in the LR table drives the spatial plot
    # directly (this replaces the separate dropdown + "Show on Map" button).
    observeEvent(input$lr_simple_table_rows_selected, {
      req(lr_solo_results(), lr_solo_score_matrix())
      row <- input$lr_simple_table_rows_selected
      if (length(row) == 0) return()
      plot_lr_pair(lr_solo_results()$LR_Pair[row])
    })

    # ── Download ────────────────────────────────────────────────────────────
    output$dl_lr_solo <- downloadHandler(
      filename = function() paste0("LR_scores_", input$lr_solo_target, "_",
                                    format(Sys.time(), "%Y%m%d"), ".csv"),
      content  = function(file) {
        req(lr_solo_results())
        write.csv(lr_solo_results(), file, row.names = FALSE)
      }
    )



    # ── Helper: load and filter L-R database ──────────────────────────────────────
    get_lr_pairs <- reactive({
      lrpair <- lr_network()
      if (is.null(lrpair)) return(NULL)
      lrpair
    })

    # # ── Main computation ───────────────────────────────────────────────────────────
    # observeEvent(input$run_lr, {

    #   # ── Guards ─────────────────────────────────────────────────────────────────

    #   req(seurat_obj)


    #   if (is.null(deconv_results())) {
    #     showNotification("Run RCTD first — L-R correlation needs cell-type proportions.", type = "error")
    #     return()
    #   }

    #   grp   <- input$lr_group
    #   spots <- if (grp == "group1") group1_spots() else group2_spots()

    #   if (length(spots) == 0) {
    #     showNotification(paste("No spots saved to", grp), type = "error")
    #     return()
    #   }

    #   rctd_props <- deconv_results()[[grp]]
    #   if (is.null(rctd_props)) {
    #     showNotification(paste("No RCTD results for", grp, "— run RCTD for this group first."),
    #                     type = "error")
    #     return()
    #   }

    #   lrpair <- get_lr_pairs()
    #   if (is.null(lrpair) || nrow(lrpair) == 0) {
    #     showNotification("No L-R pairs loaded.", type = "error")
    #     return()
    #   }

    #   # ── Align spots (RCTD ∩ spatial) ───────────────────────────────────────────
    #   shared_spots <- intersect(spots, rownames(rctd_props))
    #   if (length(shared_spots) < 5) {
    #     showNotification("Too few spots with both spatial and RCTD data (need ≥ 5).", type = "error")
    #     return()
    #   }

    #   showNotification("Running L-R colocalization...", id = "lr_running", duration = NULL)

    #   withProgress(message = "L-R Colocalization", value = 0, {
    #     tryCatch({

    #       # ── Subset Seurat to group spots ────────────────────────────────────────

    #        # adjust name if needed
    #       seurat_lr <- seurat_obj
    #       spatial_assay <- if ("Spatial" %in% names(seurat_lr@assays)) {
    #         "Spatial"
    #       } else if ("SCT" %in% names(seurat_lr@assays)) {
    #         "SCT"
    #       } else {
    #         DefaultAssay(seurat_lr)
    #       }

    #       DefaultAssay(seurat_lr) <- spatial_assay
    #       seurat_sub_lr <- subset(seurat_lr, cells = shared_spots)
    #       counts_mat   <- seurat_sub_lr@assays[[spatial_assay]]@counts
    #       avail_genes  <- rownames(counts_mat)
    #       n_spots      <- ncol(counts_mat)
    #       spot_names   <- colnames(counts_mat)




    #       # ── Spatial coordinates ─────────────────────────────────────────────────
    #       coords     <- GetTissueCoordinates(seurat_sub_lr)
    #       cc         <- colnames(coords)
    #       row_col    <- if ("imagerow" %in% cc) "imagerow" else cc[1]
    #       col_col    <- if ("imagecol" %in% cc) "imagecol" else cc[2]
    #       coords_mat <- as.matrix(coords[, c(row_col, col_col)])
    #       ## keep only within region spots so the neirghborhood won't be outsie the region
    #       coords_mat <- coords_mat[spot_names, , drop = FALSE]

    #       # ── Filter L-R pairs to genes present in this ROI ───────────────────────
    #       lrpair_avail <- dplyr::filter(lrpair,
    #                                     from %in% avail_genes,
    #                                     to   %in% avail_genes)

    #       incProgress(0.1, detail = paste(nrow(lrpair_avail), "valid L-R pairs found"))

    #       if (nrow(lrpair_avail) == 0) {
    #         removeNotification(id = "lr_running")
    #         showNotification("No L-R pairs found where both genes are expressed in this ROI.",
    #                         type = "warning")
    #         return()
    #       }

    #       k_nbrs <- 12

    #       # ── Pre-compute KNN once ─────────────────────────────────────────────────
    #       # (reused for every ligand — big speedup vs computing per pair)
    #       knn_nb <- spdep::knn2nb(spdep::knearneigh(coords_mat, k = k_nbrs))

    #       # ── Cache spatial lag for all unique ligands ─────────────────────────────
    #       unique_L <- unique(lrpair_avail$from)
    #       incProgress(0.15, detail = paste("Spatial lag for", length(unique_L), "ligands"))

    #       L_lag_cache <- list()
    #       for (gene in unique_L) {
    #         expr <- as.numeric(counts_mat[gene, ])
    #         lag  <- numeric(n_spots)
    #         for (i in seq_len(n_spots)) {
    #           nb_idx  <- knn_nb[[i]]
    #           if (length(nb_idx) > 0) {
    #             # w       <- ifelse(seq_along(nb_idx) <= 6, 1, 0.5) ## this is too strict cut-off
    #             ## to make the neighborhood region better, closer neighbors means higher weight and farther neighbors mean lower weight. This is a Gaussian kernel (smooth decay)
    #             dists <- sqrt(rowSums((coords_mat[nb_idx, ] - coords_mat[i, ])^2))
    #             sigma <- median(dists)
    #             w <- exp(-dists^2 / (2 * sigma^2))
    #             lag[i]  <- (expr[i] + sum(expr[nb_idx] * w)) / (1 + sum(w))
    #           } else {
    #             lag[i]  <- expr[i]
    #           }
    #         }
    #         L_lag_cache[[gene]] <- lag
    #       }

    #       # ── Geometric mean score: spot × pair matrix ─────────────────────────────
    #       n_pairs    <- nrow(lrpair_avail)
    #       pair_names <- paste0(lrpair_avail$from, "_x_", lrpair_avail$to)

    #       incProgress(0.3, detail = paste("Scoring", n_pairs, "L-R pairs"))

    #       score_mat <- matrix(0, nrow = n_spots, ncol = n_pairs,
    #                           dimnames = list(spot_names, pair_names))

    #       for (j in seq_len(n_pairs)) {
    #         L_lag  <- L_lag_cache[[ lrpair_avail$from[j] ]]
    #         R_expr <- as.numeric(counts_mat[ lrpair_avail$to[j], ])
    #         score_mat[, j] <- sqrt(L_lag * R_expr)   # geometric mean
    #       }

    #       # ── Correlate spot-level scores with RCTD proportions ────────────────────
    #       incProgress(0.65, detail = "Correlating with RCTD cell types")

    #                 # ── Summarise: one row per pair ──────────────────────────────────────────
    #       # ── Per-gene correlations with cell types (sender / receiver) ─────────────────
    #       incProgress(0.75, detail = "Correlating L and R genes with cell types separately")

    #       rctd_sub   <- rctd_props[spot_names, , drop = FALSE]
    #       cell_types <- colnames(rctd_sub)

    #       # Helper: best-correlated cell type for any numeric vector
    #       best_corr_celltype <- function(vec, rctd_mat, cell_types) {
    #         if (sd(vec, na.rm = TRUE) < 1e-10) return(list(ct = NA, r = NA))
    #         cors <- sapply(cell_types, function(ct) {
    #           cv <- rctd_mat[, ct]
    #           if (sd(cv, na.rm = TRUE) < 1e-10) return(NA)
    #           cor(vec, cv, use = "complete.obs", method = "spearman")
    #         })
    #         cors[cors < 0] <- NA
    #         best <- which.max(cors)
    #         if (length(best) == 0) return(list(ct = NA, r = NA))
    #         list(ct = cell_types[best], r = round(cors[best], 3))
    #       }

    #       # ── Cache per-gene cell-type correlations (reuse across pairs) ────────────────
    #       # REPLACE WITH:
    #       message(paste(slotNames(rctd_ref_cache()), collapse = ", "))

    #       # REPLACE the entire ref_means tryCatch with:
    #       ref_means <- tryCatch({
    #         obj <- rctd_ref_cache()
    #         # Compute mean expression per cell type from reference counts
    #         ct <- obj@cell_types
    #         counts <- obj@counts
    #         cell_type_names <- levels(ct)
    #         mean_mat <- sapply(cell_type_names, function(ctype) {
    #           cells <- names(ct)[ct == ctype]
    #           Matrix::rowMeans(counts[, cells, drop = FALSE])
    #         })
    #         mean_mat  # genes × cell_types matrix
    #       }, error = function(e) {
    #         message("EWA unavailable: ", e$message)
    #         NULL
    #       })
    #       use_ewa <- input$use_ewa && !is.null(ref_means)          

    #       ewa_celltype <- function(gene, props_mat, ref_means) {
    #         if (!gene %in% rownames(ref_means)) return(list(ct = NA, r = NA))
    #         ref_expr     <- ref_means[gene, ]
    #         common_types <- intersect(colnames(props_mat), names(ref_expr))
    #         if (length(common_types) == 0) return(list(ct = NA, r = NA))
    #         W            <- as.matrix(props_mat[, common_types, drop = FALSE])
    #         mu           <- ref_expr[common_types]
    #         numerator    <- sweep(W, 2, mu, "*")
    #         denom        <- rowSums(numerator) + 1e-10
    #         attribution  <- sweep(numerator, 1, denom, "/")
    #         avg_attr     <- colMeans(attribution, na.rm = TRUE)
    #         best         <- which.max(avg_attr)
    #         list(ct = names(avg_attr)[best], r = round(avg_attr[best], 3))
    #       }

    #       L_ct_cache <- list()
    #       for (gene in unique(lrpair_avail$from)) {
    #         L_ct_cache[[gene]] <- if (use_ewa) {
    #           ewa_celltype(gene, rctd_sub, ref_means)
    #         } else {
    #           best_corr_celltype(L_lag_cache[[gene]], rctd_sub, cell_types)
    #         }
    #       }

    #       R_ct_cache <- list()
    #       for (gene in unique(lrpair_avail$to)) {
    #         R_ct_cache[[gene]] <- if (use_ewa) {
    #           ewa_celltype(gene, rctd_sub, ref_means)
    #         } else {
    #           expr <- as.numeric(counts_mat[gene, ])
    #           best_corr_celltype(expr, rctd_sub, cell_types)
    #         }
    #       }         

    #       # ── LR score correlation with cell types ─────────────────────────────────────
    #       incProgress(0.85, detail = "Building summary table")

    #       LR_ct_name <- character(n_pairs)
    #       LR_ct_corr <- numeric(n_pairs)

    #       for (j in seq_len(n_pairs)) {
    #         res <- best_corr_celltype(score_mat[, j], rctd_sub, cell_types)
    #         LR_ct_name[j] <- if (!is.na(res$ct)) res$ct else "—"
    #         LR_ct_corr[j] <- if (!is.na(res$r))  res$r  else NA
    #       }

    #       # ── Build summary table ───────────────────────────────────────────────────────
    #       mean_scores <- colMeans(score_mat, na.rm = TRUE)

    #       summary_df <- data.frame(
    #         LR_Pair      = pair_names,
    #         Ligand       = lrpair_avail$from,
    #         Receptor     = lrpair_avail$to,
    #         Mean_Score   = round(mean_scores, 4),
    #         # Sender: which cell type co-localises with ligand expression
    #         L_CellType   = sapply(lrpair_avail$from, function(g) {
    #           ct <- L_ct_cache[[g]]$ct; if (is.na(ct)) "—" else ct }),
    #         L_Corr       = sapply(lrpair_avail$from, function(g) L_ct_cache[[g]]$r),
    #         # Receiver: which cell type co-localises with receptor expression  
    #         R_CellType   = sapply(lrpair_avail$to, function(g) {
    #           ct <- R_ct_cache[[g]]$ct; if (is.na(ct)) "—" else ct }),
    #         R_Corr       = sapply(lrpair_avail$to, function(g) R_ct_cache[[g]]$r),
    #         # Combined: LR score enrichment
    #         LR_CellType  = LR_ct_name,
    #         LR_Corr      = LR_ct_corr,
    #         stringsAsFactors = FALSE
    #       ) |> dplyr::arrange(dplyr::desc(Mean_Score))



    #       lr_results(summary_df)
    #       lr_score_matrix(score_mat)



    #       incProgress(1, detail = "Done")
    #       removeNotification(id = "lr_running")

    #       top_lr_pairs <- head(summary_df$LR_Pair, 3)
    #       lr_status_msg(
    #         paste0("✓ ", nrow(summary_df), " L-R pairs analyzed | ",
    #               length(shared_spots), " spots | Group: ", grp, "\n",
    #               "Top pairs: ", paste(top_lr_pairs, collapse = ", "))
    #       )

    #       showNotification(paste0("✓ L-R analysis complete: ", nrow(summary_df), " pairs"),
    #                       type = "message")

    #     }, error = function(e) {
    #       removeNotification(id = "lr_running")
    #       showNotification(paste("L-R Error:", e$message), type = "error", duration = 15)
    #     })
    #   })
    # })

    # # ── Status text ────────────────────────────────────────────────────────────────
    # output$lr_status <- renderText({
    #   req(lr_status_msg())
    #   lr_status_msg()
    # })

    # # ── Results table ──────────────────────────────────────────────────────────────

    # output$lr_table <- DT::renderDataTable({
    #   req(lr_results())
    #   df <- head(lr_results(), 20)
    #   df$Mean_Score <- round(df$Mean_Score, 4)
    #   df$L_Corr     <- round(df$L_Corr, 3)
    #   df$R_Corr     <- round(df$R_Corr, 3)
    #   df$LR_Corr    <- round(df$LR_Corr, 3)
    #   df[, c("LR_Pair","Mean_Score","L_CellType","L_Corr",
    #         "R_CellType","R_Corr","LR_CellType","LR_Corr")]
    # }, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 20))



    # # ── Download ───────────────────────────────────────────────────────────────────
    # output$dl_lr_results <- downloadHandler(
    #   filename = function() paste0("LR_colocalization_", input$lr_group, "_",
    #                                 format(Sys.time(), "%Y%m%d"), ".csv"),
    #   content  = function(file) {
    #     req(lr_results())
    #     write.csv(lr_results(), file, row.names = FALSE)
    #   }
    # )

    # # # ── Simple 2-column LR table ───────────────────────────────────────────
    # # output$lr_simple_table <- DT::renderDataTable({
    # #   req(lr_results())
    # #   df <- lr_results()[, c("LR_Pair", "Mean_Score")]
    # #   df$Mean_Score <- round(df$Mean_Score, 4)
    # #   df
    # # }, rownames = FALSE,
    # #    options  = list(scrollX = TRUE, pageLength = 10, scrollY = "200px"))

    # # ── Update LR pair selector when results are ready ─────────────────────
    # observeEvent(lr_results(), {
    #   req(lr_results())
    #   updateSelectInput(session, "lr_pair_viz", choices = lr_results()$LR_Pair)
    # })



    # # ── Show LR score on main map (same as gene set flow) ──────────────────
    # observeEvent(input$show_lr_on_map, {
    #   req(lr_score_matrix(), input$lr_pair_viz)

    #   score_mat  <- lr_score_matrix()
    #   pair       <- input$lr_pair_viz

    #   if (!pair %in% colnames(score_mat)) {
    #     showNotification("LR pair not found.", type = "error")
    #     return()
    #   }

    #   # scores vector named by spot
    #   scores     <- score_mat[, pair]
    #   all_spots  <- rownames(seurat_obj@meta.data)

    #   # fill NA for spots not in the LR analysis
    #   full_scores <- setNames(rep(NA_real_, length(all_spots)), all_spots)
    #   full_scores[names(scores)] <- scores

    #   # ── follow gene set pattern: update current_values + color scheme ──
    #   updateSelectInput(session, "feature_type", selected = "None")  # clear gene/meta
    #   showing_gene_set(FALSE)

    #   current_values(full_scores)
    #   is_categorical(FALSE)

    #   # apply the chosen color scheme from the LR section
    #   # temporarily override color_scheme used by update_map_colors
    #   # by storing it and using isolate
    #   selected_scheme <- isolate(input$lr_color_scheme)
    #   updateSelectInput(session, "color_scheme", selected = selected_scheme)

    #   update_map_colors()

    #   showNotification(paste0("Showing LR score: ", pair), type = "message", duration = 3)
    # })


    # # ── Spatial visualization of selected LR pair score ────────────────────
    # lr_spatial_plot <- reactiveVal(NULL)


    # Clustering
    observeEvent(input$run_clustering, {
      showNotification("Running clustering...", type = "message", duration = NULL, id = "clustering_run")

      tryCatch({
        # Determine which spots to cluster (All spots / an ROI / a group).
        region_key <- input$cluster_region
        spots_to_cluster <- if (is.null(region_key) || identical(region_key, "__all__")) {
          spots_sf$spot_id
        } else {
          region_spots(region_key)
        }

        # Validate that we have spots to cluster
        if (length(spots_to_cluster) == 0) {
          stop("No spots in the selected region — pick a non-empty ROI/group, or 'All spots'.")
        }

        spatial_obj <- subset(seurat_obj, cells = spots_to_cluster)

        # Get the default assay name dynamically
        default_assay <- DefaultAssay(spatial_obj)

        # Check if data needs preprocessing (works with any assay name: Spatial, RNA, etc.)
        # FIXED — works for both Assay (v4) and Assay5 (v5)
        if (!"SCT" %in% names(spatial_obj@assays)) {
          needs_preprocessing <- tryCatch({
            # v4: check scale.data and var.features slots
            current_assay <- spatial_obj@assays[[default_assay]]
            scale_data   <- slot(current_assay, "scale.data")
            var_features <- slot(current_assay, "var.features")
            is.null(scale_data) || length(scale_data) == 0 ||
              is.null(var_features) || length(var_features) == 0
          }, error = function(e) {
            # v5 Assay5: check layers and VariableFeatures instead
            layers <- names(spatial_obj@assays[[default_assay]]@layers)
            vf     <- VariableFeatures(spatial_obj)
            !("scale.data" %in% layers) || length(vf) == 0
          })

          if (needs_preprocessing) {
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
          spot_selection = if (is.null(region_key) || identical(region_key, "__all__")) "all" else region_label(region_key),
          n_spots = length(spots_to_cluster)
        ))

        removeNotification(id = "clustering_run")

        selection_text <- if (is.null(region_key) || identical(region_key, "__all__")) "all spots" else region_label(region_key)

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

    # Keep the clustering region picker in sync with ROIs/Groups.
    observe({
      updateSelectizeInput(session, "cluster_region",
                           choices = c("All spots" = "__all__", region_choices()),
                           selected = isolate(input$cluster_region))
    })

    # Show spot count for the selected clustering region.
    output$cluster_spot_count <- renderText({
      key <- input$cluster_region
      spot_count <- if (is.null(key) || identical(key, "__all__")) nrow(spots_sf) else length(region_spots(key))
      if (spot_count == 0) {
        "⚠️ Selected region has no spots — draw & save an ROI, or pick 'All spots'."
      } else {
        paste("✓", spot_count, "spots will be clustered")
      }
    })

    output$cluster_umap <- renderPlot({
      cluster_res <- cluster_results()
      if (!is.null(cluster_res)) {

        p <- DimPlot(cluster_res$seurat, reduction = "umap", label = TRUE,
                pt.size = 0.5,
                cols = if(!is.null(cluster_colors_palette())) cluster_colors_palette() else NULL) +
          ggtitle(paste("UMAP - Resolution:", cluster_res$resolution))
        cluster_umap_rv(p)
        p
      }
    })
    output$dl_cluster <- downloadHandler(

      filename = function() {
        res <- cluster_results()
        grp <- if(is.null(res)) "unknown" else res$spot_selection
        paste0("cluster_assignments_", grp, "_res", 
              cluster_results()$resolution, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        res <- cluster_results()
        req(res)
        df <- data.frame(
          spot_id = names(res$clusters),
          cluster = as.integer(res$clusters)
        )
        df <- df[order(df$cluster), ]   # ← sort by cluster number
        write.csv(df, file, row.names = FALSE)
      }
    )    

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
          cluster_colors <- RColorBrewer::brewer.pal(min(12, max(3, n_clusters)), "Set3")[1:n_clusters]
        } else {
          cluster_colors <- rainbow(n_clusters, s = 1, v = 0.9)
        }
        names(cluster_colors) <- cluster_levels

        cluster_colors_palette(setNames(cluster_colors, cluster_levels))

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
        message("Legend colors: ", paste(cluster_colors, collapse=", "))
        leafletProxy("map") %>%
          clearGroup("spots") %>%
          clearControls() %>% 
          addCircleMarkers(
            lng = spots_sf$x, lat = spots_sf$y,
            radius = input$spot_size, stroke = TRUE, color = "#333333", weight = 0.5,
            fillColor = unname(spot_colors), fillOpacity = 1, group = "spots"
          )%>%
          
          addLegend(                                   # ← add legend after markers
              position = "bottomright",
              colors = as.character(unname(cluster_colors)),
              labels = paste("Cluster", cluster_levels),
              title = "Clusters",
              opacity = 0.8
            )            
      } else {
        leafletProxy("map") %>% clearControls()
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


    observeEvent(input$load_pathway, {
      req(input$pathway_library != "None")
      genes <- current_hallmark_library()[[input$pathway_library]]$genes
      updateTextAreaInput(session, "gene_set_input", value = paste(genes, collapse = "\n"))
      updateTextInput(session, "gene_set_name", value = gsub("Hallmark: ", "", input$pathway_library))
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
        expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
        scores <- rowMeans(expr_data, na.rm = TRUE)

      } else if (method == "addmodulescore") {
        temp_seurat <- AddModuleScore(seurat_obj, features = list(genes_in_data),
                                      name = "GeneSet", assay = DefaultAssay(seurat_obj))
        scores <- temp_seurat$GeneSet1

      } else if (method == "gsva") {
        # VECTORIZED GSVA ALTERNATIVE: Fast rank-based enrichment score
        # This uses matrix operations to rank all cells at once (10-50x faster than loop)

        # Get expression matrix
        expr_mat <- tryCatch(
          GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), layer = "data"),
          error = function(e) 
          GetAssayData(seurat_obj, assay = DefaultAssay(seurat_obj), slot = "data")
        )       

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

        showNotification(
          gsub("\n", " | ", msg),
          type = "message", duration = 10
        )

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
        gene_expr <- FetchData(seurat_obj, vars = input$gene_select, layer = "data")
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


    # Redraw the Group 1 / Group 2 overlay so it PERSISTS across map redraws.
    # Clears any existing overlay first, then re-adds it only when the
    # "Show Groups on Map" checkbox is enabled. Single source of truth shared by
    # the checkbox observer and update_map_colors() (Reviewer 1, item 2: the
    # overlay should stay visible until the user explicitly unchecks the box,
    # instead of vanishing whenever the gene / palette / spot size changes).
    # Convert stored drawn features (GeoJSON polygons from the freehand tool)
    # into a list of {lng, lat} boundary rings, in the map's spot coordinate
    # space (Reviewer 1, item 3).
    extract_roi_rings <- function(feats) {
      rings <- list()
      for (feat in feats) {
        if (isTRUE(feat$type == "Feature") &&
            !is.null(feat$geometry) && identical(feat$geometry$type, "Polygon")) {
          pc <- feat$geometry$coordinates[[1]]
          if (length(pc) >= 2) {
            rings[[length(rings) + 1]] <- list(
              lng = vapply(pc, function(p) as.numeric(p[[1]]), numeric(1)),
              lat = vapply(pc, function(p) as.numeric(p[[2]]), numeric(1))
            )
          }
        }
      }
      rings
    }

    # Draw each ROI ring as a high-contrast contour: a white casing underneath
    # with a black dashed line on top, so it reads over both dark and bright
    # expression (Reviewer 1, item 3). Fixed style, no user customization.
    add_roi_contour <- function(proxy, rings) {
      for (r in rings) {
        if (length(r$lng) >= 2) {
          proxy <- proxy %>%
            addPolylines(lng = r$lng, lat = r$lat, color = "white",
                         weight = 5, opacity = 1, group = "roi_contour") %>%
            addPolylines(lng = r$lng, lat = r$lat, color = "black",
                         weight = 2, opacity = 1, dashArray = "6,8",
                         group = "roi_contour")
        }
      }
      proxy
    }

    draw_group_overlay <- function() {
      proxy <- leafletProxy("map") %>%
        clearGroup("group_display") %>%
        clearGroup("group1_display") %>%   # legacy layer names (harmless to clear)
        clearGroup("group2_display") %>%
        clearGroup("roi_contour")

      r  <- rois()
      nm <- names(r)
      # Optional focus: show all ROIs, one ROI, or the ROIs of one group.
      focus <- input$roi_show_filter
      show_i <- if (is.null(focus) || identical(focus, "__all__")) {
        seq_along(nm)
      } else if (startsWith(focus, "roi:")) {
        which(nm == sub("^roi:", "", focus))
      } else if (startsWith(focus, "grp:")) {
        members <- groups()[[sub("^grp:", "", focus)]]$members
        which(nm %in% members)
      } else {
        seq_along(nm)
      }

      # ROI member dots ("Show ROIs on Map") — one colored layer per named ROI,
      # so every saved region is visible in its own color (Reviewer 1, item 1).
      if (isTRUE(input$show_groups) && length(nm) > 0) {
        for (i in show_i) {
          spots <- r[[nm[i]]]$spots
          if (length(spots) == 0) next
          idx <- which(spots_sf$spot_id %in% spots)
          if (length(idx) == 0) next
          col <- group_color(i)
          proxy <- proxy %>%
            addCircleMarkers(
              lng = spots_sf$x[idx], lat = spots_sf$y[idx],
              radius = 4, stroke = TRUE, color = col, weight = 2,
              opacity = if (isTRUE(input$transparent_groups)) 0.6 else 1,
              fillColor = col,
              fillOpacity = if (isTRUE(input$transparent_groups)) 0 else 0.8,
              group = "group_display"
            )
        }
      }

      # ROI boundary contours ("Show ROI contours") — every saved ROI, high-contrast,
      # re-applied on every redraw so they persist (Reviewer 1, item 3).
      if (isTRUE(input$show_roi_contours) && length(nm) > 0) {
        for (i in show_i) {
          proxy <- add_roi_contour(proxy, r[[nm[i]]]$rings)
        }
      }

      invisible(NULL)
    }

    update_map_colors <- function() {
      values <- current_values()
      spot_radius <- input$spot_size

      # NOTE (Reviewer 1, item 2): the group overlay is intentionally NOT
      # force-cleared here and the "Show Groups on Map" checkbox is NOT reset.
      # The overlay is re-applied at the end via draw_group_overlay() so it
      # persists across gene / palette / spot-size changes.

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
          radius = input$spot_size,
          stroke = TRUE, color = "#333333", weight = 0.5,
          fillColor = colors, # ← ADD unname() here
          fillOpacity = 1,
          group = "spots"
        )

      # Re-apply the group overlay on top of the freshly drawn spots so it
      # persists across this redraw (Reviewer 1, item 2).
      draw_group_overlay()
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

            breaks = if (val_range[1] == val_range[2]) 
                      val_range[1] 
                    else 
                      seq(val_range[1], val_range[2], length.out = 5),
            labels = if (val_range[1] == val_range[2]) 
                      round(val_range[1], 2) 
                    else 
                      round(seq(val_range[1], val_range[2], length.out = 5), 2)

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
    # ── ROI / Group management (two-tier model, Reviewer 1, item 1) ───────────
    # ROIs are first-class named regions; Groups are named collections of ROIs.

    # Save the currently drawn region as a named ROI.
    observeEvent(input$save_roi_btn, {
      sel <- selected_spots()
      if (length(sel) == 0) {
        showNotification("Draw a region on the map first.", type = "warning"); return()
      }
      name <- input$roi_name
      if (is.null(name) || trimws(name) == "") name <- paste0("ROI ", roi_counter() + 1)
      name <- trimws(name)
      if (!is.null(rois()[[name]])) {
        showNotification(paste0("An ROI named '", name, "' already exists."), type = "warning"); return()
      }
      save_roi(name, sel, extract_roi_rings(drawn_feats()))
      roi_counter(roi_counter() + 1)
      drawn_feats(list()); selected_spots(character(0))
      session$sendCustomMessage("clearFreehandDrawings", list())
      leafletProxy("map") %>% clearGroup("drawn") %>% clearGroup("selected")
      updateTextInput(session, "roi_name", value = paste0("ROI ", roi_counter() + 1))
      showNotification(paste0("Saved region '", name, "' (", length(sel), " spots)."), type = "message")
    })

    # Saved-ROIs list: color swatch, name, spot count, remove.
    # Saved ROIs live in a collapsible dropdown so the panel stays compact even
    # after many regions are drawn (expand to review / remove).
    output$roi_chips <- renderUI({
      r <- rois(); nm <- names(r)
      if (length(nm) == 0)
        return(tags$div(style = "font-size:12px; color:#888;", "No regions saved yet."))
      tags$details(style = "font-size:12px;",
        tags$summary(style = "cursor:pointer; font-weight:600; color:#2c3e50;",
                     sprintf("Saved ROIs (%d)", length(nm))),
        tags$div(style = "margin-top:6px; display:flex; flex-direction:column; gap:4px; max-height:140px; overflow-y:auto;",
          lapply(seq_along(nm), function(i) {
            cnt <- length(r[[nm[i]]]$spots); col <- group_color(i)
            tags$div(style = "display:flex; align-items:center; gap:6px;",
                     tags$span(style = sprintf("display:inline-block; width:12px; height:12px; border-radius:3px; background:%s;", col)),
                     tags$span(style = "flex:1;", sprintf("%s  (%d spots)", nm[i], cnt)),
                     tags$button("\u2715", title = "Remove ROI",
                                 style = "border:none; background:none; cursor:pointer; color:#c0392b; font-size:13px;",
                                 onclick = sprintf("Shiny.setInputValue('remove_roi','%s',{priority:'event'});", nm[i])))
          })
        )
      )
    })

    observeEvent(input$remove_roi, {
      r <- rois(); r[[input$remove_roi]] <- NULL; rois(r)
      g <- groups()
      for (gn in names(g)) g[[gn]]$members <- setdiff(g[[gn]]$members, input$remove_roi)
      groups(g)
      showNotification(paste0("Removed ROI '", input$remove_roi, "'."), type = "message")
    })

    # Keep the "group from ROIs" multi-select in sync with available ROI names.
    observe({
      updateSelectizeInput(session, "group_member_rois",
                           choices = roi_names(),
                           selected = isolate(input$group_member_rois))
    })

    # Keep the map focus selector in sync — offer All, any ROI, or any group.
    observe({
      updateSelectizeInput(session, "roi_show_filter",
                           choices = c("All ROIs" = "__all__", region_choices()),
                           selected = isolate(input$roi_show_filter))
    })

    # Create a named group from the selected ROIs.
    observeEvent(input$create_group_btn, {
      members <- input$group_member_rois
      if (is.null(members) || length(members) == 0) {
        showNotification("Select one or more ROIs to group.", type = "warning"); return()
      }
      gname <- input$group_name
      if (is.null(gname) || trimws(gname) == "") {
        showNotification("Name the group.", type = "warning"); return()
      }
      gname <- trimws(gname)
      g <- groups(); g[[gname]] <- list(members = members); groups(g)
      updateTextInput(session, "group_name", value = "")
      updateSelectizeInput(session, "group_member_rois", selected = character(0))
      showNotification(paste0("Created group '", gname, "' with ", length(members), " ROI(s)."), type = "message")
    })

    # Groups also live in a collapsible dropdown to keep the panel compact.
    output$group_chips <- renderUI({
      g <- groups(); nm <- names(g)
      if (length(nm) == 0) return(tags$div(style = "font-size:12px; color:#888;", "No groups yet."))
      tags$details(style = "font-size:12px;",
        tags$summary(style = "cursor:pointer; font-weight:600; color:#2c3e50;",
                     sprintf("Groups (%d)", length(nm))),
        tags$div(style = "margin-top:6px; display:flex; flex-direction:column; gap:4px; max-height:120px; overflow-y:auto;",
          lapply(seq_along(nm), function(i) {
            members <- g[[nm[i]]]$members
            tags$div(style = "display:flex; align-items:center; gap:6px;",
                     tags$span(style = "flex:1;", sprintf("%s  {%s}", nm[i], paste(members, collapse = ", "))),
                     tags$button("\u2715", title = "Remove group",
                                 style = "border:none; background:none; cursor:pointer; color:#c0392b; font-size:13px;",
                                 onclick = sprintf("Shiny.setInputValue('remove_group','%s',{priority:'event'});", nm[i])))
          })
        )
      )
    })

    observeEvent(input$remove_group, {
      g <- groups(); g[[input$remove_group]] <- NULL; groups(g)
      showNotification(paste0("Removed group '", input$remove_group, "'."), type = "message")
    })

    output$group1_info <- renderText({
      g1 <- group1_spots()
      if (length(g1) == 0) "Group 1: Empty" else paste("Group 1:", length(g1), "spots")
    })

    output$group2_info <- renderText({
      g2 <- group2_spots()
      if (length(g2) == 0) "Group 2: Empty" else paste("Group 2:", length(g2), "spots")
    })

    # Redraw the overlay whenever the checkbox, the group memberships, or the
    # transparency option change. Reading each reactive here registers the
    # dependency, then delegates to draw_group_overlay() so the toggle path and
    # the persist-on-redraw path (update_map_colors) share one implementation
    # (Reviewer 1, item 2).
    observe({
      input$show_groups
      input$show_roi_contours
      rois()                   # any change to any named ROI
      input$roi_show_filter    # focus on one ROI vs all
      groups()                 # any change to group membership
      input$transparent_groups
      draw_group_overlay()
    }, priority = 1000)


    ###1
    # DEG Analysis
    # Populate the Side A / Side B pickers from all ROIs and Groups (Reviewer 1, item 1).
    observe({
      rc <- region_choices()
      updateSelectizeInput(session, "deg_side_a", choices = rc,
                           selected = isolate(input$deg_side_a))
      updateSelectizeInput(session, "deg_side_b",
                           choices = c(rc, "Rest of tissue" = "__rest__"),
                           selected = isolate(input$deg_side_b))
    })
    deg_sideA_label <- reactiveVal("Side A")
    deg_sideB_label <- reactiveVal("Side B")

    observeEvent(input$run_deg, {
      # Force sequential IMMEDIATELY before FindMarkers is called
      future::plan("sequential")

      # ── Resolve Side A / Side B from the unified ROI/group pickers ──
      vs_rest <- identical(input$deg_side_b, "__rest__")
      g1 <- region_spots(input$deg_side_a)
      g2 <- if (vs_rest) character(0) else region_spots(input$deg_side_b)
      deg_sideA_label(region_label(input$deg_side_a))
      deg_sideB_label(if (vs_rest) "Rest" else region_label(input$deg_side_b))

      # Cluster-level DEG lives in the Clustering tab; the DEG tab compares regions.
      is_cluster_analysis <- FALSE

      if (!is_cluster_analysis) {
        if (length(g1) == 0) {
          showNotification("Side A is empty — pick a non-empty ROI/group.", type = "error"); return()
        }
        if (!vs_rest && length(g2) == 0) {
          showNotification("Side B is empty — pick a non-empty ROI/group, or 'Rest of tissue'.", type = "error"); return()
        }

        showNotification("Running DEG analysis...", type = "message", duration = NULL, id = "deg_running")

        tryCatch({
          temp_idents <- rep("Other", ncol(seurat_obj))
          names(temp_idents) <- colnames(seurat_obj)
          temp_idents[g1] <- "Group1"
          if (!vs_rest) temp_idents[g2] <- "Group2"

          cells_to_use <- if (!vs_rest) c(g1, g2) else colnames(seurat_obj)

          temp_seurat <- subset(seurat_obj, cells = cells_to_use)
          temp_seurat$temp_ident <- temp_idents[cells_to_use]
          Idents(temp_seurat) <- "temp_ident"

          deg_test   <- if (is.null(input$deg_test))   "wilcox" else input$deg_test
          deg_lfc    <- if (is.null(input$deg_logfc))  0.25     else input$deg_logfc
          deg_minpct <- if (is.null(input$deg_minpct)) 0.05     else input$deg_minpct

          markers <- Seurat::FindMarkers(
            temp_seurat, ident.1 = "Group1",
            ident.2 = if (!vs_rest) "Group2" else "Other",
            test.use = deg_test, verbose = FALSE,
            min.pct = deg_minpct, logfc.threshold = deg_lfc)




          markers$gene <- rownames(markers)
          # Use Benjamini-Hochberg (FDR) correction rather than Seurat's default
          # Bonferroni p_val_adj, so the "adjusted p (BH)" control is truthful
          # and consistent with the rest of the app (Reviewer 2, item 7).
          markers$p_val_adj <- p.adjust(markers$p_val, method = "BH")
          # Filter by BH-adjusted p and the user's |log2FC| threshold
          markers <- markers[!is.na(markers$p_val_adj) &
                            markers$p_val_adj < input$volcano_fdr &
                            abs(markers$avg_log2FC) > deg_lfc, ]
          # Rank by log2FC
          markers <- markers[order(-abs(markers$avg_log2FC)), ]          

          # ── Moran's I spatial autocorrelation ──────────────────────────────
          if (requireNamespace("spdep", quietly = TRUE)) {
            tryCatch({
              coords_full <- Seurat::GetTissueCoordinates(seurat_obj)  
              message("DEBUG group colnames: ", paste(colnames(coords_full), collapse=", "))
              message("DEBUG coord_cols found: ", paste(intersect(c("imagerow", "imagecol", "pxl_row_in_fullres", "pxl_col_in_fullres"), colnames(coords_full)), collapse=", "))           
              moran_spots <- g1   # Side A (the ROI/group being characterized)
              coords <- coords_full[rownames(coords_full) %in% moran_spots, ]

              row_col <- intersect(c("imagerow", "pxl_row_in_fullres", "y"), colnames(coords_full))[1]
              col_col <- intersect(c("imagecol", "pxl_col_in_fullres", "x"), colnames(coords_full))[1]
              coords <- as.matrix(coords[, c(row_col, col_col)])

              rownames(coords) <- rownames(coords_full)[rownames(coords_full) %in% moran_spots]

              if (nrow(coords) >= 30) {
                effective_k <- min(6, nrow(coords) - 1)
                knn_obj     <- spdep::knearneigh(coords, k = effective_k)
                listw_obj   <- spdep::nb2listw(spdep::knn2nb(knn_obj), style = "W", zero.policy = TRUE)

                spatial_assay <- if ("Spatial" %in% names(seurat_obj@assays)) "Spatial" else DefaultAssay(seurat_obj)

                # Only run Moran's I on significant DEGs (the pathologist's workflow)
                candidate_genes <- intersect(markers$gene,
                  rownames(Seurat::GetAssayData(seurat_obj, assay = spatial_assay, layer = "data")))
                                    
                if (length(candidate_genes) > 0) {
                  expr_matrix <- as.matrix(
                    Seurat::GetAssayData(seurat_obj, assay = spatial_assay, layer = "data")[candidate_genes, rownames(coords), drop = FALSE]
                  )

                  moran_results <- lapply(candidate_genes, function(gene) {
                    x <- expr_matrix[gene, ]
                    if (var(x) == 0) return(data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_))
                    tryCatch({
                      mt <- spdep::moran.test(x, listw = listw_obj,
                                              zero.policy = TRUE, alternative = "greater")
                      data.frame(gene = gene, Moran_I = as.numeric(mt$estimate[1]), Moran_pval = mt$p.value)
                    }, error = function(e) data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_))
                  })

                  moran_df            <- do.call(rbind, moran_results)
                  moran_df$Moran_padj <- p.adjust(moran_df$Moran_pval, method = "BH")

                  # Merge Moran's I back into ALL markers (non-sig genes get NA)
                  markers <- merge(markers, moran_df, by = "gene", all.x = TRUE)
                  markers$gene <- as.character(markers$gene)

                  markers$spatial_class <- dplyr::case_when(
                    is.na(markers$Moran_I)                               ~ "Not tested",
                    markers$Moran_padj < 0.05 & markers$Moran_I > 0.3   ~ "Spatially structured",
                    markers$Moran_padj < 0.05                            ~ "Weakly structured",
                    TRUE                                                  ~ "Not structured"
                  )

                  # Step 1: filter to significant DEGs only
                  # Step 2: rank by Moran's I (high to low), then by DEG p_val_adj

                  markers <- markers[order(-abs(markers$avg_log2FC)), ]
                }
              }
            }, error = function(e) {
              message("Moran's I skipped: ", e$message)
            })
          }
          # ───────────────────────────────────────────────────────────────────


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
          all_markers$gene <- as.character(all_markers$gene)
          all_markers <- all_markers[order(all_markers$cluster, all_markers$p_val_adj, -abs(all_markers$avg_log2FC)), ]

          # Reorder columns for better readability
          all_markers <- all_markers[, c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")]

          # Add Moran's I for cluster-based DEG
          if (requireNamespace("spdep", quietly = TRUE)) {
            tryCatch({
              coords_full <- Seurat::GetTissueCoordinates(seurat_obj)
              message("DEBUG group colnames: ", paste(colnames(coords_full), collapse=", "))
              message("DEBUG coord_cols found: ", paste(intersect(c("imagerow", "imagecol", "pxl_row_in_fullres", "pxl_col_in_fullres"), colnames(coords_full)), collapse=", "))
              all_spots <- names(clusters)
              row_col <- intersect(c("imagerow", "pxl_row_in_fullres", "y"), colnames(coords_full))[1]
              col_col <- intersect(c("imagecol", "pxl_col_in_fullres", "x"), colnames(coords_full))[1]
              coords <- as.matrix(coords_full[all_spots, c(row_col, col_col)])              
              
              if (nrow(coords) >= 30) {
                effective_k <- min(6, nrow(coords) - 1)
                knn_obj   <- spdep::knearneigh(coords, k = effective_k)
                listw_obj <- spdep::nb2listw(spdep::knn2nb(knn_obj), style = "W", zero.policy = TRUE)
                
                spatial_assay <- if ("Spatial" %in% names(seurat_obj@assays)) "Spatial" else DefaultAssay(seurat_obj)
                
                # Filter cluster markers same way
                all_markers <- all_markers[!is.na(all_markers$p_val_adj) &
                                            all_markers$p_val_adj < input$volcano_fdr &
                                            abs(all_markers$avg_log2FC) > 0.25, ]
                candidate_genes <- intersect(all_markers$gene,
                  rownames(Seurat::GetAssayData(seurat_obj, assay = spatial_assay, layer = "data")))                
                
                if (length(candidate_genes) > 0) {
                  expr_matrix <- as.matrix(
                    Seurat::GetAssayData(seurat_obj, assay = spatial_assay, layer = "data")[candidate_genes, all_spots, drop = FALSE]
                  )
                  moran_results <- lapply(candidate_genes, function(gene) {
                    x <- expr_matrix[gene, ]
                    if (var(x) == 0) return(data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_))
                    tryCatch({
                      mt <- spdep::moran.test(x, listw = listw_obj, zero.policy = TRUE, alternative = "greater")
                      data.frame(gene = gene, Moran_I = as.numeric(mt$estimate[1]), Moran_pval = mt$p.value)
                    }, error = function(e) data.frame(gene = gene, Moran_I = NA_real_, Moran_pval = NA_real_))
                  })
                  moran_df <- do.call(rbind, moran_results)
                  moran_df$Moran_padj <- p.adjust(moran_df$Moran_pval, method = "BH")
                  all_markers <- merge(all_markers, moran_df, by = "gene", all.x = TRUE)
                  all_markers$spatial_class <- dplyr::case_when(
                    is.na(all_markers$Moran_I)                                        ~ "Not tested",
                    all_markers$Moran_padj < 0.05 & all_markers$Moran_I > 0.3        ~ "Spatially structured",
                    all_markers$Moran_padj < 0.05                                     ~ "Weakly structured",
                    TRUE                                                               ~ "Not structured"
                  )
                  all_markers <- all_markers[order(-abs(all_markers$avg_log2FC)), ]
                }
              }
            }, error = function(e) message("Moran's I skipped for clusters: ", e$message))
          }

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

    output$deg_direction <- renderText({
      req(deg_results())
      paste0(deg_sideA_label(), " vs ", deg_sideB_label(),
             "  ·  positive log2FC = higher in ", deg_sideA_label())
    })

    output$deg_table <- renderTable({
      deg <- deg_results()
      if (!is.null(deg)) {
        top_deg <- head(deg, 50)

        # Cluster-based analysis (no Moran's I)
        if ("cluster" %in% colnames(top_deg)) {
          df <- data.frame(
            Cluster = top_deg$cluster,
            Gene    = top_deg$gene,
            log2FC  = round(top_deg$avg_log2FC, 3),
            pct.1   = round(top_deg$pct.1, 3),
            pct.2   = round(top_deg$pct.2, 3),
            p_adj   = format(top_deg$p_val_adj, scientific = TRUE, digits = 3)
          )
        } else {
          # Group-based analysis — include Moran's I if available
          df <- data.frame(
            Gene   = top_deg$gene,
            log2FC = round(top_deg$avg_log2FC, 3),
            pct.1  = round(top_deg$pct.1, 3),
            pct.2  = round(top_deg$pct.2, 3),
            p_adj  = format(top_deg$p_val_adj, scientific = TRUE, digits = 3)
          )

          if ("Moran_I" %in% colnames(top_deg)) {
            df$Moran_I    <- round(top_deg$Moran_I, 3)
            df$Moran_padj <- ifelse(is.na(top_deg$Moran_padj), "—",
                                    format(top_deg$Moran_padj, scientific = TRUE, digits = 3))
            df$Spatial    <- ifelse(is.na(top_deg$spatial_class), "Not tested",
                                    as.character(top_deg$spatial_class))
          }
        }

        df
      }
    }, rownames = FALSE)

    output$deg_volcano <- renderPlot({
      req(deg_results())
      df <- deg_results()
      fc_thresh  <- if (is.null(input$volcano_log2fc)) 0.25 else input$volcano_log2fc
      fdr_thresh <- if (is.null(input$volcano_fdr)) 0.05 else input$volcano_fdr

      df$color <- dplyr::case_when(
        df$p_val_adj < fdr_thresh & df$avg_log2FC >  fc_thresh ~ "Up",
        df$p_val_adj < fdr_thresh & df$avg_log2FC < -fc_thresh ~ "Down",
        TRUE ~ "NS"
      )

      # Replace the label line with:
      top_up   <- head(df[df$color == "Up",   ][order(-df[df$color == "Up",   ]$avg_log2FC), ], 8)
      top_down <- head(df[df$color == "Down", ][order( df[df$color == "Down", ]$avg_log2FC), ], 8)
      df$label <- ifelse(df$gene %in% c(top_up$gene, top_down$gene), df$gene, NA)

      p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = color)) +
        geom_point(aes(alpha = ifelse(color == "NS", 0.3, 0.7)), size = 2.5, show.legend = FALSE) +
        scale_alpha_identity()  +
        scale_color_manual(values = c("Up" = "#8B0000", "Down" = "#4472C4", "NS" = "grey70")) +
        geom_hline(yintercept = -log10(input$volcano_fdr), linetype = "dashed", linewidth = 0.5) +
        geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", linewidth = 0.5) +
        ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 15, box.padding = 0.5,
                                  na.rm = TRUE, color = "black") +
        theme_classic(base_size = 13) +
        labs(x = paste0("log2 Fold Change  (→ higher in ", deg_sideA_label(), ")"),
             y = "-log10(FDR)",
             title = "Differential Expression",
             subtitle = paste0(deg_sideA_label(), "  vs  ", deg_sideB_label()),
             color = "") +
        theme(legend.position = "top",
              plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5),
              legend.text = element_text(size = 10))
      deg_volcano_rv(p)
      p
    }, height = 450)


    output$deg_moran_volcano <- renderPlot({
      req(deg_results())
      df <- deg_results()
      if (!"Moran_I" %in% colnames(df)) {
        plot.new()
        text(0.5, 0.5, "Moran's I not available for this comparison.", 
            cex = 1.2, col = "grey50")
        return()
      }
      
      df <- df[!is.na(df$Moran_I) & !is.na(df$Moran_padj), ]
      
      df$spatial_class <- factor(
        df$spatial_class,
        levels = c("Not structured", "Weakly structured", "Spatially structured")
      )

      ggplot(df, aes(x = Moran_I, y = -log10(Moran_padj), color = spatial_class)) +
        geom_point(alpha = 0.7, size = 2.5) +
        scale_color_manual(values = c(
          "Spatially structured" = "#E41A1C",
          "Weakly structured"    = "#FF7F00",
          "Not structured"       = "grey60"
        )) +
        geom_vline(xintercept = 0.3, linetype = "dashed", linewidth = 0.7) +   # ← Moran threshold
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.7) + # ← padj threshold
        ggrepel::geom_text_repel(
          data = subset(df, spatial_class == "Spatially structured"),
          aes(label = gene), size = 3, max.overlaps = 15,
          color = "black", na.rm = TRUE
        ) +
        theme_classic(base_size = 13) +
        labs(x = "Moran's I", y = "-log10(adj. p-value)",
            title = "Spatially Variable DEGs", color = "") +
        theme(legend.position = "top",
              plot.title = element_text(hjust = 0.5, face = "bold"))
    }, height = 400)

    # Violin plot
    # Keep the violin Side A / Side B pickers in sync with ROIs/Groups.
    observe({
      rc <- region_choices()
      updateSelectizeInput(session, "violin_side_a", choices = rc,
                           selected = isolate(input$violin_side_a))
      updateSelectizeInput(session, "violin_side_b",
                           choices = c(rc, "Rest of tissue" = "__rest__"),
                           selected = isolate(input$violin_side_b))
    })

    observeEvent(input$plot_violin, {
      if (input$violin_feature_type == "gene") {
        feature <- input$violin_gene
      } else if (input$violin_feature_type == "metadata") {
        feature <- input$violin_metadata
      } else if (input$violin_feature_type == "pathway") {
        feature <- input$violin_pathway
      } else if (input$violin_feature_type == "cellsig") {
        feature <- input$violin_cellsig        
      } else {
        feature <- input$violin_geneset
      }      

      if (is.null(feature) || feature == "") {
        showNotification("Select a feature!", type = "warning")
        return()
      }

      v_vs_rest <- identical(input$violin_side_b, "__rest__")
      g1 <- region_spots(input$violin_side_a)
      g2 <- if (v_vs_rest) character(0) else region_spots(input$violin_side_b)
      vA <- region_label(input$violin_side_a)
      vB <- if (v_vs_rest) "Rest" else region_label(input$violin_side_b)

      if (length(g1) == 0) {
        showNotification("Side A is empty — pick a non-empty ROI/group.", type = "error"); return()
      }
      if (!v_vs_rest && length(g2) == 0) {
        showNotification("Side B is empty — pick a non-empty ROI/group, or 'Rest of tissue'.", type = "error"); return()
      }

      tryCatch({
        if (input$violin_feature_type == "gene") {
          all_values <- FetchData(seurat_obj, vars = feature, layer = "data")
          is_numeric_data <- TRUE
        } else if (input$violin_feature_type == "metadata") {
          meta_col <- seurat_obj@meta.data[[feature]]
          all_values <- data.frame(value = meta_col)
          colnames(all_values) <- feature
          rownames(all_values) <- rownames(seurat_obj@meta.data)
          is_numeric_data <- is.numeric(meta_col)
          
        } else if (input$violin_feature_type == "pathway") {
          genes <- current_hallmark_library()[[feature]]$genes
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values <- data.frame(value = rowMeans(expr_data, na.rm = TRUE))
          colnames(all_values) <- feature
          rownames(all_values) <- rownames(expr_data) 
          is_numeric_data <- TRUE
        } else if (input$violin_feature_type == "cellsig") {
          genes <- cellmarker_db[[feature]]
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values <- data.frame(value = rowMeans(expr_data, na.rm = TRUE))
          colnames(all_values) <- feature
          rownames(all_values) <- rownames(expr_data)
          is_numeric_data <- TRUE          
        } else {
          scores <- gene_set_scores()[[feature]]
          all_values <- data.frame(value = scores)
          colnames(all_values) <- feature
          rownames(all_values) <- names(scores)
          is_numeric_data <- TRUE
        }

        if (!v_vs_rest) {
          cells_to_plot <- c(g1, g2)
          group_label <- c(rep(vA, length(g1)), rep(vB, length(g2)))
          names(group_label) <- cells_to_plot
        } else {
          cells_to_plot <- spots_sf$spot_id
          group_label <- rep("Rest", length(cells_to_plot))
          names(group_label) <- cells_to_plot
          group_label[g1] <- vA
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
        attr(plot_data, "vA") <- vA
        attr(plot_data, "vB") <- vB

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
        # Title from the actual Side A / Side B labels.
        vA <- attr(plot_data, "vA"); vB <- attr(plot_data, "vB")
        comparison_text <- paste0(vA, " vs ", vB)

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

        p_violin <- ggplot(plot_data, aes(x = group, y = value, fill = group)) +
          geom_violin(trim = FALSE, alpha = 0.5) +
          geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
          labs(title = comparison_text,
              subtitle = subtitle_text,
              x = "", y = "Value") +
          theme_classic(base_size = 12) +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5))

        lv <- unique(as.character(plot_data$group))
        non_rest <- setdiff(lv, "Rest")
        fill_cols <- setNames(rep("#999999", length(lv)), lv)
        fill_cols[non_rest] <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")[seq_along(non_rest)]
        p_violin <- p_violin + scale_fill_manual(values = fill_cols)

        last_group_plot(p_violin)   # ← THIS is what was missing
        violin_rv(p_violin)         # capture for figure download (R3 #5)
        print(p_violin)
      }
    })

    # Feature comparison
    # Keep the feature-comparison region picker in sync with ROIs/Groups.
    observe({
      updateSelectizeInput(session, "compare_spots_selection",
                           choices = c("All spots" = "__all__", region_choices()),
                           selected = isolate(input$compare_spots_selection))
    })

    observeEvent(input$plot_compare, {
      # Get feature 1
      if (input$compare_type1 == "gene") {
        feature1 <- input$compare_gene1
      } else if (input$compare_type1 == "metadata") {
        feature1 <- input$compare_meta1
      } else if (input$compare_type1 == "pathway") {   # ← ADD THIS BRANCH
        feature1 <- input$compare_pathway1
      } else if (input$compare_type1 == "cellsig") {
        feature1 <- input$compare_cellsig1
      } else {
        feature1 <- input$compare_geneset1
      }


      # Get feature 2
      if (input$compare_type2 == "gene") {
        feature2 <- input$compare_gene2
      } else if (input$compare_type2 == "metadata") {
        feature2 <- input$compare_meta2
      } else if (input$compare_type2 == "pathway") {   # ← ADD THIS BRANCH
        feature2 <- input$compare_pathway2
      } else if (input$compare_type2 == "cellsig") {
        feature2 <- input$compare_cellsig2
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

      # Use spots from All spots / an ROI / a group.
      cmp_key <- input$compare_spots_selection
      if (is.null(cmp_key) || identical(cmp_key, "__all__")) {
        sel_spots <- spots_sf$spot_id
        sel_label <- "all spots"
      } else {
        sel_spots <- region_spots(cmp_key)
        sel_label <- region_label(cmp_key)
        if (length(sel_spots) == 0) {
          showNotification(paste0("'", sel_label, "' has no spots — pick a non-empty ROI/group."),
                           type = "warning", duration = 5)
          return()
        }
      }

      showNotification(paste("Using", length(sel_spots), "spots from", sel_label),
                       type = "message", duration = 3)

      tryCatch({
        # Get ALL values for feature 1
        if (input$compare_type1 == "gene") {
          all_values1 <- FetchData(seurat_obj, vars = feature1, layer = "data")[, 1]
          if (is.null(names(all_values1))) {
            names(all_values1) <- colnames(seurat_obj)
          }
          is_numeric1 <- TRUE
        } else if (input$compare_type1 == "metadata") {
          all_values1 <- seurat_obj@meta.data[[feature1]]
          names(all_values1) <- rownames(seurat_obj@meta.data)
          is_numeric1 <- is.numeric(all_values1)

        } else if (input$compare_type1 == "pathway") {
          genes <- current_hallmark_library()[[feature1]]$genes
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values1 <- setNames(rowMeans(expr_data, na.rm = TRUE), rownames(expr_data))
          is_numeric1 <- TRUE
        } else if (input$compare_type1 == "cellsig") {
          genes <- cellmarker_db[[feature1]]
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values1 <- setNames(rowMeans(expr_data, na.rm = TRUE), rownames(expr_data))
          is_numeric1 <- TRUE

        } else {
          all_values1 <- gene_set_scores()[[feature1]]
          if (is.null(names(all_values1))) names(all_values1) <- colnames(seurat_obj)
          is_numeric1 <- TRUE
        }

        # Get ALL values for feature 2
        if (input$compare_type2 == "gene") {
          all_values2 <- FetchData(seurat_obj, vars = feature2, layer = "data")[, 1]
          if (is.null(names(all_values2))) {
            names(all_values2) <- colnames(seurat_obj)
          }
          is_numeric2 <- TRUE
        } else if (input$compare_type2 == "metadata") {
          all_values2 <- seurat_obj@meta.data[[feature2]]
          names(all_values2) <- rownames(seurat_obj@meta.data)
          is_numeric2 <- is.numeric(all_values2)
        } else if (input$compare_type2 == "pathway") {
          genes <- current_hallmark_library()[[feature2]]$genes
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values2 <- setNames(rowMeans(expr_data, na.rm = TRUE), rownames(expr_data))
          is_numeric2 <- TRUE
        } else if (input$compare_type2 == "cellsig") {
          genes <- cellmarker_db[[feature2]]
          genes_in_data <- intersect(genes, rownames(seurat_obj))
          expr_data <- FetchData(seurat_obj, vars = genes_in_data, layer = "data")
          all_values2 <- setNames(rowMeans(expr_data, na.rm = TRUE), rownames(expr_data))
          is_numeric2 <- TRUE          
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
                       method = "spearman")                       
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
        attr(plot_data, "comparison_mode") <- sel_label

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
          last_group_plot(p1)
          print(p1)
          

        } else {
          # Handle non-numeric comparisons
          p <- ggplot(plot_data, aes(x = feature, y = value, fill = feature)) +
            geom_boxplot() +
            labs(title = "Feature Comparison", x = "", y = "Value") +
            theme_classic() +
            theme(legend.position = "none")
            last_group_plot(p) 
          print(p)
        }
      }
    })

    output$dl_compare_plot <- downloadHandler(
      filename = function() {
        pd <- compare_data()
        req(pd)
        f1 <- gsub("[^A-Za-z0-9_]", "", attr(pd, "feature1"))
        f2 <- gsub("[^A-Za-z0-9_]", "", attr(pd, "feature2"))
        grp <- if (!is.null(attr(pd, "comparison_mode"))) 
                attr(pd, "comparison_mode") else "comparison"
        paste0(grp, "_", f1, "_vs_", f2, ".pdf")
      },
      content = function(file) {
        pd <- compare_data()
        req(pd)
        # Rebuild the plot (same logic as renderPlot)
        feature1 <- attr(pd, "feature1")
        feature2 <- attr(pd, "feature2")
        pvalue_comp <- attr(pd, "pvalue_comp")
        stat_test   <- attr(pd, "stat_test")
        n_spots     <- attr(pd, "n_spots")

        p <- ggplot(pd, aes(x = feature, y = value, fill = feature)) +
          geom_violin(trim = FALSE, alpha = 0.5) +
          geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
          labs(title = paste(feature1, "vs", feature2),
              subtitle = if (!is.na(pvalue_comp)) {
                test_name <- if (stat_test == "ttest") "Paired t-test" else "Wilcoxon test"
                paste0(test_name, " p = ", format(pvalue_comp, scientific = TRUE, digits = 3),
                        " (n = ", n_spots, " spots)")
              } else paste0("n = ", n_spots, " spots"),
              x = "", y = "Value") +
          theme_classic(base_size = 12) +
          theme(legend.position = "none",
                plot.title = element_text(hjust = 0.5, face = "bold"),
                plot.subtitle = element_text(hjust = 0.5)) +
          scale_fill_manual(values = c("#17BECF", "#1F77B4"))

        ggsave(file, plot = p, device = "pdf", width = 6, height = 5)
      }
    )

    output$dl_compare_plot_group <- downloadHandler(

      filename = function() {
        pd <- violin_data()
        grp <- gsub(" ", "_", paste0(region_label(input$violin_side_a), "_vs_", region_label(input$violin_side_b)))

        feature <- switch(input$violin_feature_type,
          "gene"     = input$violin_gene,
          "metadata" = input$violin_metadata,
          "pathway"  = input$violin_pathway,
          "cellsig"  = input$violin_cellsig,
          "geneset"  = input$violin_geneset,
          "feature"
        )
        feature <- gsub("[^A-Za-z0-9_]", "", feature)
        paste0("group_comparison_", grp, "_", feature, ".pdf")
      },
      content = function(file) {
        p <- last_group_plot()
        req(p)
        ggsave(file, plot = p, device = "pdf", width = 6, height = 5)
      }
    )
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
      filename = function() paste0(current_sample_name(), "_DEGs_FDR",
                                  input$volcano_fdr, "_",
                                  format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
      content = function(file) {
        if (!is.null(deg_results())) write.csv(deg_results(), file, row.names = FALSE)
      }
    )
  }

  shinyApp(ui, server)
}
