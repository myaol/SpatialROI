#' Run Gene Viewer - Standalone Visualization Tool
#'
#' @description
#' Launch an interactive visualization tool for exploring gene expression and metadata
#' on spatial transcriptomics data. No ROI selection - pure visualization.
#'
#' @param seurat_obj A Seurat spatial object
#' @param sample_name Character label for the dataset (default: "Sample")
#' @param show_image Logical; whether to show the H&E image if available (default: TRUE)
#'
#' @import shiny
#' @import leaflet
#' @import Seurat
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @importFrom Seurat GetTissueCoordinates GetAssayData
#'
#' @return A Shiny application for interactive gene/metadata visualization
#' @export

run_gene_viewer <- function(seurat_obj, sample_name = "Sample", show_image = TRUE) {

  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (length(seurat_obj@images) == 0) {
    stop("No spatial image found in Seurat object")
  }

  # Extract spatial data (SAME AS draw_ROI.R)
  image_name <- names(seurat_obj@images)[1]
  coords_raw <- tryCatch({
    GetTissueCoordinates(seurat_obj, image = image_name)
  }, error = function(e) seurat_obj@images[[image_name]]@coordinates)

  coords <- data.frame(spot_id = rownames(coords_raw), stringsAsFactors = FALSE)

  # Handle different coordinate column names
  if ("imagerow" %in% colnames(coords_raw) && "imagecol" %in% colnames(coords_raw)) {
    coords$x <- coords_raw$imagecol
    coords$y <- coords_raw$imagerow
  } else if ("row" %in% colnames(coords_raw) && "col" %in% colnames(coords_raw)) {
    coords$x <- coords_raw$col
    coords$y <- coords_raw$row
  } else {
    coords$x <- coords_raw[, 1]
    coords$y <- coords_raw[, 2]
  }

  # Flip y-axis (CRITICAL for alignment)
  coords$y <- max(coords$y) - coords$y + min(coords$y)

  # Get all genes and metadata
  all_genes <- rownames(seurat_obj)
  all_metadata <- colnames(seurat_obj@meta.data)

  # Prepare H&E image (SAME AS draw_ROI.R)
  he_image_base64 <- NULL
  he_image_bounds <- NULL

  if (show_image) {
    tryCatch({
      img <- seurat_obj@images[[image_name]]@image
      coords_full <- GetTissueCoordinates(seurat_obj, image = image_name)

      # Use pixel coordinates if available, otherwise use spot coordinates
      if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
        pixel_x <- coords_full$pxl_col_in_fullres
        pixel_y <- coords_full$pxl_row_in_fullres
      } else {
        # Use the same coordinates as spots
        pixel_x <- coords$x
        pixel_y <- coords$y
      }

      px_min <- max(1, floor(min(pixel_x, na.rm = TRUE)))
      px_max <- min(dim(img)[2], ceiling(max(pixel_x, na.rm = TRUE)))
      py_min <- max(1, floor(min(pixel_y, na.rm = TRUE)))
      py_max <- min(dim(img)[1], ceiling(max(pixel_y, na.rm = TRUE)))

      he_crop <- img[py_min:py_max, px_min:px_max, , drop = FALSE]

      tmpfile <- tempfile(fileext = ".png")
      png::writePNG(he_crop, target = tmpfile)
      he_image_base64 <- paste0("data:image/png;base64,", base64enc::base64encode(tmpfile))
      unlink(tmpfile)

      he_image_bounds <- list(south = py_max, west = px_min, north = py_min, east = px_max)
    }, error = function(e) {
      message("Could not load H&E image: ", e$message)
    })
  }

  # UI
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
        body {
          margin: 0;
          padding: 0;
          overflow: hidden;
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        }

        /* Top Header */
        .top-header {
          position: fixed;
          top: 0;
          left: 0;
          right: 0;
          height: 70px;
          background: linear-gradient(135deg, #0072B5 0%, #E18727 100%);
          color: white;
          display: flex;
          align-items: center;
          padding: 0 30px;
          z-index: 2000;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }
        .top-header h1 {
          margin: 0;
          font-size: 26px;
          font-weight: bold;
        }
        .top-header .subtitle {
          margin-left: 15px;
          font-size: 14px;
          opacity: 0.9;
        }

        /* Main Layout */
        .main-container {
          display: flex;
          height: 100vh;
          margin-top: 70px;
        }

        /* Left Control Panel */
        .control-panel {
          width: 350px;
          background: white;
          box-shadow: 2px 0 10px rgba(0,0,0,0.1);
          overflow-y: auto;
          padding: 20px;
        }

        /* Map Container */
        .map-container {
          flex: 1;
          position: relative;
          background: #f5f5f5;
        }

        /* Control Sections */
        .control-section {
          background: #f8f9fa;
          padding: 15px;
          border-radius: 8px;
          margin-bottom: 15px;
          border-left: 4px solid #0072B5;
        }
        .control-section h3 {
          margin: 0 0 15px 0;
          color: #2c3e50;
          font-size: 16px;
          font-weight: bold;
        }

        /* Map Controls */
        .map-controls {
          position: absolute;
          top: 20px;
          right: 20px;
          z-index: 1000;
          background: white;
          padding: 15px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          min-width: 220px;
        }
        .map-controls details {
          margin-bottom: 10px;
        }
        .map-controls summary {
          font-weight: bold;
          font-size: 14px;
          cursor: pointer;
          margin-bottom: 8px;
          color: #2c3e50;
        }
        .map-controls details[open] summary {
          color: #0072B5;
        }

        /* Spot Info */
        .spot-info {
          position: absolute;
          top: 20px;
          left: 20px;
          z-index: 1000;
          background: white;
          padding: 12px 18px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          font-weight: bold;
          font-size: 14px;
        }

        /* Color Legend */
        .legend-container {
          background: white;
          padding: 15px;
          border-radius: 8px;
          margin-top: 15px;
        }
        .legend-container h4 {
          margin: 0 0 10px 0;
          color: #2c3e50;
          font-size: 14px;
        }

        /* Buttons */
        .btn-custom {
          width: 100%;
          padding: 10px;
          font-size: 14px;
          font-weight: bold;
          border-radius: 6px;
          border: none;
          cursor: pointer;
          transition: all 0.3s;
        }
        .btn-visualize {
          background: #0072B5;
          color: white;
        }
        .btn-visualize:hover {
          background: #005a8f;
          transform: translateY(-2px);
        }
        .btn-clear {
          background: #e74c3c;
          color: white;
          margin-top: 10px;
        }
        .btn-clear:hover {
          background: #c0392b;
        }

        /* Info Box */
        .info-box {
          background: #e8f4f8;
          border-left: 4px solid #0072B5;
          padding: 12px;
          border-radius: 4px;
          margin-bottom: 15px;
          font-size: 13px;
          line-height: 1.5;
        }

        /* Select Input Styling */
        .form-control {
          border-radius: 6px;
          border: 1px solid #ddd;
        }
        .form-control:focus {
          border-color: #0072B5;
          box-shadow: 0 0 0 0.2rem rgba(0, 114, 181, 0.25);
        }

        /* Responsive */
        @media (max-width: 992px) {
          .control-panel {
            width: 300px;
          }
        }
        @media (max-width: 768px) {
          .main-container {
            flex-direction: column;
          }
          .control-panel {
            width: 100%;
            max-height: 40vh;
          }
          .map-container {
            height: 60vh;
          }
        }
      "))
    ),

    # Top Header
    tags$div(class = "top-header",
             tags$h1("🔬 SpatialScope - Gene Viewer"),
             tags$span(class = "subtitle", paste("●", sample_name))
    ),

    # Main Container
    tags$div(class = "main-container",
             # Left Control Panel
             tags$div(class = "control-panel",
                      # Instructions
                      tags$div(class = "info-box",
                               HTML("<b>📍 How to Use:</b><br>
                                    1. Select feature type (Gene or Metadata)<br>
                                    2. Choose your feature of interest<br>
                                    3. Click 'Visualize' to display on map<br>
                                    4. Adjust colors and spot size as needed")
                      ),

                      # Feature Selection
                      tags$div(class = "control-section",
                               tags$h3("🎯 Feature Selection"),
                               selectInput("feature_type", "Feature Type:",
                                           choices = c("None" = "none",
                                                       "Gene Expression" = "gene",
                                                       "Metadata" = "metadata"),
                                           selected = "none"),

                               conditionalPanel(
                                 condition = "input.feature_type == 'gene'",
                                 selectInput("gene_select", "Select Gene:",
                                             choices = c("", all_genes),
                                             selected = "")
                               ),

                               conditionalPanel(
                                 condition = "input.feature_type == 'metadata'",
                                 selectInput("meta_select", "Select Metadata:",
                                             choices = c("", all_metadata),
                                             selected = "")
                               ),

                               conditionalPanel(
                                 condition = "input.feature_type != 'none'",
                                 actionButton("visualize_btn", "🎨 Visualize on Map",
                                              class = "btn-custom btn-visualize"),
                                 actionButton("clear_viz_btn", "🗑️ Clear Visualization",
                                              class = "btn-custom btn-clear")
                               )
                      ),

                      # Color Scheme
                      conditionalPanel(
                        condition = "input.feature_type != 'none'",
                        tags$div(class = "control-section",
                                 tags$h3("🎨 Color Scheme"),
                                 selectInput("color_palette", "Color Palette:",
                                             choices = c("Grey to Red" = "greyred",
                                                         "Blue-White-Red" = "bwr",
                                                         "Viridis" = "viridis",
                                                         "Rainbow" = "rainbow"),
                                             selected = "greyred"),

                                 tags$div(class = "legend-container",
                                          tags$h4("Color Legend"),
                                          plotOutput("color_legend", height = "120px")
                                 )
                        )
                      ),

                      # Statistics
                      conditionalPanel(
                        condition = "output.feature_visualized",
                        tags$div(class = "control-section",
                                 tags$h3("📊 Statistics"),
                                 verbatimTextOutput("feature_stats")
                        )
                      )
             ),

             # Map Container
             tags$div(class = "map-container",
                      leafletOutput("map", height = "100%"),

                      # Spot Info
                      tags$div(class = "spot-info",
                               textOutput("spot_info", inline = TRUE)
                      ),

                      # Map Controls
                      tags$div(class = "map-controls",
                               tags$details(
                                 list(open = "open"),
                                 tags$summary("▼ Spot Size"),
                                 sliderInput("spot_size", NULL,
                                             min = 1, max = 15, value = 4,
                                             step = 0.5, width = "180px")
                               ),
                               tags$details(
                                 list(open = "open"),
                                 tags$summary("▼ H&E Opacity"),
                                 sliderInput("image_opacity", NULL,
                                             min = 0, max = 1, value = 0.6,
                                             step = 0.05, width = "180px")
                               )
                      )
             )
    )
  )

  # Server
  server <- function(input, output, session) {

    rv <- reactiveValues(
      current_feature = NULL,
      feature_values = NULL,
      feature_visualized = FALSE
    )

    # Initialize map (SAME CRS AS draw_ROI.R)
    output$map <- renderLeaflet({
      x_range <- range(coords$x)
      y_range <- range(coords$y)
      x_buffer <- diff(x_range) * 0.1
      y_buffer <- diff(y_range) * 0.1

      m <- leaflet(options = leafletOptions(
        crs = leafletCRS(crsClass = "L.CRS.Simple"),
        zoomControl = TRUE
      )) %>%
        fitBounds(
          lng1 = x_range[1] - x_buffer, lat1 = y_range[1] - y_buffer,
          lng2 = x_range[2] + x_buffer, lat2 = y_range[2] + y_buffer
        )

      # Add H&E image if available
      if (!is.null(he_image_base64) && !is.null(he_image_bounds)) {
        js_code <- paste0("
          function(el, x) {
            var map = this;
            var bounds = [[", he_image_bounds$south, ", ", he_image_bounds$west, "],
                          [", he_image_bounds$north, ", ", he_image_bounds$east, "]];
            window.heImageOverlay = L.imageOverlay('", he_image_base64, "', bounds, {
              opacity: 0.6,
              interactive: false
            });
            window.heImageOverlay.addTo(map);
            window.heImageOverlay.bringToBack();

            Shiny.addCustomMessageHandler('updateHEOpacity', function(message) {
              if (window.heImageOverlay) {
                window.heImageOverlay.setOpacity(message.opacity);
              }
            });
          }
        ")
        m <- m %>% htmlwidgets::onRender(js_code)
      }

      # Add default spots
      m %>%
        addCircleMarkers(
          lng = coords$x,
          lat = coords$y,
          radius = 4,
          fillColor = "#CCCCCC",
          fillOpacity = 0.7,
          stroke = TRUE,
          color = "#333333",
          weight = 1,
          layerId = coords$spot_id,
          group = "spots"
        )
    })

    # Update H&E opacity
    observe({
      if (!is.null(he_image_base64)) {
        session$sendCustomMessage("updateHEOpacity", list(opacity = input$image_opacity))
      }
    }) %>% bindEvent(input$image_opacity)

    # Visualize feature
    observeEvent(input$visualize_btn, {
      req(input$feature_type != "none")

      if (input$feature_type == "gene") {
        req(input$gene_select != "")
        feature_name <- input$gene_select
        feature_values <- GetAssayData(seurat_obj, slot = "data")[feature_name, ]
      } else if (input$feature_type == "metadata") {
        req(input$meta_select != "")
        feature_name <- input$meta_select
        feature_values <- seurat_obj@meta.data[[feature_name]]
      }

      rv$current_feature <- feature_name
      rv$feature_values <- feature_values
      rv$feature_visualized <- TRUE

      # Get color palette
      if (is.numeric(feature_values)) {
        cols <- switch(input$color_palette,
                       "greyred" = c("grey", "red"),
                       "bwr" = c("blue", "white", "red"),
                       "viridis" = viridis::viridis(100),
                       "rainbow" = rainbow(100))

        # Normalize values
        val_range <- range(feature_values, na.rm = TRUE)
        normalized <- (feature_values - val_range[1]) / (val_range[2] - val_range[1])
        normalized[is.na(normalized)] <- 0

        col_func <- colorRampPalette(cols)
        spot_colors <- col_func(100)[as.integer(normalized * 99) + 1]
      } else {
        # Categorical
        unique_vals <- unique(feature_values)
        cols <- rainbow(length(unique_vals))
        spot_colors <- cols[match(feature_values, unique_vals)]
      }

      # Update map
      leafletProxy("map") %>%
        clearGroup("spots") %>%
        addCircleMarkers(
          lng = coords$x,
          lat = coords$y,
          radius = input$spot_size,
          fillColor = spot_colors,
          fillOpacity = 0.8,
          stroke = TRUE,
          color = "#333333",
          weight = 1,
          layerId = coords$spot_id,
          group = "spots"
        )

      showNotification(paste("✅ Visualizing:", feature_name), type = "message", duration = 2)
    })

    # Clear visualization
    observeEvent(input$clear_viz_btn, {
      rv$current_feature <- NULL
      rv$feature_values <- NULL
      rv$feature_visualized <- FALSE

      leafletProxy("map") %>%
        clearGroup("spots") %>%
        addCircleMarkers(
          lng = coords$x,
          lat = coords$y,
          radius = input$spot_size,
          fillColor = "#CCCCCC",
          fillOpacity = 0.7,
          stroke = TRUE,
          color = "#333333",
          weight = 1,
          layerId = coords$spot_id,
          group = "spots"
        )

      showNotification("🗑️ Visualization cleared", type = "warning", duration = 2)
    })

    # Update spot size
    observe({
      if (!is.null(rv$feature_values)) {
        # Re-trigger visualization with new spot size
        isolate({
          if (is.numeric(rv$feature_values)) {
            cols <- switch(input$color_palette,
                           "greyred" = c("grey", "red"),
                           "bwr" = c("blue", "white", "red"),
                           "viridis" = viridis::viridis(100),
                           "rainbow" = rainbow(100))

            val_range <- range(rv$feature_values, na.rm = TRUE)
            normalized <- (rv$feature_values - val_range[1]) / (val_range[2] - val_range[1])
            normalized[is.na(normalized)] <- 0

            col_func <- colorRampPalette(cols)
            spot_colors <- col_func(100)[as.integer(normalized * 99) + 1]
          } else {
            unique_vals <- unique(rv$feature_values)
            cols <- rainbow(length(unique_vals))
            spot_colors <- cols[match(rv$feature_values, unique_vals)]
          }

          leafletProxy("map") %>%
            clearGroup("spots") %>%
            addCircleMarkers(
              lng = coords$x,
              lat = coords$y,
              radius = input$spot_size,
              fillColor = spot_colors,
              fillOpacity = 0.8,
              stroke = TRUE,
              color = "#333333",
              weight = 1,
              layerId = coords$spot_id,
              group = "spots"
            )
        })
      }
    }) %>% bindEvent(input$spot_size)

    # Spot info
    output$spot_info <- renderText({
      paste("Total Spots:", nrow(coords))
    })

    # Color legend
    output$color_legend <- renderPlot({
      req(rv$feature_values)

      if (is.numeric(rv$feature_values)) {
        cols <- switch(input$color_palette,
                       "greyred" = c("grey", "red"),
                       "bwr" = c("blue", "white", "red"),
                       "viridis" = viridis::viridis(100),
                       "rainbow" = rainbow(100))

        val_range <- range(rv$feature_values, na.rm = TRUE)
        x <- seq(val_range[1], val_range[2], length.out = 100)

        par(mar = c(3, 1, 2, 1))
        image(matrix(x, ncol = 1), col = colorRampPalette(cols)(100),
              axes = FALSE, main = rv$current_feature)
        axis(1, at = c(0, 0.5, 1),
             labels = round(c(val_range[1], mean(val_range), val_range[2]), 2),
             cex.axis = 0.8)
      }
    })

    # Feature statistics
    output$feature_stats <- renderText({
      req(rv$feature_values)

      if (is.numeric(rv$feature_values)) {
        paste0(
          "Feature: ", rv$current_feature, "\n",
          "Min: ", round(min(rv$feature_values, na.rm = TRUE), 3), "\n",
          "Max: ", round(max(rv$feature_values, na.rm = TRUE), 3), "\n",
          "Mean: ", round(mean(rv$feature_values, na.rm = TRUE), 3), "\n",
          "Median: ", round(median(rv$feature_values, na.rm = TRUE), 3)
        )
      } else {
        paste0(
          "Feature: ", rv$current_feature, "\n",
          "Type: Categorical\n",
          "Unique values: ", length(unique(rv$feature_values))
        )
      }
    })

    output$feature_visualized <- reactive({
      rv$feature_visualized
    })
    outputOptions(output, "feature_visualized", suspendWhenHidden = FALSE)
  }

  # Run app
  shinyApp(ui, server)
}
