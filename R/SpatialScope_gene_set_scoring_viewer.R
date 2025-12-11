#' Run Gene Set Score Viewer
#'
#' @description
#' Launch an interactive visualization tool for gene set scoring on spatial transcriptomics data.
#' Includes pre-defined cell type signatures from CellMarker 2.0 database.
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
#' @importFrom Seurat GetTissueCoordinates GetAssayData AddModuleScore DefaultAssay
#'
#' @return A Shiny application for interactive gene set score visualization
#' @export

run_gene_set_viewer <- function(seurat_obj, sample_name = "Sample", show_image = TRUE) {

  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  if (length(seurat_obj@images) == 0) {
    stop("No spatial image found in Seurat object")
  }

  # Print CellMarker citation to console
  message("\n", paste(rep("=", 80), collapse = ""))
  message("Gene Set Signatures from CellMarker 2.0 Database")
  message(paste(rep("=", 80), collapse = ""))
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

  # Load gene set libraries
  temp_data <- prepare_seurat_data(seurat_obj, sample_name, show_image = FALSE)
  signature_library_human <- temp_data$signature_library_human
  signature_library_mouse <- temp_data$signature_library_mouse

  # Extract spatial data (SAME AS draw_ROI.R and gene_viewer)
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

  # Prepare H&E image
  he_image_base64 <- NULL
  he_image_bounds <- NULL

  if (show_image) {
    tryCatch({
      img <- seurat_obj@images[[image_name]]@image
      coords_full <- GetTissueCoordinates(seurat_obj, image = image_name)

      if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
        pixel_x <- coords_full$pxl_col_in_fullres
        pixel_y <- coords_full$pxl_row_in_fullres
      } else {
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

  # UI (same as before - keeping your exact CSS and layout)
  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
        body {
          margin: 0;
          padding: 0;
          overflow: hidden;
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        }
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
        .main-container {
          display: flex;
          height: 100vh;
          margin-top: 70px;
        }
        .control-panel {
          width: 400px;
          background: white;
          box-shadow: 2px 0 10px rgba(0,0,0,0.1);
          overflow-y: auto;
          padding: 20px;
        }
        .map-container {
          flex: 1;
          position: relative;
          background: #f5f5f5;
        }
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
        .citation-box {
          background: #e8f4f8;
          border-left: 4px solid #0072B5;
          padding: 12px;
          border-radius: 4px;
          margin-bottom: 15px;
          font-size: 12px;
          line-height: 1.5;
        }
        .citation-title {
          font-weight: bold;
          color: #0072B5;
          font-size: 13px;
          margin-bottom: 5px;
        }
        .citation-links {
          margin-top: 8px;
          display: flex;
          gap: 10px;
        }
        .citation-links a {
          color: #0072B5;
          text-decoration: none;
          font-size: 11px;
        }
        .citation-links a:hover {
          text-decoration: underline;
        }
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
        .btn-custom {
          width: 100%;
          padding: 12px;
          font-size: 14px;
          font-weight: bold;
          border-radius: 6px;
          border: none;
          cursor: pointer;
          transition: all 0.3s;
        }
        .btn-calculate {
          background: #0072B5;
          color: white;
        }
        .btn-calculate:hover {
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
        .form-control {
          border-radius: 6px;
          border: 1px solid #ddd;
        }
        .form-control:focus {
          border-color: #0072B5;
          box-shadow: 0 0 0 0.2rem rgba(0, 114, 181, 0.25);
        }
        @media (max-width: 992px) {
          .control-panel {
            width: 350px;
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

    tags$div(class = "top-header",
             tags$h1("🧬 SpatialScope - Gene Set Analyzer"),
             tags$span(class = "subtitle", paste("●", sample_name))
    ),

    tags$div(class = "main-container",
             tags$div(class = "control-panel",
                      tags$div(class = "citation-box",
                               tags$div(class = "citation-title", "📚 Cell Marker Database"),
                               "Pre-defined signatures are curated from CellMarker 2.0, a manually curated database of ",
                               tags$b("26,915 cell markers"), " across ", tags$b("2,578 cell types"), " and ", tags$b("656 tissues."),
                               tags$br(),
                               tags$br(),
                               tags$b("Citation:"), " Hu C, Li T, Xu Y, et al. ",
                               tags$i("Nucleic Acids Res."), " 2023;51(D1):D870-D876.",
                               tags$div(class = "citation-links",
                                        tags$a(href = "https://academic.oup.com/nar/article/51/D1/D870/6775381",
                                               target = "_blank", "📄 Read Paper"),
                                        tags$a(href = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/",
                                               target = "_blank", "🌐 Visit Database"),
                                        tags$a(href = "http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html",
                                               target = "_blank", "⬇️ Download Data")
                               )
                      ),

                      tags$div(class = "control-section",
                               tags$h3("🔬 Species Selection"),
                               selectInput("species_select", "Select Species:",
                                           choices = c("Human" = "human", "Mouse" = "mouse"),
                                           selected = "human"),
                               tags$p(style = "font-size: 12px; color: #7f8c8d; margin: 5px 0 0 0;",
                                      "💡 Gene symbols will be matched to selected species")
                      ),

                      tags$div(class = "control-section",
                               tags$h3("📋 Signature Selection"),
                               selectInput("signature_select", "Pre-defined Signatures:",
                                           choices = names(signature_library_human),
                                           selected = "Custom"),
                               tags$p(style = "font-size: 12px; color: #7f8c8d; margin: 5px 0 0 0;",
                                      "💡 Select a cell type signature or enter custom genes below")
                      ),

                      tags$div(class = "control-section",
                               tags$h3("✍️ Gene Input"),
                               textAreaInput("gene_input", NULL,
                                             placeholder = "Enter genes (one per line or comma-separated):\nCD3D\nCD3E\nCD8A",
                                             height = "150px"),
                               textInput("gene_set_name", "Gene Set Name:", value = "GeneSet1")
                      ),

                      tags$div(class = "control-section",
                               tags$h3("⚙️ Parameters"),
                               selectInput("score_method", "Scoring Method:",
                                           choices = c("Mean Expression" = "mean",
                                                       "AddModuleScore" = "addmodulescore",
                                                       "GSVA (rank-based)" = "gsva"),
                                           selected = "mean"),
                               tags$p(style = "font-size: 11px; color: #7f8c8d; margin: -8px 0 10px 0;",
                                      HTML("<b>Mean</b>: Fastest, simple average<br>
                                            <b>AddModuleScore</b>: Fast, with control features<br>
                                            <b>GSVA</b>: Rank-based enrichment (may take 10-30s)")),

                               selectInput("color_palette", "Color Palette:",
                                           choices = c("Grey to Red" = "greyred",
                                                       "Blue-White-Red" = "bwr",
                                                       "Viridis" = "viridis",
                                                       "Magma" = "magma"),
                                           selected = "greyred"),

                               actionButton("calculate_btn", "🎯 Calculate & Visualize",
                                            class = "btn-custom btn-calculate"),
                               actionButton("clear_btn", "🗑️ Clear Visualization",
                                            class = "btn-custom btn-clear")
                      ),

                      conditionalPanel(
                        condition = "output.score_calculated",
                        tags$div(class = "control-section",
                                 tags$h3("📊 Score Statistics"),
                                 verbatimTextOutput("score_stats")
                        ),

                        tags$div(class = "legend-container",
                                 tags$h4("Color Legend"),
                                 plotOutput("color_legend", height = "120px")
                        )
                      )
             ),

             tags$div(class = "map-container",
                      leafletOutput("map", height = "100%"),

                      tags$div(class = "spot-info",
                               textOutput("spot_info", inline = TRUE)
                      ),

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
      current_scores = NULL,
      genes_used = NULL,
      score_calculated = FALSE
    )

    current_library <- reactive({
      if (input$species_select == "human") {
        signature_library_human
      } else {
        signature_library_mouse
      }
    })

    observe({
      updateSelectInput(session, "signature_select",
                        choices = names(current_library()),
                        selected = "Custom")
    })

    observeEvent(input$signature_select, {
      if (input$signature_select != "Custom") {
        sig_data <- current_library()[[input$signature_select]]
        genes <- if (is.list(sig_data) && "genes" %in% names(sig_data)) {
          sig_data$genes
        } else {
          sig_data
        }
        updateTextAreaInput(session, "gene_input", value = paste(genes, collapse = "\n"))
        updateTextInput(session, "gene_set_name", value = input$signature_select)
      }
    })

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

    observe({
      if (!is.null(he_image_base64)) {
        session$sendCustomMessage("updateHEOpacity", list(opacity = input$image_opacity))
      }
    }) %>% bindEvent(input$image_opacity)

    # FIXED: Calculate gene set scores
    observeEvent(input$calculate_btn, {
      gene_text <- input$gene_input
      if (is.null(gene_text) || gene_text == "") {
        showNotification("⚠️ Please enter gene names", type = "warning", duration = 3)
        return()
      }

      genes <- unlist(strsplit(gene_text, "[,\n\t ]+"))
      genes <- genes[genes != ""]
      genes <- unique(genes)

      if (length(genes) == 0) {
        showNotification("⚠️ No valid genes found", type = "warning", duration = 3)
        return()
      }

      all_genes <- rownames(seurat_obj)
      genes_found <- genes[genes %in% all_genes]
      genes_missing <- genes[!genes %in% all_genes]

      if (length(genes_found) == 0) {
        showNotification(paste("❌ None of the genes found in dataset"),
                         type = "error", duration = 5)
        return()
      }

      if (length(genes_missing) > 0) {
        showNotification(paste("⚠️", length(genes_missing), "genes not found:",
                               paste(head(genes_missing, 3), collapse = ", ")),
                         type = "warning", duration = 5)
      }

      showNotification(paste("✅ Using", length(genes_found), "genes"),
                       type = "message", duration = 3)

      withProgress(message = 'Calculating scores...', value = 0, {
        if (input$score_method == "mean") {
          expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_found, , drop = FALSE]
          scores <- colMeans(as.matrix(expr_data))

        } else if (input$score_method == "addmodulescore") {
          temp_obj <- AddModuleScore(seurat_obj, features = list(genes_found),
                                     name = "GeneSet", assay = DefaultAssay(seurat_obj))
          scores <- temp_obj$GeneSet1

        } else if (input$score_method == "gsva") {
          # VECTORIZED GSVA ALTERNATIVE: Fast rank-based enrichment score
          expr_mat <- GetAssayData(seurat_obj, slot = "data", assay = DefaultAssay(seurat_obj))

          # VECTORIZED RANKING: Rank all genes for all cells simultaneously
          rank_mat <- tryCatch({
            apply(expr_mat, 2, rank, ties.method = "average")
          }, error = function(e) {
            stop("Calculation interrupted")
          })

          # Get ranks for gene set genes across all cells
          gene_set_ranks <- rank_mat[genes_found, , drop = FALSE]

          # Calculate mean rank of gene set for each cell (vectorized)
          mean_rank_geneset <- colMeans(gene_set_ranks)

          # Calculate overall mean rank for each cell
          mean_rank_all <- colMeans(rank_mat)

          # Calculate enrichment scores (vectorized across all cells)
          # Normalize to [-1, 1] scale like GSVA
          scores <- (mean_rank_geneset - mean_rank_all) / (nrow(expr_mat) / 2)
        }

        rv$current_scores <- scores
        rv$genes_used <- genes_found
        rv$score_calculated <- TRUE
      })

      cols <- switch(input$color_palette,
                     "greyred" = c("grey", "red"),
                     "bwr" = c("blue", "white", "red"),
                     "viridis" = viridis::viridis(100),
                     "magma" = viridis::magma(100))

      val_range <- range(scores, na.rm = TRUE)
      normalized <- (scores - val_range[1]) / (val_range[2] - val_range[1])
      normalized[is.na(normalized)] <- 0

      col_func <- colorRampPalette(cols)
      spot_colors <- col_func(100)[as.integer(normalized * 99) + 1]

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

      showNotification(paste("✅ Visualizing:", input$gene_set_name),
                       type = "message", duration = 3)
    })

    observeEvent(input$clear_btn, {
      rv$current_scores <- NULL
      rv$genes_used <- NULL
      rv$score_calculated <- FALSE

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

    observe({
      if (!is.null(rv$current_scores)) {
        isolate({
          cols <- switch(input$color_palette,
                         "greyred" = c("grey", "red"),
                         "bwr" = c("blue", "white", "red"),
                         "viridis" = viridis::viridis(100),
                         "magma" = viridis::magma(100))

          val_range <- range(rv$current_scores, na.rm = TRUE)
          normalized <- (rv$current_scores - val_range[1]) / (val_range[2] - val_range[1])
          normalized[is.na(normalized)] <- 0

          col_func <- colorRampPalette(cols)
          spot_colors <- col_func(100)[as.integer(normalized * 99) + 1]

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

    output$spot_info <- renderText({
      if (rv$score_calculated) {
        paste(input$gene_set_name, "- Total Spots:", nrow(coords))
      } else {
        paste("Total Spots:", nrow(coords))
      }
    })

    output$color_legend <- renderPlot({
      req(rv$current_scores)

      cols <- switch(input$color_palette,
                     "greyred" = c("grey", "red"),
                     "bwr" = c("blue", "white", "red"),
                     "viridis" = viridis::viridis(100),
                     "magma" = viridis::magma(100))

      val_range <- range(rv$current_scores, na.rm = TRUE)
      x <- seq(val_range[1], val_range[2], length.out = 100)

      par(mar = c(3, 1, 2, 1))
      image(matrix(x, ncol = 1), col = colorRampPalette(cols)(100),
            axes = FALSE, main = input$gene_set_name)
      axis(1, at = c(0, 0.5, 1),
           labels = round(c(val_range[1], mean(val_range), val_range[2]), 3),
           cex.axis = 0.8)
    })

    output$score_stats <- renderText({
      req(rv$current_scores)

      paste0(
        "Gene Set: ", input$gene_set_name, "\n",
        "Genes Used: ", length(rv$genes_used), "\n",
        "Method: ", input$score_method, "\n\n",
        "Score Statistics:\n",
        "Min: ", round(min(rv$current_scores, na.rm = TRUE), 4), "\n",
        "Max: ", round(max(rv$current_scores, na.rm = TRUE), 4), "\n",
        "Mean: ", round(mean(rv$current_scores, na.rm = TRUE), 4), "\n",
        "Median: ", round(median(rv$current_scores, na.rm = TRUE), 4)
      )
    })

    output$score_calculated <- reactive({
      rv$score_calculated
    })
    outputOptions(output, "score_calculated", suspendWhenHidden = FALSE)
  }

  shinyApp(ui, server)
}
