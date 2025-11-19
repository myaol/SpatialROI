#' Interactive ROI Selection on Spatial Image
#'
#' Launch a minimal interactive map to manually draw a region of interest (ROI)
#' on a Seurat spatial object. Returns the spot IDs contained in the ROI.
#'
#' @param seurat_obj A Seurat spatial object with spatial coordinates and H&E image.
#' @param sample_name Optional character label for display.
#' @param show_image Logical; whether to show the histology image (default TRUE).
#' @return A character vector of spot IDs contained within the drawn ROI.
#' @import shiny
#' @import leaflet
#' @import sf
#' @importFrom Seurat GetTissueCoordinates
#' @import shinyjs
#' @export
draw_ROI <- function(seurat_obj, sample_name = "Sample1", show_image = TRUE) {

  # Check required packages
  required_pkgs <- c("sf", "leaflet", "shiny", "sp")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Required packages missing: ", paste(missing_pkgs, collapse = ", "))
  }

  # ---- Validate input ----
  if (!inherits(seurat_obj, "Seurat")) stop("Input must be a Seurat object.")
  if (length(seurat_obj@images) == 0) stop("No spatial image found in Seurat object.")

  # ---- Extract spatial data ----
  image_name <- names(seurat_obj@images)[1]
  coords <- tryCatch({
    Seurat::GetTissueCoordinates(seurat_obj, image = image_name)
  }, error = function(e) seurat_obj@images[[image_name]]@coordinates)

  spots_df <- data.frame(spot_id = rownames(coords), stringsAsFactors = FALSE)
  if ("imagerow" %in% colnames(coords) && "imagecol" %in% colnames(coords)) {
    spots_df$x <- coords$imagecol
    spots_df$y <- coords$imagerow
  } else if ("row" %in% colnames(coords) && "col" %in% colnames(coords)) {
    spots_df$x <- coords$col
    spots_df$y <- coords$row
  } else {
    spots_df$x <- coords[, 1]
    spots_df$y <- coords[, 2]
  }

  # flip y
  spots_df$y <- max(spots_df$y) - spots_df$y + min(spots_df$y)

  # Store x and y before converting to sf
  x_vals <- spots_df$x
  y_vals <- spots_df$y

  # Convert to sf
  spots_sf <- sf::st_as_sf(spots_df, coords = c("x", "y"), crs = NA)

  # Add x and y back as regular columns
  spots_sf$x <- x_vals
  spots_sf$y <- y_vals

  # ---- Prepare H&E image ----
  he_image_base64 <- NULL
  he_image_bounds <- NULL

  if (show_image) {
    tryCatch({
      img <- seurat_obj@images[[image_name]]@image
      coords_full <- Seurat::GetTissueCoordinates(seurat_obj, image = image_name)

      if ("pxl_col_in_fullres" %in% colnames(coords_full)) {
        pixel_x <- coords_full$pxl_col_in_fullres
        pixel_y <- coords_full$pxl_row_in_fullres
      } else {
        pixel_x <- spots_df$x
        pixel_y <- spots_df$y
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
      message("Could not extract H&E image: ", e$message)
    })
  }

  # ---- UI ----
  ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(
      shiny::tags$style(shiny::HTML("
        body { margin: 0; padding: 10px; }
        h3 { color: #0072B5; margin-bottom: 10px; }
        .instruction {
          background: #e8f4f8;
          border-left: 4px solid #0072B5;
          padding: 10px;
          margin: 10px 0;
          border-radius: 4px;
        }
        .map-wrapper {
          position: relative;
          height: 70vh;
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
        }
        .spot-count-display {
          position: absolute;
          top: 20px;
          left: 50%;
          transform: translateX(-50%);
          z-index: 1000;
          background: white;
          padding: 12px 18px;
          border-radius: 8px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.2);
          font-weight: bold;
          font-size: 15px;
        }
      "))
    ),
    shiny::tags$h3(paste("Draw ROI -", sample_name)),
    shiny::tags$div(class = "instruction",
                    shiny::tags$p(style = "margin: 0; color: #333;",
                                  shiny::HTML("<strong>Instructions:</strong> Look for the ✏️ pencil button in the <strong>top-left corner</strong> of the map. Click it and draw your region of interest.")
                    )
    ),
    shiny::tags$div(class = "map-wrapper",
                    leaflet::leafletOutput("map", height = "100%"),

                    # Map controls (top right)
                    shiny::tags$div(class = "map-controls",
                                    shiny::tags$details(
                                      list(open = "open"),
                                      shiny::tags$summary(
                                        shiny::tags$span("▼ Spot Size")
                                      ),
                                      shiny::sliderInput("spot_size", NULL,
                                                         min = 1, max = 10, value = 3,
                                                         step = 0.5, width = "180px")
                                    ),
                                    shiny::tags$details(
                                      list(open = "open"),
                                      shiny::tags$summary(
                                        shiny::tags$span("▼ H&E Opacity")
                                      ),
                                      shiny::sliderInput("image_opacity", NULL,
                                                         min = 0, max = 1, value = 0.6,
                                                         step = 0.05, width = "180px")
                                    )
                    ),

                    # Spot count (top center)
                    shiny::tags$div(class = "spot-count-display",
                                    shiny::textOutput("spot_count_display", inline = TRUE))
    ),
    shiny::br(),
    shiny::fluidRow(
      shiny::column(6, shiny::actionButton("clear", "🗑️ Clear Selection",
                                           class = "btn btn-danger btn-block")),
      shiny::column(6, shiny::actionButton("done", "✅ Done & Return Spots",
                                           class = "btn btn-success btn-block"))
    ),
    shiny::br(),
    shiny::verbatimTextOutput("spot_count")
  )

  # ---- Server ----
  server <- function(input, output, session) {
    rv <- shiny::reactiveValues(selected = character(0))

    output$map <- leaflet::renderLeaflet({
      x_range <- range(spots_sf$x)
      y_range <- range(spots_sf$y)
      x_buffer <- diff(x_range) * 0.1
      y_buffer <- diff(y_range) * 0.1

      m <- leaflet::leaflet(options = leaflet::leafletOptions(
        crs = leaflet::leafletCRS(crsClass = "L.CRS.Simple"),
        zoomControl = TRUE
      )) |>
        leaflet::fitBounds(
          lng1 = x_range[1] - x_buffer, lat1 = y_range[1] - y_buffer,
          lng2 = x_range[2] + x_buffer, lat2 = y_range[2] + y_buffer
        ) |>
        leaflet::addCircleMarkers(
          lng = spots_sf$x,
          lat = spots_sf$y,
          radius = 3, fillColor = "lightblue", fillOpacity = 0.7,
          stroke = TRUE, color = "black", weight = 0.5,
          group = "spots"
        )

      # COMBINE both H&E image AND freehand drawing in ONE onRender call
      js_code <- paste0("
    function(el, x) {
      var map = this;

      // Store map reference globally
      window.leafletMap = map;

      // Add H&E image if available
      ", if (!is.null(he_image_base64)) {
        paste0("
        var bounds = [[", he_image_bounds$south, ", ", he_image_bounds$west, "],
                      [", he_image_bounds$north, ", ", he_image_bounds$east, "]];
        window.heImageOverlay = L.imageOverlay('", he_image_base64, "', bounds, {
          opacity: 0.6,
          interactive: false
        });
        window.heImageOverlay.addTo(map);
        window.heImageOverlay.bringToBack();
        ")
      } else {
        ""
      }, "

      // Initialize drawn items layer
      if (!map.drawnItems) {
        map.drawnItems = new L.FeatureGroup();
        map.addLayer(map.drawnItems);
      }

      // Create custom freehand draw control
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
        Shiny.setInputValue('draw_new', feature, {priority: 'event'});

        isDrawing = false;
        freehandPoints = [];
        map.dragging.enable();
        map.getContainer().style.cursor = '';
      });

      Shiny.addCustomMessageHandler('clearDrawings', function(message) {
        console.log('Clear drawings called');
        if (map.drawnItems) {
          map.drawnItems.clearLayers();
          console.log('Cleared drawn items');
        } else {
          console.log('No drawnItems found');
        }
      });

      // Custom message handler for updating H&E opacity with unique name
      if (!window.heOpacityHandlerRegistered) {
        Shiny.addCustomMessageHandler('updateHEOpacity_roi', function(message) {
          console.log('Opacity update received:', message.opacity);
          if (window.heImageOverlay) {
            window.heImageOverlay.setOpacity(message.opacity);
            console.log('Opacity set to:', message.opacity);
          } else {
            console.log('heImageOverlay not found');
          }
        });
        window.heOpacityHandlerRegistered = true;
      }

      // Store map globally for access from R
      window.roiMap = map;
    }
  ")

      m |> htmlwidgets::onRender(js_code)
    })

    # Update spot size dynamically
    shiny::observe({
      leaflet::leafletProxy("map") |>
        leaflet::clearGroup("spots") |>
        leaflet::addCircleMarkers(
          lng = spots_sf$x,
          lat = spots_sf$y,
          radius = input$spot_size,
          fillColor = "lightblue",
          fillOpacity = 0.7,
          stroke = TRUE,
          color = "black",
          weight = 0.5,
          group = "spots"
        )
    }) |> shiny::bindEvent(input$spot_size)

    # Update H&E opacity dynamically
    shiny::observe({
      if (!is.null(he_image_base64)) {
        session$sendCustomMessage("updateHEOpacity_roi", list(opacity = input$image_opacity))
      }
    }) |> shiny::bindEvent(input$image_opacity, ignoreInit = FALSE)

    # Handle new drawings
    shiny::observeEvent(input$draw_new, {
      coords <- input$draw_new$geometry$coordinates[[1]]
      poly_x <- sapply(coords, function(p) p[1])
      poly_y <- sapply(coords, function(p) p[2])

      inside <- sp::point.in.polygon(spots_sf$x, spots_sf$y, poly_x, poly_y) > 0
      new_spots <- spots_df$spot_id[inside]
      rv$selected <- unique(c(rv$selected, new_spots))

      shiny::showNotification(paste("Added", length(new_spots), "spots"), type = "message", duration = 2)
    })

    # Display selected spots in yellow
    shiny::observe({
      sel_spots <- rv$selected

      leaflet::leafletProxy("map") |>
        leaflet::clearGroup("selected")

      if (length(sel_spots) > 0) {
        sel_indices <- which(spots_sf$spot_id %in% sel_spots)
        leaflet::leafletProxy("map") |>
          leaflet::addCircleMarkers(
            lng = spots_sf$x[sel_indices],
            lat = spots_sf$y[sel_indices],
            radius = input$spot_size + 1,
            fillColor = "yellow",
            fillOpacity = 1,
            stroke = TRUE,
            color = "red",
            weight = 2,
            group = "selected"
          )
      }
    })

    # Display spot count at top center
    output$spot_count_display <- shiny::renderText({
      paste("Selected:", length(rv$selected), "spots")
    })

    # Display spot count below map
    output$spot_count <- shiny::renderText({
      paste("Total selected spots:", length(rv$selected))
    })

    shiny::observeEvent(input$clear, {
      rv$selected <- character(0)
      session$sendCustomMessage("clearDrawings", list())

      # Also clear the selected spots visualization
      leaflet::leafletProxy("map") |>
        leaflet::clearGroup("selected")

      shiny::showNotification("Selection cleared", type = "message", duration = 2)
    })

    shiny::observeEvent(input$done, {
      if (length(rv$selected) == 0) {
        shiny::showNotification("No spots selected! Draw a region first.",
                                type = "warning", duration = 3)
      } else {
        shiny::stopApp(rv$selected)
      }
    })
  }

  # ---- Launch ----
  shiny::runApp(shiny::shinyApp(ui, server))
}
