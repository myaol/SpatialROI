#' Launch SpatialROI
#'
#' @description
#' Launch the Shiny interface for SpatialROI. Users can either:
#' - Call `SpatialROI::run_SpatialROI()` to open the app and upload their own data interactively, or
#' - Call `SpatialROI::run_SpatialROI("path/to/data.rds")` or `SpatialROI::run_SpatialROI("demo")` to launch the app preloaded with a dataset.
#'
#' @param data_path Optional path to a user data file (.rds) or "demo" to load example data.
#' @param ... Additional arguments passed to `shiny::runApp()`.
#'
#' @examples
#' \dontrun{
#' SpatialROI::run_SpatialROI()
#' SpatialROI::run_SpatialROI("demo")
#' SpatialROI::run_SpatialROI("~/Downloads/mydata.rds")
#' }
#'
#' @export
run_SpatialROI <- function(data_path = NULL, ...) {
  # Apply before the initial upload screen is created. A reverse proxy may still
  # impose a lower limit on a hosted deployment.
  options(shiny.maxRequestSize = 500 * 1024^2)
  options(SpatialROI.data_path = data_path)
  if (identical(data_path, "demo")) {
    demo_file <- system.file("extdata", "example_visium.rds", package = "SpatialROI")
    if (!file.exists(demo_file) || demo_file == "")
      stop("Demo data not found. Please ensure example_visium.rds exists in inst/extdata.")
    options(SpatialROI.data_path = demo_file)
  }
  app_dir <- system.file("app", package = "SpatialROI")
  if (app_dir == "") stop("Could not find app directory. Try reinstalling SpatialROI.", call. = FALSE)
  shiny::runApp(app_dir, display.mode = "normal", launch.browser = TRUE, ...)
}
