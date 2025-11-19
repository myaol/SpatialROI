#' Launch SpatialScope
#'
#' @description
#' Launch the Shiny interface for SpatialScope. Users can either:
#' - Call `SpatialScope::run_SpatialScope()` to open the app and upload their own data interactively, or
#' - Call `SpatialScope::run_SpatialScope("path/to/data.rds")` or `SpatialScope::run_SpatialScope("demo")` to launch the app preloaded with a dataset.
#'
#' @param data_path Optional path to a user data file (.rds, .csv, etc.) or "demo" to load example data.
#' @param ... Additional arguments passed to `shiny::runApp()`.
#'
#' @examples
#' \dontrun{
#' SpatialScope::run_SpatialScope()
#' SpatialScope::run_SpatialScope("demo")
#' SpatialScope::run_SpatialScope("~/Downloads/mydata.rds")
#' }
#'
#' @export
run_SpatialScope <- function(data_path = NULL, ...) {
  # Save user-specified data path so the Shiny app can read it
  options(SpatialScope.data_path = data_path)

  # Handle "demo" special case
  if (identical(data_path, "demo")) {
    demo_file <- system.file("extdata", "example_visium.rds", package = "SpatialScope")
    if (!file.exists(demo_file) || demo_file == "") {
      stop("Demo data not found. Please ensure example_visium.rds exists in inst/extdata.")
    }
    options(SpatialScope.data_path = demo_file)
  }

  # Locate the app directory inside the package
  appDir <- system.file("app", package = "SpatialScope")
  if (appDir == "") {
    stop("Could not find app directory. Try reinstalling SpatialScope.", call. = FALSE)
  }

  # Launch the app
  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE, ...)
}
