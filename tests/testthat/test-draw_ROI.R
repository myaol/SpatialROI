test_that("draw_ROI function is available", {
  expect_true(is.function(SpatialScope::draw_ROI))
})

test_that("demo data file exists", {
  demo_file <- system.file("extdata", "example_visium.rds", package = "SpatialScope")
  expect_true(nchar(demo_file) > 0)
})
