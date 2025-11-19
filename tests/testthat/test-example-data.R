# tests/testthat/test-example-data.R
test_that("Example data loads", {
  expect_true(file.exists(system.file("app", package = "SpatialScope")))
})
