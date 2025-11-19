test_that("SpatialScope app can launch", {
  expect_true(is.function(SpatialScope::run_spatial_selector))
})
