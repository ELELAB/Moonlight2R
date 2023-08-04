# Test DMA function

# Get data from DMA function
data(dataDMA)

# Test that output of DMA is as expected
test_that("output of DMA is a list with 2 elements", {
  expect_type(dataDMA, "list")
  expect_named(dataDMA, c("TSG", "OCG"))
  expect_true(all(sapply(dataDMA, is.character)))
})










