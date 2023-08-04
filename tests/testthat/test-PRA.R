# Test PRA function

# Run example of PRA
data(dataURA)
dataDual <- PRA(dataURA = dataURA,
                BPname = c("apoptosis", "proliferation of cells"),
                thres.role = 0)

# Test that output of PRA is as expected
test_that("output of PRA is a list with 2 elements", {
  expect_type(dataDual, "list")
  expect_named(dataDual, c("TSG", "OCG"))
  expect_named(dataDual$TSG)
  expect_named(dataDual$OCG)
  expect_true(all(sapply(dataDual, is.numeric)))
})










