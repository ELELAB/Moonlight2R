# Test PRA function

# Run example of PRA
data(dataURA)
dataPRA_test <- PRA(dataURA = dataURA,
                    BPname = c("apoptosis", "proliferation of cells"),
                    thres.role = 0)

# Load example data of PRA serving as reference point
data(dataPRA)

# Test that output of PRA is as expected
test_that("output of PRA is a list with 2 elements", {
  expect_type(dataPRA_test, "list")
  expect_named(dataPRA_test, c("TSG", "OCG"))
  expect_named(dataPRA_test$TSG)
  expect_named(dataPRA_test$OCG)
  expect_true(all(sapply(dataPRA_test, is.numeric)))
})

# Test that output of PRA is as expected compared to reference / example data
test_that("PRA output is identical to reference point", {
  expect_equal(dataPRA_test, dataPRA)
})
