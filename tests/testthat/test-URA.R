# Test URA function

# Run example of URA
data(DEGsmatrix)
data(DiseaseList)
data(EAGenes)
dataDEGs <- DEGsmatrix
data(dataGRN)

# Load example data of URA serving as reference point
data(dataURA)

dataURA_test <- URA(dataGRN = dataGRN,
                    DEGsmatrix = dataDEGs,
                    BPname = c("apoptosis",
                               "proliferation of cells"))

# Test that output of URA is as expected
test_that("output of URA is a double with 2 columns", {
  expect_type(dataURA_test, "double")
  expect_equal(ncol(dataURA_test), 2)
  expect_equal(colnames(dataURA_test), c("apoptosis",
                                    "proliferation of cells"))
  expect_true(all(sapply(dataURA_test, is.numeric)))
})

# Test that output of URA is as expected compared to reference / example data
test_that("URA output is identical to reference point", {
  expect_equal(dataURA_test, dataURA)
})

dataURA_test <- URA(dataGRN = dataGRN,
                    DEGsmatrix = dataDEGs,
                    BPname = c("apoptosis",
                               "proliferation of cells"),
                    nCores = 2)

# Test that output of URA is as expected
test_that("output of URA is a double with 2 columns", {
  expect_type(dataURA_test, "double")
  expect_equal(ncol(dataURA_test), 2)
  expect_equal(colnames(dataURA_test), c("apoptosis",
                                    "proliferation of cells"))
  expect_true(all(sapply(dataURA_test, is.numeric)))
})

# Test that output of URA is as expected compared to reference / example data
test_that("URA output is identical to reference point", {
  expect_equal(dataURA_test, dataURA)
})
