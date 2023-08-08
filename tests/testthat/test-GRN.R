# Test GRN function

# Run example of GRN
data('DEGsmatrix')
data('dataFilt')
dataDEGs <- DEGsmatrix
set.seed(317)
dataGRN_test <- GRN(TFs = sample(rownames(dataDEGs), 100),
                    DEGsmatrix = dataDEGs,
                    DiffGenes = TRUE,
                    normCounts = dataFilt,
                    nGenesPerm = 5,
                    nBoot = 5)

# Load example data of GRN serving as reference point
data(dataGRN)

# Test that output of GRN is as expected

# Test that output of GRN is a list of length 2
test_that("output of GRN is a list", {
  expect_type(dataGRN_test, "list")
  expect_type(dataGRN_test$miTFGenes, "double")
  expect_type(dataGRN_test$maxmi, "double")
  expect_length(dataGRN_test, 2)
})

# Test that names of elements in GRN are correct
test_that("names of elements in GRN output are correct", {
  expect_equal(names(dataGRN_test), c("miTFGenes", "maxmi"))
})

# Test that output of GRN is as expected compared to reference / example data
test_that("GRN output is identical to reference point", {
  expect_equal(dataGRN_test, dataGRN, tolerance = 0.1)
})








