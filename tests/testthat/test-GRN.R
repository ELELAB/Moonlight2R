# Test GRN function

# Run example of GRN
data('DEGsmatrix')
data('dataFilt')
dataDEGs <- DEGsmatrix
dataGRN <- GRN(TFs = rownames(dataDEGs)[1:10],
               DEGsmatrix = dataDEGs,
               DiffGenes = TRUE,
               normCounts = dataFilt,
               nGenesPerm = 5,
               nBoot = 5)

# Test that output of GRN is as expected

# Test that output of GRN is a list of length 2
test_that("output of GRN is a list", {
  expect_type(dataGRN, "list")
  expect_type(dataGRN$miTFGenes, "double")
  expect_type(dataGRN$maxmi, "double")
  expect_length(dataGRN, 2)
})

# Test that names of elements in GRN are correct
test_that("names of elements in GRN output are correct", {
  expect_equal(names(dataGRN), c("miTFGenes", "maxmi"))
})








