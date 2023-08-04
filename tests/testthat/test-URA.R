# Test URA function

# Run example of URA
data(DEGsmatrix)
dataDEGs <- DEGsmatrix
dataGRN <- GRN(TFs = rownames(dataDEGs)[1:10],
               DEGsmatrix = dataDEGs,
               DiffGenes = TRUE,
               normCounts = dataFilt,
               nGenesPerm = 5,
               nBoot = 5)
dataURA <- URA(dataGRN = dataGRN,
               DEGsmatrix = dataDEGs,
               BPname = c("apoptosis",
                          "proliferation of cells"))

# Test that output of URA is as expected
test_that("output of URA is a double with 2 columns", {
  expect_type(dataURA, "double")
  expect_equal(ncol(dataURA), 2)
  expect_equal(colnames(dataURA), c("apoptosis",
                                    "proliferation of cells"))
  expect_true(all(sapply(dataURA, is.numeric)))
})










