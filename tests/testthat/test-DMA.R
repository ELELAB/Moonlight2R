# Test DMA function

# Run example of DMA function
data(dataMAF)
data(DEGsmatrix)
data(dataPRA)
data(EncodePromoters)
data(NCG)
data(dataFilt)
data(cscape_somatic_output)
dataDMA <- DMA(dataMAF = dataMAF,
               dataDEGs = DEGsmatrix,
               dataPRA = dataPRA,
               runCscape = FALSE,
               results_folder = ".")

# Test that output of DMA is as expected
test_that("output of DMA is a list with 2 elements", {
  expect_type(dataDMA, "list")
  expect_named(dataDMA, c("TSG", "OCG"))
  expect_true(all(sapply(dataDMA, is.character)))
})










