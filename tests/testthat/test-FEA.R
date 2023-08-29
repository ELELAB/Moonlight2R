# Test FEA function

# Run example of FEA
data("DEGsmatrix")
data("DiseaseList")
data("EAGenes")
DEGsmatrix <- DEGsmatrix[1:10, ]
dataFEA_test <- FEA(DEGsmatrix = DEGsmatrix)

# Load example data of FEA serving as reference point
data(dataFEA)

# Test that output of FEA is as expected

# Test number of columns
test_that("Number of columns in FEA output are equal to 7", {
  expect_equal(dim(dataFEA_test), c(101, 7))
})

dataFEA_colnames <- c("Diseases.or.Functions.Annotation", "Moonlight.Z.score",
                      "p.Value", "FDR", "commonNg", "FunctionNg", "Molecules")

# Test correct column names
test_that("Column names in FEA output are correct", {
  expect_equal(colnames(dataFEA_test), dataFEA_colnames)
})

# Test expected class of values in output
test_that("Moonlight scores, p-values and FDR values are numeric", {
  expect_type(dataFEA_test$Diseases.or.Functions.Annotation, "character")
  expect_type(dataFEA_test$Moonlight.Z.score, "double")
  expect_type(dataFEA_test$p.Value, "double")
  expect_type(dataFEA_test$FDR, "double")
  expect_type(dataFEA_test$commonNg, "integer")
  expect_type(dataFEA_test$FunctionNg, "integer")
  expect_type(dataFEA_test$Molecules, "character")
})

# Test that output of FEA is as expected compared to reference 
test_that("FEA output is identical to reference point", {
  expect_equal(dataFEA_test, dataFEA)
})
