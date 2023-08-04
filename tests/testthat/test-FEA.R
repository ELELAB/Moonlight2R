# Test FEA function

# Run example of FEA
data("DEGsmatrix")
DEGsmatrix <- DEGsmatrix[1:10, ]
dataFEA <- FEA(DEGsmatrix = DEGsmatrix)

# Test that output of FEA is as expected

# Test number of columns
test_that("Number of columns in FEA output are equal to 7", {
  expect_equal(dim(dataFEA), c(101, 7))
})

dataFEA_colnames <- c("Diseases.or.Functions.Annotation", "Moonlight.Z.score",
                      "p.Value", "FDR", "commonNg", "FunctionNg", "Molecules")

# Test correct column names
test_that("Column names in FEA output are correct", {
  expect_equal(colnames(dataFEA), dataFEA_colnames)
})

# Test expected class of values in output
test_that("Moonlight scores, p-values and FDR values are numeric", {
  expect_type(dataFEA$Diseases.or.Functions.Annotation, "character")
  expect_type(dataFEA$Moonlight.Z.score, "double")
  expect_type(dataFEA$p.Value, "double")
  expect_type(dataFEA$FDR, "double")
  expect_type(dataFEA$commonNg, "integer")
  expect_type(dataFEA$FunctionNg, "integer")
  expect_type(dataFEA$Molecules, "character")
})








