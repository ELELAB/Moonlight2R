# Test FEA function

# Run example of FEA
data("DEGsmatrix")
data("DiseaseList")
data("EAGenes")
DEGsmatrix_ora <- DEGsmatrix[1:10, ]
dataFEA_test_ora <- FEA(DiffMatrix = DEGsmatrix_ora, method = "ora")
dataFEA_test_fgsea <- FEA(DiffMatrix = DEGsmatrix, method= "fgsea", seed = 123)
# Load example data of FEA serving as reference points
data(dataFEA)
data(dataFEA_fgsea)

# Test that output of FEA using ORA-based method is as expected 

# Test number of columns
test_that("Number of columns in FEA output are equal to 7", {
  expect_equal(dim(dataFEA_test_ora), c(101, 7))
})

dataFEA_colnames <- c("Diseases.or.Functions.Annotation", "Moonlight.Z.score",
                      "p.Value", "FDR", "commonNg", "FunctionNg", "Molecules")

# Test correct column names
test_that("Column names in FEA output are correct", {
  expect_equal(colnames(dataFEA_test_ora), dataFEA_colnames)
})

# Test expected class of values in output
test_that("Moonlight scores, p-values and FDR values are numeric", {
  expect_type(dataFEA_test_ora$Diseases.or.Functions.Annotation, "character")
  expect_type(dataFEA_test_ora$Moonlight.Z.score, "double")
  expect_type(dataFEA_test_ora$p.Value, "double")
  expect_type(dataFEA_test_ora$FDR, "double")
  expect_type(dataFEA_test_ora$commonNg, "integer")
  expect_type(dataFEA_test_ora$FunctionNg, "integer")
  expect_type(dataFEA_test_ora$Molecules, "character")
})

# Test that output of FEA is as expected compared to reference 
test_that("FEA output is identical to reference point", {
  expect_equal(dataFEA_test_ora, dataFEA)
})


# Test that output of FEA using fgsea method is as expected

# Test number of columns
test_that("Number of columns in FEA output are equal to 8", {
  expect_equal(dim(dataFEA_test_fgsea), c(101, 9))
})

dataFEA_colnames <- c("Diseases.or.Functions.Annotation","Moonlight.Z.score", "p.value",
                      "padj", "ES", "NES", "commonNg", "FunctionNg", "Molecules")

# Test correct column names
test_that("Column names in FEA output are correct", {
  expect_equal(colnames(dataFEA_test_fgsea), dataFEA_colnames)
})

# Test expected class of values in output
test_that("Moonlight scores, p-values and FDR values are numeric", {
  expect_type(dataFEA_test_fgsea$Diseases.or.Functions.Annotation, "character")
  expect_type(dataFEA_test_fgsea$Moonlight.Z.score, "double")
  expect_type(dataFEA_test_fgsea$p.value, "double")
  expect_type(dataFEA_test_fgsea$padj, "double")
  expect_type(dataFEA_test_fgsea$ES, "double")
  expect_type(dataFEA_test_fgsea$NES, "double")
  expect_type(dataFEA_test_fgsea$commonNg, "integer")
  expect_type(dataFEA_test_fgsea$FunctionNg, "integer")
  expect_type(dataFEA_test_fgsea$Molecules, "character")
  expect_true(all(sapply(dataFEA_test_fgsea$Molecules, is.character)))

})

# Test that output of FEA is as expected compared to reference 
test_that("FEA output is similar enough to reference point", {
  expect_equal(dataFEA_fgsea$Moonlight.Z.score, dataFEA_test_fgsea$Moonlight.Z.score)
  expect_equal(dataFEA_fgsea$NES, dataFEA_test_fgsea$NES, tolerance = 1e-6)
  expect_equal(dataFEA_fgsea$ES, dataFEA_test_fgsea$ES, tolerance = 1e-6)
})

# Check using tolerance of half order of magnitude
test_that("FEA p-values are qualitatively similar", {
expect_equal(
  -log10(pmax(dataFEA_test_fgsea$p.value, .Machine$double.xmin)),
  -log10(pmax(dataFEA_fgsea$p.value, .Machine$double.xmin)),
  tolerance = 0.5
)})