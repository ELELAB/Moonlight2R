# Test FEA function

# Run example of FEA
data("DEGsmatrix")
data("DiseaseList")
data("EAGenes")
DEGsmatrix_ora <- DEGsmatrix[1:10, ]
dataFEA_test_ora <- FEA(DiffMatrix = DEGsmatrix_ora, method="ora")
set.seed(123)
dataFEA_test_fgsea <- FEA(DiffMatrix = DEGsmatrix, method="fgsea")
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
  expect_equal(dim(dataFEA_test_fgsea), c(101, 8))
})

dataFEA_colnames <- c("Diseases.or.Functions.Annotation","Moonlight.Z.score", "p.value",
                      "padj", "ES", "NES", "size", "leadingEdge")

# Test correct column names
test_that("Column names in FEA output are correct", {
  expect_equal(colnames(dataFEA_test_fgsea), dataFEA_colnames)
})

# Test expected class of values in output
test_that("Moonlight scores, p-values and FDR values are numeric", {
  expect_type(dataFEA_test_fgsea$Diseases.or.Functions.Annotation, "character")
  expect_type(dataFEA_test_fgsea$p.value, "double")
  expect_type(dataFEA_test_fgsea$padj, "double")
  expect_type(dataFEA_test_fgsea$ES, "double")
  expect_type(dataFEA_test_fgsea$NES, "double")
  expect_type(dataFEA_test_fgsea$size, "integer")
  expect_type(dataFEA_test_fgsea$leadingEdge, "list")
  expect_true(all(sapply(dataFEA_test_fgsea$leadingEdge, is.character)))

})

# Test that output of FEA is as expected compared to reference 
test_that("FEA output is identical to reference point", {
  
  test_res <- dataFEA_test_fgsea[order(dataFEA_test_fgsea$Diseases.or.Functions.Annotation), ]
  ref_res <- dataFEA_fgsea[order(dataFEA_fgsea$Diseases.or.Functions.Annotation), ]
  
  expect_equal(test_res$Diseases.or.Functions.Annotation, ref_res$Diseases.or.Functions.Annotation)
  expect_equal(test_res$size, ref_res$size)
  
  # allow strict tolerance for ES 
  expect_equal(test_res$ES, ref_res$ES, tolerance = 1e-5)
  # allow more loose tolerance for NES
  expect_equal(test_res$NES, ref_res$NES, tolerance = 1e-1)

  # allow more loose tolerance for permutation p-values 
  expect_equal(test_res$p.value, ref_res$p.value, tolerance = 1e-1)
  expect_equal(test_res$padj, ref_res$padj, tolerance = 1e-1)

})
