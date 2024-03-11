# Test GMA function

# Run example of GMA function
data("dataMethyl")
data("dataFilt")
data("dataPRA")
data("DEGsmatrix")
data("LUAD_sample_anno")
data("NCG")
data("EncodePromoters")
data("MetEvidenceDriver")
pattern <- "^(.{4}-.{2}-.{4}-.{2}).*"
colnames(dataFilt) <- sub(pattern, "\\1", colnames(dataFilt))
dataGMA_test <- GMA(dataMET = dataMethyl, dataEXP = dataFilt,
                    dataPRA = dataPRA, dataDEGs = DEGsmatrix,
                    sample_info = LUAD_sample_anno, met_platform = "HM450",
                    prevalence_filter = NULL,
                    output_dir = "./GMAresults", cores = 1, roadmap.epigenome.ids = "E096",
                    roadmap.epigenome.groups = NULL)

# Load example data of GMA serving as reference point
data(dataGMA)

# Test that output of GMA is as expected
test_that("output of GMA is a list with 2 elements", {
  expect_type(dataGMA_test, "list")
  expect_named(dataGMA_test, c("TSG", "OCG"))
  expect_true(all(sapply(dataGMA_test, is.character)))
})

# Test that output of GMA is as expected compared to reference / example data
test_that("GMA output is identical to reference point", {
  expect_equal(dataGMA_test, dataGMA)
})
