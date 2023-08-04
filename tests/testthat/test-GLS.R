# Test GLS function

# Run example of GLS
data(dataDMA)
genes_query <- Reduce(c, dataDMA)
dataGLS <- GLS(genes = genes_query)

# Test that output of GLS is as expected

# Test number of columns
test_that("Number of columns in GLS output are equal to 8", {
  expect_equal(ncol(dataGLS), 8)
})

dataGLS_colnames <- c("pmid", "gene", "doi", "title", "abstract", "year",
                      "keywords", "pubmed_count")

# Test correct column names
test_that("Column names in GLS output are correct", {
  expect_equal(colnames(dataGLS), dataGLS_colnames)
})

# Test expected class of values in output
test_that("pmid, gene, doi, title, abstract, year, and keywords are characters", {
  expect_type(dataGLS$pmid, "character")
  expect_type(dataGLS$gene, "character")
  expect_type(dataGLS$doi, "character")
  expect_type(dataGLS$title, "character")
  expect_type(dataGLS$abstract, "character")
  expect_type(dataGLS$year, "character")
})

test_that("pubmed_count is double", {
  expect_type(dataGLS$pubmed_count, "double")
})
