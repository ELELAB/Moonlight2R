# Test GLS function

# Run example of GLS
data(dataDMA)
genes_query <- Reduce(c, dataDMA)
dataGLS_test <- GLS(genes = genes_query,
                    query_string = "AND cancer AND driver AND '1980/01/01'[Date - Publication] : '2023/01/01'[Date - Publication]")

# Load example data of GLS serving as reference point
data(dataGLS)

# Test that output of GLS is as expected

# Test number of columns
test_that("Number of columns in GLS output are equal to 8", {
  expect_equal(ncol(dataGLS_test), 8)
})

dataGLS_colnames <- c("pmid", "gene", "doi", "title", "abstract", "year",
                      "keywords", "pubmed_count")

# Test correct column names
test_that("Column names in GLS output are correct", {
  expect_equal(colnames(dataGLS_test), dataGLS_colnames)
})

# Test expected class of values in output
test_that("pmid, gene, doi, title, abstract, year, and keywords are characters", {
  expect_type(dataGLS_test$pmid, "character")
  expect_type(dataGLS_test$gene, "character")
  expect_type(dataGLS_test$doi, "character")
  expect_type(dataGLS_test$title, "character")
  expect_type(dataGLS_test$abstract, "character")
  expect_type(dataGLS_test$year, "character")
})

test_that("pubmed_count is double", {
  expect_type(dataGLS_test$pubmed_count, "double")
})

# Test that output of GLS is as expected compared to reference / example data
test_that("GLS output is identical to reference point", {
  expect_equal(dataGLS_test, dataGLS)
})
