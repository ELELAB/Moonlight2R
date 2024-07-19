# Test TFinfluence function

# Run example of TFinfluence function
specific_protein <- loadMAVISp(mavispDB = 'data/mavisp_db', 
                               mode = 'simple',
                               proteins_of_interest = c('RUNX1'))

all_proteins <- loadMAVISp(mavispDB = 'data/mavisp_db', 
                           mode = 'simple')

ensemble <- loadMAVISp(mavispDB = 'data/mavisp_db', 
                       mode = 'ensemble',
                       ensemble = 'cabsflex')

# Test that output is as expected
test_that('Output of loadMAVISp simple mode with one protein specified is a list of one tibble', {
    expect_type(specific_protein, 'list')
    expect_length(specific_protein, 1)
    expect_true(all(sapply(specific_protein, is_tibble())))
})

test_that('Output of loadMAVISp simple mode is a list of tibbles', {
    expect_type(specific_protein, 'list')
    expect_true(all(sapply(specific_protein, is_tibble())))
})

test_that('Output of loadMAVISp ensemble mode is a list of tibbles', {
    expect_type(specific_protein, 'list')
    expect_true(all(sapply(specific_protein, is_tibble())))
})
