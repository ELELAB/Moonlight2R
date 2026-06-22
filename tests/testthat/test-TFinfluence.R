# Test TFinfluence function

# Run example of TFinfluence function
data(dataPRA)
data(DEGsmatrix)
data(dataTRRUST)
data(dataMAF)
data(dataMAVISp)


reference_TFresults <- data.frame(
    Target = c("PDE2A", "GCAT", "ARMC9", "FAM155A"),
    Moonlight_gene_z_score = c(1.37642904406401, 0.810977174901446, 0.367466020609254, 1.19863281130568),
    Moonlight_Oncogenic_Mediator = c("TSG", "TSG", "TSG", "OCG"),
    logFC_target = c(-3.36278468754307, 1.0407914607103, 1.03363210942098, 1.97338960700756),
    TF = as.character(c(NA, NA, NA, NA)),
    InteractionType = as.character(c(NA, NA, NA, NA)),
    PMID = as.character(c(NA, NA, NA, NA)),
    tf_mutation = as.character(c(NA, NA, NA, NA)),
    stab_class = as.character(c(NA, NA, NA, NA)),
    in_MAVISp = c(FALSE,FALSE,FALSE,FALSE),
    mutation_available = c(FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
)

TFresults <- TFinfluence(dataTRRUST = dataTRRUST,
            dataMAF = dataMAF,
            dataDEGs = DEGsmatrix,
            dataPRA = dataPRA,
            dataMAVISp = dataMAVISp,
            stabClassMAVISp = 'rasp')

character_columns <- c('Target','Moonlight_Oncogenic_Mediator','TF','InteractionType','PMID','tf_mutation','stab_class')
numeric_columns <- c('Moonlight_gene_z_score','logFC_target')
logical_columns <- c('in_MAVISp','mutation_available')
all_columns <- c('Target', 'Moonlight_gene_z_score', 'Moonlight_Oncogenic_Mediator', 'logFC_target', 'TF', 'InteractionType', 'PMID', 'tf_mutation', 'stab_class','in_MAVISp','mutation_available')


# Test that output is as expected
test_that('Output of TFinfluence is a tibble', {
    expect_named(TFresults, all_columns)
    expect_true(all(sapply(TFresults[character_columns], is.character)))
    expect_true(all(sapply(TFresults[numeric_columns], is.numeric)))
    expect_true(all(sapply(TFresults[logical_columns], is.logical)))
})

# Test that the content of the output dataframe is as expected
test_that("Output of TFinfluence does not correspond to reference", {
    expect_equal(TFresults, reference_TFresults,ignore_attr = TRUE)
})



