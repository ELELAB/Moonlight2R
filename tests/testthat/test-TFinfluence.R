# Test TFinfluence function

# Run example of TFinfluence function
data(dataPRA)
data(DEGsmatrix)
data(dataTRRUST)
data(dataMAF)
data(dataMAVISp)


reference_TFresults <- data.frame(
    Target = c("CHST1", "GEN1", "TAP1", "TAP1"),
    Moonlight_gene_z_score = c(1.06066017177982, 1.30096115353815, 0.8262079533974, 0.8262079533974),
    Moonlight_Oncogenic_Mediator = c("TSG", "OCG", "OCG", "OCG"),
    logFC_target = c(1.94917640749552, 1.36435629696216, 1.36486453774995, 1.36486453774995),
    TF = c(NA, NA, "IRF1", "IRF2"),
    InteractionType = c(NA, NA, "Activation", "Activation"),
    PMID = c(NA, NA, "18694960", "15778351"),
    tf_mutation = as.character(c(NA, NA, NA, NA)),
    stab_class = as.character(c(NA, NA, NA, NA)),
    in_MAVISp = c(FALSE,FALSE,FALSE,FALSE),
    mutation_available = c(FALSE,FALSE,FALSE,FALSE),
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



