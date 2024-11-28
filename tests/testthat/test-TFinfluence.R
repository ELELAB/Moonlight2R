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
    stringsAsFactors = FALSE # Ensure non-character columns are not treated as factors
)

TFresults <- TFinfluence(dataTRRUST = dataTRRUST,
            dataMAF = dataMAF,
            dataDEGs = DEGsmatrix,
            dataPRA = dataPRA,
            dataMAVISp = dataMAVISp,
            stabClassMAVISp = 'rosetta')

character_columns <- c('Target','Moonlight_Oncogenic_Mediator','TF','InteractionType','PMID','tf_mutation','stab_class')
numeric_columns <- c('Moonlight_gene_z_score','logFC_target')
all_columns <- c('Target', 'Moonlight_gene_z_score', 'Moonlight_Oncogenic_Mediator', 'logFC_target', 'TF', 'InteractionType', 'PMID', 'tf_mutation', 'stab_class')


# Test that output columns are as expected
test_that("Output of TFinfluence is a tibble", {
    expect_named(TFresults, all_columns)
    expect_true(all(sapply(TFresults[character_columns], is.character)))
    expect_true(all(sapply(TFresults[numeric_columns], is.numeric)))
})

# Test that the content of the output dataframe is as expected
test_that("Output of TFinfluence does not correspond to reference", {
    expect_equivalent(TFresults, reference_TFresults)
})



