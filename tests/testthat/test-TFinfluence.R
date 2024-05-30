# Test TFinfluence function

# Run example of TFinfluence function
data(dataPRA)
data(DEGsmatrix)
data(dataTRRUST)
data(dataMAF)
data(dataMAVISp)

TFresults <- TFinfluence(dataTRRUST = dataTRRUST,
            dataMAF = dataMAF,
            dataDEGs = dataDEGs,
            dataPRA = dataPRA,
            dataMAVISp = dataMAVISp,
            stabClassMAVISp = 'rosetta')

character_columns <- c('Target','Moonlight_Oncogenic_Mediator','TF','InteractionType','PMID','tf_mutation','stab_class')
numeric_columns <- c('Moonlight_gene_z_score','logFC_target')

# Test that output is as expected
test_that('Output of TFinfluence is a tibble', {
    expect_named(TFresults, c(character_columns, numeric_columns))
    expect_true(all(sapply(TFresults[character_columns], is.character)))
    expect_true(all(sapply(TFresults[numeric_columns], is.numeric)))
})


