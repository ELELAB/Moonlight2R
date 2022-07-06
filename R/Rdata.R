#' Output example from function the Pattern Recognition Analysis
#'
#'@description The predicted TSGs and OCGs and their moonlight gene z-score.
#'The PRA() were run with expert-based approach with apoptosis and proliferation of cells. 
#'
#'@docType data
#'@usage data(dataPRA)
#'@name dataPRA
#'@aliases dataPRA
#'@return A list of two.
#'@format A list of two. 
#'

#'
"dataPRA"

#' Mutation data from TCGA LUAD
#'
#'@description An examplary MAF file from TCGA on lung cancer LUAD. It contains 3000 randomly selected mutations. 
#'
#'@docType data
#'@usage data(dataMAF)
#'@name dataMAF
#'@aliases dataMAF
#'@return A 100x141 matrix.
#'@format A 100x141 matrix. 
#'
#'
"dataMAF"

#' Output example from the function Driver Mutation Analysis
#'
#'@description The predicted driver genes, which have at least one driver mutation. 
#'
#'@docType data
#'@usage data(dataMAF)
#'@name dataMAF
#'@aliases dataMAF
#'@return A list of two, containing 23 tumor-suppressors and 39 oncogenes.
#'@format A list of two. 
#'
#'
"dataDMA"

#' Cscape-somatic annotations of TCGA-LUAD 
#'@description Output from DMA. 
#'This contains the cscape-somatic annotations for all differentially expressed genes 
#'
#'@docType data
#'@usage data(cscape_somatic_output)
#'@name cscape_somatic_output
#'@aliases cscape_somatic_output
#'@return A 645x7 matrix.
#'@format A 645x7 matrix. 
#'
#'
"cscape_somatic_output"

#' Differentially expressed genes's Mutations 
#'@description Output from DMA. 
#'This contains the differentially expressed genes's mutations and all annotations
#'generated in DMA() on the TCGA-LUAD project. 
#'
#'@docType data
#'@usage data(DEG_Mutations_Annotations)
#'@name DEG_Mutations_Annotations
#'@aliases DEG_Mutations_Annotations
#'@return A 3561x173 matrix.
#'@format A 3561x173 matrix. 
#'
#'
"DEG_Mutations_Annotations"

#' Oncogenic Mediators Mutation Summary  
#'@description Output from DMA. 
#'This contains the oncogenic mediator from the TCGA-LUAD project, and their mutation
#'assesments summarised based on CSCape-somatic and Level of Consequence. 
#'
#'@docType data
#'@usage data(Oncogenic_mediators_mutation_summary)
#'@name Oncogenic_mediators_mutation_summary
#'@aliases Oncogenic_mediators_mutation_summary
#'@return A 1105x15 matrix.
#'@format A 1105x15 matrix. 
#'
#'
"Oncogenic_mediators_mutation_summary"