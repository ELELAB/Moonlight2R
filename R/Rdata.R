#' Promoters
#'
#'@description Experimentially verified promoter sites by J. Michael Cherry, Stanford.
#' Downloaded from the ENCODE identifier ENCSR294YNI.
#' It contains chromosome, start and end sites of promoters.
#'
#'@docType data
#'@usage data(EncodePromoters)
#'@name EncodePromoters
#'@aliases EncodePromoters
#'@return A 84738x6 table
#'@format A tibble with no columnnames or rownames.
#'\enumerate{
#'\item The first column is chromosome eg. chr1
#'\item The second column is start position eg. 10451
#'\item The third column is end position eg. 10563
#'}
#'
#'@source \url{https://www.encodeproject.org/}
#'@references  ENCODE identifier: ENCSR294YNI
#'
#'Luo Y, Hitz BC, Gabdank I, Hilton JA, Kagda MS, Lam B, Myers Z, Sud P, Jou J,
#'Lin K, Baymuradov UK, Graham K, Litton C, Miyasato SR, Strattan JS, Jolanki O,
#'Lee JW, Tanaka FY, Adenekan P, O'Neill E, Cherry JM.
#'New developments on the Encyclopedia of DNA Elements (ENCODE) data portal.
#'Nucleic Acids Res. 2020 Jan 8;48(D1):D882-D889. doi: 10.1093/nar/gkz1062.
#'PMID: 31713622; PMCID: PMC7061942.
#'
"EncodePromoters"

#' Level of Consequence: Protein
#'
#'@description  A dataset binary dataset describing if a mutation of a certain class and
#' type possibly have an effect on protein structure or function.
#'@details The values are binary: 0 no effect is possible, 1 an effect is possible.
#'
#'See supplementary material for details.
#'
#'@docType data
#'@usage data(LOC_protein)
#'@name LOC_protein
#'@aliases LOC_protein
#'@return A 18x7 table
#'
#'@format A 18x7 table
#'@references  paper
#'
"LOC_protein"

#' Level of Consequence: Translation
#'
#'@description  A dataset describing if a mutation of a certain class and
#' type possibly have an effect on peptide level.
#'@details The values are binary: 0 no effect is possible, 1 an effect is possible.
#'
#'See supplementary material for details.
#'
#'@docType data
#'@usage data(LOC_translation)
#'@name LOC_translation
#'@aliases LOC_translation
#'@return A 18x7 table
#'
#'@format A 18x7 table
#'@references  paper
#'
"LOC_translation"


#' Level of Consequence: Transcription
#'
#'@description  A dataset describing if a mutation of a certain class and
#' type possibly have an effect on transcript level.
#'@details The values are binary: 0 no effect is possible, 1 an effect is possible.
#'
#'See supplementary material for details.
#'
#'@docType data
#'@usage data(LOC_transcription)
#'@name LOC_transcription
#'@aliases LOC_transcription
#'@return A 18x7 table
#'
#'@format A 18x7 table
#'@references paper
#'
"LOC_transcription"

#' Network of Cancer Genes 7.0
#'
#'@description  A dataset retrived from Network of Cancer Genes 7.0
#'@details The NCG_driver is reported as a OCG or TSG when at least one of three
#' three databases have documented it. These are cosmic gene census (cgc),
#' vogelstein et al. 2013 or saito et al. 2020. The NCG_driver is reported as a
#' candidate, when literature support the gene as a cancer driver.
#'@docType data
#'@usage data(NCG)
#'@name NCG
#'@aliases NCG
#'@return A 3347x7 table
#'
#'@format The format have been rearranged from the original.
#' <symbold>|<NCG_driver>|<NCG_cgc_annotation>|<NCG_vogelstein_annotation>|
#' <NCG_saito_annotation>|<NCG_pubmed_id>
#'@source \url{http://ncg.kcl.ac.uk/}
#'@references  Comparative assessment of genes driving cancer and somatic
#'evolution in non-cancer tissues: an update of the Network of Cancer Genes (NCG)
#'resource.
#'Dressler L., Bortolomeazzi M., Keddar M.R., Misetic H., Sartini G.,
#'Acha-Sagredo A., Montorsi L., Wijewardhane N., Repana D., Nulsen J.,
#'Goldman J., Pollit M., Davis P., Strange A., Ambrose K. and Ciccarelli F.D.
#'
"NCG"

#' Output example from function Pattern Recognition Analysis
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