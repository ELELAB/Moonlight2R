#' Moonlight2R
#'
#' @description Moonlight2R is a package designed for the identification of cancer driver genes. Please see the
#' documentation on our Bioconductor page for more details: https://www.bioconductor.org/packages/release/bioc/html/MoonlightR.html
#'
#' If you experience issues with the package, please open an Issue on our GitHub repository: https://github.com/ELELAB/Moonlight2R
#'
#' If you use this package in your research, please cite this paper: https://doi.org/10.1038/s41467-019-13803-0
#'
#' @name Moonlight2R
#' @docType package
NULL

#' Differentially expressed genes
#'
#'@description  A matrix containing differentially expressed genes between lung cancer and normal
#' samples found using TCGA-LUAD data and TCGAbiolinks. 
#'@details The matrix contains the differentially expressed genes in rows and log2 fold change
#' and FDR values in columns.
#'
#'@docType data
#'@usage data(DEGsmatrix)
#'@name DEGsmatrix
#'@aliases DEGsmatrix
#'@return A 3390x5 matrix
#'
#'@format A 3390x5 matrix
#'
"DEGsmatrix"


#' Cancer-related biological processes
#'
#'@description  A dataset containing information about 101 cancer-related biological processes. 
#'@details The dataset contains a list of the 101 biological processes which includes genes
#' playing a role in each biological processes including literature findings of the genes' function
#' in the biological processes.
#'
#'@docType data
#'@usage data(DiseaseList)
#'@name DiseaseList
#'@aliases DiseaseList
#'@return A list of 101 elements
#'
#'@format A list of 101 elements
#'
"DiseaseList"


#' Information about genes
#'
#'@description  A matrix containing information about 20038 genes 
#' including their gene description, location and family
#'@details The matrix contains the genes in rows and description, 
#' location and family in columns.
#'
#'@docType data
#'@usage data(EAGenes)
#'@name EAGenes
#'@aliases EAGenes
#'@return A 20038x5 matrix
#'
#'@format A 20038x5 matrix
#'
"EAGenes"

#' Information on GEO and TCGA data
#'
#'@description  A matrix that provides the GEO dataset
#'   matched to one of 18 TCGA cancer types
#'@details The matrix contains the cancer types in rows and 
#' information about sample type from both TCGA and GEO in columns.
#'
#'@docType data
#'@usage data(GEO_TCGAtab)
#'@name GEO_TCGAtab
#'@aliases GEO_TCGAtab
#'@return A 18x12 matrix
#'
#'@format A 18x12 matrix
#'
"GEO_TCGAtab"

#' Gene expression data from TCGA-LUAD
#'
#'@description  A matrix that provides processed gene expression 
#' data (obtained from RNA seq) from the TCGA-LUAD project
#'@details The matrix contains the genes in rows and samples in 
#' columns. The data has been downloaded and processed using 
#' TCGAbiolinks.
#'
#'@docType data
#'@usage data(dataFilt)
#'@name dataFilt
#'@aliases dataFilt
#'@return A 3000x20 matrix
#'
#'@format A 3000x20 matrix
#'
"dataFilt"

#' Gene regulatory network
#'
#'@description  The output of the GRN function which finds connections 
#' between genes.
#'@details The input to the GRN is the differentially expressed genes and #' the gene expression data.
#'
#'@docType data
#'@usage data(dataGRN)
#'@name dataGRN
#'@aliases dataGRN
#'@return A list of 2 elements where the first element is a 2x613 matrix #' and the second element is a vector of length 2
#'
#'@format A list of 2 elements where the first element is a 2x613 matrix 
#' and the second element is a vector of length 2
#'
"dataGRN"

#' Upstream regulator analysis
#'
#'@description  The output of the URA function which carries out the 
#' upstream regulator analysis
#'@details The input to URA is the output of GRN and a list of biological 
#' processes and the differentially expressed genes
#'
#'@docType data
#'@usage data(dataURA)
#'@name dataURA
#'@aliases dataURA
#'@return A 100x2 matrix
#'
#'@format A 100x2 matrix
#'
"dataURA"


#' Information of known cancer driver genes from COSMIC
#'
#'@description  A list of known cancer driver genes from COSMIC
#'@details The list contains two elements: a vector of known tumor #' suppressors and a vector of known oncogenes 
#'
#'@docType data
#'@usage data(knownDriverGenes)
#'@name knownDriverGenes
#'@aliases knownDriverGenes
#'@return A list containing two elements where the first element is a 
#' character vector of 55 and the second element is a character vector of #' 84
#'
#'@format A list containing two elements where the first element is a 
#' character vector of 55 and the second element is a character vector of #' 84
#'
"knownDriverGenes"

#' List of oncogenic mediators of 5 TCGA cancer types
#'
#'@description  A list of oncogenic mediators of 5 TCGA cancer types: 
#' BLCA, BRCA, LUAD, READ and STAD
#'@details Each element in the list contains differentially expressed 
#' genes and output from the URA and PRA functions
#'
#'@docType data
#'@usage data(listMoonlight)
#'@name listMoonlight
#'@aliases listMoonlight
#'@return A list containing 5 elements where each element contains 
#' differentially expressed genes and output from the URA and PRA 
#' functions of 5 TCGA cancer types
#'
#'@format A list containing 5 elements where each element contains 
#' differentially expressed genes and output from the URA and PRA 
#' functions of 5 TCGA cancer types
#'
"listMoonlight"

#' Information of growing/blocking characteristics of 101 biological
#' processes
#'
#'@description  A matrix with biological processes in rows and the cancer #' growing or blocking effect of the process in columns
#'@details For each biological processes the cancer growing/blocking 
#' effect is indicated 
#'
#'@docType data
#'@usage data(tabGrowBlock)
#'@name tabGrowBlock
#'@aliases tabGrowBlock
#'@return A 101x3 matrix
#'
#'@format A 101x3 matrix
#'
"tabGrowBlock"


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
#'@description The predicted TSGs and OCGs and their moonlight gene z-score based on the small sample TCGA-LUAD data.
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
#'@description An examplary MAF file from TCGA on lung cancer LUAD. It contains 500 randomly selected mutations. 
#'
#'@docType data
#'@usage data(dataMAF)
#'@name dataMAF
#'@aliases dataMAF
#'@return A 500x141 matrix.
#'@format A 500x141 matrix. 
#'
#'
"dataMAF"

#' Output example from the function Driver Mutation Analysis
#'
#'@description The predicted driver genes, which have at least one driver mutation. 
#'
#'@docType data
#'@usage data(dataDMA)
#'@name dataDMA
#'@aliases dataDMA
#'@return A list of two, containing 0 tumor-suppressor and 1 oncogene.
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
#'@return A 12x15 matrix.
#'@format A 12x15 matrix. 
#'
#'
"Oncogenic_mediators_mutation_summary"
