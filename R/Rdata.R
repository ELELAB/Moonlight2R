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

#' Functional enrichment analysis
#'
#'@description  The output of the FEA function which does enrichment analysis
#'@details The input to the FEA is the differentially expressed genes.
#'
#'@docType data
#'@usage data(dataFEA)
#'@name dataFEA
#'@aliases dataFEA
#'@return A dataframe of dimension 101x7
#'
#'@format A dataframe of dimension 101x7
#'
"dataFEA"

#' Gene regulatory network
#'
#'@description  The output of the GRN function which finds connections 
#' between genes.
#'@details The input to the GRN is the differentially expressed genes and the gene expression data.
#'
#'@docType data
#'@usage data(dataGRN)
#'@name dataGRN
#'@aliases dataGRN
#'@return A list of 2 elements where the first element is a 23x613 matrix and the second element is a vector of length 23
#'
#'@format A list of 2 elements where the first element is a 23x613 matrix 
#' and the second element is a vector of length 23
#'
"dataGRN"

#' Gene regulatory network
#'
#'@description  The output of the GRN function which finds connections 
#' between genes where the noise is set to 0 for testing reproducibility purposes.
#'@details The input to the GRN is the differentially expressed genes and the gene expression data.
#'
#'@docType data
#'@usage data(dataGRN_no_noise)
#'@name dataGRN_no_noise
#'@aliases dataGRN_no_noise
#'@return A list of 2 elements where the first element is a 23x613 matrix and the second element is a vector of length 23
#'
#'@format A list of 2 elements where the first element is a 23x613 matrix 
#' and the second element is a vector of length 23
#'
"dataGRN_no_noise"

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
#'@return A 23x2 matrix
#'
#'@format A 23x2 matrix
#'
"dataURA"

#' Upstream regulator analysis
#'
#'@description  The output of the URA function which carries out the 
#' upstream regulator analysis
#'@details This URA data is used to showcase some of the visualization
#' functions
#'
#'@docType data
#'@usage data(dataURA_plot)
#'@name dataURA_plot
#'@aliases dataURA_plot
#'@return A 12x2 matrix
#'
#'@format A 12x2 matrix
#'
"dataURA_plot"

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

#' Literature search of driver genes
#'
#'@description  A tibble containing results of literature search where predicted
#' driver genes stored in dataDMA were queried for their role as drivers in
#' PubMed 
#'@details The tibble contains PubMed IDs, doi, title, abstract, year of publication,
#' keywords, and total number of publications for the genes.
#'
#'@docType data
#'@usage data(dataGLS)
#'@name dataGLS
#'@aliases dataGLS
#'@return A 13x8 tibble.
#'
#'@format A 13x8 tibble.
#'
"dataGLS"


#' Methylation evidence table to define driver genes
#'
#'@description  A tibble containing combinations of methylation states
#' of probes used to define driver genes 
#'@details The tibble contains a value of 1 if a probe is found that
#' is either hypo-, hyper-, dualmethylated or not methylated. This is
#' compared with Moonlight's predictions of role of oncogenic mediators
#' to define driver genes based on methylation evidence. 
#'
#'@docType data
#'@usage data(MetEvidenceDriver)
#'@name MetEvidenceDriver
#'@aliases MetEvidenceDriver
#'@return A 30x6 tibble.
#'
#'@format A 30x6 tibble.
#'
"MetEvidenceDriver"


#' Methylation data matrix from TCGA-LUAD project
#'
#'@description  A data matrix containing methylation data from TCGA-LUAD
#' where CpG probes are in rows and samples are in columns. 
#'@details The CpG probes are in rows and samples are in columns.
#'
#'@docType data
#'@usage data(dataMethyl)
#'@name dataMethyl
#'@aliases dataMethyl
#'@return A 73x27 matrix.
#'
#'@format A 73x27 matrix.
#'
"dataMethyl"


#' Gene expression data from TCGA-LUAD
#'
#'@description  A matrix that provides processed gene expression 
#' data (obtained from RNA seq) from the TCGA-LUAD project
#'@details The matrix contains the genes in rows and samples in 
#' columns. The data has been downloaded and processed using 
#' TCGAbiolinks. dataFiltCol is identical to dataFilt with the
#' exception that in dataFiltCol the sample barcodes have been 
#' shortened to contain only patient barcodes.
#'
#'@docType data
#'@usage data(dataFiltCol)
#'@name dataFiltCol
#'@aliases dataFiltCol
#'@return A 3000x20 matrix
#'
#'@format A 3000x20 matrix
#'
"dataFiltCol"


#' Sample annotations of TCGA-LUAD project 
#'
#'@description  A matrix that annotates LUAD samples as either
#' cancer or normal 
#'@details The matrix contains two columns: "primary" which
#' contains patient barcodes of TCGA-LUAD and "sample.type"
#' which denotes if the sample is either a "Cancer" or 
#' "Normal" sample
#'
#'@docType data
#'@usage data(LUAD_sample_anno)
#'@name LUAD_sample_anno
#'@aliases LUAD_sample_anno
#'@return A 23x2 matrix
#'
#'@format A 23x2 matrix
#'
"LUAD_sample_anno"


#' Output example from GMA function 
#'
#'@description  The predicted driver genes based on methylation
#' evidence 
#'@details The data is a list of two elements where each 
#' element represents predicted oncogenes and tumor suppressors
#'
#'@docType data
#'@usage data(dataGMA)
#'@name dataGMA
#'@aliases dataGMA
#'@return A list of length two
#'
#'@format A list of length two
#'
"dataGMA"


#' Output example from GMA function 
#'
#'@description  The object, a list, that was returned from
#' running the EpiMix function and is one of the outputs from the 
#' GMA function. 
#'@details The data is a list of nine elements which is outputted
#' from the EpiMix function  
#'
#'@docType data
#'@usage data(EpiMix_Results_Regular)
#'@name EpiMix_Results_Regular
#'@aliases EpiMix_Results_Regular
#'@return A list of length nine
#'
#'@format A list of length nine
#'
"EpiMix_Results_Regular"


#' Output example from GMA function 
#'
#'@description  Output file from running the GMA function which
#' is a summary of the oncogenic mediators and their sum of
#' methylated CpG probes together with the evidence level of
#' their role as driver gene. 
#'@details The data is a table where each row is an oncogenic
#' mediator and the columns represent the predicted driver role
#' and the sum of hypo-, hyper-, and dualmethylated CpG sites. 
#'
#'@docType data
#'@usage data(Oncogenic_mediators_methylation_summary)
#'@name Oncogenic_mediators_methylation_summary
#'@aliases Oncogenic_mediators_methylation_summary
#'@return A 25x7 tibble
#'
#'@format A 25x7 tibble
#'
"Oncogenic_mediators_methylation_summary"


#' Output example from GMA function 
#'
#'@description  Output file from running GMA function which
#' is a summary of DEGs and associated CpG probes 
#'@details The data is a table where each row is a CpG probe
#' in a DEG. Various annotations such as start/end site of
#' CpG probe, promoter/enhancer annotations, NCG annotations
#' are included in the table.
#'
#'@docType data
#'@usage data(DEG_Methylation_Annotations)
#'@name DEG_Methylation_Annotations
#'@aliases DEG_Methylation_Annotations
#'@return A 3435x35 tibble
#'
#'@format A 3435x35 tibble
#'
"DEG_Methylation_Annotations"

