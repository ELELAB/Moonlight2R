#' @title moonlight pipeline
#' @description moonlight is a tool for identification of cancer driver genes. 
#' This function wraps the different steps of the complete analysis workflow.
#' @param BPname biological processes to use, if NULL: all processes will be used in analysis, RF for candidate; if not NULL the candidates for these processes will be determined (no learning)
#' @param dataDEGs table of differentially expressed genes
#' @param dataFilt matrix of gene expression data with genes in rows and samples in columns
#' @param Genelist Genelist
#' @param kNearest kNearest
#' @param nTF nTF
#' @param nGenesPerm nGenesPerm
#' @param nBoot nBoot
#' @param DiffGenes DiffGenes
#' @param thres.role thres.role
#' @param dataMAF A MAF file rda object for DMA
#' @param path_cscape_coding A character string to path of CScape-somatic coding file
#' @param path_cscape_noncoding A character string to path of CScape-somatic non-coding file
#' @export
#' @return table with cancer driver genes TSG and OCG.
#' @examples
#' \donttest{
#' drivers <- moonlight(dataDEGs = DEGsmatrix, 
#' dataFilt = dataFilt, 
#' BPname = c("apoptosis", "proliferation of cells"),
#' dataMAF = dataMAF, 
#' path_cscape_coding = "css_coding.vcf.gz", 
#' path_cscape_noncoding = "css_noncoding.vcf.gz")
#' }

moonlight <- function(dataDEGs, dataFilt, 
                      BPname = NULL, 
                      Genelist= NULL,
                      kNearest = 3, nGenesPerm = 2000, DiffGenes = FALSE,
                      nBoot = 400, nTF = NULL,thres.role = 0, 
                      dataMAF, 
                      path_cscape_coding, path_cscape_noncoding){

        # Check user input
  
        if (.row_names_info(dataDEGs) < 0) {
		stop("Row names were generated automatically. The input DEG table needs to have
			the gene names as rownames. Double check that genes are rownames.")
	}
  
	if (is.null(dim(dataFilt))) {
		stop("The expression data must be non-empty with genes in rows and samples in columns")
	}
  
	if (!is.null(BPname) && all(BPname %in% names(DiseaseList)) == FALSE) {
		stop("BPname should be NULL or a character vector containing one or more BP(s) 
			among possible BPs stored in the DiseaseList object.")
	}
  
	if (!is.null(Genelist) && !is.character(Genelist)) {
		stop("Genelist must be NULL or a character vector containing gene names")
	}
  
	if (!is.numeric(kNearest) | !is.numeric(nGenesPerm) | !is.numeric(nBoot) | !is.numeric(thres.role)) {
		stop("kNearest, nGenesPerm, nBoot, and thres.role must be numeric values")
	}
  
	if (!is.logical(DiffGenes)) {
		stop("DiffGenes must be either TRUE or FALSE")
	}
  
	if (is.null(dim(dataMAF))) {
		stop("The mutation data must be a non-empty table")
	}
  
	if (!is.character(path_cscape_coding) | !is.character(path_cscape_noncoding)) {
		stop("Paths to cscape coding and non-coding files must be character vectors")
	}
 
    	res <- NULL

        ### functional enrichment analysis -> carried out in URA, not necessary hear
        # print("-----------------------------------------")
        # print("Functional enrichment analysis")
        # print("-----------------------------------------")
        	
        # dataFEA <- FEA(BPname=BPname, DEGsmatrix = dataDEGs)

        ### gene regulatory network
        print("-----------------------------------------")
        print("Gene regulatory network")
        print("-----------------------------------------")

        #### parameter nTF for testing purposes
        if(is.null(nTF)){
            nTF <- nrow(dataDEGs)
        }

        if(is.null(Genelist)){
            Genelist <- rownames(dataDEGs)[seq.int(nTF)]
        }
        dataGRN <- GRN(TFs = Genelist, normCounts = dataFilt,
                       DEGsmatrix = dataDEGs,DiffGenes = FALSE,
                       nGenesPerm = nGenesPerm, kNearest = kNearest, nBoot = nBoot)

        ### upstream regulator analysis
        print("-----------------------------------------")
        print("Upstream regulator analysis")
        print("-----------------------------------------")

        dataURA <- URA(dataGRN = dataGRN, DEGsmatrix = dataDEGs, BPname = BPname)

        ### get TSG/OCG candidates using PRA
        print("-----------------------------------------")
        print("Get candidates")
        print("-----------------------------------------")
        listCandidates <- PRA(dataURA = dataURA, BPname = BPname, thres.role = thres.role)
        
        ### get TSG/OCG driver genes using DMA
        print("-----------------------------------------")
        print("Get driver genes")
        print("-----------------------------------------")
        driverGenes <- DMA(dataMAF = dataMAF, dataDEGs = dataDEGs, dataPRA = listCandidates, 
                           runCscape = TRUE, coding_file = path_cscape_coding, noncoding_file = path_cscape_noncoding, 
                           results_folder = "./DMAresults")
        

        res <- list(dataDEGs = dataDEGs,
                dataURA = dataURA, 
                listCandidates = listCandidates,
                driverGenes = driverGenes)

    return(res)
}

utils::globalVariables(c("DiseaseList"))
