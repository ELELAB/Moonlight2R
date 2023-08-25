#' LPA
#'
#' This function carries out the literature phenotype analysis (LPA)
#' @param dataDEGs is output from DEA
#' @param BP is biological process
#' @param BPlist is list of genes annotated in BP
#' @import RISmed
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#' @return table with number of pubmed that affects, increase or decrase genes annotated in BP
#' @examples
#' data('DEGsmatrix')
#' data('DiseaseList')
#' BPselected <- c("apoptosis")
#' BPannotations <- DiseaseList[[match(BPselected, names(DiseaseList))]]$ID
LPA <- function(dataDEGs,
                BP,
                BPlist) {
  
  # Check user input
  
  if (.row_names_info(dataDEGs) < 0) {
    stop("Row names were generated automatically. The input DEG table needs to
         have the gene names as rownames. Double check that genes are rownames.")
  }
  
  if ("logFC" %in% colnames(dataDEGs) == FALSE) {
    stop("The input DEG table must contain a column called logFC.")
  }
  
  if (all(BP %in% names(DiseaseList)) == FALSE) {
    stop("BPname should be a character vector containing one or more BP(s)
         among possible BPs stored in the DiseaseList object.")
  }
  
  if (!is.character(BPlist)) {
    stop("BPlist must be a character vector of genes")
  }
  
  BPgenesDEGs <- intersect(BPlist, rownames(dataDEGs))
  
  dataDEGsBP <- dataDEGs[BPgenesDEGs, ]
  
  DiseaseMN <- matrix(0, nrow(dataDEGsBP), 7)
  
  colnames(DiseaseMN) <- c("ID",
                           "Genes.in.dataset",
                           "Prediction..based.on.expression.direction.",
                           "Exp.Log.Ratio",
                           "Findings",
                           "PubmedDecreases",
                           "PubmedIncreases")
  
  DiseaseMN <- as.data.frame(DiseaseMN)
  DiseaseMN$Prediction..based.on.expression.direction. <- rep("Decreased", nrow(dataDEGsBP))
  DiseaseMN$ID <- BPgenesDEGs
  DiseaseMN$Genes.in.dataset <- BPgenesDEGs
  DiseaseMN$Exp.Log.Ratio <- dataDEGs[BPgenesDEGs, "logFC"]
  rownames(DiseaseMN) <- DiseaseMN$ID
  
  pb <- txtProgressBar(min = 0, max = nrow(DiseaseMN), style = 3)
  
  for (i in seq.int(nrow(DiseaseMN))) {
    
    setTxtProgressBar(pb, i)
    
    curG <- DiseaseMN$ID[i]
    
    keypubmed <- paste(curG, BP, "decreases")
    res <- EUtilsSummary(keypubmed, type = "esearch", db = "pubmed")
    fetch <- EUtilsGet(res, type = "efetch", db = "pubmed")
    RecordsDec <- length(fetch@PMID)
    DiseaseMN[curG, "PubmedDecreases"] <- RecordsDec
    
    keypubmed <- paste(curG, BP, "increases")
    res <- EUtilsSummary(keypubmed, type = "esearch", db = "pubmed")
    fetch <- EUtilsGet(res, type = "efetch", db = "pubmed")
    RecordsInc <- length(fetch@PMID)
    DiseaseMN[curG, "PubmedIncreases"] <- RecordsInc
    
    if (DiseaseMN[curG, "PubmedDecreases"] == DiseaseMN[curG, "PubmedIncreases"]) {
      DiseaseMN[curG, "Findings"] <- paste0("Affects (", RecordsInc, ")")
    }
    
    if (DiseaseMN[curG, "PubmedDecreases"] < DiseaseMN[curG, "PubmedIncreases"]) {
      DiseaseMN[curG, "Findings"] <- paste0("Increases (", RecordsInc, ")")
    }
    
    if (DiseaseMN[curG, "PubmedDecreases"] > DiseaseMN[curG, "PubmedIncreases"]) {
      DiseaseMN[curG, "Findings"] <- paste0("Decreases (", RecordsDec, ")")
    }
    
  }
  
  close(pb)
  
  return(DiseaseMN)
  
}

utils::globalVariables(c("DiseaseList"))

