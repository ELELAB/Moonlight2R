#' FEA
#'
#' This function carries out the functional enrichment analysis (FEA).
#' The analysis can be performed using as input a data matrix from Differential Expression Analysis of genes or Differential Abundance Analysis of proteins.
#' The analysis can be performed using the ORA (Over Representation Analysis)-based method which implements Fisher's Exact Test or the pre-ranked FCS(Functional Class Scoring)-based method
#' which implements fgsea (Fast Gene Set Enrichment Analysis).
#' @param BPname BPname biological process such as "proliferation of cells", "ALL" (default) if FEA should be carried out for all 101 biological processes
#' @param DiffMatrix output from DEA such as dataDEGs or from DAA such as differential abundace matrix
#' @param method method to perform Functional Enirchment analysis. Must be either "ora" (Over-Representation Analysis) or "fgsea" (Fast Gene Set Enrichment Analysis)
#' @param seed internal seed before the 'fgsea' function is called to ensure reproducibility of the permutation-based p-values. Default is 123
#' @importFrom stats p.adjust
#' @importFrom fgsea fgsea
#' @return matrix from FEA
#' @export
#' @examples
#' data(DEGsmatrix)
#' data(DiseaseList)
#' data(EAGenes)
#' DEGsmatrix <- DEGsmatrix[seq.int(2), ]
#' dataFEA <- FEA(DiffMatrix = DEGsmatrix, BPname = "apoptosis")
FEA <- function(BPname = NULL,
                DiffMatrix,
                method = "ora",
                seed = 123) {

  # List of variable names
  variables_to_check <- c("DiseaseList", "EAGenes")

  # Check and load variables if they do not exist
  for (variable_name in variables_to_check) {
    if (! variable_name %in% names(.GlobalEnv)) {
      data(list=c(variable_name))
    }
  }

  DiseaseList <- get("DiseaseList")
  EAGenes <- get("EAGenes")

  # Check user input

  if ("logFC" %in% colnames(DiffMatrix) == FALSE) {
    stop("The input data table must contain a column called logFC.")
  }

  if (.row_names_info(DiffMatrix) < 0) {
    stop("Row names were generated automatically. The input data table needs to
    have the gene names as rownames. Double check that genes are rownames.")
  }

  if (!is(BPname, "NULL") && all(BPname %in% names(DiseaseList)) == FALSE) {
    stop("BPname should be NULL or a character vector containing one or more
    BP(s) among possible BPs stored in the DiseaseList object.")
  }
  if (!method %in% c("ora", "fgsea")) {
    stop("Invalid method name. Must be either ora or fgsea")

  }

  if (is(BPname, "NULL")) {
    lf2 <- names(DiseaseList)
  } else {
    lf2 <- BPname
  }
  TableDiseasesNew <- NULL

  pb <- txtProgressBar(min = 0, max = length(DiseaseList), style = 3)
 
   
  if (method == "fgsea") {

    pvals <- pmax(DiffMatrix$PVal, .Machine$double.xmin)
    rankings <- sign(DiffMatrix$logFC)*(-log10(pvals))
    names(rankings)  <- rownames(DiffMatrix)
    # sorting strategy to avoid different sorting in case of tie ranks
    rankings <- rankings[order(-rankings, names(rankings))]
    
    
    pathwayNamesList <- list()
    for (k in seq_along(lf2)) {

      selected_diseases <- as.data.frame(DiseaseList[[which(names(DiseaseList) == lf2[k])]])
      name <- lf2[k]
      values <- selected_diseases$ID
      pathwayNamesList[[name]] <- values
    }

    pathway_sizes <- sapply(DiseaseList, function(pathway) length(unique(pathway$ID)))
    max_pathway_size <- max(pathway_sizes)
    
    set.seed(seed)
    fgseaRes <- fgsea(pathways = pathwayNamesList, stats = rankings, scoreType = 'std', minSize = 1, maxSize = max_pathway_size, nproc = 1)
    
    bp_score_collection <- list()
    for (k in seq_along(lf2)) {
      
      setTxtProgressBar(pb, k)

      res <- as.data.frame(matrix(0, nrow = 1, ncol = 2,
                                  dimnames = list(1, c("pathway",
                                                     "Moonlight.Z.score"))))
      
      selected_diseases <- as.data.frame(DiseaseList[[which(names(DiseaseList) == lf2[k])]])
      
      res$pathway <- lf2[k]
      
      Zscore <- .compute_moonlight_zscore(selected_diseases)
      
      res$Moonlight.Z.score <- Zscore
      
      bp_score_collection[[k]] <- res
    
    }

    close(pb)

    bp_score_collection_merged <- bind_rows(bp_score_collection)

    TableDiseasesNew <- fgseaRes %>% left_join(bp_score_collection_merged, by = "pathway")
    TableDiseasesNew <- TableDiseasesNew %>% select(-log2err)
    TableDiseasesNew <- TableDiseasesNew %>% 
      rename("Diseases.or.Functions.Annotation" = pathway, "p.value" = pval) %>% 
      select("Diseases.or.Functions.Annotation", "Moonlight.Z.score","p.value", "padj", "ES", "NES", "size", "leadingEdge")

  
  }
  else if (method == "ora") {

    for (k in seq.int(lf2)) {

      setTxtProgressBar(pb, k)

      res <- as.data.frame(matrix(0, nrow = 1, ncol = 7,
                                  dimnames = list(1, c("Diseases.or.Functions.Annotation",
                                                      "p.Value", "Moonlight.Z.score",
                                                      "commonNg",
                                                      "FunctionNg",
                                                      "Delta",
                                                      "Molecules"))))

      GeneList <- data.frame(PROBE_ID = rownames(DiffMatrix),
                            logFC = DiffMatrix$logFC,
                            stringsAsFactors = FALSE)


      rownames(GeneList) <- GeneList$PROBE_ID

      selected_diseases <- as.data.frame(DiseaseList[[which(names(DiseaseList) == lf2[k])]])

      selected_diseases$ID <- selected_diseases$Genes.in.dataset

      res$commonNg <- length(intersect(GeneList$PROBE_ID, selected_diseases$ID))

      res$FunctionNg <- nrow(selected_diseases)

      res$Diseases.or.Functions.Annotation <- lf2[k]

      allgene <- unique(c(as.character(unique(EAGenes[, "ID"])),
                        GeneList$PROBE_ID,
                        selected_diseases$ID))

      seta <- allgene %in% GeneList$PROBE_ID
      setb <- allgene %in% selected_diseases$ID

      if (res$commonNg > 1) {
        ft <- fisher.test(seta, setb)
        FisherpvalueTF <- ft$p.value
        res$p.Value <- FisherpvalueTF
      } else {
        res$p.Value <- 1
      }

      GeneList <- GeneList[GeneList$PROBE_ID %in% selected_diseases$ID, ]

      selected_diseases <- selected_diseases[selected_diseases$ID %in% GeneList[, "PROBE_ID"], ]
      selected_diseases[, "Exp.Log.Ratio"] <- gsub(",", ".", selected_diseases[, "Exp.Log.Ratio"])
      selected_diseases[, "Exp.Log.Ratio"] <- as.numeric(selected_diseases[, "Exp.Log.Ratio"])

      rownames(selected_diseases) <- selected_diseases$ID

      res$Molecules <- paste0(GeneList$PROBE_ID, collapse = ",")

      for (idx in seq.int(nrow(selected_diseases))) {

        currTR <- selected_diseases$ID[idx]

        if (length(grep("Increases", selected_diseases[currTR, "Findings"])) == 1) {

          if (sign(GeneList[currTR, "logFC"]) > 0) {
            selected_diseases[currTR, "Prediction..based.on.expression.direction."] <- "Increased"
          } else if (sign(GeneList[currTR, "logFC"]) < 0) {
            selected_diseases[currTR, "Prediction..based.on.expression.direction."] <- "Decreased"
          }
        }

        if (length(grep("Decreases", selected_diseases[currTR, "Findings"])) == 1) {
          if (sign(GeneList[currTR, "logFC"]) < 0) {
            selected_diseases[currTR, "Prediction..based.on.expression.direction."] <- "Increased"
          } else if (sign(GeneList[currTR, "logFC"]) > 0) {
            selected_diseases[currTR, "Prediction..based.on.expression.direction."] <- "Decreased"
          }
        }
      }

      Zscore <- .compute_moonlight_zscore(selected_diseases)

      res$Moonlight.Z.score <- Zscore
      TableDiseasesNew <- rbind(TableDiseasesNew, res)
    }

    close(pb)

    TableDiseasesNew <- cbind(TableDiseasesNew, FDR = p.adjust(TableDiseasesNew$p.Value,
                                                             method = "fdr"))
    TableDiseasesNew <- subset(TableDiseasesNew, select = c("Diseases.or.Functions.Annotation",
                                                          "Moonlight.Z.score",
                                                          "p.Value",
                                                          "FDR",
                                                          "commonNg",
                                                          "FunctionNg",
                                                          "Molecules"))
    
  }

  return(TableDiseasesNew)
}

#' compute_moonlight_zscore
#' 
#' This internal function compute the Moonlight Z-score for each enriched pathway 
#' extracted from the Functional Enrichment Analysis
#' @keywords internal
#' @param selected_diseases A dataframe containing pathway's genes expression direction and changes 
#' @return Moonlight Z-score A numeric score that measures how well the direction of change agrees between observed expression changes and literature findings
#' @noRd
.compute_moonlight_zscore <- function(selected_diseases) {
  
  Correlation <- matrix(0, nrow(selected_diseases), 1)

  selected_diseases <- cbind(selected_diseases, Correlation)

  if (length(grep("Decreases", selected_diseases$Findings)) != 0) {

        selected_diseases[grep("Decreases", selected_diseases$Findings), "Findings"] <- -1
        selected_diseases[grep("Increases", selected_diseases$Findings), "Findings"] <- 1
        selected_diseases[grep("Affects", selected_diseases$Findings), "Findings"] <- 0
        selected_diseases[, "Findings"] <- as.numeric(selected_diseases[, "Findings"])

        selected_diseases[, "Exp.Log.Ratio"] <- gsub(",", ".", selected_diseases[, "Exp.Log.Ratio"])
        selected_diseases[, "Exp.Log.Ratio"] <- as.numeric(selected_diseases[, "Exp.Log.Ratio"])
        
        PredictionIncreased <- which(sign(selected_diseases$Exp.Log.Ratio) == selected_diseases$Findings)
        PredictionDecreased <- which(sign(selected_diseases$Exp.Log.Ratio) != selected_diseases$Findings)
        PredictionAffected <- which(sign(selected_diseases$Findings) == 0)

        selected_diseases[PredictionIncreased, "Correlation"] <- 1
        selected_diseases[PredictionDecreased, "Correlation"] <- -1
        selected_diseases[PredictionAffected, "Correlation"] <- 0

        Zscore <- sum(selected_diseases$Correlation) / sqrt( length(PredictionIncreased) + length(PredictionDecreased))
      } else {
        Zscore <- 0
      }

      return(Zscore)

}