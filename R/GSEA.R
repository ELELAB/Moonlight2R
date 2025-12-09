#' GSEA
#'
#' This function carries out the GSEA enrichment analysis.
#' @param DEGsmatrix DEGsmatrix output from DEA such as dataDEGs
#' @param top is the number of top BP to plot
#' @param plot if TRUE return a GSEA's plot
#' @import org.Hs.eg.db
#' @importFrom grDevices dev.list
#' @importFrom grDevices graphics.off
#' @importFrom clusterProfiler bitr
#' @importFrom DOSE gseDO
#' @importFrom DOSE simplot
#' @return return GSEA result
#' @export
#' @examples
#' data("DEGsmatrix")
#' DEGsmatrix_example <- DEGsmatrix[1:2,]
#' dataFEA <- GSEA(DEGsmatrix = DEGsmatrix_example)
GSEA <- function(DEGsmatrix,
                 top,
                 plot = FALSE) {

  # Check user input

  if (.row_names_info(DEGsmatrix) < 0) {
    stop("Row names were generated automatically. The input DEG table needs to
         have the gene names as rownames. Double check that genes are rownames.")
  }

  dataDEGsnew <- cbind(mRNA = rownames(DEGsmatrix), DEGsmatrix)

  eg <- as.data.frame(bitr(dataDEGsnew$mRNA,
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = "org.Hs.eg.db"))
  eg <- eg[!duplicated(eg$SYMBOL), ]

  dataDEGsFiltLevel <- dataDEGsnew
  dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$SYMBOL, ]
  dataDEGsFiltLevel <- dataDEGsFiltLevel[order(dataDEGsFiltLevel$mRNA, decreasing = FALSE), ]
  eg <- eg[order(eg$SYMBOL, decreasing = FALSE), ]

  dataDEGsFiltLevel$GeneID <- eg$ENTREZID

  dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))

  genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
  names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID

  genelistDEGs_sort <- sort(genelistDEGs, decreasing = TRUE)

  y <- gseDO(genelistDEGs_sort,
             nPerm = 100,
             minGSSize = 120,
             pvalueCutoff = 0.2,
             pAdjustMethod = "BH",
             verbose = FALSE)

  res <- as.matrix(summary(y))

  if (plot == TRUE) {
    topID <- res[1, 1]
    pdf("GSEAplot.pdf")
    DOSE::simplot(y, geneSetID = topID)

    if (!(is.null(dev.list()["RStudioGD"]))) {
      graphics.off()
    }
  }

  return(res)

}
