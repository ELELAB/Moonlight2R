#' @title plotURA: Upstream regulatory analysis heatmap plot
#' @description
#' This function visualizes the URA in a heatmap
#' @param dataURA output URA function
#' @param additionalFilename figure name
#' @importFrom gplots heatmap.2
#' @importFrom gplots redblue
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @return heatmap
#' @export
#' @examples
#' data(dataURA)
#' data(DiseaseList)
#' data(tabGrowBlock)
#' data(knownDriverGenes)
#' dataDual <- PRA(dataURA = dataURA,
#' BPname = c("apoptosis","proliferation of cells"),
#' thres.role = 0)
#' TSGs_genes <- names(dataDual$TSG)
#' OCGs_genes <- names(dataDual$OCG)
#' plotURA(dataURA = dataURA[c(TSGs_genes, OCGs_genes),],additionalFilename = "_example")
plotURA <- function(dataURA,
                    additionalFilename = "URAplot") {

  # Check user input

  if (!is(dataURA, "matrix")) {
    stop("dataURA must be a numeric matrix")
  }

  if (nrow(dataURA) < 2 | ncol(dataURA) < 2) {
    stop("dataURA must have at least two rows and two columns")
  }

  if (!is(additionalFilename, "NULL") & !is(additionalFilename, "character")) {
    stop("additionalFilename must be either NULL or a character vector 
containing part of the filename of plot")
  }

  if (nrow(dataURA) > 70) {
    cexRow <- 0.1 + 1 / (2 * nrow(dataURA)) * (log10(nrow(dataURA)))
  } else {
    cexRow <- 0.2 + 1 / (10 * log10(nrow(dataURA)))
  }

  if (!is(additionalFilename, "NULL")) {
    pdf(file = paste0("plotURA", additionalFilename, ".pdf"))
  }

  par(oma = c(6, 4, 4, 2))

  heatmap.2(dataURA,
            trace = "none",
            col = rev(gplots::redblue(128)),
            Colv = TRUE,
            dendrogram = "row",
            notecex = 10 / nrow(dataURA),
            cexCol = 0.2 + 1 / (3 * log10(ncol(dataURA))),
            cexRow = cexRow,
            cellnote = signif(dataURA, 4),
            notecol = "black",
            key = TRUE,
            keysize = 2)

  if (!is(additionalFilename, "NULL")) {
    dev.off()
  }
}

