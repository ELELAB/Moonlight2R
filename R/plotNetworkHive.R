#' @title plotNetworkHive: Hive network plot
#' @description
#' This function visualizes the GRN as a hive plot
#' @param dataGRN output GRN function
#' @param namesGenes list TSG and OCG to define axes
#' @param thres threshold of edges to be included
#' @param additionalFilename additionalFilename
#' @import HiveR
#' @importFrom graphics plot.new
#' @importFrom grid gpar
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @export
#' @return no results Hive plot is executed
#' @examples
#' data(knownDriverGenes)
#' data(dataGRN)
#' plotNetworkHive(dataGRN = dataGRN, namesGenes = knownDriverGenes, thres = 0.55)
plotNetworkHive <- function(dataGRN,
                            namesGenes,
                            thres,
                            additionalFilename = NULL) {

  # Check user input

  if (!is(namesGenes, "list")) {
    stop("namesGenes must be a list containing genes")
  }

  if (!is(thres, "numeric")) {
    stop("thres must be numeric defining threshold of edges to be included")
  }

  if (!is(additionalFilename, "NULL") & !is(additionalFilename, "character")) {
    stop("additionalFilename must be either NULL or a character vector
containing part of the filename of plot")
  }

  if (!is(dataGRN, "list")) {
    stop("dataGRN must be a list")
  }

  names.genes.all <- intersect(as.character(unique(c(unlist(namesGenes),
                                                     rownames(dataGRN[[1]])))),
                               colnames(dataGRN[[1]]))

  tmp <- dataGRN[[1]][, (names.genes.all)]
  tmp[tmp < thres] <- 0

  genes.missing <- setdiff(names.genes.all, rownames(dataGRN[[1]]))
  tmp <- rbind(tmp, matrix(0,
                           nrow = length(genes.missing),
                           ncol = length(names.genes.all)))
  rownames(tmp) <- c(rownames(dataGRN[[1]]), genes.missing)

  tmp <- tmp[names.genes.all, names.genes.all]
  diag(tmp) <- 0

  myadj <- adj2HPD(M = tmp, axis.cols = "lightgray")
  myadj$nodes$axis <- as.integer(rep(1, length(names.genes.all)) +
                                   c(2 * as.numeric(names.genes.all %in% namesGenes$TSG) +
                                     as.numeric(names.genes.all %in% namesGenes$OCG)))

  n.axis <- table(myadj$nodes$axis)
  names(n.axis) <- paste0("a.", names(n.axis))

  mycols <- c("darkgrey", "darkgreen", "goldenrod")

  myadj$nodes$color <- mycols[myadj$nodes$axis]
  myadj$nodes$size <- 0.1

  for (i in seq.int(nrow(myadj$nodes))) {
    myadj$nodes$radius[i] <- n.axis[paste0("a.", myadj$nodes$axis[i])]
    n.axis[paste0("a.", myadj$nodes$axis[i])] <- n.axis[paste0("a.", myadj$nodes$axis[i])] - 1
  }

  ind.ocg <- which(rownames(tmp)[myadj$edges$id1] %in% namesGenes$OCG)
  myadj$edges$color[ind.ocg] <- "darkgreen"

  ind.ocg <- which(rownames(tmp)[myadj$edges$id2] %in% namesGenes$OCG)
  myadj$edges$color[ind.ocg] <- "darkgreen"

  ind.tsg <- which(rownames(tmp)[myadj$edges$id1] %in% namesGenes$TSG)
  myadj$edges$color[ind.tsg] <- "goldenrod"

  ind.tsg <- which(rownames(tmp)[myadj$edges$id2] %in% namesGenes$TSG)
  myadj$edges$color[ind.tsg] <- "goldenrod"

  if (!is(additionalFilename, "NULL")) {
    pdf(paste0("networkHive", additionalFilename, ".pdf"))
  }

  HiveR::plotHive(myadj,
                  axLabs = c("remaining TFs", "OCG", "TSG"),
                  bkgnd = "white",
                  anNode.gpar = gpar(fontsize = 10, col = "black", lwd = 0.5))

  if (!is(additionalFilename, "NULL")) {
    dev.off()
  }
}

utils::globalVariables(c("is"))

