#' @title Generate network
#' @description This function carries out the gene regulatory network inference using parmigene
#' @param TFs a vector of genes.
#' @param DEGsmatrix DEGsmatrix output from DEA such as dataDEGs
#' @param DiffGenes if TRUE consider only diff.expr genes in GRN
#' @param normCounts is a matrix of gene expression with genes in rows and samples in columns.
#' @param kNearest the number of nearest neighbors to consider to estimate the mutual information.
#' Must be less than the number of columns of normCounts.
#' @param nGenesPerm nGenesPerm
#' @param nBoot nBoot
#' @param noise_mi noise in knnmi.cross function. Default is 1e-12.
#' @importFrom parmigene knnmi.cross
#' @export
#' @return an adjacent matrix
#' @examples
#' data('DEGsmatrix')
#' data('dataFilt')
#' dataGRN <- GRN(TFs = sample(rownames(DEGsmatrix), 30),
#' DEGsmatrix = DEGsmatrix,
#' DiffGenes = TRUE,
#' normCounts = dataFilt,
#' nGenesPerm = 2,
#' nBoot = 2)
GRN <- function(TFs,
                DEGsmatrix,
                DiffGenes = FALSE,
                normCounts,
                kNearest = 3,
                nGenesPerm = 2000,
                nBoot = 400,
                noise_mi = 1e-12) {

  # Check user input

  if (!is(TFs, "character") | length(TFs) == 0) {
    stop("TFs must be a non-empty character vector containing gene names")
  }

  if (.row_names_info(DEGsmatrix) < 0) {
    stop("Row names were generated automatically. The input DEG table needs to
have the gene names as rownames. Double check that genes are rowname")
  }

  if (!is(DiffGenes, "logical")) {
    stop("DiffGenes must be either TRUE or FALSE")
  }

  if (is(dim(normCounts), "NULL")) {
    stop("The expression data must be non-empty with genes in rows and samples 
in columns")
  }

  if (!is(kNearest, "numeric") | !is(nGenesPerm, "numeric") | !is(nBoot, "numeric") | !is(noise_mi, "numeric")) {
    stop("kNearest, nGenesPerm, nBoot and noise_mi must be numeric values")
  }

  normCountsA <- normCounts
  normCountsB <- normCounts

  if (DiffGenes == TRUE) {
    commonGenes <- intersect(rownames(DEGsmatrix), rownames(normCountsB))
    normCountsB <- normCountsB[commonGenes,]

  } else {
    normCountsB <- normCountsA
  }

  MRcandidates <- intersect(rownames(normCountsA), TFs)

  # Mutual information between TF and genes
  sampleNames <- colnames(normCounts)
  geneNames <- rownames(normCounts)

  miTFGenes <- knnmi.cross(normCountsA[MRcandidates, ],
                           normCountsB,
                           k = kNearest,
                           noise = noise_mi)

  # Threshold with bootstrap
  tfListCancer <- TFs
  tfListCancer <- intersect(tfListCancer, rownames(normCountsA))

  maxmi <- rep(0, length(tfListCancer))

  Cancer_null_distr <- matrix(0, length(tfListCancer), nBoot)
  rownames(Cancer_null_distr) <- tfListCancer

  for (i in seq.int(nBoot)) {

    SampleS <- sample(seq.int(ncol(normCountsA)))
    g <- sample(seq.int(nrow(normCountsA)), nGenesPerm)
    mi <- knnmi.cross(normCountsA[tfListCancer, ],
                      normCountsA[g, SampleS],
                      k = kNearest,
                      noise = noise_mi)

    maxmiCurr <- apply(mi, 1, max)
    Cancer_null_distr[,i] <- maxmiCurr
    index <- maxmi < maxmiCurr
    maxmi[index] <- maxmiCurr[index]

  }

  names(maxmi) <- rownames(Cancer_null_distr)

  return(list(miTFGenes = miTFGenes, maxmi = maxmi))

}
