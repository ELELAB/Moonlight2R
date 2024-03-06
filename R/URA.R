#' @title URA Upstream Regulator Analysis
#' @description
#' This function carries out the upstream regulator analysis
#' @param dataGRN output GNR function
#' @param DEGsmatrix output DPA function
#' @param BPname biological processes
#' @param nCores number of cores to use
#' @importFrom stats fisher.test
#' @import doParallel
#' @import foreach
#' @export
#' @return an adjacent matrix
#' @examples
#' data(DEGsmatrix)
#' dataDEGs <- DEGsmatrix
#' data(dataGRN)
#' data(DiseaseList)
#' data(EAGenes)
#' dataURA <- URA(dataGRN = dataGRN,
#' DEGsmatrix = dataDEGs,
#' BPname = c("apoptosis",
#' "proliferation of cells"))
URA <- function(dataGRN,
                DEGsmatrix,
                BPname,
                nCores = 1) {

  # List of variable names
  variables_to_check <- c("DiseaseList", "EAGenes")

  # Check and load variables if they do not exist
  for (variable_name in variables_to_check) {
    if (! variable_name %in% names(.GlobalEnv)) {
      data(list=c(variable_name))
    }
  }

  DiseaseList <- get("DiseaseList")

  # Check user input
  if (.row_names_info(DEGsmatrix) < 0) {
    stop("Row names were generated automatically. The input DEG table needs to
have the gene names as rownames. Double check that genes are rownames.")
  }

  if (!is(BPname, "NULL") && all(BPname %in% names(DiseaseList)) == FALSE) {
    stop("BPname should be NULL or a character vector containing one or more
BP(s) among possible BPs stored in the DiseaseList object.")
  }

  doParallel::registerDoParallel(cores = nCores)

  if (is(BPname, "NULL")) {
    BPname <- names(DiseaseList)
  }

  tRlist <- rownames(dataGRN$miTFGenes)
  pb <- txtProgressBar(min = 0, max = length(tRlist), style = 3)
  j <- NULL

  TableDiseases <- foreach(j = seq.int(tRlist), .combine = "rbind", .packages="foreach") %dopar% {

    currentTF <- as.character(tRlist[j])
    currentTF_regulon <- names(which(dataGRN$miTFGenes[currentTF, ] > as.numeric(dataGRN$maxmi[currentTF])))
    currentTF_regulon <- as.matrix(currentTF_regulon)
    DEGsregulon <- intersect(rownames(DEGsmatrix), currentTF_regulon)

    if (length(DEGsregulon) > 2) {
      tabFEA <- FEA(BPname = BPname, DEGsmatrix = DEGsmatrix[DEGsregulon, ])
      return(tabFEA$Moonlight.Z.score)
    } else {
      return(rep(0, length(BPname)))
    }
  }

  dimnames(TableDiseases) <- list(tRlist, BPname)

  stopImplicitCluster()

  close(pb)

  return(TableDiseases)

}
