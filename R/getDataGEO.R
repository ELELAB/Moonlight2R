#' getDataGEO
#'
#' This function retrieves and prepares GEO data
#' @param GEOobject GEOobject
#' @param platform platform
#' @param TCGAtumor tumor name
#' @importFrom GEOquery getGEO
#' @importMethodsFrom Biobase fvarLabels fvarLabels<-
#' @export
#' @return return GEO gset
#' @examples
#' data(GEO_TCGAtab)
#' dataGEO <-  getDataGEO(GEOobject = "GSE15641", platform = "GPL96")
getDataGEO <- function(GEOobject = "GSE39004",
                       platform = "GPL6244",
                       TCGAtumor = NULL) {

  if (! "GEO_TCGAtab" %in% names(.GlobalEnv)) {
    data(list=c("GEO_TCGAtab"))
  }

  GEO_TCGAtab <- get("GEO_TCGAtab")

  # Check user input

  if (!is(GEOobject, "character")) {
    stop("GEOobject must be a character vector")
  }

  if (!is(platform, "character")) {
    stop("Platform must be a character vector")
  }

  if (!is(TCGAtumor, "NULL") & !is(TCGAtumor, "character")) {
    stop("TCGA tumor must either be NULL or a character vector")
  }

  GEO_TCGAtab <- get("GEO_TCGAtab")

  if (length(TCGAtumor) != 0) {
    GEOobject <- GEO_TCGAtab[GEO_TCGAtab$Cancer ==  TCGAtumor,"Dataset"]
    platform <- GEO_TCGAtab[GEO_TCGAtab$Cancer ==  TCGAtumor,"Platform"]
  }

  gset <- getGEO(GEOobject, GSEMatrix = TRUE, AnnotGPL = TRUE)

  if (length(gset) > 1) {
    idx <- grep(platform, attr(gset, "names"))
  } else {
    idx <- 1
  }

  gset <- gset[[idx]]

  fvarLabels(gset) <- make.names(fvarLabels(gset))

  return(gset)

}

