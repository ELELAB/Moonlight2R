#' EpiMix_getInfiniumAnnotation
#' This function gets Infinium probe annotation from the
#' sesameData library. This function is from the EpiMix package
#' https://bioconductor.org/packages/release/bioc/html/EpiMix.html.
#' Zheng Y, Jun J, Gevaert O (2023). EpiMix: 
#' EpiMix: an integrative tool for the population-level analysis of 
#' DNA methylation. R package version 1.1.2.
#' @param plat A character string representing the methylation platform 
#' which can either be HM27, HM450 or EPIC
#' @param genome A character string representing the genome build version
#' which can either be hg19 or hg38
#' @import ExperimentHub
#' @return Probe annotations
#' @keywords internal
EpiMix_getInfiniumAnnotation <- function(plat = "EPIC", genome = "hg38") {
  hubID <- NULL
  if (tolower(genome) == "hg19" & toupper(plat) == "HM27")
    hubID <- "EH3672"
  if (tolower(genome) == "hg38" & toupper(plat) == "HM27")
    hubID <- "EH3673"
  if (tolower(genome) == "hg19" & toupper(plat) == "HM450")
    hubID <- "EH3674"
  if (tolower(genome) == "hg38" & toupper(plat) == "HM450")
    hubID <- "EH3675"
  if (tolower(genome) == "hg19" & toupper(plat) == "EPIC")
    hubID <- "EH3670"
  if (tolower(genome) == "hg38" & toupper(plat) == "EPIC")
    hubID <- "EH3671"
  ProbeAnnotation <- ExperimentHub::ExperimentHub()[[hubID]]
  return(ProbeAnnotation)
}



