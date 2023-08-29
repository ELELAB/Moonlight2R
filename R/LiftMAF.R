#' LiftMAF
#'
#' This function lifts a MAF file to a different genomic build.
#' @param Infile A tibble of MAF.
#' @param Current_Build A charcter string, either \code{GRCh38} or \code{GRCh37}
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer import.chain
#' @importFrom AnnotationHub AnnotationHub
#' @return MAF tibble with positions lifted to another build
#' @export
#' @examples
#' data(dataMAF)
#' LiftMAF(dataMAF, Current_Build = 'GRCh38')
LiftMAF <- function(Infile,
                    Current_Build) {

  # The input file is assumed to be maf_tibble, this file is then lifted to
  # either 38 or 37 and return as a tibble

  ah <- AnnotationHub()

  if (Current_Build == "GRCh38") {
    chain <- ah[['AH14108']]   #UCSC  hg38ToHg19.over.chain.gz
  } else if (Current_Build == "GRCh37") {
    chain <- ah[['AH14150']]   #UCSC  hg19ToHg38.over.chain.gz
  } else {
    print("Error: Build must be either GRCh38 or GRCh37")
    return
  }

  # Change to Grange format
  infile_GRange <- makeGRangesFromDataFrame(Infile,
                                            start.field = "Start_Position",
                                            end.field = "End_Position",
                                            seqnames.field = "Chromosome",
                                            keep.extra.columns = TRUE)

  # Do liftover
  infile_GRange_lifted <- liftOver(x = infile_GRange, chain = chain)

  # Recreate maf column names and tibble format
  outfile_tibble_lifted <- as_tibble(infile_GRange_lifted) %>%
    dplyr::rename(Start_Position = start,
                  End_Position = end,
                  Chromosome = seqnames,
                  Strand = strand) %>%
    dplyr::select(-c(group, group_name, width))

  return(outfile_tibble_lifted)
}

utils::globalVariables(c("liftOver", "start", "end", "seqnames", "strand",
                         "group", "group_name", "width"))
