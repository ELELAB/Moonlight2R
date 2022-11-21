#' LiftMAF
#'
#' This function lifts a MAF file to a different genomic build.
#' @param Infile A tibble of MAF. 
#' @param Current_Build A charcter string, either \code{GRCh38} or \code{GRCh37}.
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rtracklayer import.chain liftOver
#' @return MAF tibble with positions lifted to another build 
#' @export
#' @examples
#' 
#' LiftMAF(Infile, Current_Build = 'GRCh38')


LiftMAF <- function(Infile, Current_Build){
  #The input file is assumed to be maf_tibble, this file is then lifted to 
  # either 38 or 37 and return as a tibble
  flag <- FALSE
  if(Current_Build == 'GRCh38'){
    chainBuild <- "hg38ToHg19.over.chain"
    #Make chain from 38 to 19
    path = system.file(package="liftOver", "extdata", chainBuild)
    flag <- TRUE
  } else if(Current_Build == "GRCh37"){
    chainBuild <- "hg19ToHg38.over.chain"
    path <- system.file("data", "hg19ToHg38.over.chain", 
                        package = "Moonlight2R", mustWork = TRUE)
    flag <- TRUE
  } else {
    print("Error: Build must be either GRCh38 or GRCh37")
  }

  if(flag == TRUE){
    #Change to Grange format
    infile_GRange <- makeGRangesFromDataFrame(Infile, 
                                              start.field = "Start_Position", 
                                              end.field = "End_Position",
                                              seqnames.field = "Chromosome",
                                              keep.extra.columns = TRUE)
    #Import chain
    chain <- import.chain(path)
    
    #Do liftover
    infile_GRange_lifted <- liftOver(x = infile_GRange, chain = chain)
    
    #recreate maf column names and tibble format
    outfile_tibble_lifted <- as_tibble(infile_GRange_lifted) %>% 
      dplyr::rename(Start_Position = start,
                    End_Position = end, 
                    Chromosome = seqnames,
                    Strand = strand) %>% 
      dplyr::select(-c(group, group_name, width))
    
    return(outfile_tibble_lifted)
  }
}
