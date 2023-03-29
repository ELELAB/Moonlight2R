#' tabix_func
#'
#' This function retrives the individial score for a SNP
#' @param Ranges The position
#' @param Reference_Allele The reference nucleotide
#' @param Mutant The mutant nucleotide
#' @param file_coding cscape_table with coding scores
#' @param file_noncoding cscape_table with noncoding scores
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite nest unnest
#' @importFrom stringr str_split str_replace_all
#' @importFrom seqminer tabix.read.table
#' @importFrom tibble tibble
#'
#' @return returns the score
#' @keywords internal
#' @examples
#' \dontrun{ 
#' data <- tabix_func(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding)
#'}


tabix_func <- function(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding){
  
  flag <- FALSE
  bases <- c("A", "T","G", "C")
  chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13",
                   "14", "15", "16", "17", "18", "19", "20", "21", "22") #X and Y not avail
  chromosome <- str_split(Ranges, pattern = ':',simplify = TRUE) %>% .[[1]]
  remark <- 'std'
  type <- ''
  
  reference <- as.character(Reference_Allele)
  mutant <- as.character(Mutant)
  # Check Reference, Mutant and Chromosome are correct format and possible 
  if( !(reference %in% bases) | !(mutant %in% bases)){
    remark <- "Error: Reference and mutant nucleotides must both be ACGT"
    score <- NA
  } else if( !(chromosome %in% chromosomes)){
    remark <- "Error: Unexpected chromosome (X and Y not available)"
    score <- NA
    # If all okay, continue:
  } else{
    #Look for score in coding region
    x <- as_tibble(tabix.read.table(tabixFile = file_coding,
                                    tabixRange = Ranges)) %>%
      mutate(across(where(is.logical),as.character)) %>%
      mutate(across(.cols = everything(),
                    .fns =~ str_replace_all(string =., pattern = "TRUE", replacement = "T")))
    
    # Has the data been found
    if(dim(x)[1] != 0){
      flag <- TRUE
      type <- 'Coding'
    } else {
      # If no annotation is found in this position try the noncoding file
      x <- as_tibble(tabix.read.table(tabixFile = file_noncoding,
                                      tabixRange = Ranges)) %>%
        mutate(across(where(is.logical),as.character)) %>%
        mutate(across(.cols = everything(),
                      .fns =~ str_replace_all(string =., pattern = "TRUE", replacement = "T")))
      
      
      if(dim(x)[1] != 0){
        flag <- TRUE
        type <- 'Noncoding'
      } else{
        # if there still is no data found at this position
        flag <- FALSE
        remark <- 'Error: Unexcepted position'
        score <- NA
        print('here 3')
      }
    }
    
    # If the position has been found
    if(flag == TRUE){
      if (reference == mutant){
        remark <- 'Error: Reference and mutant must be different'
        score <- NA
      }
      else if (reference == x$V3[1]) {
        score <- x %>% filter(V4 == mutant) %>% dplyr::select(V5) %>% pull()
        remark <- confidence(score, type)
      } else {
        remark <- 'Error: Unexpected reference'
        score <- NA
      }
    }}
  
  if (type == 'Coding'){
    data <- tibble(Coding_score = score, Remark = remark)
  } else if (type == 'Noncoding'){
    data <- tibble(Noncoding_score = score, Remark = remark)
  } else {
    data <- tibble(Remark = remark) #tibble(Coding_score = '', Noncoding_score = '', Remark = remark)
  }
  return(data)
}
