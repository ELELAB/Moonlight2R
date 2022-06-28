#' RunCscape_somatic
#'
#' This function retrive cscape-scores to SNPs
#' @param input Input matching cscape input
#' @param coding_file cscape_table with coding scores
#' @param noncoding_file cscape_table with noncoding scores
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr unite nest unnest
#'
#'
#' @return returns a tibble with a score and remark for each SNP
#' @export
#' @examples
#' 
#' cscape_out <- RunCscape_somatic(input, coding_file, noncoding_file)


RunCscape_somatic <- function(input, coding_file, noncoding_file){
  cscape_in <- input %>%
    unite(col = Ranges, c("Chr", "Start_Position"), sep = ":", remove = FALSE) %>%
    unite(col = Ranges, c("Ranges","Start_Position"), sep = "-", remove = FALSE)
  
  #When using Pmap important to call the list and varibales the same name
  cscape_annot <- cscape_in %>% 
    mutate(file_coding = coding_file,
           file_noncoding = noncoding_file) %>%
    dplyr::select(-Start_Position) %>% 
    mutate(Mydata = pmap(.l= list(Ranges, Reference_Allele, Mutant, file_coding, file_noncoding),
                         .f = tabix_func)) %>% 
    unnest(cols = Mydata)
  
  # Make standard Cscape output format
  cscape_out <- cscape_annot %>% 
    separate(col = Ranges, into = c('Chr', 'Position', NA)) %>% 
    dplyr::rename(Reference = Reference_Allele) %>% 
    dplyr::select(-c(file_coding, file_noncoding)) %>% 
    relocate(Remark, .after = last_col()) %>% 
    mutate(across(.fns = ~guess_parser)) 
  
  return(cscape_out)
}
