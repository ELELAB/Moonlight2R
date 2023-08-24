#' MAFtoCscape
#'
#' This function extracts columns from a MAF tibble to fit CScape input format
#' @param MAF tibble of MAF 
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @importFrom tidyr separate
#' @return tibble of cscape-somatic input
#' @export
#' @examples
#' print(data(dataMAF))
#' MAFtoCscape(dataMAF)

MAFtoCscape <- function(MAF){
  cscape <- MAF %>% 
    filter(Variant_Type == 'SNP') %>% 
    dplyr::select(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2) %>% 
    mutate(Mutant = case_when(Reference_Allele == Tumor_Seq_Allele1 ~ Tumor_Seq_Allele2,
                              Reference_Allele == Tumor_Seq_Allele2 ~ Tumor_Seq_Allele1)) %>% 
    separate(Chromosome, into = c(NA, "Chr"), sep = 3) %>% 
    dplyr::select(Chr, Start_Position, Reference_Allele, Mutant)
  return(cscape)
}

utils::globalVariables(c("Variant_Type", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Chr", "Mutant"))
