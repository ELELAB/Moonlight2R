#' PRAtoTibble
#'
#' This function changes the PRA output to tibble format
#' @param dataPRA RDA object (list of two) from PRA
#' @import dplyr  
#' @importFrom magrittr "%>%"
#'
#' @return tibble with drivers
#' @export
#' @examples
#' 
#' PRAtoTibble(dataPRA)



PRAtoTibble <- function(dataPRA){

  # Wrangle Data:
  TSG <- as_tibble(dataPRA$TSG, rownames = NA) %>% 
    rownames_to_column(var = "Hugo_Symbol") %>%  
    mutate(Hugo_Symbol = str_trim(Hugo_Symbol, side = 'both'),
           Moonlight_Oncogenic_Mediator = "TSG") %>%
    dplyr::rename(Moonlight_gene_z_score = value)
  
  OCG <- as_tibble(dataPRA$OCG, rownames = NA) %>% 
    rownames_to_column(var = "Hugo_Symbol") %>% 
    mutate(Hugo_Symbol = str_trim(Hugo_Symbol, side = 'both'),
           Moonlight_Oncogenic_Mediator = "OCG") %>% 
    dplyr::rename(Moonlight_gene_z_score = value)
  
  drivers <- full_join(TSG, OCG, by = c("Hugo_Symbol", 
                                        "Moonlight_gene_z_score", 
                                        "Moonlight_Oncogenic_Mediator") ) 
  
  return(drivers)
}
