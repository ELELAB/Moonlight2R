#' PRAtoTibble
#'
#' This function changes the PRA output to tibble format
#' @param dataPRA RDA object (list of two) from PRA
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_trim
#'
#' @return tibble with drivers
#' @export
#' @examples
#' data('dataPRA')
#' PRAtoTibble(dataPRA)
PRAtoTibble <- function(dataPRA) {

  # Check user input

  if (all(names(dataPRA) %in% c("TSG", "OCG")) == FALSE) {
    stop("The two list elements in PRA data must be named TSG and OCG")
  }

  # Wrangle data

  if (!is(dataPRA$TSG, "NULL")) {
    TSG <- as_tibble(dataPRA$TSG, rownames = NA) %>%
      rownames_to_column(var = "Hugo_Symbol") %>%
      mutate(Hugo_Symbol = str_trim(Hugo_Symbol, side = "both"),
             Moonlight_Oncogenic_Mediator = "TSG") %>%
      dplyr::rename(Moonlight_gene_z_score = value)
  } else {
    TSG <- NULL
  }

  if (!is(dataPRA$OCG, "NULL")) {
    OCG <- as_tibble(dataPRA$OCG, rownames = NA) %>%
      rownames_to_column(var = "Hugo_Symbol") %>%
      mutate(Hugo_Symbol = str_trim(Hugo_Symbol, side = "both"),
             Moonlight_Oncogenic_Mediator = "OCG") %>%
      dplyr::rename(Moonlight_gene_z_score = value)
  } else {
    OCG <- NULL
  }

  if (!is(TSG, "NULL") & !is(OCG, "NULL")) {
    drivers <- full_join(TSG, OCG, by = c("Hugo_Symbol",
                                          "Moonlight_gene_z_score",
                                          "Moonlight_Oncogenic_Mediator"))
  } else if (is(TSG, "NULL")) {
    drivers <- OCG
  } else if (is(OCG, "NULL")) {
    drivers <- TSG
  }

  return(drivers)

}

utils::globalVariables(c("Hugo_Symbol", "value"))
