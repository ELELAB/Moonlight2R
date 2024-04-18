#' loadMAVISp
#' 
#' This function loads the MAVISp database from the directory specified by the user. 
#' 
#' @param mavispDB path to the MAVISp database as a string
#' @param proteins_of_interest vector containing specific proteins of interest in HUGO format
#' @param mode string determining whether to use simple or ensemble mode of mavisp. Default is simple mode.
#' Takes values: 
#' \itemize{
#' \item simple
#' \item ensemble
#'} 
#' @importFrom stringr str_c str_split_i
#' @importFrom purrr map
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble_col
#' @importFrom dplyr pull filter
#' @importFrom rlang set_names
#' @importFrom withr with_options
#' @return returns a list of tibbles each containing the MAVISp entry of one protein
#' @export
#' @examples
#' 
#' mavisp_data <- loadMAVISp(mavispDB = "/data/raw_data/computational_data/mavisp_database/biorxiv_v4_17112023",
#'           proteins_of_interest = c('NQO1','TP53'),
#'           mode = 'ensemble')
#' 
#' mavisp_data <- loadMAVISp(mavispDB = "/data/raw_data/computational_data/mavisp_database/biorxiv_v4_17112023") 

loadMAVISp <- function(mavispDB = NULL,
                       proteins_of_interest = NULL,
                       mode = 'simple'){
    # Look in simple mode index.csv if the protein is in the database
    if (file.exists(str_c(mavispDB,'/dataset_info.csv')) == FALSE){
        stop('MAVISp database file not found at the provided path')
    } else if (mode == 'simple'){
        table_location <- str_c(mavispDB,'/simple_mode/dataset_tables/')
    } else if (mode == 'ensemble'){
        table_location <- str_c(mavispDB,'/ensemble_mode/dataset_tables/')
    } else {
        stop('Mode not specified correctly. Takes values "simple" or "ensemble"')
    }

    # Load data for proteins of interest or all proteins
    if (is.null(proteins_of_interest)){
        rawFiles <- list.files(table_location,
                               full.names = TRUE)

    } else {
        proteins_of_interest <- str_c(proteins_of_interest,'-')
        rawFiles <- list.files(table_location,
                               full.names = TRUE) |>
                               as_tibble_col(column_name = 'filepath') |>
                    filter(grepl(paste(proteins_of_interest, 
                                        collapse = '|'),
                                        filepath)) |>
                    pull(filepath)      
    }

    mavispData <- rawFiles |>
        set_names(str_split_i(basename(rawFiles), '-', 1)) |>
        # Supress non-fatal warnings
        map(function(x) withr::with_options(
                        list(rlib_name_repair_verbosity = 'quiet'),
                        suppressWarnings(
                        classes = 'vroom_parse_issue',
                        read_csv(file = x,
                                progress = FALSE,
                                show_col_types = FALSE))))

    return(mavispData)
}
