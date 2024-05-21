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
#' @param simulations string specifying the simulation to use. The simulations possible can be found in the MAVISp database.
#' @importFrom stringr str_c str_split_i
#' @importFrom purrr map
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble_col
#' @importFrom dplyr pull filter
#' @importFrom rlang set_names
#' @importFrom withr with_options
#' @return returns a list of tibbles each containing the MAVISp entry of one protein. Each entry will contain an extra column specifying what data
#' the stability classification is based on.
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
                       mode = 'simple',
                       simulation = 'md'){
    # Look in simple mode index.csv if the protein is in the database
    if (file.exists(str_c(mavispDB,'/dataset_info.csv')) == FALSE){
        stop('MAVISp database file not found at the provided path, or the database_info.csv file is missing.')
    } else if (mode == 'simple'){
        table_location <- str_c(mavispDB,'/simple_mode/dataset_tables/')
    } else if (mode == 'ensemble'){
        if (length(simulation) > 1){
            stop('Only one simulation can be specified for ensemble mode.')
        }
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

    # Double check there are files
    if (length(rawFiles) == 0){
        stop('No MAVISp files matching the criteria were found.')
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
    
    # Combine stability results based on user specification
    if (mode == 'ensemble' & length(simulation) == 1){
        mavispData <- mavispData |>
            map(function(x) rename(x, 'Stability classification, (Rosetta, FoldX)' = matches(paste0('Stability classification, [A-Za-z0-9, ]*\\(Rosetta, FoldX\\)( \\[',simulation,'\\])?')),
                                    'Stability classification, (RaSP, FoldX)' = matches(paste0('Stability classification, [A-Za-z0-9, ]*\\(RaSP, FoldX\\)( \\[',simulation,'\\])?'))) |>
                            mutate(stab_class_ros_source = paste0('ensemble_mode_',simulation),
                                   stab_class_rasp_source = paste0('ensemble_mode_',simulation)))
    } else if (mode == 'simple'){
        mavispData <- mavispData |>
            map(function(x) mutate(stab_class_data_type = 'simple_mode'))
    }


    return(mavispData)
}
