#' loadMAVISp
#' 
#' This function loads the MAVISp database from the directory specified by the user. 
#' 
#' @param mavispDB path to the MAVISp database
#' @param proteins_of_interest vector containing specific proteins of interest in HUGO format
#' @importFrom stringr str_c str_split_i
#' @importFrom purrr map
#' @importFrom readr read_csv
#' @importFrom tibble as_tibble_col
#' @importFrom dplyr pull filter
#' @importFrom rlang set_names
#' @return returns a list of tibbles each containing the MAVISp entry of one protein
#' @export
#' @examples
#' 
#' mavisp_data <- loadMAVISp(mavispDB = "/data/raw_data/computational_data/mavisp_database/biorxiv_v4_17112023",
#'           proteins_of_interest = c('NQO1','TP53'))
#' 
#' mavisp_data <- loadMAVISp(mavispDB = "/data/raw_data/computational_data/mavisp_database/biorxiv_v4_17112023") 

loadMAVISp <- function(mavispDB = NULL,
                       proteins_of_interest = NULL){
    # Look in simple mode index.csv if the protein is in the database
    if (file.exists(str_c(mavispDB,'/simple_mode/index.csv')) == FALSE){
        stop("MAVISp database file not found at the provided path or does not contain 'simple mode'")
    } else {
        table_location <- str_c(mavispDB,'/simple_mode/dataset_tables/')
    }

    # Load data for proteins of interest or all proteins
    if (is.null(proteins_of_interest)){
        rawFiles <- list.files(table_location,
                               full.names = TRUE)

        mavispData <- rawFiles |>
                        set_names(str_split_i(basename(rawFiles), '-', 1)) |>
                        # Supress non-fatal warnings
                        map(function(x) suppressWarnings(
                                        classes = 'vroom_parse_issue',
                                        read_csv(file = x,
                                                progress = FALSE,
                                                show_col_types = FALSE)))
        
        return(mavispData)

    } else {
        proteins_of_interest <- str_c(proteins_of_interest,'-')
        rawFiles <- list.files(table_location,
                               full.names = TRUE) |>
                               as_tibble_col(column_name = 'filepath')
        
        filtered_files <- rawFiles |>
                            filter(grepl(paste(proteins_of_interest, 
                                               collapse = '|'),
                                         filepath)) |>
                            pull(filepath)

        mavispData <- filtered_files |>
                        set_names(str_split_i(basename(filtered_files), '-', 1)) |>
                        # Supress non-fatal warnings
                        map(function(x) suppressWarnings(
                                        classes = 'vroom_parse_issue',
                                        read_csv(file = x,
                                                progress = FALSE,
                                                show_col_types = FALSE)))
        
        return(mavispData)
    }

}


