

# Input
# Path to mavisp database
# List of proteins to look for (HUGO)
suppressPackageStartupMessages(library('tidyverse'))


# Return tibble with data

loadMAVISp <- function(proteins_of_interest = NULL,
                       mavispDB =   NULL){
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
        print(rawFiles)
        mavispData <- rawFiles |>
                        set_names(str_split_i(basename(rawFiles), '-', 1)) |>
                        map(function(x) read_csv(file = x,
                                                progress = FALSE,
                                                show_col_types = FALSE))

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
                        map(function(x) read_csv(file = x,
                                                progress = FALSE,
                                                show_col_types = FALSE))
    }

}

loadMAVISp(proteins_of_interest = c('NQO1','TP53'))


# loadMAVISp(mavispDB = "/data/raw_data/computational_data/mavisp_database/biorxiv_v4_17112023")
