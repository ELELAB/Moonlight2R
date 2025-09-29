#' GLS
#' This function carries out gene literature search.
#' @param genes A character string containing the genes
#' to search in PubMed database
#' @param query_string A character string containing words
#' in query to follow the gene of interest. Default is
#' "AND cancer AND driver" resulting in a final query of
#' "Gene AND cancer AND driver". Standard PubMed syntax
#' can be used in the query. For example Boolean operators
#' AND, OR, NOT can be applied and tags such as [AU],
#' [TITLE/ABSTRACT], [Affiliation] can be used.
#' @import easyPubMed
#' @import tibble
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @return A tibble containing results of literature search
#' where PubMed was queried for information of input genes.
#' Each row in the tibble contains a PubMed ID matching the
#' query, doi, title, abstract, year of publication, mesh_terms,
#' and total number of PubMed publications, resulting in a
#' total of eight columns.
#' @export
#' @examples
#' genes_query <- "BRCA1"
#' dataGLS <- GLS(genes = genes_query,
#'                query_string = "AND cancer AND driver")
GLS <- function(genes,
                query_string = "AND cancer AND driver") {

  # Check user input

  if (!is(genes, "character")) {
    stop("Genes must be a character vector containing gene names to search
in PubMed")
  }

  if (!is(query_string, "character")) {
    stop("The query string must be a character vector")
  }

  # Initialize empty tibble to store results
  pubmed_mining <- tibble()

  # For each gene x in input, search PubMed based on specified
  # query

  search_pubmed <- function(x) {

    pubmed_query <- paste(x, query_string)

    # Search and retrieve results from PubMed
    gene_pubmed <- epm_query(pubmed_query)

    # Retrieve number of publications
    count_pubmed <- gene_pubmed@meta$exp_count %>%
      as.numeric()

    # If query matches any pubmed records
    if (count_pubmed > 0) {

      # Fetch data of PubMed records searched via above query
      top_results <- epm_fetch(gene_pubmed)

      # Extract information from PubMed records into a an easyPubMed object
      record_info <- epm_parse(top_results, max_authors = 1)
      # Extract information from PubMed records into a table
      record_info_data <- get_epm_data(record_info)

      # Select only PubMed id, doi, title, abstract, year, and mesh terms of
      # PubMed records
      record_info_wrangled <- record_info_data %>%
        as_tibble() %>%
        dplyr::select(c(pmid, doi, title, abstract, year, mesh_terms)) %>%
        mutate(gene = x, pubmed_count = count_pubmed) %>%
        dplyr::relocate(gene, .after = pmid)

      # Bind table to table containing results from previous gene(s)
      pubmed_mining <- pubmed_mining %>%
        bind_rows(record_info_wrangled)

      # If no records of query is found in PubMed
    } else {

      # Create tibble of one row of gene that did not have any PubMed results
      no_results_tbl <- tibble(pmid = NA, gene = x, doi = NA,
                               title = NA, abstract = NA, year = NA,
                               mesh_terms = NA, pubmed_count = count_pubmed)

      # Bind tibble of gene without PubMed information to table containing
      # results of previous gene(s)
      pubmed_mining <- pubmed_mining %>%
        bind_rows(no_results_tbl)

    }

  }

  literature_search <- map(genes, search_pubmed) %>%
    bind_rows()

  return(literature_search)

}

utils::globalVariables(c("pmid", "doi", "title", "abstract", "year",
                         "mesh_terms", "gene"))
