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
#' @param max_records An integer containing the maximum
#' number of records to be fetched from PubMed. 
#' @import easyPubMed
#' @import tibble
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @return A tibble containing results of literature search
#' where PubMed was queried for information of input genes. 
#' Each row in the tibble contains a PubMed ID matching the
#' query, doi, title, abstract, year of publication, keywords, 
#' and total number of PubMed publications, resulting in a 
#' total of eight columns.
#' @export
#' @examples 
#' data("dataDMA") 
#' genes_query <- Reduce(c, dataDMA)
#' query <- 
#' "AND cancer AND driver AND '1980/01/01'[Date - Publication] : '2023/01/01'[Date - Publication]"
#' dataGLS <- GLS(genes = genes_query,
#'                query_string = query)

utils::globalVariables(c("pmid", "doi", "title", "abstract", "year", "keywords", "gene"))

GLS <- function(genes, 
                query_string = "AND cancer AND driver",
		max_records = 20) {

  # Check user input
  
  if (!is.character(genes)) {
    stop("Genes must be a character vector containing gene names to search
         in PubMed")
  }
  
  if (!is.character(query_string)) {
    stop("The query string must be a character vector")
  }
  
  if (!is.numeric(max_records)) {
    stop("The maximum number of records to retrieve must be numeric")
  }
  
  # Initialize empty tibble to store results
  pubmed_mining <- tibble()
  
  # For each gene x in input, search PubMed based on specified 
  # query
  literature_search <- map(genes, function(x) {
    
    pubmed_query <- paste(x, query_string) 
    
    # Search and retrieve results from PubMed
    gene_pubmed <- get_pubmed_ids(pubmed_query)
    
    # Retrieve number of publications
    count_pubmed <- gene_pubmed$Count %>% 
      as.numeric()
    
    # If query matches any pubmed records 
    if (count_pubmed > 0) {
      
      # Fetch data of PubMed records searched via above query
      top_results <- fetch_pubmed_data(gene_pubmed, 
                                       retstart = 0, 
                                       retmax = max_records)
  
      # Extract information from PubMed records into a table
      record_info <- table_articles_byAuth(top_results, 
                                           included_authors = "first", 
                                           max_chars = -1, 
                                           getKeywords = TRUE) 
      
      # Select only PubMed id, doi, title, abstract, year, and keywords of
      # PubMed records
      record_info_wrangled <- record_info %>% 
        as_tibble() %>% 
        dplyr::select(c(pmid, doi, title, abstract, year, keywords)) %>% 
        mutate(gene = x,
               pubmed_count = count_pubmed) %>%
        dplyr::relocate(gene, 
                        .after = pmid)
   
      # Bind table to table containing results from previous gene(s)
      pubmed_mining <- pubmed_mining %>% 
        bind_rows(record_info_wrangled)
      
      # If no records of query is found in PubMed
    } else {
      
      # Create tibble of one row of gene that did not have any PubMed results
      no_results_tbl <- tibble(pmid = NA,
                               gene = x,
                               doi = NA,
                               title = NA,
                               abstract = NA,
                               year = NA,
                               keywords = NA,
                               pubmed_count = count_pubmed)
      
      # Bind tibble of gene without PubMed information to table containing
      # results of previous gene(s)
      pubmed_mining <- pubmed_mining %>% 
        bind_rows(no_results_tbl)
      
    }
    
  }) %>% 
    bind_rows()
  
  return(literature_search)
  
}
