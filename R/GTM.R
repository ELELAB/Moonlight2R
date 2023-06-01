#' GTM
#' This function carries out gene text mining.
#' @param genes A character string containing the genes
#' to perform text mining using PubMed database
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
#' @import stringr
#' @import tidyr
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @return A tibble containing results of text mining where
#' PubMed was queried for information of input genes. Each
#' row in the tibble contains a sentence from an abstract
#' that contains either of the words supplied in the query.
#' Each row also contains corresponding PubMed ID, gene, 
#' total number of PubMed publications, doi, title, full
#' abstract and year of publications, resulting in a total
#' of eight columns.
#' @export
#' @examples 
#' data("dataDMA") 
#' genes_query <- Reduce(c, dataDMA)
#' dataGTM <- GTM(genes = genes_query) 
GTM <- function(genes, 
                query_string = "AND cancer AND driver",
		max_records = 20) {
  
  # Initialize empty tibble to store results
  pubmed_mining <- tibble()
  
  # For each gene x in input, search PubMed based on specified 
  # query
  text_mining <- map(genes, function(x) {
    
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
                                           getKeywords = FALSE) 
      
      # Select only PubMed id, doi, title, abstract and year of
      # PubMed records
      record_info_wrangled <- record_info %>% 
        as_tibble() %>% 
        dplyr::select(c(pmid, doi, title, abstract, year)) %>% 
        mutate(gene = x) %>% 
        dplyr::rename(abstract_full = abstract)
      
      # Split pubmed_query into words and remove presence of tags
      pubmed_query_split <- str_remove_all(pubmed_query, "\\[.*?\\]") %>% 
        str_replace_all(., "\\b(?: AND | OR | NOT )\\b", "|")
      
      # Select all abstracts from tibble containing all information 
      all_abstracts <- record_info_wrangled %>% 
        dplyr::select(abstract_full)
      
      # For each abstract
      abstract_query <- map(all_abstracts$abstract_full, function(abstract) {
        
        # Split abstract into sentences 
        abstract_split <- str_split(string = abstract, 
                                    pattern = "\\. ") %>% 
          unlist %>% 
          as.list()
        
        # Extract those sentences in abstract that contain either of the words
        # supplied in the query
        abstract_split <- abstract_split[str_detect(string = abstract_split, 
                                                    pattern = pubmed_query_split)] %>% 
          as.list()
        
      })
      names(abstract_query) <- record_info_wrangled$pmid
      
      # Combine abstract_query with record_info_wrangled 
      merged_data <- record_info_wrangled %>%
        mutate(abstract_sentences = map(pmid, ~ abstract_query[[.]])) %>% 
        unnest(cols = abstract_sentences) %>% 
        mutate(abstract_sentences = unlist(abstract_sentences)) %>% 
        dplyr::relocate(gene,
                        abstract_sentences, 
                        .after = pmid)
      
      # Add number of total publications found in PubMed to table
      merged_data <- merged_data %>% 
        mutate(pubmed_count = count_pubmed) %>% 
        dplyr::relocate(pubmed_count,
                        .before = doi)
      
      # Bind table to table containing results from previous gene(s)
      pubmed_mining <- pubmed_mining %>% 
        bind_rows(merged_data)
      
      # If no records of query is found in PubMed
    } else {
      
      # Create tibble of one row of gene that did not have any PubMed results
      no_results_tbl <- tibble(pmid = NA,
                               gene = x,
                               abstract_sentences = NA,
                               pubmed_count = count_pubmed,
                               doi = NA,
                               title = NA,
                               abstract_full = NA,
                               year = NA)
      
      # Bind tibble of gene without PubMed information to table containing
      # results of previous gene(s)
      pubmed_mining <- pubmed_mining %>% 
        bind_rows(no_results_tbl)
      
    }
    
  }) %>% 
    bind_rows()
  
  return(text_mining)
  
}
