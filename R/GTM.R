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
#' @import easyPubMed
#' @import tibble
#' @importFrom magrittr "%>%"
#' @importFrom purrr map
#' @importFrom dplyr rename
#' @return A tibble containing results of text mining where
#' PubMed was queried for information of input genes. Each
#' row in the tibble contains a PubMed ID for a given gene.
#' If more than 20 results from the specified query is obtained, 
#' only the first 20 are included. The tibble contains four
#' columns: Hugo_Symbol which contains the given gene, 
#' pubmed_query which contains the query that was searched,
#' pubmed_count which is the number of total publications 
#' that the query resulted in and pubmed_ids which contains the
#' PubMed IDs. 
#' @export
#' @examples 
#' data("dataDMA") 
#' genes_query <- Reduce(c, dataDMA)
#' dataGTM <- GTM(genes = genes_query) 
GTM <- function(genes, 
                query_string = "AND cancer AND driver") {
  
  # Initialize empty tibble to store results
  pubmed_mining <- tibble(x = character(), 
                          pubmed_query = character(),
                          pubmed_count = numeric(),
                          pubmed_ids = character())
  
  # For each gene x in input, search PubMed based on specified 
  # query
  text_mining <- map(genes, function(x) {
    
    pubmed_query <- paste(x, query_string) 
    
    # Search and retrieve results from PubMed
    gene_pubmed <- get_pubmed_ids(pubmed_query)
    
    # Obtain number of publications
    pubmed_count <- gene_pubmed$Count
    pubmed_ids <- unlist(gene_pubmed$IdList)
    
    # Combine results of gene x to tibble
    tmp_tibble <- tibble(x, pubmed_query, pubmed_count, pubmed_ids) %>% 
      mutate(pubmed_count = as.numeric(pubmed_count))
    pubmed_mining <- pubmed_mining %>% 
      bind_rows(tmp_tibble)
    
  }) %>% 
    bind_rows() %>% 
    dplyr::rename(Hugo_Symbol = x)
  
  return(text_mining)
  
}
