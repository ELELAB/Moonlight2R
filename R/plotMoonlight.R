#' plotMoonlight 
#' 
#' This function creates a heatmap of Moonlight gene z-scores for selected genes.
#'
#' @param DEG_Mutations_Annotations A tibble, output file from DMA. 
#' @param Oncogenic_mediators_mutation_summary A tibble, output file from DMA.
#' @param dataURA Output URA function.
#' @param gene_type A character string either \code{"mediators"} or \code{"drivers"}. 
#' \itemize{
#' \item If \code{NULL} defaults to \code{"drivers"}. 
#' \item \code{"mediators"} will show the oncogenic mediators with the highst number of mutations regardless of driver/passenger classification.
#' \item \code{"drivers"} will show the driver genes with the highest number of driver mutations.
#' }
#' @param n Numeric. The top number of genes to be plotted. If \code{NULL} defaults to 50.
#' @param genelist A vector of strings containing Hugo Symbols of genes. 
#' Overwrites \code{gene_type} argument.
#' @param additionalFilename A character string. Adds prefix or filepath to the filename of the pdf.
#'
#' @import dplyr  
#' @importFrom magrittr "%>%"
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom tidyr drop_na pivot_longer replace_na
#' @importFrom tibble rownames_to_column
#'
#' @return No return value. Moonlight scores are plotted for selected genes.
#' @export
#'
#' @examples
#' plotMoonlight(DEG_Mutations_Annotations, 
#'               Oncogenic_mediators_mutation_summary, 
#'               dataURA, gene_type = "drivers", n = 50,
#'               additionalFilename = "path/myplot_")
#' plotMoonlight(DEG_Mutations_Annotiontions, 
#'               Oncogenic_Mediators_mutation_summary, 
#'               dataURA, 
#'               genelist = c("BRCA1", "BRCA2", "GATA3", "RB1"))

plotMoonlight <- function(DEG_Mutations_Annotations, 
                          Oncogenic_mediators_mutation_summary,
                          dataURA,
                          gene_type = "drivers",
                          n = 50, 
                          genelist = c(),
                          additionalFilename = ""){
  
  # The differentially expressed genes, that are annotated as TSG/OCG
  DEGs <- DEG_Mutations_Annotations %>% 
    select(Hugo_Symbol, Moonlight_gene_z_score, logFC) %>% 
    unique() %>% 
    drop_na(Moonlight_gene_z_score)
  
  # restructure URA to tibble
  ura <- as_tibble(dataURA, rownames = NA) %>% 
    rownames_to_column(var = "Genes")
  
  ura_wrangled <- ura %>%  
    pivot_longer(cols = !c('Genes'), 
                 names_to = 'Biological_Process',
                 values_to = 'Moonlight_score') %>% 
    right_join(Oncogenic_mediators_mutation_summary, 
               by = c("Genes" = "Hugo_Symbol")) %>% 
    right_join(DEGs, by = c("Genes" = "Hugo_Symbol")) %>%
    replace_na(list(CScape_Driver = 0, 
                    CScape_Passenger = 0, 
                    CScape_Unclassified = 0)) #%>% 
  
  # Type of plot:
  if (length(genelist) > 0 ){
    ura_wrangled <- ura_wrangled %>% 
      filter(Genes %in% genelist)
    n <- Summary_wrangled %>% summarise() %>% count() %>% pull()
    if (n <= 1){
      stop("The genelist must contain at least one OCG and one TSG.")
    }
    
    } else if (gene_type == "mediators"){
      ura_wrangled <- ura_wrangled %>% 
        slice_max(Total_Mutations, n = n, with_ties = FALSE)
        
    } else{
    ura_wrangled <- ura_wrangled %>% 
      slice_max(CScape_Driver, n = n, with_ties = FALSE)
      
    }
  
  # Variable for color scaling in legend
  max_driver <- ura_wrangled %>% arrange(desc(CScape_Driver)) %>% 
    select(CScape_Driver) %>%  head(1) %>% pull
  
  # Plot Heatmap
  bp_heatmap <- heatmap(ura_wrangled,
                        .row = Biological_Process,
                        .column = Genes,
                        .value = Moonlight_score,
                        scale = "none",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "complete",
                        cluster_rows = FALSE) %>%
    add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>%
    add_tile(logFC, palette = c("chartreuse4","firebrick3")) %>%
    add_tile(CScape_Driver, palette = colorRamp2(c(0,max_driver), c("white", "dodgerblue3"))) %>%
    add_bar(Total_Mutations) 
  
  save_pdf(bp_heatmap, height = 15, width = 35, units = "cm",
           filename = paste(additionalFilename,"moonlight_heatmap.pdf", sep =""))
}
