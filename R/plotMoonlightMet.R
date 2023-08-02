#' plotMoonlightMet
#' This function visualizes the effect of genes on biological processes
#' and total number of hypo/hyper/dual methylated CpG sites in genes.
#' @param DEG_Methylation_Annotations A tibble which is outputted
#' from the EDA function.
#' @param Oncogenic_mediators_methylation_summary A tibble which
#' is outputted from the EDA function.
#' @param dataURA Output of the URA function: a table containing
#' the effect of oncogenic mediators on biological processes. This
#' effect is quantified through Moonlight Gene Z-scores.
#' @param genes A character string containing the genes to be 
#' visualized.
#' @param additionalFilename A character string that can be used
#' to add a prefix or filepath to the filename of the pdf 
#' visualizing the heatmap. Default is an empty string. 
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' @importFrom magrittr "%>%"
#' @return No value is returned. Visualization in form of a heatmap
#' is saved.
#' @export
#' @examples 
#' data("DEG_Methylation_Annotations")
#' data("Oncogenic_mediators_methylation_summary")
#' data("dataURA")
#' genes <- c("ABCG2", "ABHD6", "ACADL", "ACAN", "ACE2", "ACSS2",
#'            "ACTG2", "ADAM19", "ADCK5", "ADHFE1", "ADRB1", "AFAP1L1")
#' plotMoonlightMet(DEG_Methylation_Annotations = DEG_Methylation_Annotations, 
#'                  Oncogenic_mediators_methylation_summary = Oncogenic_mediators_methylation_summary, 
#'                  dataURA = dataURA, 
#'                  genes = genes,
#'                  additionalFilename = "./EDAresults/")
plotMoonlightMet <- function(DEG_Methylation_Annotations,
                             Oncogenic_mediators_methylation_summary,
                             dataURA,
                             genes,
                             additionalFilename = "") { 
  
  # Get logFC of DEGs
  DEGs_logFC <- DEG_Methylation_Annotations %>% 
    dplyr::select(Hugo_Symbol, 
                  logFC) %>% 
    distinct()
  
  # Convert URA data to long format needed for visualization
  ura <- as_tibble(dataURA, 
                   rownames = NA) %>% 
    rownames_to_column(var = "Genes") %>% 
    pivot_longer(cols = !c("Genes"), 
                 names_to = "Biological_process", 
                 values_to = "Moonlight_gene_z_score") 
  
  # Count total number of hypo/hyper/dual methylated CpG sites 
  # for each oncogenic mediator
  combined_data <- Oncogenic_mediators_methylation_summary %>% 
    dplyr::select(-c(contains("NCG"),
                     No_methyl_num)) %>% 
    mutate(Total_methyl_num = rowSums(dplyr::select(.,
                                                    contains("methyl_num")), 
                                      na.rm = TRUE)) %>% 
    #join with ura data
    left_join(x = .,
              y = ura,
              by = c("Hugo_Symbol" = "Genes")) %>% 
    #join with DEA data
    left_join(x = .,
              y = DEGs_logFC,
              by = "Hugo_Symbol") %>% 
    #filter data to contain only genes supplied in input
    filter(Hugo_Symbol %in% genes)
  
  # Get number of Moonlight_Oncogenic_Mediator groups in subsetted data
  driver_groups <- combined_data$Moonlight_Oncogenic_Mediator %>% 
    unique %>% 
    length
  
  # If data contains both oncogenic mediator gene types
  if (driver_groups == 2) {
    
    # Create heatmap of genes grouped by oncogenic mediator types
    met_heatmap <- combined_data %>% 
      group_by(Moonlight_Oncogenic_Mediator) %>% 
      heatmap(.data = .,
              .row = Biological_process,
              .column = Hugo_Symbol, 
              .value = Moonlight_gene_z_score,
              scale = "none",
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              cluster_rows = FALSE) %>%
      add_tile(Moonlight_Oncogenic_Mediator, 
               palette = c("goldenrod2", "dodgerblue3")) %>% 
      add_tile(logFC, 
               palette = c("chartreuse4",  "firebrick3")) %>% 
      add_bar(Total_methyl_num) 
    save_pdf(met_heatmap, 
             height = 15, 
             width = 35, 
             units = "cm", 
             filename = paste(additionalFilename,
                              "moonlight_heatmap_met.pdf", 
                              sep = ""))
    
    # If data only contains one oncogenic mediator gene type
  } else {
    
    # Create heatmap of genes 
    met_heatmap <- combined_data %>% 
      heatmap(.data = .,
              .row = Biological_process,
              .column = Hugo_Symbol, 
              .value = Moonlight_gene_z_score,
              scale = "none",
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              cluster_rows = FALSE) %>% 
      add_tile(logFC, 
               palette = c("chartreuse4",  "firebrick3")) %>% 
      add_bar(Total_methyl_num) 
    save_pdf(met_heatmap, 
             height = 15, 
             width = 35, 
             units = "cm", 
             filename = paste(additionalFilename,
                              "moonlight_heatmap_met.pdf", 
                              sep = ""))
    
  }
  
}


