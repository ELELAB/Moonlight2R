#' plotHeatmap
#'
#' This function creates a unclustered heatmap from the inputted data tibble 
#' and saves it
#' 
#' @param df a tibble 
#' 
#' @import dplyr 
#' @importFrom magrittr "%>%"
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' 
#' @return The name of the alphabeatically first gene in the tibble
#' @export
#'

plotHeatmap <- function(df){

  # Check user input
  if (c("Moonlight_Oncogenic_Mediator") %in% colnames(Oncogenic_mediators_mutation_summary) == FALSE) {
    stop("Moonlight_Oncogenic_Mediator must be a column name in Oncogenic_mediators_mutation_summary")
  }

  first_gene <- pull(df[1,1])
  
  driver_mut_heatmap <-  df %>% 
    group_by(Moonlight_Oncogenic_Mediator) %>% 
    heatmap(.row = Mutation_type,
            .column = Hugo_Symbol,
            .value = Count,
            scale = "none",
            cluster_rows = FALSE, cluster_columns = FALSE, 
            show_column_dend = FALSE, 
            show_column_names = TRUE,
            palette_value = c("white", "blue", "darkblue"),
            column_title = paste("Heatmap Driver annotation by CScape-Somatic, from gene:",
                                 first_gene,
                                 "\n Hugo_Symbol")) %>% 
    add_tile(Moonlight_Oncogenic_Mediator, 
             palette = c("goldenrod2", "dodgerblue3")) %>% 
    add_bar(Total_Mutations) %>% 
    add_tile(logFC, palette = c("chartreuse4","firebrick3"))
  
  #Save plot
  save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
           filename = paste("heatmaps/heatmap_",first_gene , ".pdf", sep =""))
  return(first_gene)
}  

utils::globalVariables(c("Mutation_type", "Hugo_Symbol", "Count",
"Total_Mutations", "logFC", "Oncogenic_mediators_mutation_summary"))
