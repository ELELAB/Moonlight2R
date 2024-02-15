#' plotGMA
#' This function plots results of the Gene Methylation Analysis.
#' It visualizes the number of hypo/hyper/dual methylated CpG
#' sites in oncogenic mediators or in a user-supplied gene list. 
#' The results are visualized either in a single heatmap or split 
#' into different ones which is specified in the function's three 
#' modes: split, complete and genelist.
#' @param DEG_Methylation_Annotations A tibble which is outputted
#' from the GMA function.
#' @param Oncogenic_mediators_methylation_summary A tibble which
#' is outputted from the GMA function.
#' @param type A character string which can either be split,
#' complete or genelist. If type is set to split, the entire
#' dataset is split into groups of 40 genes and individual
#' heatmaps of groups each containing 40 genes will be created and
#' subsequently merged into one pdf where each page in the pdf
#' is an individual heatmap. The genes will be sorted alphabetically. 
#' If type is set to complete, a single heatmap is created where
#' the number of differentially methylated CpG sites are shown for
#' all oncogenic mediators. If type is set to genelist, a single
#' heatmap will be created for genes supplied by the user in the
#' genelist parameter. Default is split. 
#' @param genelist A character string containing HUGO symbols
#' of genes to be visualized in a single heatmap. Default is NULL.
#' @param additionalFilename A character string that can be used
#' to add a prefix or filepath to the filename of the pdf 
#' visualizing the heatmap. Default is an empty string. 
#' @import dplyr
#' @import tidyr
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' @importFrom qpdf pdf_combine
#' @importFrom purrr map
#' @importFrom magrittr "%>%"
#' @return No value is returned. Visualizations in form of heatmaps
#' are saved.
#' @export
#' @examples 
#' data("DEG_Methylation_Annotations")
#' data("Oncogenic_mediators_methylation_summary")
#' genes <- c("ACAN", "ACE2", "ADAM19", "AFAP1L1")
#' plotGMA(DEG_Methylation_Annotations = DEG_Methylation_Annotations, 
#'         Oncogenic_mediators_methylation_summary = Oncogenic_mediators_methylation_summary, 
#'         type = "genelist", genelist = genes, 
#'         additionalFilename = "./GMAresults/")
plotGMA <- function(DEG_Methylation_Annotations,
                    Oncogenic_mediators_methylation_summary,
                    type = "split",
                    genelist = NULL,
                    additionalFilename = "") {

  # Check type input
  if (!type %in% c("split", "complete", "genelist")) {
    stop("type must be a character string of either split, complete or genelist")
  }
  
  # Check genelist input
  if (type == "genelist") {
    if (class(genelist) != "character") {
      stop("genelist must be a character string")
    }
    if (length(genelist) <= 0) {
      stop("genelist must not be empty")
    }
  }
  
  # Get logFC of DEGs
  DEGs_logFC <- DEG_Methylation_Annotations %>% 
    dplyr::select(Hugo_Symbol, 
                  logFC) %>% 
    distinct()
  
  # Reformat Oncogenic_mediators_methylation_summary data to long format
  # needed for visualization in heatmap
  summary_heatmap <- Oncogenic_mediators_methylation_summary %>% 
    pivot_longer(data = ., 
                 cols = contains("methyl_num"),
                 names_to = "Met_type",
                 values_to = "Met_count") %>% 
    dplyr::select(-contains("NCG")) %>% 
    replace_na(list(Met_count = 0)) %>% 
    left_join(x = .,
              y = DEGs_logFC,
              by = "Hugo_Symbol")
  
  ## Genelist mode
  
  # If genelist mode is chosen
  if (type == "genelist") {
    
    # Subset data to contain only supplied genes in genelist
    summary_heatmap <- summary_heatmap %>% 
      filter(Hugo_Symbol %in% genelist)
    
    # Get number of Moonlight_Oncogenic_Mediator groups in subsetted data
    driver_groups <- summary_heatmap$Moonlight_Oncogenic_Mediator %>% 
      unique %>% 
      length

    met_heatmap <- summary_heatmap
    
    # If genelist contain both oncogenic mediator gene types
    if (driver_groups == 2) {
       met_heatmap <- met_heatmap %>% group_by(Moonlight_Oncogenic_Mediator)
    }

   # Create heatmap of genes in genelist
   met_heatmap <- met_heatmap %>% heatmap(.data = .,
                                           .row = Met_type,
                                           .column = Hugo_Symbol, 
                                           .value = Met_count,
                                           cluster_rows = FALSE,
                                           cluster_columns = TRUE,
                                           show_column_dend = FALSE,
                                           show_column_names = TRUE,
                                           scale = "none",
                                           palette_value = c("white", "blue", "darkblue"),
                                           column_title = "Heatmap of methylation state by EpiMix") %>% 
      add_tile(logFC, 
               palette = c("chartreuse4",  "firebrick3"))

    # If genelist contain both oncogenic mediator gene types
    # add gene type annotation to heatmap
    if (driver_groups == 2) {
      met_heatmap <- met_heatmap %>% 
        add_tile(Moonlight_Oncogenic_Mediator, 
                 palette = c("goldenrod2", "dodgerblue3")) 
    }

    # Save heatmap
    save_pdf(met_heatmap, 
             height = 15, 
             width = 35, 
             units = "cm", 
             filename = paste(additionalFilename,
                              "heatmap_genelist_met.pdf", 
                              sep = ""))
 
  }
  
  ## Complete mode
  else if (type == "complete") {
    
    # Heatmap of all genes grouped by oncogenic mediator types
    met_heatmap <- summary_heatmap %>% 
      group_by(Moonlight_Oncogenic_Mediator) %>% 
      heatmap(.data = .,
              .row = Met_type,
              .column = Hugo_Symbol, 
              .value = Met_count,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              show_column_names = TRUE,
              scale = "none",
              palette_value = c("white", "blue", "darkblue"),
              column_title = "Heatmap of methylation state by EpiMix") %>% 
      add_tile(Moonlight_Oncogenic_Mediator, 
               palette = c("goldenrod2", "dodgerblue3")) %>% 
      add_tile(logFC, 
               palette = c("chartreuse4",  "firebrick3")) 
      save_pdf(met_heatmap, 
               height = 15, 
               width = 35, 
               units = "cm", 
               filename = paste(additionalFilename,
                                "heatmap_complete_met.pdf", 
                                sep = ""))
    
  }
  
  ## Split mode
  else if(type == "split") {
    
    if (dir.exists(paste(tempdir(), "/heatmaps", sep = ""))) {
      stop("'heatmaps' folder already exits. Please remove it.")
    }
    else {
      print("A temporary heatmaps folder is created. It will be removed before finalising.")
      dir.create(path = paste(tempdir(), "/heatmaps", sep = ""), 
                 showWarnings = TRUE, 
                 recursive = TRUE)
    }
    
    # Create vector of group assignments to be used for assigning genes
    # a group. Each group will be of length 160 corresponding to 
    # 4*40 = 160 rows for 40 genes as each gene has 4 rows. Last group
    # will not necessarily be of length 160 depending on number of
    # rows in data.
    split_vector <- rep(seq(1, 
                            ceiling(nrow(summary_heatmap) / 160), 
                            by = 1), 
                        each = 160)
    split_vector <- split_vector[1:nrow(summary_heatmap)]
    
    # Group data into groups of 40 genes
    # where each nested tibble will be a group of 40 genes
    grouped_data <- summary_heatmap %>% 
      arrange(Hugo_Symbol) %>% 
      ungroup() %>% 
      mutate(gr = split_vector) %>% 
      group_by(gr) %>% 
      tidyr::nest(data = -gr)
    
    print("Start creating indivual heatmaps, this can take a while")
    
    # Create heatmaps for each nested tibble in grouped_data
    # where each heatmap contains 40 genes grouped by oncogenic mediator type
    # Save resulting heatmap as a pdf
    split_heatmaps <- map(1:nrow(grouped_data), 
                          function(x) {
      
      met_heatmap <- grouped_data[[2]][[x]] %>% 
        group_by(Moonlight_Oncogenic_Mediator) %>% 
        heatmap(.data = .,
                .row = Met_type,
                .column = Hugo_Symbol, 
                .value = Met_count,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_column_dend = FALSE,
                show_column_names = TRUE,
                scale = "none",
                palette_value = c("white", "blue", "darkblue"),
                column_title = "Heatmap of methylation state by EpiMix") %>% 
        add_tile(Moonlight_Oncogenic_Mediator, 
                 palette = c("goldenrod2", "dodgerblue3")) %>% 
        add_tile(logFC, 
                 palette = c("chartreuse4",  "firebrick3"))  
        save_pdf(met_heatmap, 
                 height = 15, 
                 width = 35, 
                 units = "cm", 
                 filename = paste(tempdir(), "/heatmaps/heatmap_", 
                                  x, 
                                  ".pdf", 
                                  sep = ""))
      
    })
    
    print("Indivial heatmaps finished, start combining")
    
    # Get names of saved pdfs containing individual heatmaps
    pdf_names <- map(1:nrow(grouped_data), 
                     function(x) {
      pdf_names <- paste(tempdir(), "/heatmaps/heatmap_", 
                         x, 
                         ".pdf", 
                         sep="")
    }) %>% 
      unlist 
    
    # Combine individual heatmaps into a single pdf
    pdf_combine(input = pdf_names, 
                output = paste(additionalFilename,
                               "heatmaps_split_met.pdf", 
                               sep = ""))
    
    # Remove temporary directory
    unlink(paste(tempdir(), "/heatmaps", sep = ""), 
           recursive = TRUE)
    
  }
}
