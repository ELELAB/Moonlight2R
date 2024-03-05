#' plotDMA
#'
#' This function creates one or more heatmaps on the output from DMA.
#' It visualises the CScape-Somatic annotations per oncogenic mediator either
#' in a single heatmap or split into several different ones.
#' It is also possible to provide a personalised genelist to visualise.
#'
#' @param DEG_Mutations_Annotations A tibble, output file from DMA.
#' @param Oncogenic_mediators_mutation_summary A tibble, output file from DMA.
#' @param type A character string. It can take the values \code{"split"} or \code{"complete"}.
#' If both type and genelist are \code{NULL}, the function will default to \code{"split"}.
#' \itemize{
#' \item \code{"split"} will split the entire dataset into sections of 40 genes and create individual plots.
#' These plots will be merged into one pdf. The genes will be sorted alphabeatically.
#' \item \code{"complete"} will create one plot, though it will not be possible to see the individual gene names.
#' The heatmap will be clustered hierarchically.
#' }
#' @param genelist A character vector containing HUGO symbols.
#' A single heatmap will be created with only these genes.
#' The heatmap will be hierarchically clustered. This will overwrite \code{type}.
#'
#' @param additionalFilename A character string. Adds prefix or filepath to the filename of the pdf.
#'
#' @import dplyr
#' @import tidyr
#' @importFrom magrittr "%>%"
#' @import ComplexHeatmap
#' @import tidyHeatmap
#' @importFrom qpdf pdf_combine
#' @importFrom purrr map
#'
#' @return No return value. DMA results are plotted.
#' @export
#'
#' @examples
#' data('DEG_Mutations_Annotations')
#' data('Oncogenic_mediators_mutation_summary')
#' plotDMA(DEG_Mutations_Annotations,
#'         Oncogenic_mediators_mutation_summary,
#'         genelist = c("ACSS2", "AFAP1L1"), additionalFilename = "myplots_")
plotDMA <- function(DEG_Mutations_Annotations,
                    Oncogenic_mediators_mutation_summary,
                    type = "split",
                    genelist = c(),
                    additionalFilename = "") {

  # Check user input

  if (all(c("Hugo_Symbol", "logFC") %in% colnames(DEG_Mutations_Annotations)) == FALSE) {
    stop("DEG_Mutations_Annotations must contain Hugo_Symbol and logFC as column names")
  }

  if (all(c("CScape_Driver", "CScape_Passenger", "CScape_Unclassified", "Moonlight_Oncogenic_Mediator") %in% colnames(Oncogenic_mediators_mutation_summary)) == FALSE) {
    stop("Oncogenic_mediators_mutation_summary must contain CScape_Driver,
CScape_Passenger, CScape_Unclassified, and Moonlight_Oncogenic_Mediator as column names")
  }

  if (!is(type, "NULL") && (type %in% c("split", "complete")) == FALSE) {
    stop("Type must either be NULL, split or complete")
  }

  if (!is(genelist, "NULL") & !is(genelist, "character")) {
    stop("Genelist must be either NULL or a character vector containing gene names")
  }

  if (!is(additionalFilename, "character")) {
    stop("additionalFilename must be a character vector adding a prefix or filepath to the filename of the pdf")
  }

  # Modify input
  DEGs <- DEG_Mutations_Annotations %>%
    dplyr::select(Hugo_Symbol, logFC) %>%
    distinct()

  Summary_wrangled <- Oncogenic_mediators_mutation_summary %>%
    pivot_longer(cols = c(CScape_Driver, CScape_Passenger, CScape_Unclassified),
                 names_to = "Mutation_type",
                 values_to = "Count") %>%
    replace_na(list(Count = 0)) %>%
    left_join(DEGs)

  if (length(genelist) > 0) {

    Summary_wrangled <- Summary_wrangled %>%
      filter(Hugo_Symbol %in% genelist) %>%
      group_by(Moonlight_Oncogenic_Mediator)

    # Check that there are both OCGs and TGS
    n <- Summary_wrangled %>%
      summarise() %>%
      count() %>%
      pull()

    if (n > 1) {

      driver_mut_heatmap <- Summary_wrangled %>%
        heatmap(.row = Mutation_type,
                .column = Hugo_Symbol,
                .value = Count,
                scale = "none",
                cluster_rows = FALSE,
                cluster_columns = TRUE,
                show_column_dend = FALSE,
                show_column_names = TRUE,
                palette_value = c("white", "blue", "darkblue"),
                column_title = paste("Heatmap Driver annotation by CScape-Somatic",
                                     "\n Hugo_Symbol")) %>%
        add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>%
        add_bar(Total_Mutations) %>%
        add_tile(logFC, palette = c("chartreuse4", "firebrick3"))

      # Save plot
      save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
               filename = paste(additionalFilename, "heatmap_genelist.pdf",
               sep = ""))
    } else {
      stop("The genelist must contain at least one OCG and one TSG")
    }
  } else if (type == "complete") {

    driver_mut_heatmap <- Summary_wrangled %>%
      group_by(Moonlight_Oncogenic_Mediator) %>%
      heatmap(.row = Mutation_type,
              .column = Hugo_Symbol,
              .value = Count,
              scale = "none",
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              show_column_names = TRUE,
              palette_value = c("white", "blue", "darkblue"),
              column_title = paste("Heatmap Driver annotation by CScape-Somatic",
                                   "\n Hugo_Symbol")) %>%
      add_tile(Moonlight_Oncogenic_Mediator, palette = c("goldenrod2", "dodgerblue3")) %>%
      add_bar(Total_Mutations) %>%
      add_tile(logFC, palette = c("chartreuse4", "firebrick3"))

    # Save plot
    save_pdf(driver_mut_heatmap, height = 15, width = 35, units = "cm",
             filename = paste(additionalFilename, "heatmap_complete.pdf", sep = ""))

  } else if (type == "split") {

    # Create temporary heatmap folder
    if (dir.exists("./heatmaps")) {
      stop("'heatmaps' folder already exits. Please remove it.")
    } else {
      print("A temporary heatmaps folder is created. It will be removed before finalising.")
      dir.create(path = "./heatmaps", showWarnings = TRUE, recursive = TRUE)
    }

    # Make vector with groups (40 genes in each plot x 3 (driver, pas, unclas) per gene) = 120
    split_vector <- rep(seq(1, ceiling(nrow(Summary_wrangled) / 120), by = 1), each = 120)
    split_vector <- split_vector[seq.int(nrow(Summary_wrangled))]

    grouped_data <- Summary_wrangled %>%
      arrange(Hugo_Symbol) %>%
      ungroup() %>%
      mutate(gr = split_vector) %>%
      group_by(gr) %>%
      tidyr::nest(data = -gr)

    print("Start creating indivual heatmaps, this can take a while")
    test<- grouped_data %>%
      ungroup %>%
      mutate(First_gene = map(grouped_data$data, ~plotHeatmap(.))) %>%
      mutate(heatmaps = "heatmaps/heatmap_", pdf = ".pdf") %>%
      unite(heatmaps, c(heatmaps, First_gene, pdf), sep = "")
    print("Indivial heatmaps finished, start combining")

    # Stable pdf into one pdf
    pdfs <- test %>%
      pull(heatmaps)
    pdf_combine(input = pdfs, output = paste(additionalFilename, "heatmaps_split.pdf", sep = ""))

    # Remove heatmap-folder with redudant pdfs
    unlink("heatmaps", recursive = TRUE)

  }
}

utils::globalVariables(c("Hugo_Symbol", "logFC", "CScape_Driver", "CScape_Passenger",
                         "CScape_Unclassified", "Moonlight_Oncogenic_Mediator",
                         "Mutation_type", "Count", "Total_Mutations", "gr",
                         "heatmaps", "First_gene"))

