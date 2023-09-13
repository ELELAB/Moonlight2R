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
#' @param BPlist A vector of strings. Selection of the biological processes to visualise.
#'  If left empty defaults to every BP provided in the URA file.
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
#' data(DEG_Mutations_Annotations)
#' data(Oncogenic_mediators_mutation_summary)
#' data(dataURA_plot)
#' plotMoonlight(DEG_Mutations_Annotations,
#'               Oncogenic_mediators_mutation_summary,
#'               dataURA_plot, gene_type = "drivers", n = 10,
#'               additionalFilename = "myplot_")
plotMoonlight <- function(DEG_Mutations_Annotations,
                          Oncogenic_mediators_mutation_summary,
                          dataURA,
                          gene_type = "drivers",
                          n = 50,
                          genelist = c(),
                          BPlist = c(),
                          additionalFilename = "") {

  # Check user input

  if (all(c("Hugo_Symbol", "Moonlight_gene_z_score", "logFC") %in% colnames(DEG_Mutations_Annotations)) == FALSE) {
    stop("DEG_Mutations_Annotations must contain Hugo_Symbol,
Moonlight_gene_z_score, and logFC as column names")
  }

  if (c("Moonlight_Oncogenic_Mediator") %in% colnames(Oncogenic_mediators_mutation_summary) == FALSE) {
    stop("Oncogenic_mediators_mutation_summary must contain Moonlight_Oncogenic_Mediator
as column name")
  }

  if (is.null(dim(dataURA))) {
    stop("The URA data must be non-empty with genes in rows and BPs in columns")
  }

  if (!is.null(gene_type) && (gene_type %in% c("mediators", "drivers")) == FALSE) {
    stop("Gene type must either be NULL, mediators or drivers")
  }

  if (!is.null(genelist) & !is.character(genelist)) {
    stop("Genelist must be either NULL or a character vector containing gene names")
  }

  if (!is.null(BPlist) & !is.character(BPlist)) {
    stop("BPlist must be either NULL or a character vector containing BP names")
  }

  if (!is.character(additionalFilename)) {
    stop("additionalFilename must be a character vector adding a prefix or
filepath to the filename of the pdf")
  }

  # The differentially expressed genes that are annotated as TSG/OCG
  DEGs <- DEG_Mutations_Annotations %>%
    select(Hugo_Symbol, Moonlight_gene_z_score, logFC) %>%
    unique() %>%
    drop_na(Moonlight_gene_z_score)

  # Restructure URA to tibble
  ura <- as_tibble(dataURA, rownames = NA) %>%
    rownames_to_column(var = "Genes")

  ura_wrangled <- ura %>%
    pivot_longer(cols = !c("Genes"),
                 names_to = "Biological_Process",
                 values_to = "Moonlight_score") %>%
    right_join(Oncogenic_mediators_mutation_summary,
               by = c("Genes" = "Hugo_Symbol")) %>%
    right_join(DEGs, by = c("Genes" = "Hugo_Symbol")) %>%
    replace_na(list(CScape_Driver = 0,
                    CScape_Passenger = 0,
                    CScape_Unclassified = 0))

  # Biological processes
  n_BP <- ura_wrangled %>%
    group_by(Biological_Process) %>%
    summarise() %>%
    count() %>%
    pull()

  if (length(BPlist) > 0) {
    ura_wrangled <- ura_wrangled %>%
      filter(Biological_Process %in% BPlist)
    n_BP <- length(BPlist)
  }

  n <- n * n_BP

  # Type of plot:
  if (length(genelist) > 0) {
    ura_wrangled <- ura_wrangled %>%
      filter(Genes %in% genelist)
  } else if (gene_type == "mediators") {
    ura_wrangled <- ura_wrangled %>%
      slice_max(Total_Mutations, n = n, with_ties = FALSE)
  } else {
    ura_wrangled <- ura_wrangled %>%
      slice_max(CScape_Driver, n = n, with_ties = FALSE)
  }

  n_types <- ura_wrangled %>%
    group_by(Moonlight_Oncogenic_Mediator) %>%
    summarise() %>%
    count() %>%
    pull()

  if (n_types <= 1) {
    stop("Must contain at least one OCG and one TSG. Increase n or add gene to genelist")
  }

  # Variable for color scaling in legend
  max_driver <- ura_wrangled %>%
    arrange(desc(CScape_Driver)) %>%
    select(CScape_Driver) %>%
    head(1) %>%
    pull

  # Plot Heatmap
  bp_heatmap <- heatmap(ura_wrangled,
                        .row = Biological_Process,
                        .column = Genes,
                        .value = Moonlight_score,
                        scale = "none",
                        clustering_distance_columns = "euclidean",
                        clustering_method_columns = "complete",
                        cluster_rows = FALSE) %>%
    add_tile(Moonlight_Oncogenic_Mediator,
             palette = c("goldenrod2", "dodgerblue3")) %>%
    add_tile(logFC, palette = c("chartreuse4", "firebrick3")) %>%
    add_tile(CScape_Driver, palette = colorRamp2(c(0, max_driver), c("white", "dodgerblue3"))) %>%
    add_bar(Total_Mutations)

  save_pdf(bp_heatmap, height = 15, width = 35, units = "cm",
           filename = paste(additionalFilename, "moonlight_heatmap.pdf", sep = ""))
}

utils::globalVariables(c("Hugo_Symbol", "Moonlight_gene_z_score", "logFC",
                         "Biological_Process", "Genes", "Total_Mutations",
                         "CScape_Driver", "Moonlight_Oncogenic_Mediator",
                         "head", "Moonlight_score"))
