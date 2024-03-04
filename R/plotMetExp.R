#' plotMetExp
#' This function visualizes results of EpiMix. 
#' @param EpiMix_results The object, a list, that was returned from
#' running the EpiMix function and is one of the outputs from the 
#' GMA function.
#' @param probe_name A character string containing the name of the
#' CpG probe that will be plotted. 
#' @param dataMET A data matrix containing the methylation data
#' where the CpG probes are in the rows and samples are in the
#' columns
#' @param dataEXP A data matrix containing the gene expression
#' data where the genes are in the rows and the samples are in 
#' the columns
#' @param gene_of_interest A character string containing the 
#' name of the gene that will be plotted. 
#' @param additionalFilename A character string that can be used
#' to add a prefix or filepath to the filename of the pdf 
#' visualizing the heatmap. Default is an empty string. 
#' @importFrom EpiMix EpiMix_PlotModel
#' @importFrom ggplot2 ggsave
#' @return No value is returned. Visualizations are saved. 
#' @export
#' @examples 
#' data("EpiMix_Results_Regular")
#' data("dataMethyl")
#' data("dataFilt")
#' pattern <- "^(.{4}-.{2}-.{4}-.{2}).*"
#' colnames(dataFilt) <- sub(pattern, "\\1", colnames(dataFilt))
#' plotMetExp(EpiMix_results = EpiMix_Results_Regular,
#'            probe_name = "cg03035704",
#'            dataMET = dataMethyl,
#'            dataEXP = dataFilt,
#'            gene_of_interest = "ACVRL1",
#'            additionalFilename = "./GMAresults/")
plotMetExp <- function(EpiMix_results,
                       probe_name,
                       dataMET,
                       dataEXP,
                       gene_of_interest,
                       additionalFilename = "") {
  
  # Run EpiMix_PlotModel function
  EpiMix_plots <- EpiMix_PlotModel(EpiMixResults = EpiMix_results, 
                                   Probe = probe_name, 
                                   methylation.data = dataMET, 
                                   gene.expression.data = dataEXP, 
                                   GeneName = gene_of_interest)
  
  # Save mixture model plot
  ggsave(plot = EpiMix_plots$MixtureModelPlot, 
         width = 8,
         height = 8,
         filename = paste(additionalFilename,
                          "MixtureModelPlot.pdf", 
                          sep = ""),
         path = ".")
  
  # Save violin plot
  ggsave(plot = EpiMix_plots$ViolinPlot, 
         width = 8,
         height = 8,
         filename = paste(additionalFilename,
                          "ViolinPlot.pdf", 
                          sep = ""),
         path = ".")
  
  # Save correlation plot
  ggsave(plot = EpiMix_plots$CorrelationPlot, 
         width = 8,
         height = 8,
         filename = paste(additionalFilename,
                          "CorrelationPlot.pdf", 
                          sep = ""),
         path = ".")
}

