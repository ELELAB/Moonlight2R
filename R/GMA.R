#' GMA
#' This function carries out Gene Methylation Analysis
#' @param dataMET A data matrix containing the methylation data
#' where the CpG probes are in the rows and samples are in the
#' columns
#' @param dataEXP A data matrix containing the gene expression
#' data where the genes are in the rows and the samples are in 
#' the columns
#' @param dataPRA A table containing the output of the PRA 
#' function
#' @param dataDEGs A table containing the output of a DEA
#' where gene names are rownames
#' @param sample_info A table containing information on the
#' samples. This table needs to contain two columns called
#' primary and sample.type. The primary column contains 
#' sample names which should be the same as the column names
#' in dataMET. The sample.type column indicates for each
#' sample if it is a Cancer or Normal sample.
#' @param met_platform A character string representing the
#' microarray type that was used to collect the methylation
#' data. This can either be HM27, HM450 or EPIC. Default
#' is HM450. 
#' @param prevalence_filter An float or NULL representing if
#' a prevalence filter should be applied or not. Default is NULL,
#' meaning a prevalence filter will not be applied. If a float
#' is specified, a prevalence filter will be applied where 
#' methylation states of probes will be altered depending on the
#' threshold of prevalence supplied as prevalence_filter. For 
#' example, if prevalence_filter = 20, it means that if the prevalence 
#' of the hyper- or hypomethylated CpG probe exceeds 20, the methylation 
#' state will be unchanged but if the prevalence is lower than 20 the 
#' methylation state will be changed to NA, meaning no methylation state
#' was detected. In case of dual methylated probes, the methylation 
#' state will stay dual if both the prevalence of hyper- and 
#' hypomethylations exceed 20, but if only one of the prevalences exceed 
#' 20 the dual state will be changed to the state exceeding 20. If none 
#' of the prevalences exceed 20, the dual state will be changed to NA.  
#' @param output_dir Path to where the results will be stored.
#' If this directory does not exist, it will be created by the 
#' function. Default is ./GMAresults. 
#' @param cores Number of cores to be used. Default is 1. 
#' @param roadmap.epigenome.ids A character string representing
#' the epigenome ID that will be used to select enhancers. 
#' Since enhancers are tissue-specific, the tissue type needs 
#' to be specified in EpiMix. The enhancers are found from the
#' RoadmapEpigenome project and the IDs can be found from
#' Figure 2 in the publication with doi: 10.1038/nature14248.
#' Default is NULL.
#' @param roadmap.epigenome.groups A character string representing
#' the epigenome group that will be used to select enhancers. 
#' Details are provided above with the roadmap.epigenome.ids 
#' parameter. Default is NULL.
#' @import EpiMix
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import readr
#' @importFrom fuzzyjoin genome_left_join
#' @importFrom stringr str_c
#' @importFrom magrittr "%>%"
#' @importFrom BiocGenerics as.data.frame
#' @return List of two elements, containing predicted oncogenes
#' and tumor suppressors. Additionally, various output files are 
#' saved in the specified output directory: 
#' DEG_Methylation_Annotations.rda, Oncogenic_mediators_methylation_summary.rda,
#' EpiMix_Results_Enhancer.rds, EpiMix_Results_Regular.rds, 
#' FunctionalPairs_Enhancer.csv, FunctionalPairs_Regular.csv, 
#' FunctionalProbes_Regular.rds
#' @export
#' @examples 
#' data("dataMethyl")
#' data("dataFiltCol")
#' data("dataPRA")
#' data("DEGsmatrix")
#' data("LUAD_sample_anno")
#' data("NCG")
#' data("EncodePromoters")
#' data("MetEvidenceDriver")
#' dataGMA <- GMA(dataMET = dataMethyl, dataEXP = dataFiltCol, 
#' dataPRA = dataPRA, dataDEGs = DEGsmatrix, 
#' sample_info = LUAD_sample_anno, met_platform = "HM450",
#' prevalence_filter = NULL,
#' output_dir = "./GMAresults", cores = 1, roadmap.epigenome.ids = "E096", 
#' roadmap.epigenome.groups = NULL)
GMA <- function(dataMET, 
                dataEXP, 
                dataPRA, 
                dataDEGs, 
                sample_info, 
                met_platform = "HM450", 
                prevalence_filter = NULL,
                output_dir = "./GMAresults",
                cores = 1, 
                roadmap.epigenome.ids = NULL, 
                roadmap.epigenome.groups = NULL) {
  
  if (dir.exists(output_dir)) {
    print("Output folder already exits")
  }
  else {
    dir.create(path = output_dir, 
               showWarnings = TRUE,
               recursive = TRUE)
  }
  
  
  ## EpiMix ---------------------------------------------------------------
  
  ## Run EpiMix regular mode
  EpiMix_regular <- EpiMix(methylation.data = dataMET,
                           gene.expression.data = dataEXP,
                           sample.info = sample_info,
                           mode = "Regular",
                           group.1 = "Cancer",
                           group.2 = "Normal",
                           met.platform = met_platform,
                           correlation = "negative",
                           OutputRoot = output_dir,
                           cores = cores)
  
  ## Run EpiMix enhancer mode
  EpiMix_enhancer <- EpiMix(methylation.data = dataMET,
                            gene.expression.data = dataEXP,
                            sample.info = sample_info,
                            mode = "Enhancer",
                            roadmap.epigenome.ids = roadmap.epigenome.ids,
                            roadmap.epigenome.groups = roadmap.epigenome.groups,
                            group.1 = "Cancer",
                            group.2 = "Normal",
                            met.platform = met_platform,
                            correlation = "negative",
                            OutputRoot = output_dir,
                            cores = cores)
  
  # Convert functional pairs results from EpiMix to tibble
  EpiMix_regular_FP <- as_tibble(EpiMix_regular$FunctionalPairs) %>% 
    dplyr::rename("Hugo_Symbol" = "Gene")
  EpiMix_enhancer_FP <- as_tibble(EpiMix_enhancer$FunctionalPairs) %>% 
    dplyr::rename("Hugo_Symbol" = "Gene")
  
  # Add enhancer prefix to columns in EpiMix_enhancer_FP
  # Create column that specifies that identified probes are in enhancers 
  EpiMix_enhancer_FP <- EpiMix_enhancer_FP %>% 
    dplyr::rename_with(.cols = -1, 
                       .fn = ~str_c(., "_enhancer")) %>% 
    dplyr::mutate(Probe_in_enhancer = "Enhancer") %>% 
    dplyr::relocate(Probe_in_enhancer, 
                    .after = Probe_enhancer) 
  
  
  ## Initialize data --------------------------------------------------------
  
  # Convert PRA data to tibble
  PRA_tbl <- PRAtoTibble(dataPRA = dataPRA)
  
  # Get NCG data
  NCG <- get("NCG")
  
  
  ## DEG_Methylation_Annotations --------------------------------------------
  
  ## Intersect DEA results with EpiMix results
  ## Create DEG_Methylation_Annotations
  
  # Convert DEA data to tibble
  dataDEGs <- dataDEGs %>% 
    tibble::rownames_to_column(var = "Hugo_Symbol") %>% 
    as_tibble()
  
  # Bind EpiMix regular and enhancer functional-pairs results
  # Join with DEA data and with PRA data containing information 
  # about oncogenic mediators
  # Add NCG annotations
  DEG_Methylation_Annotations <- EpiMix_regular_FP %>% 
    bind_rows(EpiMix_enhancer_FP) %>% 
    right_join(x = .,
               y = dataDEGs,
               by = "Hugo_Symbol") %>% 
    left_join(x = .,
              y = PRA_tbl,
              by = "Hugo_Symbol") %>%
    left_join(x = .,
              y = NCG,
              by = c("Hugo_Symbol" = "symbol")) %>%
    dplyr::select(-c("Fold change of gene expression",
                     "Fold change of gene expression_enhancer"))
  
  ## Get chromosome, start and end locations of probes
  
  # Retrieve annotations of probes and select only chromosome, start and
  # end positions
  probe_annotation <- EpiMix_getInfiniumAnnotation(plat = met_platform, 
                                                   genome = "hg38") %>% 
    BiocGenerics::as.data.frame() %>% 
    as_tibble(rownames = "Probe") %>% 
    dplyr::select(Probe,
                  seqnames, 
                  start,
                  end) %>% 
    dplyr::rename(Chromosome = seqnames)
  
  # Add probe annotations (chromosome, start and end positions) 
  # to DEG_Methylation_Annotations 
  DEG_Methylation_Annotations <- DEG_Methylation_Annotations %>% 
    left_join(x = .,
              y = probe_annotation,
              by = "Probe") %>% 
    left_join(x = .,
              y = probe_annotation,
              by = c("Probe_enhancer" = "Probe")) %>% 
    unite(start, 
          start.x, 
          start.y, 
          na.rm = TRUE) %>% 
    unite(end, 
          end.x, 
          end.y, 
          na.rm = TRUE) %>% 
    unite(chromosome, 
          Chromosome.x, 
          Chromosome.y, 
          na.rm = TRUE) %>% 
    mutate(start = na_if(start, 
                         ""),
           end = na_if(end, 
                       ""),
           chromosome = na_if(chromosome, 
                              "")) %>% 
    dplyr::relocate(chromosome, 
                    start, 
                    end, 
                    .after = "Hugo_Symbol") %>% 
    mutate(start = as.numeric(start),
           end = as.numeric(end))
  
  ## Add Encode promoter annotations to DEG_Methylation_Annotations
  
  # Get promoter annotations and join with DEG_Methylation_Annotations
  promoters <- get("EncodePromoters") %>% 
    mutate(Annotation = "Promoter") %>% 
    dplyr::select(X1, 
                  X2, 
                  X3, 
                  Annotation) %>% 
    dplyr::rename(chromosome_annot = "X1", 
                  promoter_start = "X2", 
                  promoter_end = "X3")
  DEG_Methylation_Annotations <- genome_left_join(x = DEG_Methylation_Annotations,
                                                  y = promoters,
                                                  by = c(chromosome = "chromosome_annot", 
                                                         start = "promoter_start", 
                                                         end = "promoter_end")) %>% 
    dplyr::select(-c("chromosome_annot")) %>% 
    dplyr::rename(Probe_in_promoter = Annotation) %>% 
    dplyr::relocate(contains("promoter"),
                    .before = NCG_driver)
  
  
  ## Oncogenic_mediators_methylation_summary --------------------------------
  
  # Get MetEvidenceDriver
  MetEvidenceDriver <- get("MetEvidenceDriver")
  
  ## Intersect oncogenic mediator results with EpiMix functional pairs results
  ## Create Oncogenic_mediators_methylation_summary
  
  # If prevalence filter is not selected
  if (is.null(prevalence_filter)) {
    
    # Left join PRA data and functional pairs results from EpiMix
    # Count number of CpG probes that are hyper/hypo/dual methylated for 
    # each oncogenic mediator
    # Add oncogenic mediator role (putative TSG/OCG) for each oncogenic mediator
    Oncogenic_mediators_methylation_summary <- PRA_tbl %>% 
      left_join(x = ., 
                y = EpiMix_regular_FP, 
                by = "Hugo_Symbol") %>% 
      group_by(Hugo_Symbol,
               State) %>% 
      dplyr::count() %>% 
      pivot_wider(., 
                  id_cols = "Hugo_Symbol",
                  names_from = "State",
                  values_from = "n") %>% 
      dplyr::rename("No" = "NA") %>% 
      full_join(x = .,
                y = PRA_tbl,
                by = "Hugo_Symbol") %>% 
      dplyr::select(-"Moonlight_gene_z_score") %>% 
      dplyr::relocate(.data = ., 
                      "Moonlight_Oncogenic_Mediator",
                      .after = "Hugo_Symbol") %>% 
      replace(is.na(.), 0) %>% 
      ungroup
    
    # If prevalence filter is selected 
  } else {
    
    # Left join PRA data and functional pairs results from EpiMix
    # Update methylation state of probes based on prevalence
    # Count number of CpG probes that are hyper/hypo/dual methylated for 
    # each oncogenic mediator
    # Add oncogenic mediator role (putative TSG/OCG) for each oncogenic mediator
    Oncogenic_mediators_methylation_summary <- PRA_tbl %>% 
      left_join(x = ., 
                y = EpiMix_regular_FP, 
                by = "Hugo_Symbol") %>% 
      dplyr::select(Hugo_Symbol, Moonlight_Oncogenic_Mediator,
                    Probe, State, contains("Prevalence")) %>% 
      # Update methylation state of probes based on prevalence of
      # methylation in patients
      dplyr::mutate(State_upd = case_when(State == "Hyper" & 
                                            `Prevalence of hyper (%)` > prevalence_filter ~ "Hyper",
                                          State == "Hypo" & 
                                            `Prevalence of hypo (%)` > prevalence_filter ~ "Hypo",
                                          State == "Dual" &
                                            `Prevalence of hyper (%)` > prevalence_filter &
                                            `Prevalence of hypo (%)` > prevalence_filter ~ "Dual",
                                          State == "Dual" &
                                            `Prevalence of hyper (%)` < prevalence_filter &
                                            `Prevalence of hypo (%)` > prevalence_filter ~ "Hypo",
                                          State == "Dual" &
                                            `Prevalence of hyper (%)` > prevalence_filter &
                                            `Prevalence of hypo (%)` < prevalence_filter ~ "Hyper")) %>% 
      group_by(Hugo_Symbol,
               State_upd) %>% 
      dplyr::count() %>% 
      pivot_wider(., 
                  id_cols = "Hugo_Symbol",
                  names_from = "State_upd",
                  values_from = "n") %>% 
      dplyr::rename("No" = "NA") %>% 
      full_join(x = .,
                y = PRA_tbl,
                by = "Hugo_Symbol") %>% 
      dplyr::select(-"Moonlight_gene_z_score") %>% 
      dplyr::relocate(.data = ., 
                      "Moonlight_Oncogenic_Mediator",
                      .after = "Hugo_Symbol") %>% 
      replace(is.na(.), 0) %>% 
      ungroup
    
  }
  
  # If Hypo, Hyper, Dual or No columns do not exist, create them containing all 0s
  if (!"Hypo" %in% colnames(Oncogenic_mediators_methylation_summary)) {
    # Create the Hypo column and fill it with 0s
    Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>%
      mutate(Hypo = 0)
  }
  if (!"Hyper" %in% colnames(Oncogenic_mediators_methylation_summary)) {
    # Create the Hyper column and fill it with 0s
    Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>%
      mutate(Hyper = 0)
  }
  if (!"Dual" %in% colnames(Oncogenic_mediators_methylation_summary)) {
    # Create the Dual column and fill it with 0s
    Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>%
      mutate(Dual = 0)
  }
  if (!"No" %in% colnames(Oncogenic_mediators_methylation_summary)) {
    # Create the No column and fill it with 0s
    Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>%
      mutate(No = 0)
  }
  
  # Temporary copy of Oncogenic_mediators_methylation_summary where all values
  # greater than 1 are set to 1 to allow for overlap with MetEvidenceDriver 
  # table
  Oncogenic_mediators_methylation_summary_tmp <- Oncogenic_mediators_methylation_summary %>% 
    mutate(across(where(is.integer), ~ if_else(. > 0, 1, .))) %>% 
    replace(is.na(.), 0)
  
  # Left join Oncogenic_mediators_methylation_summary_tmp with MetEvidenceDriver
  Oncogenic_mediators_methylation_summary_tmp <- Oncogenic_mediators_methylation_summary_tmp %>% 
    left_join(MetEvidenceDriver, 
              by = c("Hypo", "Hyper", "No", "Dual", "Moonlight_Oncogenic_Mediator"))
  
  # Add Evidence column to original Oncogenic_mediators_methylation_summary
  Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>% 
    mutate(Evidence = Oncogenic_mediators_methylation_summary_tmp$Evidence)
  
  # Rename methylation state columns
  Oncogenic_mediators_methylation_summary <- Oncogenic_mediators_methylation_summary %>% 
    dplyr::rename("No_methyl_num" = "No",
                  "Hypo_methyl_num" = "Hypo",
                  "Hyper_methyl_num" = "Hyper",
                  "Dual_methyl_num" = "Dual")

  # Print number of oncogenic mediators with different level of evidence
  count_evidence <- map(Oncogenic_mediators_methylation_summary$Evidence %>% 
                          unique, 
                        function(x) { 
                          num_genes <- Oncogenic_mediators_methylation_summary %>% 
                            dplyr::filter(Evidence == x) %>% 
                            distinct(Hugo_Symbol) %>% 
                            pull %>% 
                            length 
                          print(paste(num_genes, 
                                      "oncogenic mediator(s) out of", 
                                      Oncogenic_mediators_methylation_summary %>% 
                                        distinct(Hugo_Symbol) %>% 
                                        pull %>% 
                                        length,
                                      "were found in the", 
                                      x, 
                                      "evidence category")) })
  
  
  ## Get driver genes --------------------------------
  
  # Get driver genes which are those oncogenic mediators that have a 
  # significant association between expression and methylation and
  # where pattern of presumed driver role matches with methylation
  # pattern
  TSG <- Oncogenic_mediators_methylation_summary %>% 
    dplyr::filter(Evidence == "Agreement",
                  Moonlight_Oncogenic_Mediator == "TSG") %>% 
    pull(Hugo_Symbol)
  OCG <- Oncogenic_mediators_methylation_summary %>% 
    dplyr::filter(Evidence == "Agreement",
                  Moonlight_Oncogenic_Mediator == "OCG") %>% 
    pull(Hugo_Symbol)
  
  
  ## Save data ------------------------------------------------------------
  
  # Save DEG_Methylation_Annotations as rda file
  save(DEG_Methylation_Annotations,
       file = paste(output_dir, 
                    "DEG_Methylation_Annotations.rda",
                    sep = "/"))
  
  # Save Oncogenic_mediators_methylation_summary as rda file
  save(Oncogenic_mediators_methylation_summary,
       file = paste(output_dir, 
                    "Oncogenic_mediators_methylation_summary.rda",
                    sep = "/"))
  
  return(list(TSG = TSG, 
              OCG = OCG))
  
  
}

