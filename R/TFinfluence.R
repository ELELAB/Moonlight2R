#' TFinfluence
#' 
#' This function finds mutations in the transcription factors (TF) of the DEGs that have a TF in the database.
#' 
#' @param dataTRRUST Tibble containing gene-TF pairs and type of interaction
#' Must contain the following columns:
#' \itemize{
#' \item TF (HUGO symbol of TF)
#' \item Target (HUGO symbol of target gene of TF)
#' \item InteractionType (The effect the TF has on the target gene (Activation, Repression))
#' }
#' @param dataMAF Tibble containing mutation info from MAF
#' Must at least contain the following columns:
#' \itemize{
#' \item Hugo_Symbol eg. BRCA1
#' \item HGVSp_Short eg. p.V83F
#' }
#' @param dataDEGs Output DEA function. Tibble containing differentially expressed genes.
#' Must contain the following columns
#' \itemize{
#' \item GENE (HUGO symbol of DEG)
#' \item logFC (The log fold change of DEG)
#'}
#' @param dataPRA Output PRA function. List of TSG and OCG
#' Must contain the following elements
#' \itemize{
#' \item TSG (The TSGs identified by moonlight along with the moonlight score)
#' \item OCG (The OCGs identified by moonlight along with the moonlight score)
#'}
#' @param dataMAVISp Output loadMAVISp function. List of tibbles, one for each protein.
#' The tibbles must contain at least
#' \itemize{
#' \item (First column must contain the mutation (e.g. A54W). It is assumed that the column name is empty)
#'}
#' Then the tibbles must contain columns mathing either of the following names(or they will be excluded) 
#' \itemize{
#' \item (Stability classification, [A-Za-z0-9]+, \\(Rosetta, FoldX\\)) (Values: Stabilizing, Neutral, Destabilizing, Uncertain)
#' \item (Stability classification, [A-Za-z0-9]+, \\(RaSP, FoldX\\)) (Values: Stabilizing, Neutral, Destabilizing, Uncertain)
#'}
#' @param stabClassMAVISp The protocol to use for mutation stability classification. Default is rosetta.
#' Accepts one of the following strings
#' \itemize{
#' \item rasp (uses the FoldX/Rosetta protocol)
#' \item rasp (uses the FoldX/RaSP protocol)
#' }
#' 
#' @import dplyr
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @importFrom purrr map keep
#' @importFrom data.table rbindlist
#' @importFrom stringr str_extract
#' 
#' @return returns a tibble containing:
#' \itemize{
#' \item GENE
#' \item TF 
#' \item InteractionType 
#' \item PMID
#' \item tf_mutation
#' \item stab_class (the effect on stability as classified by MAVISp)
#' }
#'
#' @examples
#'
#' TFinfluence(dataTRRUST = dataTRRUST,
#'            dataMAF = dataMAF,
#'            dataDEGs = dataDEGs,
#'            dataPRA = dataPRA,
#'            dataMAVISp = dataMAVISp)
#'
#' @export

TFinfluence <- function(dataPRA,
                        dataDEGs,
                        dataTRRUST,
                        dataMAF,
                        dataMAVISp,
                        stabClassMAVISp = 'rosetta'){ 
    # Control user input -------------
    # dataPRA
    if (all(names(dataPRA) %in% c("TSG", "OCG")) == FALSE) {
        stop("The two list elements in PRA data must be named TSG and OCG")
    }

    # dataDEGs
    if (is.null(dim(dataDEGs)) | nrow(dataDEGs) == 0) {
        stop("The DEG data must be a non-empty table")
    }

    # dataTRRUST
    if (is.null(dim(dataTRRUST)) | nrow(dataTRRUST) == 0) {
        stop("The transcription factor data must be a non-empty table")
    }

    trrust_columns <- c('TF', 'Target', 'InteractionType')

    if (all(trrust_columns %in% names(dataTRRUST)) == FALSE) {
        stop("TRRUST dataframe does not contain the correct columns")
    }

    # dataMAF
    if (is.null(dim(dataMAF)) | nrow(dataMAF) == 0) {
        stop("The mutation data must be a non-empty table")
    }

    maf_columns <- c("Hugo_Symbol",
                "HGVSp_Short")

    if (all(maf_columns %in% names(dataMAF)) == FALSE) {
        stop("MAF file does not contain the correct columns")
    }

    # stabClassMAVISp
    if (!(stabClassMAVISp %in% c('rasp', 'rosetta'))){
        stop("The protocol for stability classification is not specified correctly. Accepts strings 'rasp' or 'rosetta'")
    } else {
        stabClassMAVISp <- paste0('stab_class_', stabClassMAVISp)
    }

    # Load data --------------------------------
    drivers <- PRAtoTibble(dataPRA)

    # Read maf and add ID number to each mutation
    dataMAFFiltered <- dataMAF |> 
        select(c("Hugo_Symbol",
                "HGVSp_Short")) |>
        mutate(HGVSp_Short = str_extract(HGVSp_Short, pattern = "[A-Z]\\d+[A-Z]")) |>
        rename(mutation = HGVSp_Short) |>
        drop_na()

    # Filter MAVISp data 
    # Keep only the stability classification
    dataMAVISpFiltered <- dataMAVISp |>
                    map(function(x) rename(x, 'stab_class_rosetta' = matches('Stability classification, [A-Za-z0-9, ]*\\(Rosetta, FoldX\\)'),
                                              'stab_class_rasp' = matches('Stability classification, [A-Za-z0-9, ]*\\(RaSP, FoldX\\)'))) |>
                    keep(function(x) stabClassMAVISp %in% colnames(x)) |>
                    map(function(x) rename(x, 
                                           'mutation' = 1,
                                           'stab_class' = stabClassMAVISp) |>
                                    select(mutation, stab_class)) |>
                    rbindlist(idcol = 'protein') |>
                    as_tibble()

    # Analysis -------------------
    # Convert rownames to column for DEGs
    dataDEGs <- dataDEGs |>
        rownames_to_column(var = "GENE") |>
        select(GENE, logFC) 

    # Join drivers with expression
    drivers_expr <- drivers |>
        left_join(dataDEGs, 
                  by = join_by(Hugo_Symbol == GENE))

    # Join drivers with TF
    drivers_expr_tf <- drivers_expr |>
        left_join(dataTRRUST,
                  by = join_by(Hugo_Symbol == Target)) |>
        rename('logFC_target' = logFC,
               'Target' = Hugo_Symbol)

    # Map mutation file to TF
    drivers_TF_mut <- drivers_expr_tf |>
        left_join(dataMAFFiltered,
                  by = join_by(TF == Hugo_Symbol),
                  relationship = "many-to-many") |>
        rename('tf_mutation' = mutation)

    # Match TF-mut with mavisp to see the effect
    drivers_mut_mavisp <- drivers_TF_mut |>
        left_join(dataMAVISpFiltered, 
                  by = join_by(TF == protein, tf_mutation == mutation))
    
    # Check biological implications
    # Activation -> destabilising mutation -> decrease
    # Repression -> destabilising mutation -> increase
    filtered_effect <- drivers_mut_mavisp |>
        filter((InteractionType == 'Activation' & stab_class == 'Destabilizing' & logFC_target < 0) |
                (InteractionType == 'Repression' & stab_class == 'Destabilizing' & logFC_target > 0) |
                (stab_class == 'Uncertain') |
                (stab_class == 'Neutral'))

    return(drivers_mut_mavisp)
}

