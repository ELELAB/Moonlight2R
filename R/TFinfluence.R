#' TFinfluence
#' 
#' This function finds mutations in the transcription factors (TF) of the DEGs that have a TF in the database.
#' 
#' @param dataTRRUST RDA object of gene-TF pairs
#' Must contain the following columns:
#' \itemize{
#' \item TF (HUGO symbol of TF)
#' \item GENE (HUGO symbol of gene)
#' }
#' @param dataMAF A MAF file RDA object
#' The MAF file must at least contain the following columns:
#' \itemize{
#' \item Hugo_Symbol eg. BRCA1
#' \item HGVSp_Short eg. p.V83F
#' }
#' @param dataDEGs Output DEA function
#' @param dataMAVISp Output loadMAVISp function
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
#' @export



TFinfluence <- function(dataTRRUST,
                         dataMAF,
                         dataDEGs,
                         dataMAVISp,
                         dataTFexpr = FALSE){ 
    # Control user input -------------
    # dataTRRUST
    if (is.null(dim(dataTRRUST))) {
        stop("The transcription factor data must be a non-empty table")
    }

    trrust_columns <- c('TF', 'Target', 'InteractionType')

    if (all(trrust_columns %in% names(dataTRRUST)) == FALSE) {
        stop("TRRUST dataframe does not contain the correct columns")
    }

    # dataMAF
    if (is.null(dim(dataMAF))) {
        stop("The mutation data must be a non-empty table")
    }

    maf_columns <- c("Hugo_Symbol",
                "HGVSp_Short")

    if (all(maf_columns %in% names(dataMAF)) == FALSE) {
        stop("MAF file does not contain the correct columns")
    }

    # dataDEGs
    if (is.null(dim(dataDEGs))) {
        stop("The DEG data must be a non-empty table")
    }

    # dataTFexpr
    if (!(dataTFexpr == FALSE) & is.null(dim(dataTFexpr))) {
        stop("The TF expression data must be a non-empty table")
    }
    

    # Load data --------------------------------
    # Read maf and add ID number to each mutation
    dataMAFFiltered <- dataMAF |> 
        select(c("Hugo_Symbol",
                "HGVSp_Short")) |>
        mutate(HGVSp_Short = str_extract(HGVSp_Short, pattern = "[A-Z]\\d+[A-Z]")) |>
        rename(mutation = HGVSp_Short)

    # Filter MAVISp data 
    dataMAVISpFiltered <- mavispFiltering(dataMAVISp)


    # Analysis -------------------
    # convert rownames to column for DEGs
    dataDEGs <- dataDEGs |>
        rownames_to_column(var = "GENE") |>
        select(GENE, logFC) 

    # join DEGs with TF
    DEG_TF <- dataDEGs |>
        left_join(dataTRRUST,
                  by = join_by(GENE == Target)) |>
        drop_na() |>
        rename('logFC_gene' = logFC)

    # Map TF to mutation file
    mut_DEG_TF <- DEG_TF |>
        left_join(dataMAFFiltered,
                  by = join_by(TF == Hugo_Symbol),
                  relationship = "many-to-many") |>
        rename('tf_mutation' = mutation)

    # Match TF-mut with mavisp to see the effect
    mavisp_mut_DEG_TF <- mut_DEG_TF |>
        left_join(dataMAVISpFiltered, 
                  by = join_by(TF == protein, tf_mutation == mutation))

    
    
    # Check biological implications
    # Activation -> destabilising mutation -> decrease
    # Repression -> destabilising mutation -> increase
    filtered_effect <- mavisp_mut_DEG_TF |>
        filter((InteractionType == 'Activation' & stab_class == 'Destabilizing' & logFC_gene < 0) |
                (InteractionType == 'Repression' & stab_class == 'Destabilizing' & logFC_gene > 0) |
                (stab_class == 'Uncertain'))

    if (!(dataTFexpr == FALSE)){
        filtered_effect <- filtered_effect |>
            left_join(dataTFexpr, 
                  by = join_by(TF == GENE))
    }
    
    return(mavisp_mut_DEG_TF)
}

mavispFiltering <- function(dataMAVISp){
    # Keep only the stability classification
    dataMAVISpFiltered <- dataMAVISp |>
                    map(function(x) rename(x, 'stab_class' = matches('(Stability classification, [A-Za-z0-9]+, \\(Rosetta, FoldX\\))'))) |>
                    keep(function(x) 'stab_class' %in% colnames(x)) |>
                    map(function(x) select(x, 1, stab_class) |>
                                    rename('mutation' = 1)) |>
                    rbindlist(idcol = 'protein') |>
                    as_tibble()
    return(dataMAVISpFiltered)
}
