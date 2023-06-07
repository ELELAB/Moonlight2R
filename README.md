
Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

Repository associated to the publications:

Interpreting pathways to discover cancer driver genes with Moonlight. Colaprico A, Olsen C, Bailey MH, Odom GJ, Terkelsen T, Silva TC, Olsen AV, Cantini L, Zinovyev A, Barillot E, Noushmehr H, Bertoli G, Castiglioni I, Cava C, Bontempi G, Chen XS, Papaleo E. Nat Commun. 2020 Jan 3;11(1):69. doi: 10.1038/s41467-019-13803-0., PMID: 31900418

An Automatized Workflow to Study Mechanistic Indicators for Driver Gene Prediction with Moonlight2 Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico, Catharina Olsen, Matteo Tiberti, Elena Papaleo*, bioRxiv 2022. doi: 10.1101/2022.11.18.517066

contacts for repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Matteo Tiberti: tiberti-at-cancer.dk

## Introduction

This repository contains a new release fo our MoonlightR R package, called
Moonlight2R, which implements the new and improved Moonlight2 workflow. Moonlight2R
contains a number of differences and improvements respect to the original MoonlightR
package. For a full overview of what the package can do, please see and cite when
appropriate:

> Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico,
> Catharina Olsen, Matteo Tiberti, Elena Papaleo*
> An Automatized Workflow to Study Mechanistic Indicators for Driver Gene 
> Prediction with Moonlight2, bioRxiv 2022. doi: 10.1101/2022.11.18.517066

### Changes respect to MoonlightR

- The following user-facing functions were added to Moonlight2R:
    - `DMA` and `plotDMA`, for the new driver mutation analysis layer
    - `GLS`, for automated gene literature search for the most interesting
    identified driver genes
    - `plotMoonlight` function to generate a heatmap of Moonlight gene z-scores
    for selected genes
    - new helper functions: `LiftMAF`, `MAFtoCscape`, `PRAtoTibble`, `RunCscape_somatic`, 
`confidence`, `plotHeatmap`, and `tabix_func`

- Deprecation and removal of the following functions, which were available in 
MoonlightR: `DPA` and `getDataTCGA`. This means that users now needs to supply
their own differential expression analysis results (see vignettes) to run
predictions with Moonlight2R

- The `moonlight` function, which implements the whole Moonlight pipeline in a
single function,  was updated accordingly to the new changes in Moonlight2R

- Example data has been updated which follows changes in the functions. 
Specifically, the following data files have been added to Moonlight2R: 
`DEG_Mutations_Annotations.rda`, `EncodePromoters.rda`, `LOC_protein.rda`, `LOC_transcription.rda`, 
`LOC_translation.rda`, `NCG.rda`, `Oncogenic_mediators_mutation_summary.rda`, `cscape_somatic_output.rda`, 
`dataDMA.rda`, `dataGLS.rda`, and `dataMAF.rda`. The following data files were not included in 
Moonlight2R: `GDCprojects.rda` and `geneInfo.rda` as these files were connected to the deleted functions. 

- the vignette has been updated in light of changes implemented in Moonlight2R

### Installation from BioConductor

To install this package from the BioConductor repositories, start R (version "4.3")
and enter:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Moonlight2R")
```

### Installation from GitHub
```R
devtools::install_github(repo = "ELELAB/Moonlight2R")
```

#### Installation from GitHub with accompanying vignette

You need the BiocStyle Bioconductor package to install Moonlight2R with the vignette.
This package can be installed like:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BiocStyle")
```

To install Moonlight2R with its vignette:
```R
devtools::install_github(repo = "ELELAB/Moonlight2R", build_vignettes = TRUE)
```

To view the vignette:
```R
library(Moonlight2R)
```

```R
vignette( "Moonlight2R", package="Moonlight2R")
```

### Installation of other requirements

To make the best out of Moonlight2R the user should download the [CSCapeSomatic
pre-computed scores](http://cscape-somatic.biocompute.org.uk/#download). More
detailed instructions are available in the [INSTALL file](./INSTALL)