
Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

Repository associated to the publications:

Interpreting pathways to discover cancer driver genes with Moonlight. Colaprico A, Olsen C, Bailey MH, Odom GJ, Terkelsen T, Silva TC, Olsen AV, Cantini L, Zinovyev A, Barillot E, Noushmehr H, Bertoli G, Castiglioni I, Cava C, Bontempi G, Chen XS, Papaleo E. Nat Commun. 2020 Jan 3;11(1):69. doi: 10.1038/s41467-019-13803-0., PMID: 31900418

An Automatized Workflow to Study Mechanistic Indicators for Driver Gene Prediction with Moonlight2 Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico, Catharina Olsen, Matteo Tiberti, Elena Papaleo*, bioRxiv 2022. doi: 10.1101/2022.11.18.517066

contacts for repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Matteo Tiberti: tiberti-at-cancer.dk

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


### Changes between MoonlightR and Moonlight2R

The following main functions were added to Moonlight2R to reflect its new features: 
`DMA.R`, `GLS.R`, `plotDMA.R`, and `plotMoonlight.R` which are associated with the following 
new helper functions: `LiftMAF.R`, `MAFtoCscape.R`, `PRAtoTibble.R`, `RunCscape_somatic.R`, 
`confidence.R`, `plotHeatmap.R`, and `tabix_func.R`. 

The following functions were not included in Moonlight2R: `DPA.R` and `getDataTCGA.R`. 
They were removed to make Moonlight2R independent of TCGAbiolinks. 

The `moonlight.R` function, which wraps the whole Moonlight pipeline in a single function, 
was updated to reflect the new implementations in Moonlight2R. 

Accordingly, example data has been updated which follows changes in the functions. 
Specifically, the following data files have been added to Moonlight2R: 
`DEG_Mutations_Annotations.rda`, `EncodePromoters.rda`, `LOC_protein.rda`, `LOC_transcription.rda`, 
`LOC_translation.rda`, `NCG.rda`, `Oncogenic_mediators_mutation_summary.rda`, `cscape_somatic_output.rda`, 
`dataDMA.rda`, `dataGLS.rda`, and `dataMAF.rda`. The following data files were not included in 
Moonlight2R: `GDCprojects.rda` and `geneInfo.rda` as these files were connected to the deleted functions. 

Finally, the vignette has been updated in light of changes implemented in Moonlight2R. 
