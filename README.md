
Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen

Repository associated to the publications:

Interpreting pathways to discover cancer driver genes with Moonlight. Colaprico A, Olsen C, Bailey MH, Odom GJ, Terkelsen T, Silva TC, Olsen AV, Cantini L, Zinovyev A, Barillot E, Noushmehr H, Bertoli G, Castiglioni I, Cava C, Bontempi G, Chen XS, Papaleo E. Nat Commun. 2020 Jan 3;11(1):69. doi: 10.1038/s41467-019-13803-0., PMID: 31900418

An Automatized Workflow to Study Mechanistic Indicators for Driver Gene Prediction with Moonlight2 Astrid Saksager, Mona Nourbakhsh, Nikola Tom, Xi Steven Chen, Antonio Colaprico, Catharina Olsen, Matteo Tiberti, Elena Papaleo*, bioRxiv 2022. doi: 10.1101/2022.11.18.517066 





contacts for repository: Elena Papaleo, elpap-at-dtu.dk, elenap-at-cancer.dk; Matteo Tiberti: tiberti-at-cancer.dk


# Identify oncogenes and tumor suppressor genes from genomic data.

### Installation from GitHub ###
```R
devtools::install_github(repo = "ELELAB/Moonlight2R")
```

#### Installation from GitHub with accompanying vignette ####

You need the BiocStyle Bioconductor package to install Moonlight2R with the vignette. This package can be installed like:  
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
vignette( "Moonlight2R", package="Moonlight2R")
```


