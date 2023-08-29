# Moonlight2R 0.99.10

## Summary

* updated code style in all functions 

# Moonlight2R 0.99.9

## Summary

* fixed hardcoded plot title in `plotFEA`
* updated reference in docs to published paper
* fixed excessively long example line in GLS function
* moved `globalVariables` calls after function definitions
* turned `LazyData` to `false` and issued connected fixes

# Moonlight2R 0.99.8

## Summary

* replaced number sequences generation (e.g. using `seq(n)` instead of `1:n`) in all functions

# Moonlight2R 0.99.7

## Summary

* added bindings for global variables in majority of functions

# Moonlight2R 0.99.6

## Summary

* switched `\dontrun` to `\donttest` in some examples

* fixed vignette to have fewer `eval=FALSE` chunks

# Moonlight2R 0.99.5

## Summary

* added several checks for correctness of main function arguments

# Moonlight2R 0.99.4

## Summary

* changed instances of `sapply` to `vapply`

# Moonlight2R 0.99.3

## Summary

* added tests with testthat

* updated following example data: dataFEA, dataGRN, dataURA, dataPRA and cscape_somatic_output

* added following example data: dataURA_plot, dataGRN_no_noise 

# Moonlight2R 0.99.2

## Summary

* added `GLS` (Gene Literature Search) function

* fixed library problems with vignettes

# Moonlight2R 0.99.1

## Summary

* removed package documentation from Rdata.R

# Moonlight2R 0.99.0

## Summary

* first release of Moonlight2R.

## New features (added or significantly changed respect to MoonlightR)

* `DMA`			Driver mutation analysis

* `plotDMA`		Creates one or more heatmap of the output from DMA

* `plotMoonlight`	Creates heatmap of Moonlight Gene Z-scores for selected genes

* `EncodePromoters`	Experimentially verified promoter sites

* `LOC_protein`		Level of consequence protein

* `LOC_translation`	Level of consequence translation

* `LOC_transcription`	Level of consequence transcription

* `NCG`			Network of Cancer Genes 7.0

* `dataPRA`		output from PRA function

* `dataMAF`		Mutation data from TCGA-LUAD

* `dataDMA`		Output from DMA function

* `cscape_somatic_output`	Cscape-somatic annotations of TCGA-LUAD 

* `DEG_Mutations_Annotations`	Differentially expressed genes's Mutations 

* `Oncogenic_mediators_mutation_summary`	Oncogenic Mediators Mutation Summary  

* `moonlight`		Function to run moonlight pipeline

