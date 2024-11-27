# Test loadMAVISp function

# reference datasets

specific_protein_ref <- as_tibble(data.frame(
  Mutation = c("R49C", "R49H", "R49P", "M51I", "M51L"),
  "Stability (FoldX5, alphafold, kcal/mol)" = c(0.49477, 0.60357, 1.84988, 0.64349, 0.01443),
  "Stability (Rosetta Cartddg2020, alphafold, kcal/mol)" = c(0.32483, 0.35782, 0.03844, 0.26667, 0.36531),
  "Stability (RaSP, alphafold, kcal/mol)" = c(6.24453672193644, 6.22116346712404, 6.24439788955118, 4.34939570830432, 4.19197330424889),
  "Stability classification, alphafold, (Rosetta, FoldX)" = c("Neutral", "Neutral", "Neutral", "Neutral", "Neutral"),
  "Stability classification, alphafold, (RaSP, FoldX)" = c("Uncertain", "Uncertain", "Uncertain", "Uncertain", "Uncertain"),
  "Local Int. (Binding with CBFB_AFmulti, heterodimer, FoldX5, kcal/mol)" = c(0.12688, -0.13046, -0.28584, 0.04618, 0.1536),
  "Local Int. (Binding with CBFB_AFmulti, heterodimer, Rosetta Talaris 2014, kcal/mol)" = c(-0.02707, 0.11159, -0.02728, -0.09484, -0.00663),
  "Local Int. classification (CBFB_AFmulti)" = c("Neutral", "Neutral", "Neutral", "Neutral", "Neutral"),
  "Relative Side Chain Solvent Accessibility in wild-type" = c(58.2, 58.2, 58.2, 21.8, 21.8),
  "gnomAD genome allele frequency" = as.numeric(c(NA, NA, NA, NA, NA)),
  "gnomAD exome allele frequency" = as.numeric(c(NA, NA, NA, NA, NA)),
  "REVEL score" = as.numeric(c(NA, NA, NA, NA, NA)),
  "Mutation sources" = c("mutations_clinvar", "mutations_clinvar", "mutations_clinvar", "mutations_clinvar", "mutations_clinvar"),
  PTMs = as.character(c(NA, NA, NA, NA, NA)),
  "is site part of phospho-SLiM" = c(FALSE, FALSE, FALSE, TRUE, TRUE),
  "PTM residue SASA (%)" = c(58.2, 58.2, 58.2, 21.8, 21.8),
  "Change in stability with PTM (FoldX5, kcal/mol)" = as.numeric(c(NA, NA, NA, NA, NA)),
  "Change in binding with PTM (FoldX5, kcal/mol)" = as.logical(c(NA, NA, NA, NA, NA)),
  "Change in binding with mutation (FoldX5, kcal/mol)" = as.logical(c(NA, NA, NA, NA, NA)),
  "PTM effect in regulation" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in stability" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in function" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "AlloSigma2 mutation type" = c("DOWN", "DOWN", "DOWN", "-", "-"),
  "AlloSigma2 predicted consequence - active sites" = c("neutral", "neutral", "neutral", "uncertain", "uncertain"),
  "AlloSigma2 predicted consequence - cofactor sites" = c("neutral", "neutral", "neutral", "uncertain", "uncertain"),
  "AlloSigma2 predicted consequence - pockets and interfaces" = c("stabilizing", "stabilizing", "stabilizing", "uncertain", "uncertain"),
  "ClinVar Variation ID" = as.character(c(576717, 1421576, 2422003, 1061802, 239046)),
  "ClinVar Interpretation" = c("Uncertain significance", "Uncertain significance", "Uncertain significance", "Uncertain significance", "Uncertain significance"),
  "ClinVar Review Status" = as.character(c(1, 1, 1, 1, 1)),
  "AlphaFold2 model pLDDT score" = c(73.63, 73.63, 73.63, 83.75, 83.75),
  "AlphaFold2 model secondary structure" = c("-", "-", "-", "H", "H"),
  "DeMaSk delta fitness" = c(-0.3905, -0.2693, -0.4559, -0.1432, -0.0892),
  "DeMaSk Shannon entropy" = c(0.2854, 0.2854, 0.2854, 1.0975, 1.0975),
  "DeMaSk log2 variant frequency" = c(-21.8538, -5.7587, -8.5659, -7.5924, -3.4632),
  "DeMaSk predicted consequence" = c("loss_of_function", "loss_of_function", "loss_of_function", "loss_of_function", "loss_of_function"),
  "GEMME Score" = c(-4.59051627158367, -2.0166488770451, -5.21762245227481, -0.430978418892934, -1.08466102827133),
  "EVE score" = c(2.56798164871032e-05, 1.09938059686617e-05, 0.0053544417229454, 2.01767420794855e-09, 2.28711094549718e-14),
  "EVE classification (25% Uncertain)" = c("Benign", "Benign", "Benign", "Benign", "Benign"),
  "AlphaMissense pathogenicity score" = c(0.5407, 0.239, 0.7822, 0.6019, 0.2257),
  "AlphaMissense classification" = c("ambiguous", "benign", "pathogenic", "pathogenic", "benign"),
  References = c("https://doi.org/10.1101/2022.10.22.513328", "https://doi.org/10.1101/2022.10.22.513328", "https://doi.org/10.1101/2022.10.22.513328", "https://doi.org/10.1101/2022.10.22.513328", "https://doi.org/10.1101/2022.10.22.513328"),
  stab_class_data_type = c("simple_mode", "simple_mode", "simple_mode", "simple_mode", "simple_mode"),
  stringsAsFactors = FALSE, check.names=FALSE
), .name_repair=function(x) {return(x)})

ensemble_ref <- as_tibble(data <- data.frame(
  "Mutation" = c("W91C", "W91R", "W91S", "P92A", "P92H"),
  "Stability (FoldX5, kcal/mol) [md]" = c(1.09645, 0.95517, 1.32043, 1.18621, 1.46241),
  "Stability (Rosetta Cartddg2020, kcal/mol) [md]" = c(1.36905, 1.37777666666667, 1.11122333333333, 0.227326666666667, 0.140363333333333),
  "Stability (RaSP, kcal/mol) [md]" = c(0.948174598045667, 0.270774913774021, 0.298933573850963, 0.233867585586445, -0.028912863696222),
  "Stability classification, (Rosetta, FoldX)1" = c("Neutral", "Neutral", "Neutral", "Neutral", "Neutral"),
  "Stability classification, (RaSP, FoldX)1" = c("Neutral", "Neutral", "Neutral", "Neutral", "Neutral"),
  "Stability (FoldX5, kcal/mol) [cabsflex]" = c(-0.00274, -0.12478, 0.01862, 1.39309, 3.21532),
  "Stability (Rosetta Cartddg2020, kcal/mol) [cabsflex]" = c(1.76043333333333, 1.43288, 1.66825666666667, 0.268706666666667, 1.03208666666667),
  "Stability (RaSP, kcal/mol) [cabsflex]" = c(1.0890734956992, 0.308181035152004, 0.509681924125949, 0.426519134963358, -0.0911043495475189),
  "Stability classification, (Rosetta, FoldX)2" = c("Neutral", "Neutral", "Neutral", "Neutral", "Uncertain"),
  "Stability classification, (RaSP, FoldX)2" = c("Neutral", "Neutral", "Neutral", "Neutral", "Uncertain"),
  "Relative Side Chain Solvent Accessibility in wild-type (average) [md]" = c(70.277, 70.277, 70.277, 50.692, 50.692),
  "Relative Side Chain Solvent Accessibility in wild-type (standard deviation) [md]" = c(22.26, 22.26, 22.26, 19.319, 19.319),
  "Relative Side Chain Solvent Accessibility in wild-type (average) [cabsflex]" = c(66.47, 66.47, 66.47, 37.285, 37.285),
  "Relative Side Chain Solvent Accessibility in wild-type (standard deviation) [cabsflex]" = c(25.088, 25.088, 25.088, 20.61, 20.61),
  "gnomAD genome allele frequency" = as.numeric(c(NA, NA, NA, NA, NA)),
  "gnomAD exome allele frequency" = as.numeric(c(NA, NA, NA, NA, NA)),
  "REVEL score" = c(NA, NA, "0.485", "0.509", NA),
  "Mutation sources" = c("somatic_muts", "clinvar", "cBioPortal", "COSMIC,somatic_muts", "COSMIC"),
  "PTMs [md]" = as.character(c(NA, NA, NA, NA, NA)),
  "is site part of phospho-SLiM [md]" = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  "PTM residue SASA (%) [md]" = c(70.277, 70.277, 70.277, 50.692, 50.692),
  "Change in stability with PTM (FoldX5, kcal/mol) [md]" = as.numeric(c(NA, NA, NA, NA, NA)),
  "Change in binding with PTM (FoldX5, kcal/mol) [md]" = c(NA, NA, NA, NA, NA),
  "Change in binding with mutation (FoldX5, kcal/mol) [md]" = c(NA, NA, NA, NA, NA),
  "PTM effect in regulation [md]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in stability [md]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in function [md]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTMs [cabsflex]" = as.character(c(NA, NA, NA, NA, NA)),
  "is site part of phospho-SLiM [cabsflex]" = c(TRUE, TRUE, TRUE, TRUE, TRUE),
  "PTM residue SASA (%) [cabsflex]" = c(66.47, 66.47, 66.47, 37.285, 37.285),
  "Change in stability with PTM (FoldX5, kcal/mol) [cabsflex]" = as.numeric(c(NA, NA, NA, NA, NA)),
  "Change in binding with PTM (FoldX5, kcal/mol) [cabsflex]" = c(NA, NA, NA, NA, NA),
  "Change in binding with mutation (FoldX5, kcal/mol) [cabsflex]" = c(NA, NA, NA, NA, NA),
  "PTM effect in regulation [cabsflex]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in stability [cabsflex]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "PTM effect in function [cabsflex]" = c("neutral", "neutral", "neutral", "neutral", "neutral"),
  "ClinVar Variation ID" = c(NA, "845798", NA, NA, NA),
  "ClinVar Interpretation" = c(NA, "Uncertain significance", NA, NA, NA),
  "ClinVar Review Status" = c(NA, "1", NA, NA, NA),
  "DeMaSk delta fitness" = c(-0.2748, -0.3228, -0.248, -0.1124, -0.2307),
  "DeMaSk Shannon entropy" = c(2.4818, 2.4818, 2.4818, 1.3832, 1.3832),
  "DeMaSk log2 variant frequency" = c(-20.417, -20.417, -2.9594, -4.9818, -20.5915),
  "DeMaSk predicted consequence" = c("loss_of_function", "loss_of_function", "loss_of_function", "loss_of_function", "loss_of_function"),
  "GEMME Score" = c(-0.527962776148948, -3.32436611299892, -1.25226234139329, -1.1402449252468, -3.34022868100625),
  "AlphaMissense pathogenicity score" = c(0.7778, 0.9346, 0.447, 0.1419, 0.3065),
  "AlphaMissense classification" = c("pathogenic", "pathogenic", "ambiguous", "benign", "benign"),
  "References" = c("10.1101/2022.10.22.513328", "10.1101/2022.10.22.513328; 37085483", "10.1101/2022.10.22.513328; 37085483", "10.1101/2022.10.22.513328; 37085483", "10.1101/2022.10.22.513328; 37085483"),
  "stab_class_ros_source" = c("ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex"),
  "stab_class_rasp_source" = c("ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex", "ensemble_mode_cabsflex"),
  stringsAsFactors = FALSE, check.names=FALSE
), .name_repair=function(x) {return(x)})

# gathering test data
mavisp_db_location <- system.file('extdata', 'mavisp_db', package='Moonlight2R')

specific_protein <- loadMAVISp(mavispDB = mavisp_db_location,
                               mode = 'simple',
                               proteins_of_interest = c('RUNX1'))

all_proteins <- loadMAVISp(mavispDB = mavisp_db_location,
                           mode = 'simple')

ensemble <- loadMAVISp(mavispDB = mavisp_db_location,
                       mode = 'ensemble',
                       ensemble = 'cabsflex')

# Test that output is as expected
test_that('Output of loadMAVISp simple mode with one protein specified is a list of one tibble', {
    expect_type(specific_protein, 'list')
    expect_length(specific_protein, 1)
    expect_true(all(sapply(specific_protein, is_tibble)))
})

test_that('Output of loadMAVISp on RUNX1 (single protein) is as expected', {
    expect_equal(specific_protein$RUNX1[1:5, ],
                 specific_protein_ref,
                 ignore_attr=TRUE)
})

test_that('Output of loadMAVISp simple mode is a list of tibbles', {
    expect_type(specific_protein, 'list')
    expect_length(all_proteins, 2)
    expect_true(all(sapply(specific_protein, is_tibble)))
})

test_that('Output of loadMAVISP on RUNX1 (all proteins) is as expected', {
    expect_equal(all_proteins$RUNX1[1:5, ],
                 specific_protein_ref,
                 ignore_attr=TRUE)
})

test_that('Output of loadMAVISp ensemble mode is a list of tibbles', {
    expect_type(ensemble, 'list')
    expect_true(all(sapply(ensemble, is_tibble)))
})

test_that('Output of loadMAVISP on TP53 (ensemble mode) is as expected', {
    expect_equal(ensemble$TP53[1:5, ],
                 ensemble_ref,
                 ignore_attr=TRUE)
})
