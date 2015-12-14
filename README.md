# RettProject

1. Download this repository and download the "Data", "Annot" and "Results" folders from the Voineagu Lab website. 
2. Save all folders in the same directory. 
3. The codes in this repository are now executable. 

## geneCounts.R ######################################
This script generates gene-level quantification of RNASeq reads. 

## geneCounts.pbs ####################################
Job script if one wants to run geneCounts.R on HPC.

## Deconvolution.R ################################
This file uses deconvolution methods to estimate the proportions of neurons in the samples.

## DEAnalysis.R ##################################
This file performs differential expression analysis on RNASeq data.

## validation.R #################################
This file uses micro-array data to calculate the validation rate of the RNA-Seq differentially expressed genes - 74%. 

## ontology.R ##################################
Gene ontology analysis using goseq.
