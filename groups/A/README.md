This folder contains all the code for running the module-based LDA classifier on 
data sets in the Synapse repository of the CRC Subtyping Consortium. The gene modules 
and the trained classifier are available from Synapse.

## groupA_subtyping_pipeline.R
---
The subtyping pipeline. This script automatically downloads all the files from Synapse, 
performs the subtyping on each one of them and uploads the results to Synapse. 

## synapseHelper.R
---
Help functions for downloading data from Synapse and mapping probe set IDs or 
gene symbols to Entrez Gene IDs. 

## moduleScores.R
---
The function for computing the module (meta-gene) scores for the samples in a given 
data set.

## MetaGenesSubtyping.R
---
The function for training or applying a module-based LDA classifier. 

## run_gene_filtering_pipeline.R
---
The script which runs the probe selection/gene filtering/dataset normalization/merging pipeline

## get_datasets_from_synapse.R
---
script to load datasets from synapse

# SageDatasetAnnotator.R
---
script to annotate datasets and select reference genes by mad

# probes.cor.selection.R
script to select probes by computing cor(cor)

# normalize.datasets.R
Script to normalize datasets with ComBat using MSI-like signature as covariate

# gene.filtering.R
Script to filter genes by quantile range (0.05, 0.95), and correlations

# simple.dataset.merger.R
Script to merge datasets and write the merged matrix

