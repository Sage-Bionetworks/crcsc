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


