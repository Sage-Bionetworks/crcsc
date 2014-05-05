library(rGithubClient)
crcRepo <- getRepo("sib-bcf/crcsc")

work.dir <- "/export/scratch/paolo/SageCRCsubtyping/analysis/work/test"
setwd(work.dir)

# get helper functions
sourceRepoFile(crcRepo, "groups/A/pipeline/synapseHelper.R")

# Load datasets from synapse
sourceRepoFile(crcRepo, "groups/A/datasets_merge/get_datasets_from_synapse.R")

# Annotate datasets and select reference genes by mad
sourceRepoFile(crcRepo, "groups/A/datasets_merge/SageDatasetAnnotator.R")

# Select probes by computing cor(cor)
sourceRepoFile(crcRepo, "groups/A/datasets_merge/probes.cor.selection.R")

# normalize datasets with ComBat using MSI-like signature as covariate
sourceRepoFile(crcRepo, "groups/A/datasets_merge/normalize.datasets.R")

# filter genes by quantile range (0.05, 0.95), and correlations
sourceRepoFile(crcRepo, "groups/A/datasets_merge/gene.filtering.R")

# merge datasets and write the merged matrix
sourceRepoFile(crcRepo, "groups/A/datasets_merge/simple.dataset.merger.R")