All code for running the iNMF subtyping ([full text](http://dx.doi.org/10.1186/1755-8794-5-66)) on data sets in the Synapse repository of the CRC Subtyping Consortium. 

synapse_helper.r
-----------------------
Helper functions for downloading data from Synapse and mapping of IDs

inmf.r
--------
Contains all functions necessary for applying the iNMF subtyping on any data set. The iNMF signatures can be obtained from the Synapse repository or the supplementary material published along with the article. 

1. methods for mapping signatures to data set-specific IDs

2. methods for running the hierarchical iNMF subtyping in a fully automated fashion

3. methods for running a bootstrapped version of the subtyping; bootstraps are run by subsampling samples in the data set and give an estimation of the subtyping probability for each sample

subtyping.r
---------------
Script that automatically downloads all data sets, subtypes them and stores the sample/subtype probability matrix in Synapse. The script also takes care of provenance recording.

plotting.r
-----------
Some functions for creating diagnostic plots for the subtyping.

1. probabilityPlot: Visualizes the probability of assigning a sample to a specific subtype. Samples (x-axis) are assigned to the subtype with the highest probability and are ordered with descending probability within each subtype.

2. facetedProbabilityPlot: Same information as probabilityPlot but plots each subtype in a separate facet

3. marginPlot: Margin of assignment to a subtype. The margin is defined as difference between the two largest probabilities. The larger the margin, the higher the subtyping certainty. Samples are ordered as in the probability plots.

4. coclustering: Coclustering heatmap. Visualizes how often two samples were assigned to the same subtype.

5. forestPlot: For each data set, visualizes the mean or median probability with which samples are assigned to each subtype. Additionally, shows the mean or median probability for each subtype for all samples that were not assigned to this subtype. Error bars represent standard deviations. The third option is to show the mean margin for samples assigned to each subtype for each data set. Data sets are ordered by descending number of samples, which is also reflected in the size of the plotting symbol.

6. heatmap: Heatmap showing expression values of signature genes. Subtype assignment is depicted using colored bars.

plot_subtype_results.r:
------------------------------
Generate all plots (except for heatmaps) for all results of subtyping.r. 
