# Generate plots for subtyping results. Subtyping has to be run first using
# "subtyping.r".
# 
# Author: schlandi
###############################################################################

library(rGithubClient)

# GitHib repository
crcRepo = getRepo("andreas-schlicker/crcsc")

# Plotting functions
sourceRepoFile(crcRepo, "groups/F/pipeline/plotting.r")
plottingFunctions = getPermlink(crcRepo, "groups/F/pipeline/plotting.r")

# Plotting of all results
plotting.cols = c("gray10", "gray40", "gray70", "springgreen4", "springgreen2")
names(plotting.cols) = c("1.1", "1.2", "1.3", "2.1", "2.2")
plotting.shapes = 15:19
for (x in names(allResults)) {
	png(paste(x, "probability.png", sep="_"), width=6000, height=3000, res=300)
	print(probabilityPlot(allResults[[x]][[1]], colors=plotting.cols))
	dev.off()
	
	png(paste(x, "probability_faceted.png", sep="_"), width=6000, height=3000, res=300)
	print(facetedProbabilityPlot(allResults[[x]][[1]], colors=plotting.cols))
	dev.off()
	
	png(paste(x, "margins_margin.png", sep="_"), width=6000, height=3000, res=300)
	print(marginPlot(allResults[[x]][[1]], colors=plotting.cols, shapes=plotting.shapes))
	dev.off()
	
	png(paste(x, "cooccurrence.png", sep="_"), width=3000, height=3000, res=300)
	coclusteringPlot(allResults[[x]][[2]])
	dev.off()
}

png("forest_plot_mean.png", width=10000, height=8000, res=300)
forestPlot(allResults, colors=plotting.cols, xaxis="statistic", statistic="mean")
dev.off()

png("forest_plot_median.png", width=10000, height=8000, res=300)
forestPlot(allResults, colors=plotting.cols, xaxis="statistic", statistic="median")
dev.off()

png("forest_plot_margin.png", width=10000, height=8000, res=300)
forestPlot(allResults, colors=plotting.cols, xaxis="statistic", statistic="margin")
dev.off()
