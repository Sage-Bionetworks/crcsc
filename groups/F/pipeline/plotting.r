# Functions for plotting subtyping results.
# 
# Author: Andreas Schlicker
###############################################################################

##' Calculates the margin of a vector. The margin is defined as 
##' largest element - second largest element. If the vector contains
##' only one element, this is returned.
##' @param vec the data vector
##' @param na.rm ignored, for compatibility reasons
##' @return the margin
margin = function(vec, na.rm=TRUE) {
	x = sort(vec[!is.na(vec)], decreasing=TRUE)
	
	if (length(x) < 2) {
		return(x[1])
	}
	
	x[1] - x[2]
}

##' Converts a subtyping matrix into the data frame used for plotting
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @return the plotting data frame
##' @author Andreas Schlicker
pMat2PlotDf = function(matrix) {
	data.frame(pvalue=as.vector(matrix), # Concatenate the probability by subtype
			   Subtype=rep(colnames(matrix), each=nrow(matrix)), # Add the subtype label
			   sample.name=rep(rownames(matrix), times=ncol(matrix))) # Add the sample label
}

##' Converts a subtyping matrix into a maring data frame. The margin
##' is defined as probability difference between the two most likely
##' subtypes. 
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @return the plotting data frame
##' @author Andreas Schlicker
pMat2MarginDf = function(matrix) {
	data.frame(margin=apply(matrix, 1, function(x) { margin(x) }),
			   Subtype=apply(matrix, 1, function(x) { names(x)[which(x == max(x))][1] }),
			   sample.name=rownames(matrix))
}

##' Convert a subtyping matrix into a data frame for a forest plot.
##' @param matrix the subtyping matrix
##' @param name the data set name
##' @param statistic name of a function used to calculate the statistic
##' @return the plotting data frame
##' @author Andreas Schlicker
pMat2ForestDf = function(matrix, name, statistic="mean") {
	# Assign samples to subtypes
	subtypes = apply(matrix, 1, function(x) { names(x)[which(x == max(x))][1] })
	
	# Which samples do or do not belong to each subtype?
	sample2subtype = list()
	for (n in colnames(matrix)) {
		sample2subtype[[n]] = names(subtypes)[subtypes == n]
		# If we're calculating the margin, we're not interested in samples that are not 
		# assigned to a subtype
		if (statistic != "margin") {
			sample2subtype[[paste("not.", n, sep="")]] = names(subtypes)[subtypes != n]
		}
	}
	
	# Calculate statistics
	means = double(length(sample2subtype))
	stddev = double(length(sample2subtype))
	for (i in 1:length(means)) {
		if (statistic == "margin") {
			# The margin is a per sample statistic, so get the mean value for a subtype
			if (length(sample2subtype[[i]]) == 1) {
				means[i] = get(statistic)(matrix[sample2subtype[[i]], ], na.rm=TRUE)
			} else if (length(sample2subtype[[i]]) > 1) {
				means[i] = mean(apply(matrix[sample2subtype[[i]], ], 1, get(statistic), na.rm=TRUE))
			}
		} else {
			means[i] = get(statistic)(matrix[sample2subtype[[i]], ceiling(i/2)], na.rm=TRUE)
		}
		names(means)[i] = names(sample2subtype)[i]
		
		#stddev[i] = means[i] / nrow(matrix)
		stddev[i] = sd(matrix[sample2subtype[[i]], ceiling(i/2)], na.rm=TRUE)
		names(stddev)[i] = names(sample2subtype)[i]
	}
	
	# Combine
	df = data.frame(Dataset=name,
			 	   	Subtype=names(means),
			   		statistic=means,
			   		ymin=sapply(names(means), function(x) { max(0, means[x] - stddev[x]) }),
			   		ymax=sapply(names(means), function(x) { min(1, means[x] + stddev[x]) }),
			   		Dataset.size=nrow(matrix))
	df[, "Subtype"] = factor(df[, "Subtype"], levels=names(sample2subtype))
	
	df
}

##' Create a sample ordering. Each sample is assigned to the subtype with the
##' highest probability. Within each subtype, the samples are ordered by decreasing
##' probability for this subtype.
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @return vector with ordered sample names
##' @author Andreas Schlicker
orderSamples = function(matrix) {
	# Get the subtype with the highest probability for each sample
	sample.subtype = apply(matrix, 1, function(x) { names(x)[which(x == max(x))][1] })
	
	# Final ordered vector of sample names
	sample.order = character(length(sample.subtype))
	# Current index we're at
	l = 1
	# Going through the subtypes
	for (k in colnames(matrix)) {
		# Samples of this subtype
		temp.samp = names(sample.subtype)[sample.subtype == k]
		n.temp.samp = length(temp.samp)

		if (n.temp.samp > 1){
			# Sort according to decreasing probability and get the sample names
			sample.order[l:(l+n.temp.samp - 1)] = names(matrix[temp.samp, k])[order(matrix[temp.samp, k], decreasing=TRUE)]
			names(sample.order)[l:(l+n.temp.samp - 1)] = k
			l = l + n.temp.samp
		} else if (n.temp.samp == 1) {
			# If we only have one sample, R doesn't preserve the names
			sample.order[l:(l+n.temp.samp - 1)] = temp.samp
			names(sample.order)[l:(l+n.temp.samp - 1)] = k
			l = l + n.temp.samp
		}
	}
	
	sample.order
}

##' Plots the probabilities for each sample and subtype as stacked bar plot.
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @param colors vector with colors to use; default: NULL
##' @param line.color color for the line dividing subtypes; default: black
##' @return the ggplot2 object
##' @author Andreas Schlicker
probabilityPlot = function(matrix, colors=NULL, line.color="red") {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2PlotDf(matrix)
	sample.order = orderSamples(matrix)
	subtype.boundaries = cumsum(table(names(sample.order))) + 0.5
	
	p = ggplot(plotting.df, aes(x=sample.name, y=pvalue, fill=Subtype)) +
		geom_bar(stat="identity") + 
		geom_vline(xintercept=subtype.boundaries[-1*(length(subtype.boundaries))], linetype="longdash", color=line.color) + 
		scale_x_discrete(limits=sample.order) +
		xlab("Samples") + 
		ylab("Subtyping probability") +
		theme(axis.ticks=element_blank(), 
			  axis.text.x=element_blank(), 
			  axis.title.x=element_text(size=20, face="bold"),
			  axis.text.y=element_text(size=20, face="bold"),
			  axis.title.y=element_text(size=20, face="bold"))

	if (!is.null(colors)) {
		p = p + scale_fill_manual(values=colors)
	}
	
	p
}

##' Plots the probabilities for the subtyping. The probabilities of each
##' sample to belong to a subtype are visualized in a separate facet.
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @param ncol number of facets to plot per column, default: 1
##' @param colors vector with colors to use; default: NULL
##' @param line.color color for the line dividing subtypes; default: black
##' @return the ggplot2 object
##' @author Andreas Schlicker
facetedProbabilityPlot = function(matrix, ncol=1, colors=NULL, line.color="red") {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2PlotDf(matrix)
	sample.order = orderSamples(matrix)
	subtype.boundaries = cumsum(table(names(sample.order))) + 0.5
	
	p = ggplot(plotting.df, aes(x=sample.name, y=pvalue, fill=Subtype)) +
		geom_bar(stat="identity", position=position_dodge()) +
		facet_wrap(~ Subtype, ncol=ncol) + 
		geom_vline(xintercept=subtype.boundaries[-1*(length(subtype.boundaries))], linetype="longdash", color=line.color) +
		guides(fill=FALSE) +
		scale_x_discrete(limits=sample.order) +
		scale_y_continuous(limits = c(0, 1)) +
		xlab("Samples") + 
		ylab("Subtyping probability") +
		theme(axis.ticks=element_blank(), 
			  axis.text.x=element_blank(), 
			  axis.title.x=element_text(size=20, face="bold"),
			  axis.text.y=element_text(size=20, face="bold"),
			  axis.title.y=element_text(size=20, face="bold"),
			  strip.text=element_text(size=20, face="bold"))
	
	if (!is.null(colors)) {
		p = p + scale_fill_manual(values=colors)
	}
	
	p
}

##' Visualize the margin between the two subtypes with the highest probabilities.
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @param colors vector with colors to use; default: NULL
##' @param shapes vector with shapes to use; default: NULL
##' @param line.color color for the line dividing subtypes; default: black
##' @return the ggplot2 object
##' @author Andreas Schlicker
marginPlot = function(matrix, colors=NULL, shapes=NULL, line.color="red") {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2MarginDf(matrix)
	sample.order = orderSamples(matrix)
	subtype.boundaries = cumsum(table(names(sample.order))) + 0.5
	
	p = ggplot(plotting.df, aes(x=sample.name, y=margin, color=Subtype, group=Subtype, shape=Subtype)) +
		geom_point(size=5) + 
		scale_x_discrete(limits=sample.order) +
		scale_y_continuous(limits = c(0, 1)) + 
		geom_vline(xintercept=subtype.boundaries[-1*(length(subtype.boundaries))], linetype="longdash", color=line.color) +
		xlab("Samples") + 
		ylab("Subtyping margin") +
		theme(axis.ticks=element_blank(), 
				axis.text.x=element_blank(), 
				axis.title.x=element_text(size=20, face="bold"),
				axis.text.y=element_text(size=20, face="bold"),
				axis.title.y=element_text(size=20, face="bold"))

	if (!is.null(colors)) { 
		p = p + scale_color_manual(values=colors)
	}
	if (!is.null(shapes)) { 
		p = p + scale_shape_manual(values=shapes)
	}
	
	p
}

##' Generate coclustering heatmap showing how often two samples end up in the
##' same subtype. Rows and columns are clustered using Ward's method.
##' @param matrix cooccurrence matrix. 
##' @param normalize boolean; if TRUE, normalize each row by the diagonal element
##' @return heatmap
##' @author Andreas Schlicker
coclusteringPlot = function(matrix, normalize=FALSE) {
	require(gplots) || stop("Can't load package \"gplots\"")
	
	if (normalize) {
		matrix = matrix / diag(matrix)
	}
	
	heatmap.2(matrix, scale="none", trace="none",
			  Colv=as.dendrogram(hclust(as.dist(1-matrix), method="ward")),
			  Rowv=as.dendrogram(hclust(as.dist(1-matrix), method="ward")),
			  col=colorpanel(49, low="white", high="royalblue4"),
			  breaks=seq(0, 1, length.out=50),
			  density.info="none",
			  labCol="",
			  labRow=""
	)
}

##' Generates a forest plot for the subtyping results.
##' For each data set and subtype, the plot contains two entries, one for
##' all samples assigned to that subtype and one for all samples not assigned
##' to the subtype. 
##' @param results named list with all results from the subtyping
##' @param colors vector with colors to use; default: NULL
##' @param xaxis either "dataset" or "statistic", denoting which is plotted on the xaxis
##' @param statistic "mean", "median" or "margin"
##' @param combine boolean indicating whether values for sample belonging to and not 
##' belonging to a subtype should be combined in one facet
##' @return the ggplot object
##' @author Andreas Schlicker
forestPlot = function(results, colors=NULL, xaxis=c("dataset", "statistic"), statistic=c("mean", "median", "margin"), combine=TRUE) {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	xaxis = match.arg(xaxis)
	statistic = match.arg(statistic)

	plotting.df = do.call("rbind", lapply(names(results), function(x) { pMat2ForestDf(results[[x]][[1]], x, statistic) }))
	
	# Add indicator for plotting symbol
	in.subtype = rep("y", times=nrow(plotting.df))
	in.subtype[grep("not", plotting.df[, "Subtype"])] = "n"
	plotting.df = cbind(plotting.df, in.subtype=in.subtype)
	
	# Reorder data sets from largest to smallest
	plotting.df[, "Dataset"] = factor(plotting.df[, "Dataset"], 
									  levels=as.character(unique(plotting.df[unique(order(plotting.df[, "Dataset.size"])), "Dataset"])))
	
	if (combine) {
		plotting.df[, "Subtype"] = gsub("not.", "", plotting.df[, "Subtype"])
	}
	
	total = sum(unlist(lapply(results, function(x) { nrow(x[[1]]) })))
	plotting.df[, "Dataset.size"] = sapply(plotting.df[, "Dataset.size"], function(x) { (x / total) } )
	plotting.df[, "Dataset.size"] = plotting.df[, "Dataset.size"] / max(plotting.df[, "Dataset.size"]) * 3 + 4.5
	
	p = ggplot(plotting.df, aes(x=Dataset, y=statistic, ymin=ymin, ymax=ymax, color=Subtype, shape=in.subtype)) +
		geom_point(aes(size=Dataset.size)) + 
		geom_linerange(aes(size=1.5)) +
		guides(color=FALSE, size=FALSE, shape=guide_legend(override.aes=list(size=5))) +
		facet_grid(. ~ Subtype) +
		xlab("Data set") + 
		ylab(statistic) + 
		theme(axis.text.x=element_blank(),
			  axis.title.x=element_text(size=20, face="bold"),
			  axis.text.y=element_text(size=20, face="bold"),
			  axis.title.y=element_text(size=20, face="bold"),
			  strip.text=element_text(size=20, face="bold"),
			  legend.title=element_text(size=20, face="bold"),
			  legend.text=element_text(size=20, face="bold"))
	
	if (statistic == "margin") {
		p = p + guides(shape=FALSE)
	} else {
		p = p + scale_shape_manual(name="Assigned to\nsubtype",
						     	   breaks=c("n", "y"),
						   		   values=c(15, 16),
								   labels=c("no", "yes"))
	}

	if (xaxis == "statistic") {
		p = p + coord_flip() + theme(axis.text.x=element_text(size=20, face="bold"))
	} else {
		p = p + theme(axis.ticks.x=element_blank())
	}

	if (!is.null(colors)) {
		notcols = colors
		names(notcols) = paste("not.", names(notcols), sep="")
		p = p + scale_color_manual(values=c(colors, notcols))
	}
	
	p
}

##' Create the quality control heatmap. Column and row colors indicate cluster
##' membership for samples and features, respectively. 
##' @param exprs expression matrix with samples in columns and features in rows
##' @param clustering named list with sample to cluster assignment
##' @param signatures named list with feature to signature assignment
##' @return heatmap.2 for plotting
##' @author Andreas Schlicker
createHeatmap = function(exprs, clustering, signatures) {
	require(gplots) || stop("I need package \"gplots\" for doing this.")
	
	samples = c()
	for (n in sort(names(clustering))) {
		samples = c(samples, clustering[[n]])
	}
	
	features = c()
	for (n in sort(names(signatures))) {
		features = c(features, signatures[[n]])
	}
	
	bounds = quantile(exprs[unlist(signatures), unlist(clustering)], probs=c(0.05, 0.95), na.rm=TRUE)
	csd = c()
	for (i in 1:length(clustering)) {
		csd = c(csd, rep(palette()[i %% 8], times=length(clustering[[sort(names(clustering))[i]]])))
	}
	rsd = c()
	for (i in 1:length(signatures)) {
		rsd = c(rsd, rep(palette()[i %% 8], times=length(signatures[[sort(names(signatures))[i]]])))
	}
	heatmap.2(exprs[features, samples],
			trace="none",
			scale="none",
			col=colorpanel(49, low="blue", high="yellow"),
			breaks=seq(bounds[1], bounds[2], length.out=50),
			ColSideColors=csd,
			RowSideColors=rsd,
			Colv=NA,
			Rowv=NA,
			dendrogram="none")
}

