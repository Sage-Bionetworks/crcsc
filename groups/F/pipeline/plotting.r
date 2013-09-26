# Functions for plotting subtyping results.
# 
# Author: Andreas Schlicker
###############################################################################

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
	data.frame(margin=apply(matrix, 1, function(x) { sort(x, decreasing=TRUE)[1] - sort(x, decreasing=TRUE)[2] }),
			   Subtype=apply(matrix, 1, function(x) { names(x)[which(x == max(x))][1] }),
			   sample.name=rownames(matrix))
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

		if (n.temp.samp > 0){
			# Sort according to decreasing probability and get the sample names
			sample.order[l:(l+n.temp.samp - 1)] = names(matrix[temp.samp, k])[order(matrix[temp.samp, k], decreasing=TRUE)]
			l = l + n.temp.samp
		} else if (n.temp.samp == 1) {
			# If we only have one sample, R doesn't preserve the names
			sample.order[l:(l+n.temp.samp - 1)] = temp.samp
			l = l + n.temp.samp
		}
	}
	
	sample.order
}

##' Plots the probabilities for each sample and subtype as stacked bar plot.
##' @param matrix the subtyping matrix with samples in rows and subtypes
##' in columns
##' @param colors vector with colors to use; default: NULL
##' @return the ggplot2 object
##' @author Andreas Schlicker
probabilityPlot = function(matrix, colors=NULL) {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2PlotDf(matrix)
	sample.order = orderSamples(matrix)
	
	p = ggplot(plotting.df, aes(x=sample.name, y=pvalue, fill=Subtype)) +
		geom_bar(stat="identity") + 
		scale_x_discrete(limits=sample.order) +
		xlab("Samples") + 
		ylab("Subtyping probability") +
		theme(axis.ticks=element_blank(), 
			  axis.text.x=element_blank(), 
			  axis.title.x=element_text(size=20, face="bold"),
			  axis.text.y=element_text(size=16),
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
##' @return the ggplot2 object
##' @author Andreas Schlicker
facetedProbabilityPlot = function(matrix, ncol=1, colors=NULL) {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2PlotDf(matrix)
	sample.order = orderSamples(matrix)
	
	p = ggplot(plotting.df, aes(x=sample.name, y=pvalue, fill=Subtype)) +
		geom_bar(stat="identity", position=position_dodge()) +
		facet_wrap(~ Subtype, ncol=ncol) + 
		guides(fill=FALSE) +
		scale_x_discrete(limits=sample.order) +
		xlab("Samples") + 
		ylab("Subtyping probability") +
		theme(axis.ticks=element_blank(), 
			  axis.text.x=element_blank(), 
			  axis.title.x=element_text(size=20, face="bold"),
			  axis.text.y=element_text(size=16),
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
##' @return the ggplot2 object
##' @author Andreas Schlicker
marginPlot = function(matrix, colors=NULL, shapes=NULL) {
	require(ggplot2) || stop("Can't load package \"ggplot2\".")
	
	plotting.df = pMat2MarginDf(matrix)
	sample.order = orderSamples(matrix)
	
	p = ggplot(plotting.df, aes(x=sample.name, y=margin, color=Subtype, group=Subtype, shape=Subtype)) +
		geom_point(size=5) + 
		scale_x_discrete(limits=sample.order) +
		xlab("Samples") + 
		ylab("Subtyping margin") +
		theme(axis.ticks=element_blank(), 
				axis.text.x=element_blank(), 
				axis.title.x=element_text(size=20, face="bold"),
				axis.text.y=element_text(size=16),
				axis.title.y=element_text(size=20, face="bold"))

	if (!is.null(colors)) { 
		p = p + scale_color_manual(values=cols)
	}
	if (!is.null(shapes)) { 
		p = p + scale_shape_manual(values=shapes)
	}
	
	p
}

##' Generate cooccurrence heatmap showing how often two samples end up in the
##' same subtype. Rows and columns are clustered using Ward's method.
##' @param matrix cooccurrence matrix. 
##' @param normalize boolean; if TRUE, normalize each row by the diagonal element
##' @return heatmap
##' @author Andreas Schlicker
cooccurrencePlot = function(matrix, normalize=FALSE) {
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
