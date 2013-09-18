## Prepare data to plot in Euler Grid.
## Quantity or binary table with columns for samples and rows for features



#### FUNCTIONS #######


## convert a (table, matrix,) data.frame to either 0/1 or FALSE/TRUE based on whether values exceed threshold.
## to invert a threshold use: new.table <- table == FALSE
binarizeTable <- function(data, threshold=0, returnBoolean=FALSE)  {

	stopifnot(all(apply(data, 2, is.numeric)))
	if(returnBoolean)  {
		return.table <- data > threshold

	} else  {
		return.table <- apply(data, 2, FUN=function(x) ifelse(x > threshold, 1, 0))
		
	}
	return(as.data.frame(return.table))
}  


# if a binary table has multiple columns of same type these can be reduce to single columns per type following 'groupRule'
# groupRule  "any" (default), "all", "majority"
# groupDesign  a data frame with two columns "sampleName" and "groupName". All sampleName values must be unique.
# Currently, binaryTable must be a data.frame (to add columns by name).
reduceGroups <- function(binaryTable, groupDesign, groupRule="any", discardNonGroupSamples=FALSE) {
	# All sampleName values must be unique.
	stopifnot(length(unique(groupDesign[,"sampleName"])) == nrow(groupDesign))

	groups <- unique(groupDesign[,"groupName"])
	nonGroupSamples <- setdiff(colnames(binaryTable),groupDesign[,"sampleName"] )	# get list of samples NOT in groups
	#print(groups)
	#print(nonGroupSamples)	

	if(discardNonGroupSamples)  {
		reducedDat <- data.frame()
	} else {
		reducedDat <- binaryTable[,nonGroupSamples]
	}

	for(thisGroup in groups)  {
		theseSamples <- groupDesign[groupDesign[,"groupName"]==thisGroup,"sampleName"]
		sumVector <- rowSums(binaryTable[,theseSamples])
		
		reducedDat[,thisGroup] <- switch(groupRule,
					any =  ifelse(sumVector > 0 , 1, 0),
					all =  ifelse(sumVector == length(theseSamples) , 1, 0),
					majority =  ifelse(sumVector >= length(theseSamples)/2 , 1, 0)  )

	}
	
	return(reducedDat)
}



### Takes a binary (1/0) table and returns a table of unique rows with an additional column of counts for the occurrence of each row.
## Final table is sorted by count (decreasing order).
scoreCardinalities <- function(binaryTable )  {
	n.cols <- ncol(binaryTable)
	print(paste(nrow(binaryTable), "rows"))
	print(paste(n.cols, "columns. ", 2^n.cols, "possible patterns."))
	
	
	uniquePatterns <- unique(binaryTable )
	row.names(uniquePatterns) <- apply(uniquePatterns, 1, FUN=function(x) paste(x,collapse=""))
	print(paste(nrow(uniquePatterns), "patterns present"))
	
	patterns <- apply(binaryTable , 1, FUN=function(x) paste(x,collapse=""))
	patternCounts <- table(patterns)
	
	uniquePatterns[names(patternCounts),"count"] <- patternCounts
	uniquePatterns <- uniquePatterns[order(uniquePatterns$count, decreasing=T),]
	return(uniquePatterns)
}


##prettyTicks(31,n.ticks=2, inc.zero=T)
prettyTicks.old <- function(valueRange, n.ticks=3, inc.zero = FALSE)  {
	max.value <- max(valueRange)
	n.signif <- min(nchar(max.value)-1, 2)
	max.tick <- signif(max.value, n.signif)
	interval <- max.tick/n.ticks
	ticks <- seq(0, max.tick, interval)
	if(!inc.zero)	ticks <- ticks[-1]
	ticks <- signif(ticks, n.signif)
	return(ticks)
}


# returns a vector of numbers within a valueRange to use as tickmarks on an axis.
##prettyTicks(31,n.ticks=2, inc.zero=T)
# TODO improve this!
prettyTicks <- function(valueRange, n.ticks=3, inc.zero = FALSE)  {
	max.value <- max(valueRange)

	# find suitable interval
	n.signif <- 1 
	interval <- signif(max.value/n.ticks, n.signif)
	while(interval*(n.ticks-1) >= max.value) {
			n.signif <- n.signif +  1
			interval <- signif(max.value/n.ticks, n.signif)
	}
	ticks <- seq(0, by=interval, length.out=n.ticks)
	if(!inc.zero)	ticks <- ticks[-1]
	return(ticks)
}




## prototype code using built-in functions barplot() and image()
## Two separate plotting areas required and two plots would not align (something weird about 'image()')
plotEuler.old <- function(binaryGrid, counts, labels)  {
	par(mfrow=c(2,1))
	n.groups <- length(counts)
	n.samples <- ncol(binaryGrid)
	barplot(counts, space=0, beside=T)
	image(as.matrix(binaryGrid), col=c("grey","darkolivegreen4"), axes = F)
	axis(side=2, at=seq(0,1,length.out=n.samples),labels=labels, las=2)
	grid(nx=nrow(binaryGrid),ny=n.samples, col="grey20", lwd=2)
	box( col="grey20", lwd=2)

}


### main plotting function for EulerGrid
# binaryGrid	a data.frame containing only 0/1. One column per sample. 
# counts		vector of counts. Must match number of rows in binary grid (and in same order)
# labels
plotEuler <- function(binaryGrid, counts, labels, y_buffer=0.1, dropEmptySet=TRUE, dropFullSet=FALSE, dropSets='', fg.colour="darkolivegreen4",bg.colour="grey")  {

	n.samples <- ncol(binaryGrid)

	#rownames(binaryGrid) <- apply(binaryGrid)     # should already be binary chain matching the row if user used scoreCardinalities()

	# allow removal of certain sets from the table of counts (e.g. empty set , full set)
	if(dropEmptySet)  {
		dropSets <- intersect(rownames( binaryGrid),unique(c(dropSets, paste(rep(0,n.samples),collapse=""))))
	}
	if(dropFullSet)  {
		dropSets <- intersect(rownames( binaryGrid),unique(c(dropSets, paste(rep(1,n.samples),collapse=""))))
	}
	print(paste("Dropping:", paste(dropSets, collapse=",")))
	
	keepSet <- setdiff(rownames(binaryGrid), dropSets)
	keepSetIndex <- match(keepSet , rownames(binaryGrid))
	
	# remove unwanted sets.
	binaryGrid <- binaryGrid[keepSetIndex,]
	counts <- counts[keepSetIndex]
	n.counts <- length(counts)
	
	bar.bottom <- 1 + y_buffer
	max.count <- max(counts)

	# all rects in grid, specified by rows from bottom, moving left to right
	grid.x1 <- rep(seq(0,1,length.out=n.counts+1)[-(n.counts + 1)] , n.samples)
	grid.x2 <- rep(seq(0,1,length.out=n.counts+1)[-1] , n.samples)
	grid.y1 <- rep(seq(0,1,length.out=n.samples+1)[-(n.samples + 1)] , each=n.counts)
	grid.y2 <- rep(seq(0,1,length.out=n.samples+1)[-1] ,each= n.counts)

	colVector <- unlist(binaryGrid)   # concatenation by colums - to be used as rows from bottom, left to right.
	colVector <- ifelse(colVector == 1,fg.colour, bg.colour )

	# begin plotting
	plot.new()
	plot.window(xlim=c(0,1),ylim=c(0,1+bar.bottom))
	# draw the grid and add labels
	rect(grid.x1, grid.y1, grid.x2, grid.y2 , col=colVector)
	labelPosVector <- (seq(0,1,length.out=n.samples+1)[-(n.samples + 1)])  + (seq(0,1,length.out=n.samples+1)[2] /2)
	mtext(colnames(binaryGrid), side=2, at=labelPosVector, las=2)

	# draw the bargraph and add an axis
	rect(seq(0,1,length.out=n.counts+1)[-(n.counts+1)],  bar.bottom , seq(0,1,length.out=n.counts+1)[-1], (counts/max.count) + bar.bottom , col="grey")
	tickVector <- prettyTicks(range(counts), n.ticks=4,inc.zero=T)
	tickPosVector <- (tickVector/max.count) + bar.bottom
	axis(at=tickPosVector , side=2, labels =tickVector, las=2)

}





