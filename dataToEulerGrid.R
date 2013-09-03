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
plotEuler <- function(binaryGrid, counts, labels, y_buffer=0.1, fg.colour="darkolivegreen4",bg.colour="grey")  {

	n.counts <- length(counts)
	n.samples <- ncol(binaryGrid)
	
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





stopifnot(FALSE)
#####EXAMPLES ########

rawData <- data.frame(	sample.A.1=c(0.2, 0, 0, 0.4, 0.3, 0,0,0.8,0.6,0.1), 
				sample.A.2=c(0.1, 0, 0, 0.5, 0.3, 0,0,0.7,0.9,0.2),
				sample.B.1=c(0, 0.8, 0, 0.5, 0.3, 0,0,0.8,0.3,0), 
				sample.B.2=c(0, 0.9, 0, 0.6, 0.3, 0,0.1,0.9,0.8,0),
				sample.C.1=c(0.1, 0, 0, 0.2, 0.3, 0,0,0.7,0.2,0), 
				sample.C.2=c(0, 0, 0, 0.3, 0.3, 0,0,0.8,0.5,0)	)

rawDataLong <- rbind(rawData,rawData,rawData)




my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B"))

my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2", "sample.C.1", "sample.C.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B", "group.C", "group.C"))

head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign))
head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign,groupRule="all"))

scoreCardinalities(binarizeTable(rawDataLong))
rowCounts <-   scoreCardinalities(binarizeTable(rawDataLong))


plotEuler(rowCounts[,1:6], rowCounts$count, names(rowCounts[,1:6]))




stopifnot(FALSE)
############ DEVELOPMENT ############


## Prepare data to plot in Euler Grid.
## Quantity or binary table with columns for samples and rows for features


#  If non-binary, set threshold to count detection of feature in a sample
# 	Hence, all data should be pre-normalised to use same threshold

# Given that with 5 distinct samples or groups, there will be 32 combinations and with 6, 64, 
#	Would be good to set a plot limit and only group a mx number of proportion of samples 
#	Proportin could also be determined from table data (e.g. 90% of binarized data probbaly in 30% of columns).
# Give text output on %of combinations shown and 

## How to deal with replicates or sample groups?  "sll", "majority", "any"




rawData <- data.frame(	sample.A.1=c(0.2, 0, 0, 0.4, 0.3, 0,0,0.8,0.6,0.1), 
				sample.A.2=c(0.1, 0, 0, 0.5, 0.3, 0,0,0.7,0.9,0.2),
				sample.B.1=c(0, 0.8, 0, 0.5, 0.3, 0,0,0.8,0.3,0), 
				sample.B.2=c(0, 0.9, 0, 0.6, 0.3, 0,0.1,0.9,0.8,0),
				sample.C.1=c(0.1, 0, 0, 0.2, 0.3, 0,0,0.7,0.2,0), 
				sample.C.2=c(0, 0, 0, 0.3, 0.3, 0,0,0.8,0.5,0)

		)

rawDataLong <- rbind(rawData,rawData,rawData)

table(rawDataLong)



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

binarizeTable(rawDataLong)
binarizeTable(rawDataLong,threshold=0.1)
binarizeTable(rawDataLong, returnBoolean=TRUE)



# if a table has multiple columns of same type these can be reduce to a single table following 'groupRule'
# groupRule  "any" (default), "all", ("majority")
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
					all =  ifelse(sumVector >= length(theseSamples)/2 , 1, 0)  )

	}
	
	return(reducedDat)
}

my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B"))

my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2", "sample.C.1", "sample.C.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B", "group.C", "group.C"))

head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign))
head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign,groupRule="all"))

head(binarizeTable(rawDataLong))



### Takes a binary (1/0) table and returns a table of unique rows with a column of counts for the occurrence of each row.
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


# unique rows reveal all the patterns present

scoreCardinalities(binarizeTable(rawDataLong))
rowCounts <-   scoreCardinalities(binarizeTable(rawDataLong))


par(mfrow=c(2,1), mar=c(2,7,2,2))
barplot(rowCounts$count)
image(as.matrix(rowCounts[,1:6]), col=c("grey","darkolivegreen4"), axes = F)
axis(side=2, at=seq(0,1,length.out=6),labels=names(rowCounts[,1:6]), las=2)
grid(6, col="black", lwd=2)


# or as a function (with call below)

plotEuler <- function(binaryGrid, counts, labels)  {
	par(mfrow=c(2,1), mar=c(2,7,2,2))
	n.groups <- length(counts)
	n.samples <- ncol(binaryGrid)
	barplot(counts, space=0)
	image(as.matrix(binaryGrid), col=c("grey","darkolivegreen4"), axes = F)
	axis(side=2, at=seq(0,1,length.out=n.samples),labels=labels, las=2)
	grid(n.samples, col="black", lwd=2)

}




plotEuler(rowCounts[,1:6], rowCounts$count, names(rowCounts[,1:6]))



# found that barplot was not perfectly aligning with grid when using large number of combinations.
# redraw using rect? 
#rect(xleft, ybottom, xright, ytop, density = NULL, angle = 45,
 #    col = NA, border = NULL, lty = par("lty"), lwd = par("lwd"),
#     ...)

counts <- c(30, 12, 7, 6, 3, 3,3,3,2,2,2,1,1,1,1)


n.counts <- length(counts)
max.count <- max(counts)
plot.new()
rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)
image(as.matrix(rowCounts.plot), col=c("grey","darkolivegreen4"), axes = F)

# heatmap creates a mess
plot.new()
rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)
heatmap(as.matrix(rowCounts.plot), col=c("grey","darkolivegreen4"))


# shows that the barplot region is inside the image region. Don't know why.
plot.new()
image(as.matrix(rowCounts.plot), col=c("grey","darkolivegreen4"), axes = F, add=T)
rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)


############ USING rect() for both parts (on same plot) allows perfect alignment.


prettyTicks <- function(valueRange, n.ticks=3, inc.zero = FALSE)  {
	max.value <- max(valueRange)
	n.signif <- min(nchar(max.value)-1, 2)
	max.tick <- signif(max.value, n.signif)
	#return(max.tick)
	
	interval <- max.tick/n.ticks

	
	ticks <- seq(0, max.tick, interval)
	if(!inc.zero)	ticks <- ticks[-1]
	return(ticks)
}

prettyTicks(31)

prettyTicks(31,n.ticks=2, inc.zero=T)


y_buffer <- 0.1
fg.colour <- "darkolivegreen4"
bg.colour <- "grey"
rowCounts.plot <- rowCounts.plot[,-13]
n.samples <- ncol(rowCounts.plot)
bar.bottom <- 1 + y_buffer

par(mfrow=c(1,1))
plot.new()
plot.window(xlim=c(0,1),ylim=c(0,1+bar.bottom))
#rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)

# all rects in grid, specified by rows from bottom, moving left to right
grid.x1 <- rep(seq(0,1,length.out=n.counts+1)[-(n.counts + 1)] , n.samples)
grid.x2 <- rep(seq(0,1,length.out=n.counts+1)[-1] , n.samples)
grid.y1 <- rep(seq(0,1,length.out=n.samples+1)[-(n.samples + 1)] , each=n.counts)
grid.y2 <- rep(seq(0,1,length.out=n.samples+1)[-1] ,each= n.counts)

colVector <- unlist(rowCounts.plot)   # concatenation by colums - to be used as rows from bottom, left to right.
#colVector <- as.character(ifelse(colVector == 1,fg.colour, bg.colour ))
colVector <- ifelse(colVector == 1,fg.colour, bg.colour )
rect(grid.x1, grid.y1, grid.x2, grid.y2 , col=colVector)
labelPosVector <- (seq(0,1,length.out=n.samples+1)[-(n.samples + 1)])  + (seq(0,1,length.out=n.samples+1)[2] /2)
mtext(colnames(rowCounts.plot), side=2, at=labelPosVector, las=2)
rect(seq(0,1,length.out=n.counts+1)[-(n.counts+1)],  bar.bottom , seq(0,1,length.out=n.counts+1)[-1], (counts/max.count) + bar.bottom , col="grey")
tickVector <- prettyTicks(range(counts), inc.zero=T)
tickPosVector <- (tickVector/max.count) + bar.bottom
axis(at=tickPosVector , side=2, labels =tickVector, las=2)





############################

test <- barplot(counts, plot=F)	# just returns a vector of centroal heights.


plot.new(xlim=c(-1,15))


barplot(counts)
rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)


barplot(counts, space=-0.1)

rect(0:(n.counts-1),  0, 1:n.counts, counts/max.count)

#seq(0,1,length.out=n.counts)

rect(0:(n.counts-1),  0, 1:n.counts, counts/max.count)
rect(seq(0,1,length.out=n.counts+1)[-n.counts],  0, seq(0,1,length.out=n.counts+1)[-1], counts/max.count)

rect(1,0,2,20)
### trying to use code from Eulergrid.R  (Biostars).
plotEulergrid(plotTitle="Footprint__overlaps__for__multiple__cell__lines\n(FDR__0.001)", offCellColor="gray80", onCellColor="springgreen4", 
				setNames=c("GM06990","HepG2","K562","SKNSH","TH1"), setCardinalities=c(212350,233552,270586,287731,240701,93351,64049,89860,110579,62852,96806,89476,62075,64644,90129,30893,51178,53416,29083,32041,51033,28922,28279,48629,27407,22805,23548,39400,22418,21029,17172), 
				setTotal=689952, setTotalWithout=-1, 
				outputFilename="overlaps.fdr0p001.112409.png", showWholeSets=-1, ctsCardinalities=c(65897,97624,173336,150753,91965))




binaryCombinations <- function(n, r=0, values=c(0,1), resultVector="")  {
	
	newVector <- c(paste(resultVector,0,sep=""),paste(resultVector,1,sep=""))
	if(n > 1) 	{
		newVector <- binaryCombinations(n-1, r+1, resultVector=newVector)
	}

	return(newVector)
}



binaryCombinations(2)
binaryCombinations(5)


# or learn how to use built in from R!    permutations(2,5, repeats.allowed=T) -1
apply(permutations(2,5, repeats.allowed=T) -1   , 1, FUN=function(x) paste(x,  collapse=""))
