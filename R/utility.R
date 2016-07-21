#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

# source("C:/Users/dave/scripts/utility.R")

############ FUNCTIONS
## inserts "\n" into a text line at a given length spacing. Useful with textplot
utility.insertNewLines <- function(thisText,length=60)  {
	#test if string contains whitespace
	#
	splitLine <- character()
	tokens <- unlist(strsplit(thisText," "))
	while(length(tokens) > 0)  {
		thisLineIndex <- which(cumsum(nchar(tokens)+1) < (length+2) )
		thisLine <- paste(tokens[thisLineIndex],collapse=" ")
		splitLine <- paste(splitLine,thisLine,"\n",sep="")
		tokens <- tokens[-thisLineIndex]

	}
	return(splitLine)
}




## calculate expanded or contracted axis limits for a plot. Returns vector of two values.
wideLimits <- function(data, scale=1.2) {
	range <- max(data) - min(data)
	midpoint <- (max(data) + min(data) ) /2
	range <- range * scale
	minMark <- signif(midpoint - (range/2),3)
	maxMark <- signif(midpoint + (range/2),3)
	return(c(minMark,maxMark))	
}


splitStringToTable <- function(stringVector, sep=";")  {
	
	resultList <- strsplit(as.character(stringVector), sep, fixed=T)
	n.rows <- length(resultList)
	if(n.rows != length(stringVector)) stop("Failed!")
		
	splitMatrix <- matrix(unlist(resultList), nrow=n.rows, byrow=T)	
	return(as.data.frame(splitMatrix))
}
#testVector <- paste("name", 1:10, letters[1:10], sep=".")
#table.out <- splitStringToTable(testVector, sep=".")


# Add text notes to specific rows of a data.frame
# work in progress
note <- function(df, text, append=TRUE, row.ids, id.column="name", note.column="notes", append.sep=" ; ")  {

	# cannot use is.null(df['notes']) as gives error, even though is.null(df$notes) would work
	if(is.na(match(note.column, names(df))))   df[,note.column] <- ""


	# test for valid row.ids
	if(length(row.ids) < 1)  stop("no row.ids supplied")
	
	# test for multiple matches

	# is text a vector? if so, must match length of row.ids


	row.index <- match(row.ids, df[,id.column])
	if(length(row.index) < 1)  stop("No mathcing rows found")
	if(append) {
		df[row.index,note.column] <- paste(df[row.index,note.column], text, sep=append.sep)
	} else {
		df[row.index,note.column] <- text	
	}
	
	return(df)
}  
#test.df <- data.frame(name=letters[1:10], length=runif(10))
#test.df <- note(test.df, "test",  row.ids="c")



# x is a matrix or data frame with two columns.
# should getting plotting zone from graphics device?
arrangeLabels <- function(x, style="vertical", weights=NULL, xlim=c(0,1), ylim=c(0,1))  {
    x.order <- rank(x[,1], ties.method="first")
    y.order <- rank(x[,2], ties.method="first")
    
    
    
    
    y.pos <- seq(min(ylim), max(ylim), length.out=length(y.order))
    x.pos <- rep(max(xlim), length(x.order))
    
    if(!is.null(weights))  {
        y.order.wt <- y.order * (weights/sum(weights))
        y.order <- rank(y.order.wt, ties.method="first")
    }
    
    return(list(points=x, lines=cbind(x.pos[x.order], y.pos[y.order])))
}
