#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

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



