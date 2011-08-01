#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

############ FUNCTIONS

##  for use in sapply e.g. > sapply(two,findit,one) 
## returns a single integer index
findGreatestLowerValue <- function(x,vec){   #  modified from http://r.789695.n4.nabble.com/find-closest-match-between-two-vectors-td3091320.html
 	y <- vec - x
	y[y>=0] <- NA
 	if(all(is.na(y)))NA else which.max(y)
}

## returns a single integer index
findLeastLargerValue <- function(x,vec){   #  modified from http://r.789695.n4.nabble.com/find-closest-match-between-two-vectors-td3091320.html
 	y <- x - vec
	y[y>=0] <- NA
 	if(all(is.na(y)))NA else which.max(y)
}

countExactMatches <- function(x, vec)  {
	sum(!is.na(match(x, vec)))
}

## returns a single integer index
findClosestSingleValue <- function(x, vec) {
	if(countExactMatches(x,vec) > 0)  {
		return(match(x,vec)[1])				# if multiple exact matches, take the first. 
	} else {
		larger <- findLeastLargerValue(x, vec)
		smaller <- findGreatestLowerValue(x, vec)
		if(is.na(larger))  {
			return(smaller)
		} else if (is.na(smaller)) {
			return(larger)
		} else {
			y <- abs(x - c(vec[smaller],vec[larger]))
			return(c(smaller,larger)[which.min(y)])
		}	
	}
}


#findExactMatches <- function(x, vec)  {
#
#}

#could iterate over closest values to find chain of closest values. 























