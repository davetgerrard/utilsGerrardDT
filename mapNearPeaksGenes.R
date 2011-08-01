#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

## Some function to help find closest features in one set of ordered features within another set of ordered features. 
## designed for genome features. Such as TF binding sites and chip-peaks. 



############ FUNCTIONS

##  for use in sapply e.g. > sapply(two,findGreatestLowerValue,one) 
## returns a single integer index
findGreatestLowerValue <- function(x,vec){   #  modified from http://r.789695.n4.nabble.com/find-closest-match-between-two-vectors-td3091320.html
 	y <- vec - x
	y[y>=0] <- NA
 	if(all(is.na(y)))NA else which.max(y)
}

##  for use in sapply e.g. > sapply(two,findLeastLargerValue,one) 
## returns a single integer index
findLeastLargerValue <- function(x,vec){   #  modified from http://r.789695.n4.nabble.com/find-closest-match-between-two-vectors-td3091320.html
 	y <- x - vec
	y[y>=0] <- NA
 	if(all(is.na(y)))NA else which.max(y)
}

countExactMatches <- function(x, vec)  {
	sum(!is.na(match(x, vec)))
}

##  for use in sapply e.g. > sapply(two,findClosestSingleValue,one) 
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


# window_type: one of "centred" (default), "left", "right" 
# overlap_type: one of "contained", "partial"
# impossible values (e.g. beyond end of chromosomes) are NOT accounted for 
# because purpose is to return existing features.
### NON-FUNCTIONAL
findFeaturesInWindow <- function(x, vec, window_size, window_type="centred", overlap_type="contained")  {

	#findInterval  not sure if this can be used in proper vector manner.
	
	
	#findInterval(feature.starts, window)
}


# can be used to find ranged features if used on start and end vectors and then intersected.
# e.g. intersect( findPointFeaturesInWindow(x, feature_STARTS, win.size), findPointFeaturesInWindow(x, feature_ENDS,  win.size)
# feature_starts <- seq(1000, 100000, by=1000)
# findPointFeaturesInWindow(15000, feature_starts, 5000)
# intersect(findPointFeaturesInWindow(15000, feature_starts, 5000), findPointFeaturesInWindow(15000, feature_ends, 5000))
findPointFeaturesInWindow <- function(x, vec, window_size, window_type="centred")  {
	
	#window_start <- x - floor(window_size/2)
	#window_end <- window_start + window_size
	window <- makeWindow(x,window_size, window_type)
	which(findInterval(vec,window) == 1)
}

## generic window making function. Works on single points or vectors. 
# window_type: one of "centred" (default), "left", "right" 
makeWindow <- function(x, window_size, window_type="centred")  {
	window_start <- switch(window_type,
		centred = x - floor(window_size/2),
		left = x - window_size,
		right = x
	)
	window_end <- window_start + window_size
	if(length(x) ==1) {
		return(c(window_start, window_end))
	} else {
		return(cbind(window_start,window_end))
	}
}



###example

#	closestSinglesVector <- sapply(set.1.subset[,sortBy.1[2]], findClosestSingleValue, set.2.subset[,sortBy.2[2]])
	# WHAT TO DO ABOUT SHARED COLUMNS NAMES? prefix? 
#	thisLevelTable <- cbind(set.1.subset, set.2.subset[closestSinglesVector,])






















