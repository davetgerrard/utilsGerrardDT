

# returns a table (bottom-left including diagonal) denoting whether pairs within a set of features overlap.
getOverlapMatrix <- function(starts, ends, extension=0, half_matrix=FALSE)  {
	n <- length(starts)
	overlap.matrix <- matrix(NA, nrow=n, ncol=n)

	for(i in 1:length(starts))  {
		overlap.matrix[i:n,i]  <-  starts[i:n] < ends[i]
	}
		
	if(!half_matrix) {
		for(i in 1:(n-1)) {
			for(j in (i+1):n) {
				overlap.matrix[i,j] <- overlap.matrix[j,i]
			}	
		}
	}
	return(overlap.matrix)

}


recursePaths <- function(map, currentPath, index, pathList=list()) {
	max.n <- length(map)
	  if(length(setdiff(currentPath, map[[index]])) > 0 ||  identical( map[[index]] , as.integer( currentPath))) {		
		# reached an end
		pathList[[length(pathList) + 1]] <- currentPath
	} else {	
		newPath <- sort(unique(c(currentPath, index)))
		for(i in setdiff(map[[index]], currentPath))  {
			pathList <- recursePaths(map, currentPath=newPath, index=i, pathList=pathList)
		}
	}
	return(pathList)
}


minPathNumber <- function(junctionsTable, extension=0)  {
	


	# make a copy of the junctions to add extensions to act as exons
	# also remove junctions with exact same start and end.
	newJunctions <- unique(subset(junctionsTable, select=c(chr, start, end)))
	#newJunctions$start <- newJunctions$start - extension
	#newJunctions$end <- newJunctions$end + extension
	
	n <- nrow(newJunctions)		

	if(n == 1)  {		# only 1 junction
		return(1)
	}

	# create a matrix reporting which junctions overlap with each other.
	overlap.matrix <- getOverlapMatrix(newJunctions$start, newJunctions$end)

	#juncs <- 1:n
	if(sum(overlap.matrix , na.rm=T )  == n)  {	# no overlapping junctions, except with themselves.
		return(1)
	}

	compat.list <- list()
	for(i in 1:n)  {
		compat.list[[i]] <-  sort(c(i,which(!overlap.matrix[,i])))
	}
	
	# Find all the paths for each index of the map. 
	# This seems a very ugly (inefficient) way to do it.
	fullList <- list()
	for(i in 1:n)  {
		fullList <- c(fullList , unique(recursePaths(compat.list, currentPath=c(i), index=i, pathList=list())) )
	}
	fullList <- unique(fullList)

	# Find redundant sets from the list.   (e.g. if have c(1,2) and c(1,2,3), want to keep only c(1,2,3).
	redundant.index <- integer()
	for(i in 1:length(fullList))  {
		for(j in 1:length(fullList))  {
			if(i != j)  {
				if(length(setdiff(fullList[[i]], fullList[[j]])) == 0 )  {
					redundant.index <- c(redundant.index, i)
				} 
			}
		}
	}

	# Remove the redundant elements from the list of paths.
	if(length(redundant.index) > 1)  {
		finalList <- fullList[-redundant.index]
	}  else  {
		finalList <- fullList
	}
	return(length(finalList))

}




  




