###
### Given a sequence of numbers, splits into clusters separated by at least split.at

clusterSplitter <- function(vec, split.at, returnType="group")  {
    
    intervals <- vec[2:length(vec)] - vec[1:(length(vec)-1)]
    splits <- c(0,which(intervals >= split.at), length(vec))
    
    n.groups <- length(splits) - 1
    run.lengths <- splits[2:length(splits)] - splits[1:(length(splits)-1)]
    #return(which(intervals >= split.at))
    #return(run.lengths)
    #return(splits)
        
    return(rep(1:n.groups, times=run.lengths))
}


test.vec <- c(1,2,6,7,11,30, 33,34,50, 56, 59, 65, 67, 69, 89, 91)

clusterSplitter(test.vec, 10)

#get range of members in first group
range(which(clusterSplitter(test.vec, 10)==1))

#get range of members in third group
range(which(clusterSplitter(test.vec, 10)==3))


# things I tried that didn't work
#density(test.vec)
#library(classInt)
#classIntervals(test.vec)
