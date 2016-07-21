

groupStats <- function(x, groupedCols=list(names(x)), groupFun="sum")    {
  result <- data.frame()
  for(i in   1: length(groupedCols)) {
     thisGroup <- names(groupedCols)[i]
     
     if(length(groupedCols[[i]]) > 1) {
      thisResult <- apply(x[, groupedCols[[i]] ], 1, groupFun)
     } else {
       thisResult <- x[, groupedCols[[i]] ]
     }
     if(i==1) {
       result <- data.frame(v1=thisResult)
       names(result) <- thisGroup
     } else {
       result[, thisGroup] <- thisResult
     }
  }
  return(result)
}



testData <- data.frame(A=1:10, B=2:11, C=3:12, D=21:30, E=LETTERS[1:10])
testGroups <- list(group1=c("A", "C") , group2=c("B", "D"))  
  
groupStats(testData, groupedCols=testGroups)  
groupStats(testData, groupedCols=testGroups, groupFun="mean")  

rm(testData, testGroups)
