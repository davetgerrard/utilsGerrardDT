
# rescale a vector into a new range or expand/contract the range by a constant value centred on midpoint.
rescaleVector <- function(x, newRange=NULL, proportion=1)  {
  span.x <- max(x) - min(x)
  if(is.null(newRange)) {
    midpoint <- min(x) + (span.x/2)
    newRange <- c( midpoint- ((midpoint-min(x))*proportion), midpoint + ((max(x) - midpoint)*proportion) )
  } 

    span.y <- max(newRange) - min(newRange)
  x2 <-   min(newRange) +  (((x-min(x)) / span.x) * span.y) 
  return(x2)
}


#rescaleVector(x= c(1:10))   # should be unchanged

#rescaleVector(x= c(1:10), newRange=40:60) 
#rescaleVector(x= c(1:10), newRange=40:60, proportion=5)   # should ignore the proportion
#rescaleVector(x= c(1:10), proportion=5)
#range(rescaleVector(x= c(1:10), proportion=5))
#range(rescaleVector(x= c(1:10), proportion=.5))