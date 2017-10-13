
# calculate the ratio between the intersect and the union of two GR ranges. 
inurt_GR <- function(x, y,  useCanonical=TRUE, ignore.strand=TRUE, extend=0,
                     verbose=FALSE)  {
 
  stopifnot(class(x) == "GRanges")
  stopifnot(class(y) == "GRanges")
  
  if(useCanonical) {
    seqlevels(x, pruning.mode="coarse") <- grep("_", seqlevels(x), invert=T, value=T)
    seqlevels(y, pruning.mode="coarse") <- grep("_", seqlevels(y), invert=T, value=T)
  }
   
  intersectTotal <- sum(as.numeric(width(GenomicRanges::intersect(x,y, ignore.strand=ignore.strand))))
  unionTotal <- sum(as.numeric(width(GenomicRanges::reduce(c(x,y), ignore.strand=ignore.strand))))
  return(intersectTotal/unionTotal)
}