

# phi correlation from two genomic ranges objects
# requires seqinfo for whole genome or genome size
# AB  Both Groups
# A   A , not B
# B   B, not A
# C   neither A nor B
calcPhiCoef <- function(AB, A, B, C) {
  return(( AB * C - A*B )   / sqrt((AB+A)*(B+C)*(AB+B)*(A+C)))
}

#calcPhiCoef(1,0,0,1)   #1
#calcPhiCoef(0,1,1,0)  #-1

# requires vcd package for assocstats
# requires GenomicRanges
phiCorGR <- function(x, y, genomeSize=NULL, useCanonical=TRUE, ignore.strand=TRUE,
                     verbose=FALSE)  {
  # require(vcd)   # removed because phi coeff is absolute value and I need sign
  require(GenomicRanges)
  
  # test we have what we need
  stopifnot(class(x) == "GRanges")
  stopifnot(class(y) == "GRanges")
  stopifnot(sum(as.numeric(seqlengths(x))) == sum(as.numeric(seqlengths(y))) | is.numeric(genomeSize))
  
  if(useCanonical) {
    seqlevels(x, pruning.mode="coarse") <- grep("_", seqlevels(x), invert=T, value=T)
    seqlevels(y, pruning.mode="coarse") <- grep("_", seqlevels(y), invert=T, value=T)
  }
  
  if(is.null(genomeSize)) genomeSize <-  sum(as.numeric(seqlengths(x)))
  
  intersectTotal <- sum(as.numeric(width(GenomicRanges::intersect(x,y, ignore.strand=ignore.strand))))
  xOnlyTotal <- sum(as.numeric(width(GenomicRanges::setdiff(x,y, ignore.strand=ignore.strand))))
  yOnlyTotal <- sum(as.numeric(width(GenomicRanges::setdiff(y,x, ignore.strand=ignore.strand))))
  neitherTotal <- genomeSize - (intersectTotal + xOnlyTotal + yOnlyTotal)
  tablePhi <- matrix(c(intersectTotal, xOnlyTotal, yOnlyTotal, neitherTotal), nrow=2, ncol=2, dimnames=list(x=c(TRUE, FALSE), y=c(TRUE,FALSE)))
  if(verbose) print(tablePhi)
  #return(assocstats(tablePhi)$phi)
  return(calcPhiCoef(AB=intersectTotal, A=xOnlyTotal, B=yOnlyTotal, C=neitherTotal))
}

#library(BSgenome.Hsapiens.UCSC.hg38)
#x.GR <- GRanges(c("chr1", "chr2", "chrX"), IRanges(start=c(10000, 20000, 30000), width=400), seqinfo=seqinfo(Hsapiens))
#y.GR <- GRanges(c("chr1", "chr2", "chrX","chrY"), IRanges(start=c(10200, 20100, 30050, 40000), width=600), seqinfo=seqinfo(Hsapiens))

#phiCorGR(x.GR, y.GR )
#phiCorGR(x.GR, x.GR )  # should be 1

#calcPhiCoef(1,0,0,1)   #1
#calcPhiCoef(0,1,1,0)  #-1
