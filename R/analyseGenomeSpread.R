


# given a GRanges object, measure the features and make statements about spread around the genome.
# ideally, the GRanges object should have seqlengths info.

analyseGenomeSpread.GR <- function(x, show.plots=FALSE, n.plotChroms=10)  {
  stopifnot(class(x) == "GRanges")
  
  formatBig <- function(x) format(x,big.mark=",",scientific=FALSE, trim=TRUE)  # for repeated use on large numbers.
  
  allChromNames <- seqnames(seqinfo(x))
  chrom.canonical <- grep("_", allChromNames, invert=TRUE, value=TRUE)
  chrom.canonical <- chrom.canonical[order(seqlengths(x)[chrom.canonical])]
   # order by size
  genomeSize <- sum(as.numeric(seqlengths(x)))
  
  x.canonical <- x[seqnames(x) %in% chrom.canonical]
  #seqinfo(x.canonical) <- seqinfo(x)[chrom.canonical]
  seqlevels(x.canonical, pruning.mode="coarse") <- chrom.canonical   # will throw out any data on noncanonical chroms.
  genomeSize.canonical <- sum(as.numeric(seqlengths(x.canonical)))
  
  ### Drop all unused seqlevels:  - could be useful later..
  #seqlevels(gr) <- seqlevelsInUse(gr)
  
  #print(seqinfo(x))
  print(paste("Object with", formatBig(length(x)), "regions from genome", unique(genome(x))))
  print(paste("The full genome has length", formatBig(genomeSize), "over", length(allChromNames), "sequences"))  
  print(paste("The canonical genome has length", formatBig(genomeSize.canonical), "over", length(chrom.canonical), "sequences")) 
  
  n.chroms.used <- length(seqlevelsInUse(x.canonical))
  choms.notUsed <- setdiff(chrom.canonical,seqlevelsInUse(x.canonical))
  n.features <- length(x.canonical)
    
  print(paste("Your track contains", formatBig(n.features), "features over", n.chroms.used, "sequences"))
  if(length(choms.notUsed) > 0) {
    print(paste("The following chromosomes have zero features:-", paste(choms.notUsed, collapse=",")))
  } else {
    print("All canonical chroms have features")
  }
  
  coverage.canonical <- sum(width(x.canonical))
  x.reduced <- reduce(x.canonical)
  n.features.reduced <- length(x.reduced)
  coverage.reduced <- sum(width(x.reduced))
  
  print(paste("The features cover", formatBig(coverage.canonical), "bases or ", round((coverage.canonical/genomeSize.canonical)*100, digits=3)  , "% of the genome"))
  if(n.features.reduced < n.features) {
    print(paste("Some features overlap each other."))
    print(paste("The intersect of features cover", formatBig(coverage.reduced), "bases or ", round((coverage.reduced/genomeSize.canonical)*100, digits=3)  , "% of the genome"))
    x.use <- x.reduced
  } else {
    print(paste("No features overlap each other"))
    x.use <- x.canonical
  }
  
  # calc counts per chrom
  chrom.counts <- table(seqnames(x.canonical))
  max.count <- max(chrom.counts)
  max.count.chrom <- names(chrom.counts)[which.max(chrom.counts)]
  min.count <- min(chrom.counts)
  min.count.chrom <- names(chrom.counts)[which.min(chrom.counts)]
  print(paste("The sequence with most features is", max.count.chrom, "with" , formatBig(max.count)))
  print(paste("The sequence with fewest features is", min.count.chrom, "with" , formatBig(min.count)))
  
  # calc count per kilobase
  cpkb <- (chrom.counts[chrom.canonical] / seqlengths(x.canonical)) * 1000000
  max.cpkb <- max(cpkb)
  max.cpkb.chrom <- names(cpkb)[which.max(cpkb)]
  min.cpkb <- min(cpkb)
  min.cpkb.chrom <- names(cpkb)[which.min(cpkb)]
  print(paste("The sequence with most features per Mb is", max.cpkb.chrom, "with" , round(max.cpkb, digits=3)))
  print(paste("The sequence with fewest features per Mb is", min.cpkb.chrom, "with" , round(min.cpkb, digits=3)))
  
  # calc coverage per chrom.
  chrom.sums <- by(x.canonical, INDICES=as.character(seqnames(x.canonical)), FUN=function(y) sum(width(y)))
  chrom.sums.vec <- as.integer(chrom.sums)
  chrom.sums.reduced <- by(x.reduced, INDICES=as.character(seqnames(x.reduced)), FUN=function(y) sum(width(y)))
  chrom.sums.vec.reduced <- as.integer(chrom.sums.reduced)
  names(chrom.sums.vec) <- names(chrom.sums)
  names(chrom.sums.vec.reduced) <- names(chrom.sums.reduced)
  if(show.plots) {
    plot(seqlengths(x.canonical)[chrom.canonical],chrom.sums.vec[chrom.canonical] , ylim=c(0, max(chrom.sums.vec)),xlab="sequence length", ylab="coverage")
    points(seqlengths(x.reduced)[chrom.canonical],chrom.sums.vec.reduced[chrom.canonical], col="blue")
  }
  # plot coverage per chrom against seqlengths.
  # add on non-redundant coverage per chrom.
  
  print(paste("Correlation between coverage and chromosome lengths:" , round(cor(seqlengths(x.use)[chrom.canonical],chrom.sums.vec[chrom.canonical], use="complete.obs"), digits=3 )))
  
  x.mean <- mean(width(x.use))
  x.meanTrimmed <- mean(width(x.use), trim=.1)
  x.sd <- sd(width(x.use))
  x.CoV <- x.sd / x.mean
  
  # do some stats on distribution of lengths. 
  # Are all features from the same distribution?
  print(paste("Mean length:", round(x.mean, digits=3)))
  print(paste("Trimmed mean length:", x.meanTrimmed))
  print(paste("Standard deviation:", round(x.sd, digits=3)))
  print(paste("Coefficient of variation:", round(x.CoV, digits=3)))
  print(paste())
  # What about clumping within chromosomes?  
  
  print(paste())
}


# source("C:/Users/Dave/utilsGerrardDT/R/analyseGenomeSpread.R")

#library(BSgenome.Hsapiens.UCSC.hg38)
#myData <- GRanges(paste0("chr", rep(c(1:22, "X", "Y", "M", "1_GL383518v1_alt"),  length.out=1000)) , IRanges(1:1000, width=1:1000), seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))
#analyseGenomeSpread.GR(myData, show.plots=TRUE)


# ideas for good variance measures
# http://stats.stackexchange.com/questions/18661/mean-and-variance-of-a-zero-inflated-poisson-distribution
# https://stat.ethz.ch/pipermail/r-help/2007-November/146285.html
