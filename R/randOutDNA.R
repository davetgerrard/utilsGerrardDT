

randOutDNA <- function(n, stringLen=50,  dna=c("A", "C", "G", "T"))  {
 for (i in 1:n) {
  x <- paste(sample( dna, stringLen, replace=T), collapse="")
  print(x, quote=F)
 }
}

#randOutDNA(10)
#randOutDNA(10000)
