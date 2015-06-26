make.unique.all <- function(x, sep=".")  {   # slower than make.unique() but does not  leave first occurrence unchanged
  y <- as.character(x)
  n_occur <- data.frame(table(y))
  for(i in which(n_occur$Freq > 1)) {
    this.y <- n_occur$y[i]
    this.freq <- n_occur$Freq[i]
    this.index <- which(y == this.y)
    y[this.index] <- paste(this.y, 1:this.freq, sep=sep)
  }  
  return(y)
}

#make.unique.all(c("A", "B", "C", "A", "B"), sep="-")