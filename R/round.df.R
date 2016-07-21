
# applies digits to all numeric columns on all numeric columns and 
round.df <- function(df, digits=3)   {
  #make index of numeric columns. 
  # could not get apply to work.
  num.index <- unlist(lapply(df,  class))  == "numeric"
  
  df[,num.index] <- round(df[,num.index], digits=digits)
  
  return(df)
  
}


round.file <- function(file,digits=3,  sep="\t", quote=F, row.names=F, suffix=".tab") {
 df <- read.delim(file)
 new.df <- round.df(df)
 outFile <- sub(suffix, paste(".dig", digits, suffix, sep=""), file)
 stopifnot(file != outFile) 
 write.table(new.df, outFile, sep=sep, quote=quote, row.names=row.names)
}

#df <- read.delim("results/xie.jennings.rank.pca.tab")
#head(df)
#apply(df, 2, is)

#new.df <- round.df(df)

#fileName <- "results/xie.jennings.rank.pca.tab"
#round.file(fileName)
#write.table(new.df, file="results/xie.jennings.rank.pca..dig3.tab", quote=F, row.names=F, sep="\t")
