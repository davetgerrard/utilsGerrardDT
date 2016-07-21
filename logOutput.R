

# utility script to register each outfile with the script (and line/block) that created it.

# ?findLineNum  # probably not that helpful

# write out date, outputFileName, source script, getwd().


logOutput <- function(outFile, logFile="", tag=NA)  {
  if(!exists("getCurrentSCript") ) source( "C:/Users/Dave/utilsGerrardDT/getCurrentScript.R")
  script <- try(getCurrentScript())
  
  outLine <- paste(Sys.time(), outFile,  script, tag, getwd(),  sep="\t")
  cat(outLine, file=logFile, append=TRUE)
  cat("\n", file=logFile, append=TRUE)
}


#logOutput(outFile="test.out")
