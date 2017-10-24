
# ToDO breakdown summary by directory? or have a second wrapper function to do that.
# might need to be carefule how this is combined with using different suffixes. 

# read all script files from a directory (with recursion option)
# count total lines, blanks lines, commented lines.
# tabulate
# optional count occurence of certain phrases (within non comment lines).
codeCountTable <- function(dir="./", suffixes="R", verbose=FALSE)  {
  
  suffixPattern <- paste0("*.", suffixes, "$")   # suffix is last part of filename.
  
  fileVec <- dir(path=dir,pattern=suffixPattern, recursive=T, full.names = T)
  # file.info(fileVec)  # might be useful in future
  
  codeTable <- data.frame()
  for(thisFile in fileVec) {
    lines <- readLines(thisFile)
    cleanLines <- trimws(lines, "left")    # remove leading white space.
    
    
    commentLineIndex <- grep("^#", cleanLines)
    blankLineIndex <- which(cleanLines == "")
    annotCodeIndex <- setdiff(grep("#", cleanLines), commentLineIndex)
    thisRow <- data.frame(files=1,total=length(cleanLines), code=length(cleanLines)-length(c(commentLineIndex, blankLineIndex)), 
                          comment=length(commentLineIndex), blank=length(blankLineIndex), annot.code=length(annotCodeIndex),
                          file=thisFile)
    codeTable <- rbind(codeTable, thisRow)
  }
  
  if(verbose) print(colSums(codeTable[,c("files","total", "code", "comment", "blank", "annot.code")]))
  return(codeTable)
  


}
  
#x <- codeCountTable(dir="analysis/", verbose=T)

