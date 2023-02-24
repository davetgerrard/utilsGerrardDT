

makeRandomNames <- function(alphabet=c( letters, 0:9), req_len=1, string_len=10, unique=TRUE)  {
  
  #vec <- c( letters, 0:9)
  #req_len <- 10000
  #string_len <- 4
  
  finalSet <- character()
  while(length(finalSet) < req_len) {
    
    
    finalSet <- unique(c(finalSet, paste(sample(alphabet, string_len), collapse="", sep="")))
  }
  
  
  return(finalSet)
}


makeRandomNames()
makeRandomNames(req_len=10, string_len=4) 
makeRandomNames(alphabet=LETTERS, req_len=10, string_len=4) 

#length(unique(finalSet))
#length(finalSet)


# how many are required
#length(alphabet) ^ string_len
#length(alphabet) ^4



#?clip
#cat(finalSet[1:50], file="clipboard-", sep="\n")
#cat(finalSet, file="clipboard-16384", sep="\n")  # need the ੢
#cat(finalSet, file="C:/Temp/tempIds.txt"  , sep="\n")
#瘴楺眊㕡 to 
