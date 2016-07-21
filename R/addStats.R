#source("C:/Users/Dave/utilsGerrardDT/calcTau.R")   # for tau

# add a list of stats to a data.frame  
# mean, median, sd, tau, min, max

addStats.df <- function(x, columns=colnames(x), stats=c("min", "max", "sd", "tau"))  {
  
  require(utilsGerrardDT)
  check.functions <- logical()
  for(thisFunc in stats)  check.functions[thisFunc] <- is.function(get(thisFunc))   # could not get mget() to work.
  stopifnot(all(check.functions))
  
  #check if columns already exist with these names
  stopifnot(length(intersect(colnames(x), stats)) == 0)
  
  for(thisFunc in stats)  {
    x[,thisFunc] <- apply(x[,columns], 1, thisFunc)
  }
  
  return(x)
}
  
#c(min, max, sd, tau)
