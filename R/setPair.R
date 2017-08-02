setPair <- function(setA, setB) {
  
  out <- data.frame(unique.A=length(setdiff(setA,setB)),
                    unique.B=length(setdiff(setB,setA)),
                    shared=length(intersect(setA,setB)))
  return(out)
}



#setPair(c(1:5), c(3:10))