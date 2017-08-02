

arseFromElbow <- function(x, y, plot.it=TRUE, plot.elbow=TRUE, ...)  {
  
  
  # need to re-scale x to match y , and allow for uneven spacing of x values.
  span.x <- max(x) - min(x)
  span.y <- max(y) - min(y)
  x2 <-   min(y) +  (((x-min(x)) / span.x) * span.y) # x2 values are evenly spaced between min and max of y
  
  max.x2 <- max(x2)
  min.y <- min(y)
  
  br.dist <- sqrt((max.x2 - x2)^2 + (y-min.y)^2 )
  
  if(plot.it) {
    plot(x,y, xaxt='n', ...)
    segments(max(x), min(y), x[which.min(br.dist)], y[which.min(br.dist)])
    
    points(x[which.min(br.dist)], y[which.min(br.dist)] , col="red", pch=19, cex=1.5)
    
    if(plot.elbow)  text(x[which.min(br.dist)], y[which.min(br.dist)], labels=signif(y[which.min(br.dist)], digits=4), adj=c(2,-1), cex=2)
    
  }
  return(y[which.min(br.dist)])
}
