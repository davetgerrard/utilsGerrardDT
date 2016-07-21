
# makes a scatterplot of x and y and adds a density plot to the top and right
# NEED to draw your own axis afterwards. 
scatterWithDensity <- function(x, y, xlab="x", ylab="y", smoothScatter=FALSE, ...)  {
  margin.x <- c(max(x), max(x)*1.2)
  margin.y <- c(max(y), max(y)*1.2) 
  
  #layout(matrix(c(2,0,1,3), 2,2 ,byrow=T))
  if(smoothScatter)  {
    smoothScatter(x,y, xlim=c(min(x), margin.x[2]), ylim=c(min(y), margin.y[2]), frame.plot=F, xlab=xlab, ylab=ylab, ...)
  } else {
  plot(x,y, xlim=c(min(x), margin.x[2]), ylim=c(min(y), margin.y[2]), frame.plot=F, xlab=xlab, ylab=ylab, ...)
  }
  x.dens <- density(x)
  lines( x.dens$x[x.dens$x < max(x)], ((x.dens$y[x.dens$x < max(x)]/max(x.dens$y)) * (margin.y[2] - margin.y[1])) + max(y),type="l", lwd=2)
  #abline(h=margin.y[1], col="grey")
  #lines(density(x)+1)
  y.dens <- density(y)
  lines(((y.dens$y[y.dens$x < max(y)]/max(y.dens$y)) * (margin.x[2] - margin.x[1])) + max(x), y.dens$x[y.dens$x < max(y)] ,  type="l", lwd=2)
  #abline(v=margin.x[1], col="grey")
  segments(margin.x[1],max(y),margin.x[2], max(y))
  segments(margin.x[1],max(y),margin.x[1],margin.y[2])
  text(x=mean(margin.x), y=mean(margin.y), "density")
  text(x=mean(margin.x), y=mean(margin.y), "density")
  text(x=margin.x[2], y=margin.y[1], "y", pos=3)
  text(x=margin.x[1], y=margin.y[2], "x", pos=4)
  segments(max(min(x),0),max(y),max(x), max(y),  col="grey")
  segments(max(x),max(min(y),0),max(x), max(y), col="grey")
}


#scatterWithDensity(x=gene.reads.filt$log10.max.reads, y=gene.reads.filt$tau, xlab="max reads (log10)", ylab="tissue-specificity (tau)", col="grey", axes=F, cex=1)
#title(paste("n=", nrow(gene.reads.filt)))
#axis(1,at=1:6, labels=10^(1:6))
#axis(2,at=seq(0, 1, 0.2), labels=seq(0, 1, 0.2), las=2)