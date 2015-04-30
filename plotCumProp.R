plotCumProp <- function(x, line.col="blue", xlims=c(-1,1), add=FALSE, xlab="score", res=1000, ...)  {
  if((min(x) < xlims[1]) | (max(x) > xlims[2]))  warning("scores outside plotting range!")
  
  breaks <- seq(xlims[1], xlims[2], length.out=res)
  
  mean.cut <- cut(x, breaks, right=F)
  mean.freq <- table(mean.cut)
  cumfreq0 = c(0, cumsum(mean.freq)) / length(x)
  if(add)  {
    lines(breaks, cumfreq0 , xlab=xlab, ylab="cumulative proportion", type="l", lwd=2, col=line.col)
    
  } else {
    #plot.new()
    plot(breaks, cumfreq0 , xlab=xlab, ylab="cumulative proportion", type="l", lwd=2, col=line.col, xlim=xlims, ylim=c(0,1), ...)
  }
}