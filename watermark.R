
# prints text faintly in bottom right hand corner. Good with getCurrentScritp.R
watermark <- function(x, inc.date=T, side=1, line=floor(par()$mar[side]) -1, adj=1, cex=.7, col="lightgrey", ...)  {
  
  if(inc.date) { x <- paste(x, date()) }
  
  mtext(x, side=side, line=line, adj=adj, col=col, cex=cex)

}

#watermark(getCurrentSCript())