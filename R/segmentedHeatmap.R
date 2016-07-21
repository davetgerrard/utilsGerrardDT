



#source('C:/Users/Dave/utilsGerrardDT/utility.R')

segmentedHeatmap <- function(x,  groupList, group.labels=names(groupList), sample.labels = colnames(x), cex.ylabels=1, cex.xlabels=1, add.grid=TRUE)  {
  #old.par <- par()
  #par(mar=c(5, 10, 3,3))
  full.index <- match(unlist(groupList), row.names(x))
  
  group.bounds <- as.integer(cumsum(lapply(groupList, length)))
  groupBreaks <- c(0, group.bounds + .5) / max(group.bounds)
  # find the midpoints between successive edges.
  groupCentres <- groupBreaks[-1] - (( groupBreaks[-1] - groupBreaks[1:(length(groupBreaks)-1)]) / 2)
  
  (groupLabels <- as.character(unlist(lapply(group.labels, utility.insertNewLines, length=25))))
  
  
  #image(t(x[full.index,])                    )
  image(scale(t(x[full.index,])), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ))
  box()
  if(add.grid) grid(nx=ncol(x), ny=length(full.index), col="grey")
  
  par(xpd=T)
  segments(-1, y0=groupBreaks, x1=-.3, lwd=2, lty=1)
  #segments(1.07, y0=groupBreaks, x1=1.2, lwd=2, lty=1)
  
  #for(i in 1:(length(group.bounds)-1))  {
  #  
  #}
  mtext(groupLabels, side=2, at=groupCentres, las =2, line=2, cex=cex.ylabels)
  mtext(sample.labels, side=1, at=seq(0, 1, length.out=length(sample.labels)), las=2, line=1, cex=cex.xlabels)
  par(xpd=F)
}


# myData <- as.matrix(data.frame(time.1 = c(rep(3, 14), rep(0, 12)), 
#                     time.2= c(rep(2, 18), rep(1, 8)), time.3=c(rep(1, 10), rep(3, 16))))
# row.names(myData) <- LETTERS

# groups <- list(A= LETTERS[1:8], B=LETTERS[7:18], C=LETTERS[14:24])
# group.labels <- list(A="a long name for the first group", B="a long name for the second group", C=" a shorter name")                 

# segmentedHeatmap(x=myData, groupList=groups, group.labels=group.labels)
