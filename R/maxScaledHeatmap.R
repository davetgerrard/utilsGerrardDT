




maxScaledHeatmap <- function(df, minMax=0, showMaxScale=FALSE, maxScale.ints=1:5, showIntenseScale=TRUE, showLegend=TRUE, useLog10scale=FALSE,  x.adjust = .020 ,
                             legend.title="intensity", mar.vec =c(10,10,10,10), name.sides=1:4, scale.colour="forestgreen", 
                             font.row=1,font.col=1, cex.row =1, cex.col=1, row.labs=row.names(df), col.labs=colnames(df), showGrid=TRUE)  {
  
  max.vec <- apply(df, 1, max)
  exclude.index <- max.vec < minMax
  df[exclude.index,] <- 0
  #df <- scale(df, center=F, scale=apply(df, 2, max))
  
  if(useLog10scale)  max.vec <- log10(max.vec)
  
  if(nrow(df) > 0 ) {
    if(showIntenseScale)  {
      if(useLog10scale) {
        legend.ints <- 1:ceiling(max(max.vec))
      } else {
        legend.ints <- pretty(max.vec, min.n=3, n=4)
      }
      
      
      par(mar=mar.vec)
      #x.adjust <- .020  # adjustment to align grid with image(). Now a parameter
      y.adjust <- x.adjust * (ncol(df)/nrow(df))
      # more complex plot with scale bar and gap between it and hte rest of the image/grid.
      image(scale(t(df), center=F, scale=apply(df, 1, max)), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ), xlim=c(-.2,1 + x.adjust), add=F)
      image(t(as.matrix(max.vec)), zlim=c(0, max(max.vec)), add=T,col=colorRampPalette(c("white", scale.colour))( 30), x=c(-.2, -.1), y=seq(0,1, length.out=length(max.vec)))
      #abline(v=seq(0,1, length.out=ncol(df)), h=seq(0,1, length.out=nrow(df)))
      #grid()
     
      x.corners <- seq(0-x.adjust,1 + x.adjust, length.out=ncol(df)+1)
      y.corners <- seq(0-y.adjust,1 + y.adjust, length.out=nrow(df)+1)
      if(showGrid) rect(xleft=rep(x.corners[1:ncol(df)], each=nrow(df)), ybottom=rep(y.corners[1:nrow(df)], times=ncol(df)), xright=rep(x.corners[1:(ncol(df)+1)], each=nrow(df)), ytop=rep(y.corners[2:(nrow(df)+1)], times=ncol(df)), lty=3, lwd=1, border="grey")
      if(showGrid) rect(xleft=-.2, ybottom=y.corners[1:nrow(df)], xright=-.1, ytop=y.corners[2:(nrow(df)+1)], lty=3, border="grey")
      if(4 %in% name.sides) mtext(row.labs, side=4, at=seq(0,1,length.out=nrow(df)), las=2, line=1, font=font.row, cex=cex.row)
      if(1 %in% name.sides) mtext(col.labs, side=1, at=seq(0,1,length.out=ncol(df)), las=2, line=1, font=font.col, cex=cex.col)  
      if(2 %in% name.sides) mtext(row.labs, side=2, at=seq(0,1,length.out=nrow(df)), las=2, line=1, font=font.row, cex=cex.row)
      if(3 %in% name.sides) mtext(col.labs, side=3, at=seq(0,1,length.out=ncol(df)), las=2, line=1, font=font.col, cex=cex.col)  
      
      # add a legend for intensity scaling
      if(showLegend) {
        
        if(useLog10scale)  {
          
          
          par(xpd=TRUE)
          legend(-.1, max(y.corners)+y.adjust, legend=rev(10^legend.ints), fill=rev(colorRampPalette(c("white", scale.colour))( length(legend.ints)  )), title=legend.title, xjust=1, yjust=0)
          par(xpd=FALSE)
          
        } else {
          
          par(xpd=TRUE)
          legend(-.1,  max(y.corners)+y.adjust, legend=rev(legend.ints), fill=rev(colorRampPalette(c("white", scale.colour))( length(legend.ints) )), title=legend.title, xjust=1, yjust=0)
          par(xpd=FALSE)
        }
      }
    } else {
      par(mar=mar.vec)
      # scale using max value
      image(scale(t(df), center=F, scale=apply(df, 1, max)), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ))
      #image(t(scale(df, center=F, scale=TRUE)), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ))
      #image(t(df), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ))
      box()
      if(showGrid) grid(nx=ncol(df), ny=nrow(df), col="grey")
      if(4 %in% name.sides) mtext(row.labs, side=4, at=seq(0,1,length.out=nrow(df)), las=2, line=1, font=font.row, cex=cex.row)
      if(1 %in% name.sides) mtext(col.labs, side=1, at=seq(0,1,length.out=ncol(df)), las=2, line=1, font=font.col, cex=cex.col)  
      if(2 %in% name.sides) mtext(row.labs, side=2, at=seq(0,1,length.out=nrow(df)), las=2, line=1, font=font.row, cex=cex.row)
      if(3 %in% name.sides) mtext(col.labs, side=3, at=seq(0,1,length.out=ncol(df)), las=2, line=1, font=font.col, cex=cex.col)  
      #meta <- paste("TFgrid.R: read counts", normalisation, "minMax", minMax)
      #mtext(meta, side=1, at=1, line=11, las=1, cex=.8)
    } 
    if(showMaxScale)  {
      par(xpd=TRUE)
      maxScaleText <- rep("", length(maxScale.ints))
      maxScaleText[1] <- "0"
      maxScaleText[length(maxScale.ints)] <- "Max"
      legend(1.1, max(y.corners)+y.adjust, legend=rev(maxScaleText), fill=rev(colorRampPalette(c("white", "blue"))( length(maxScale.ints)  )), title="Row scale", xjust=0, yjust=0)
      par(xpd=FALSE)
    }
  } else {
    warning("No data to plot") 
  }
  
}


#data.subset <- na.omit(gene.reads.subset[TF.vec ,samples.ordered ] )
#colnames(data.subset) <- sampleLabels[colnames(data.subset)]  # need to convert colnames to standard labels

#maxScaledHeatmap(data.subset, minMax=100)
#maxScaledHeatmap(data.subset, minMax=100, , showMaxScale=TRUE,  showIntenseScale=TRUE)

