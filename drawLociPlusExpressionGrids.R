
scaleCoords <- function(x, input.range=range(x), output.range=c(0,1))  {
  input.interval <- input.range[2] - input.range[1]
  output.interval <- output.range[2] - output.range[1]
  x.prop <- (x - input.range[1]) / input.interval
  out.x <- (x.prop * output.interval)  + output.range[1]
  return(out.x)
}

#scaleCoords(c(400, 0, 45))

drawLoci <- function(x, xlim=wideLimits(c(min(x$start), max(x$end)), ...), chrom=x$chr[1],  ...)  {
  x$midpoint <- x$start + round((x$end - x$start)/2)
  
  plot.units<- 1/ (xlim[2] - xlim[1])
  plot.new()
  #segments(xlim[1]*plot.units, .5, xlim[2]*plot.units, .5)
  segments(.1, .2, .9, .2)
  segments(c(.1,.9), .18, c(.1,.9), .2)
  text(c(.1, .9), c(.1), xlim)
  text(.1, .35, chrom, cex=1.2)
  
  x$propStart <- scaleCoords(x$start, input.range=xlim, output.range=c(.1, .9))
  x$propEnd <- scaleCoords(x$end, input.range=xlim, output.range=c(.1, .9))
  x$propMidpoint <- scaleCoords(x$midpoint, input.range=xlim, output.range=c(.1, .9))  
  x$arrowStart <- ifelse(x$strand=="+", x$propStart, x$propEnd)
  x$arrowEnd <- x$propMidpoint
  
  y.heights <- rep_len(c(.4, .5, .6), nrow(x))
  
  segments(x$propStart,  y.heights,  x$propEnd , y.heights)    
  #plot arrow at midpoint of each feature different direction for +/- strands
  
  arrows(x$arrowStart, y.heights, x$arrowEnd, y.heights, length=.05)
  
  
  text(x$propMidpoint, y.heights, row.names(x), pos=3, cex=.8)                 
  
}




# the idea is to cleverly space the positions of a set of
# values into individual bins with clumping linked to the original data.
# Returns a vector of NA  with integers at some positions which index the original vector
assignToBins <- function(x, n.bins = length(x)*2)  {
  bin.limits <- seq(min(x), max(x), length.out=n.bins)
  index <- as.integer(rep(NA, n.bins))
  x.index <- order(x)  # ties are sorted on index, should be fine
  i <- 1
  
  for(j in 1:n.bins)  {
    rem <- length(x) - i
    if(rem < 0)  {break}
    k <- n.bins - j    # the number of spaces left
    x.next <- if(i >= length(x))  {x[x.index[i]]} else {x[x.index[i+1]]}
    if(k <= rem  | (x[x.index[i]] <= bin.limits[j])  | (x.next <= bin.limits[j])  )  {
      index[j]  <- x.index[i]
      i <- i + 1 
    }      
  }
  return(index)
}

#test.vec <- c(3,4,5,16,20, 3,4)
#assignToBins(test.vec, n.bins=20)
#test.vec[assignToBins(test.vec, n.bins=20)]

# draw a set of loci with matched floating columns of expressino above them.
drawLociPlusExpressionGrids <- function(x, value.table,  minMax = 0, mar.vec=c(0,0,0,0),
                                        text.edge=.2, text.cex=1, box.top = .95, box.bottom = .3, max.box.radius=.01,
                                        column.left = .3, column.right = .9,
                                        max.column.radius = .01, column.res=50, xlim=wideLimits(c(min(x$start), max(x$end)), ...), chrom=x$chr[1],
                                        chrom.line=.1, line.start=.1, line.end=.9, loci.bottom = chrom.line +.025, loci.cex =.8* text.cex,
                                        loci.top = box.bottom - .1, loci.lines =3, clever.spacing=TRUE, intensity.scale= FALSE, 
                                        scale.colour="forestgreen", scale.prop=.15,draw.connectors=FALSE, ...)  {
  
  
  if(intensity.scale)  {
    scale.top <- box.top
    box.top <- box.bottom + ((scale.top - box.bottom) * (1-scale.prop))
    scale.mid <- (scale.top + box.top)/2
  }
  
  
  df <- value.table
  max.vec <- apply(df, 1, max)
  exclude.index <- max.vec < minMax
  df[exclude.index,] <- 0
  
  # image(scale(t(df), center=F, scale=apply(df, 1, max)), axes=F, mar=c(10,7),col=colorRampPalette(c("white", "blue"))( 30 ), add=F)
  #   image(t(as.matrix(max.vec)), zlim=c(0, max(max.vec)), add=T,col=colorRampPalette(c("white", scale.colour))( 30), x=c(-.2, -.1), y=seq(0,1, length.out=length(max.vec)))
  
  # do the same using rect?
  scale.df <- scale(t(df), center=F, scale=apply(df, 1, max))
  scale.df[is.na(scale.df)] <- 0.0    # do not want NAs resulting from dividing by zero when max = 0
  
  x$midpoint <- x$start + round((x$end - x$start)/2)
  x <- x[order(x$midpoint),]   # this probably makes feature.order redundant
  # need to order the loci.table and then use this to order the expression table
  scale.df <- scale.df[,row.names(x)]
  max.vec <- apply(df[row.names(x),],1, max)
  
  feature.order <- order(x$midpoint)
  
  # how many boxes within each column
  box.n <- nrow(scale.df)
  box.mids <- seq(box.bottom, box.top, length.out=box.n)
  box.radius <- min(max.box.radius, (box.top - box.bottom) / (box.n *2))
  
  # how many columns
  # calculate where to put different columns
  # don't fill up the whole space
  column.n <- ncol(scale.df)
  column.radius <- min(max.column.radius, ((column.right - column.left) / (column.n + 1)))
  #if(clever.spacing) {
  n.bin <- max(column.n, column.res)
  which.bin <- match(1:column.n, assignToBins(x$midpoint, n.bins=n.bin))
  empty.mids <- seq(column.left, column.right, length.out=n.bin )
  column.mids <- as.numeric(na.omit(empty.mids[which.bin]))
  #} else {  #space evenly
  #  column.mids <- seq(column.left, column.right, length.out=column.n )   # not working
  #} 
  # begin plotting
  par(mar=mar.vec)  # need to have plenty of space to plot names.
  plot.new()
  # next bit seems crazy, but struggled to find way to get the actual colour values for each box.
  s.colours <- apply(colorRamp(c("white", scale.colour))( max.vec /max(max.vec))/255, 1 , FUN=function(x) rgb(x[1], x[2], x[3]))
  
  for(j in 1:length(feature.order))  {
    i <- feature.order[j]
    thisFeature <- row.names(x)[i]
    i.colours <- apply(colorRamp(c("white", "blue"))( scale.df[,i] )/255, 1 , FUN=function(x) rgb(x[1], x[2], x[3]))
    
    rect(xleft=column.mids[i]-column.radius, xright=column.mids[i]+column.radius, ybottom=box.mids-box.radius, ytop =box.mids+box.radius, col=i.colours)
    if(intensity.scale)  {
      rect(xleft=column.mids[i]-column.radius, xright=column.mids[i]+column.radius, ybottom=scale.mid-box.radius, ytop =scale.mid+box.radius, col=s.colours[i])
    }
  }
  #legend("center", row.names(scale.df), fill=i.colours)  # works well for one sample
  
  
  
  text(text.edge, box.mids, labels=row.names(scale.df), adj=1, cex=text.cex)
  if(intensity.scale)  text(text.edge, scale.mid, labels="intensity", adj=1, cex=text.cex)
  #segments(.1, box.mids, .5, box.mids)  lines from names to boxes not good.
  
  # draw the chromosome scale bar 
  segments(line.start, chrom.line, line.end, chrom.line)
  segments(c(line.start,line.end), chrom.line-.02, c(line.start,line.end), chrom.line)
  text(c(line.start, line.end), c(chrom.line-.05), prettyNum(xlim,big.mark=",", preserve.width="none"))
  text(line.start, chrom.line+.025, chrom, cex=1.2)
  
  # draw the individual loci
  x$propStart <- scaleCoords(x$start, input.range=xlim, output.range=c(line.start, line.end))
  x$propEnd <- scaleCoords(x$end, input.range=xlim, output.range=c(line.start, line.end))
  x$propMidpoint <- scaleCoords(x$midpoint, input.range=xlim, output.range=c(line.start, line.end))  
  x$arrowStart <- ifelse(x$strand=="+", x$propStart, x$propEnd)
  x$arrowEnd <- x$propMidpoint
  
  y.heights <- rep_len(seq(loci.bottom, loci.top, length.out=loci.lines), nrow(x))
  
  segments(x$propStart,  y.heights,  x$propEnd , y.heights)    
  #plot arrow at midpoint of each feature different direction for +/- strands
  
  arrows(x$arrowStart, y.heights, x$arrowEnd, y.heights, length=.05)
  
  text(x$propMidpoint, y.heights, row.names(x), pos=3, cex=loci.cex)  
  
  if(draw.connectors)  {    # draw lines from the bottom of each column to a position above the midpoint of the corresponding gene
    connect.bottom <- loci.top + ((box.bottom - loci.top)/2) 
    bevel.pos  <- connect.bottom + ((box.bottom - connect.bottom)/5)     # used to stragten up the connector at the end
    segments(x$propMidpoint, bevel.pos, column.mids,box.bottom )    
    segments(x$propMidpoint, connect.bottom, x$propMidpoint,bevel.pos )    
  }
  
}



test.drawLoci <- function()  {
  n.loci <- 10
  starts <- floor(runif(n.loci,min=500000, max=600000))
  ends <- starts + floor(rnorm(n.loci, mean=5000, sd=500))
  strands <- sample(c("+", "-") , n.loci, replace=T)
  loci.table <- data.frame(chr="chrX", start=starts, end=ends, strand=strands, row.names=LETTERS[1:n.loci])
  loci.table$tss <- ifelse(loci.table$strand == "+", loci.table$start, loci.table$end)
  
  drawLoci(loci.table)

}
# test.drawLoci()



test.drawLociPlusExpressionGrids  <- function(...)  {
  n.loci <- 10
  starts <- floor(runif(n.loci,min=500000, max=600000))
  ends <- starts + floor(rnorm(n.loci, mean=5000, sd=500))
  strands <- sample(c("+", "-") , n.loci, replace=T)
  loci.table <- data.frame(chr="chrX", start=starts, end=ends, strand=strands, row.names=LETTERS[1:n.loci])
  loci.table$tss <- ifelse(loci.table$strand == "+", loci.table$start, loci.table$end)
  
  
  samples <- rev(paste("sample", 1:10, sep="."))
  
  exp.table <- matrix(runif( length(samples)*n.loci, 0,10), nrow=n.loci, dimnames=list(row.names(loci.table),samples))
  
  drawLociPlusExpressionGrids(x= loci.table, value.table=exp.table, ...)
}
#test.drawLociPlusExpressionGrids() 
# test.drawLociPlusExpressionGrids(text.edge=.1, column.left = .15,box.bottom = .4 , column.res=100, max.column.radius = .01) 
# test.drawLociPlusExpressionGrids(text.edge=.1, column.left = .15,box.bottom = .4 , column.res=100, max.column.radius = .01) 


