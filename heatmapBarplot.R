

# create an image from a table and show row and column sums (means etc) as histograms at the sides

# optionally, scale the image to rowsums/colsums



# rowBar.min  0, NULL

heatmapBarplot <- function(mat, rowBar.min=0, colBar.min=0, grid.prop.x=.7, grid.prop.y=grid.prop.x, spacer=.05, 
                           useLog10=FALSE, row.labels=dimnames(mat)[[1]],  col.labels=dimnames(mat)[[2]], 
                           xlab="", ylab="",
                           r.vals = apply(mat, 1, sum), c.vals = apply(mat, 2, sum))  {

  
  
  grid.min.x <- .1
  grid.min.y <- .1
    
  grid.max.x <- grid.min.x + grid.prop.x
  grid.max.y <- grid.min.y + grid.prop.y
  
  cell.width <- (grid.max.x - grid.min.x) / dim(mat)[1]
  cell.height <- (grid.max.y - grid.min.y) / dim(mat)[2]
  
  grid.left <- seq(from=grid.min.x, by = cell.width, length.out = dim(mat)[1])
  grid.right <- seq(to=grid.max.x, by = cell.width, length.out = dim(mat)[1])
  grid.horMid <- grid.left + (cell.width/2)
  
  grid.bottom <- seq(from=grid.min.y, by = cell.height, length.out = dim(mat)[2])
  grid.top <- seq(to=grid.max.y, by = cell.height, length.out = dim(mat)[2])
  grid.verMid <- grid.bottom + (cell.height/2)
  
  
  if(useLog10) {
    mat.scale <- log10(mat)    # may need to deal with 0s and negs.
  }
  else {
    mat.scale <- mat
  }
  mat.prop <- mat.scale/max(mat.scale)
  
  
  plot.new()

  palette(colorRampPalette(c("blue","white", "red"))( 30 ))
  rect(xleft = rep(grid.left, each=dim(mat)[2]), 
       ybottom = rep(grid.bottom, times=dim(mat)[1]), 
       xright = rep(grid.right, each=dim(mat)[2]), 
       ytop = rep(grid.top, times=dim(mat)[1]), col = mat.prop  * length(palette()))


  # add labels
  
  text(x = grid.horMid, y= 0.05, labels=col.labels)
  text(x = mean(c(grid.min.x,grid.max.x)), y=0.02, xlab)
  text(x = 0.05, y = grid.verMid, labels=row.labels)
  text(x = 0.02, y= mean(c(grid.min.x,grid.max.x)), ylab, srt=90)
  
  # row summarising bar plot

  r.bar.min.x <- grid.max.x + spacer
  r.bar.max.x <- 1.0
  
  #r.bar.min.y <- grid.min.y 
  #r.bar.max.y <- grid.max.y
  
  
  #r.vals <- apply(mat, 1, sum)   # can change this function later as a parameter
  # calculate the 
  
  if(useLog10) r.vals <- log10(r.vals)
  
  if(is.null(rowBar.min)) {
    r.vals.scale <-  pretty(r.vals)
  } else {
    r.vals.scale <- pretty(c(rowBar.min,max(r.vals)))
  }
  # calculate the horizontal height of the row summarizing values on the scale between r.bar.min.x and r.bar.max.x
  r.vals.right <- (( (r.vals -  r.vals.scale[1])  / (max(r.vals.scale) - min(r.vals.scale)) )   *  (r.bar.max.x - r.bar.min.x)) + r.bar.min.x 
  
  rect(xleft=r.bar.min.x, ybottom=grid.bottom, xright=r.vals.right, ytop = grid.top)
  #rect(xleft=r.bar.min.x, ybottom=rev(grid.top), xright=r.vals.right, ytop = rev(grid.bottom))  # to start row at top, would need compensatory change in colours

  
  # column summarising bar plot
  
  c.bar.min.y <- grid.max.y + spacer
  c.bar.max.y <- 1.0
  
  #c.bar.min.x <- grid.min.x 
  #c.bar.max.x <- grid.max.x
  
  #c.vals <- apply(mat, 2, sum)   # can change this function later as a parameter
  if(useLog10) c.vals <- log10(c.vals)
  # calculate the 
  
 if(is.null(colBar.min)) {
   c.vals.scale <-  pretty(c.vals)
 } else {
   c.vals.scale <- pretty(c(rowBar.min,max(c.vals)))
 }
  # calculate the vertical height of the column summarizing values on the scale between c.bar.min.x and c.bar.max.x
  c.vals.top <- (( (c.vals -  c.vals.scale[1])  / (max(c.vals.scale) - min(c.vals.scale)) )   *  (c.bar.max.y - c.bar.min.y)) + c.bar.min.y 
  
  rect(xleft=grid.left, ybottom=c.bar.min.y, xright=grid.right, ytop = c.vals.top) 
  
 # add an axis for the row summary bars
  segments(x0=r.bar.min.x, y0= .05, x1=r.bar.max.x, y1=.05)
  segments(x0=c(r.bar.min.x,r.bar.max.x), y0= 0.045, x1= c(r.bar.min.x,r.bar.max.x),  y1=.05)
   text(x= c(r.bar.min.x, r.bar.max.x) , y= .05, labels=r.vals.scale[c(1,length(r.vals.scale))], pos=1)
 
 
  # add an axis for the column summary bars
  segments(x0=.05, y0= c.bar.min.y, x1= .05,  y1=c.bar.max.y)
  segments(x0=.045, y0= c(c.bar.min.y,c.bar.max.y), x1= .05,  y1=c(c.bar.min.y,c.bar.max.y))
  text(x= 0.05 , y= c(c.bar.min.y, c.bar.max.y), labels=c.vals.scale[c(1,length(c.vals.scale))], pos=2)
 
}

#par(mar=c(2,2,2,2))

#test_data <- matrix(rnorm(100, mean=100, sd=10), nrow = 10, byrow=T)
#heatmapBarplot(test_data)


#test_data <- matrix(1:100, nrow = 10, byrow=T)
#heatmapBarplot(test_data)

# TODO
# implement row and column scaling (and both simultaneously), keep barplot on real (or log10 scaled) values.
# investigate failed plots.





