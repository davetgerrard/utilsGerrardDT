

### Chain regions of genome based on similarity score (expression profiles, histone mods (not sequence))


# rules on chaining:  reciprocal = all region-pairs must pass threshold
#                     sequential = pairs are included if they pass threshold at least once
#   distance:    
#   same.chrom =TRUE  (regions must be on same chromosome.)

# regions   a set of regions in GR format (fields used: chr, start, strand)
# scores: a distance matrix 
# cut.method: how to determine the similarity cut-off to form chains.  
# cut.prop:  0.2
chain.regions <- function(regions, scores, dist.m="euclidean", hclust.m="complete", same.strand=FALSE, chain.method="reciprocal", same.chrom=TRUE, 
                          max.distance=-1, cut.prop=0.1, cut.method="global")  {
  
  # checks and further params
  if(nrow(regions) != nrow(scores))  stop("lengths of regions and scores do not match")
  cut.level <- max(dist(scores)) * cut.prop
  
  
  #d.m <- cor(scores) # calculate distance between regions
  
  #d.d <- dist(d.m)  # calculate distance
  
  groups <- cutree(hclust(dist(scores, method=dist.m), method=hclust.m), h=cut.level)
  
  return(groups)
  
}



region.distance <- function(bed.data, method="midpoint")  {
  midpoints <-  bed.data$start + floor((bed.data$end - bed.data$start) / 2)
  
}



# should the main function have to calculate scores? or 
# can save a lot of compute by only calculating distance for valid groups (eg. within chromosome). 

# 
# ### SET UP TEST DATA
# 
# t.r <- data.frame(chr=c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2"),
#                   start=c(100, 200, 300, 400, 500, 600, 700, 800),
#                   end=c(150, 250, 350, 450, 550, 650, 750, 850),
#                   strand=c('+', '+', '-', '+', '-', '+', '+', '+') )
# t.r$name <- paste(t.r$chr,"-",  t.r$start, "_" ,t.r$end, sep="")
# nrow(t.r)
# 
# pattern.1 <- c(50, 100, 10, 0, 200, 50)
# names(pattern.1) <- paste("sample", 1:6, sep=".")
# pattern.2 <- pattern.1 + 1:6     # similar array
# pattern.3 <- pattern.1 + (1:6)*10     # 
# 
# # need to rbind 8 patterns including some very similar ones.
# 
# t.s.data <- as.data.frame(rbind(pattern.1, pattern.2, pattern.3, rev(pattern.1), sort(pattern.1), rev(pattern.3), sort(pattern.3) , rev(pattern.2) ), row.names=t.r$name)
# 
# cor(t.s.data)
# 
# 
# 
# #### RUN THE CHAINING
# 
# r.chains <- chain.regions(t.r, t.s.data)





