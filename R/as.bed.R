
#http://genome.ucsc.edu/FAQ/FAQformat.html#format1.7
#               strand - Defines the strand - either '+' or '-'.
#               thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
#   thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
#   itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
#   blockCount - The number of blocks (exons) in the BED line.
#   blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
#   blockStarts - 

# give bed format column names to bed table lacking a header.
as.bed <- function(data, bed.col.names = c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",  "blockStarts") ) {
    stopifnot(is.data.frame(data))
    n.col <- ncol(data)
    names(data)[1:n.col]  <- bed.col.names[1:n.col]
    return(data)
    
    
}
                   
