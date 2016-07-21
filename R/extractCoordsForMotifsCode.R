# data 	must have: chr start end name strand 
# nameColumn	name of the column with the feature name
# add.upstream
# add.downstream
# offset 	move the coords by this amount. -ve for upstream, +ve for downstream relative to strand.
# reference
# delim
extractCoords <- function(data, nameColumn=NA, add.upstream = 500, add.downstream = 0, offset = 0,reference = "tss", delim="\t")  {
	result.chr <- data$chr
	result.name <- paste(data[,nameColumn], add.upstream, reference, add.downstream,"offset",offset,  sep=".")
	
	data$start <- ifelse(data$strand=='-', data$start - offset, data$start + offset )
	data$end <- ifelse(data$strand=='-', data$end - offset, data$end + offset )

	if(reference=="tss")  {	# coords about a specific point
		result.start <- ifelse(data$strand=='-', data$end - add.downstream, data$start - add.upstream)
		result.end <- ifelse(data$strand=='-', data$end + add.upstream, data$start + add.downstream)
	} 
	if(reference=="gene")  { # coords extended around the gene.
		#result.start <- data$start - add.upstream
		#result.end <-  data$end + add.downstream
		result.start <- ifelse(data$strand=='-', data$start - add.downstream, data$start - add.upstream)
		result.end <- ifelse(data$strand=='-', data$end + add.upstream, data$end + add.downstream)
	}



	result.start <- ifelse(result.start < 0 ,0, result.start)
	result.end <-  ifelse(result.end <= 0, 1, result.end)

	#stopifnot(result.end >= result.start)

	#return(paste(result.chr,result.start,result.end,result.name,"\n",sep=delim))
	#return(c(result.chr,result.start,result.end,result.name))
	return(data.frame(chr=result.chr,start=result.start,end=result.end,name=result.name))
}

