
# write old 
dropMeta <- function(df, oldFile, droppedFile="META_CHANGES.LOG", sep="")  {
	if(file.exists(droppedFile))  {
		df[,'META_FILE']  <- oldFile
		write.table(df, file=droppedFile, append=T)
	}  else  {
		# set up file with some annotation
		#TODO
		write("# Changelog of all .META files in this directory", file=droppedFile)
		write.table(df, file=droppedFile, append=T)
	}

}


# template (if there are other meta files with the same descriptions that could be used). TODO
metaTable <- function(baseFile, headerLinePos=1, sep=" ", templateFiles="")  {
	# test if file exists
	if(!file.exists(baseFile))  stop(paste("Cannot find data file:", baseFile))

	# read column names from file. 
	col.names <- as.character(unlist(read.delim(baseFile, nrows=1, sep="\t", header=F, skip = headerLinePos-1)[1,] ))

	# test if metadata file exists
	# if it does, load it because it likely contains 
	# descriptions that have been manually added and we don't want to lose.
	metaDataFileName <- paste(baseFile, "META", sep=".")
	if(file.exists(metaDataFileName))  {
		orig.metaData <- read.delim(metaDataFileName, sep=sep, header=T)
		row.names(orig.metaData) <- orig.metaData[,1]	# test for unique or let R fail here anyway?
		# test if all of col.names are accounted for.
		shared <- intersect(orig.metaData[,1], col.names)
		old.unique <- setdiff(orig.metaData[,1], col.names)	
		new.unique <- setdiff(col.names, orig.metaData[,1])	
		new.metaData <- orig.metaData
		if(length(old.unique) > 0)  {
			# some old names are not in new file
			# backupOldFile? Start or add to a META.CHANGELOG file.
			dropMeta(orig.metaData[old.unique,], oldFile=metaDataFileName, sep=sep)
			print(paste(length(new.unique)  ,"old annotations removed to", "META_CHANGES.LOG" ))
			new.metaData <- orig.metaData[shared,]
		} 
		if(length(new.unique) > 0) {
			# some new names are not in old file
			# keep old description that are still in new set
			print(paste(length(new.unique), "new annotations added"  ))
			newRows <- data.frame(META.NAMES = new.unique, META.DESCRIPTION=NA)
			new.metaData <- rbind(new.metaData, newRows)

		} 
	}   else  {
		# if metadata file does not exist, make a simple one and tell the user.	
		print(paste("Creating new meta data file:", metaDataFileName))
		new.metaData <- data.frame(META.NAMES = col.names, META.DESCRIPTION=NA)
	}
	# write out 
	print("Writing up-dated meta-data")
	write.table(new.metaData, file=metaDataFileName, sep=sep, quote=F, row.names=F)
}


# examples
#metaTable("blob")
#metaTable("UNEX_XJ_1000bp.lincRNA.overlaps.txt", sep="\t")
#metaTable("UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab", sep="\t")



