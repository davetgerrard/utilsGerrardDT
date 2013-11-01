
# write old 
dropMeta <- function(df, oldFile, droppedFile="META_CHANGES.LOG", sep="")  {
	df[,'META_FILE']  <- oldFile
	df[,'MOD_TIME']  <- date()
	colOrder <- c('MOD_TIME','META_FILE', setdiff(names(df), c('MOD_TIME','META_FILE')))
	df <- df[, colOrder]
	if(file.exists(droppedFile))  {		
		write.table(df, file=droppedFile, append=T, row.names=F, col.names=F, sep=sep)
	}  else  {
		# set up file with some annotation
		#TODO
		write("# Changelog of all .META files in this directory", file=droppedFile)
		write.table(df, file=droppedFile, append=T, row.names=F, col.names=F, sep=sep)
	}

}


# template (if there is another meta file with the same descriptions that could be used). TODO
metaTable <- function(baseFile, headerLinePos=1, sep=" ", templateFile="")  {
	# test if file exists
	if(!file.exists(baseFile))  stop(paste("Cannot find data file:", baseFile))


	# read column names from file. 
	#col.names <- as.character(unlist(read.delim(baseFile, nrows=1, sep=sep, header=F, skip = headerLinePos-1)[1,] ))
	col.names <- scan(baseFile, what="character", nlines = 1, skip = headerLinePos-1, sep=sep, quiet=T)
	col.names <- col.names[col.names != ""]     # don't want blank entries
 	# test if metadata file exists
	# if it does, load it because it likely contains 
	# descriptions that have been manually added and we don't want to lose.
	metaDataFileName <- paste(baseFile, "META", sep=".")

	if(templateFile != "") {
		if(!file.exists(metaDataFileName))  {
			print(paste("Using template file",templateFile ) )
			file.copy(templateFile, metaDataFileName)
		} else {
			print(paste("Template file",templateFile , "ignored because",metaDataFileName, "already exists") )
		}
	}


	if(file.exists(metaDataFileName))  {
		orig.metaData <- read.delim(metaDataFileName, sep=sep, header=T)
		row.names(orig.metaData) <- orig.metaData[,'META.NAMES']	# test for unique or let R fail here anyway?
		# test if all of col.names are accounted for.
		shared <- intersect(orig.metaData[,'META.NAMES'], col.names)
		old.unique <- setdiff(orig.metaData[,'META.NAMES'], col.names)	
		new.unique <- setdiff(col.names, orig.metaData[,'META.NAMES'])	
		new.metaData <- orig.metaData
		if(length(old.unique) > 0)  {
			# some old names are not in new file
			# backupOldFile? Start or add to a META.CHANGELOG file.
			dropMeta(orig.metaData[old.unique,], oldFile=metaDataFileName, sep=sep)
			print(paste(length(old.unique)  ,"old annotations removed to", "META_CHANGES.LOG" ))
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


stopifnot(FALSE)

# examples
#metaTable("blob")
#metaTable("UNEX_XJ_1000bp.lincRNA.overlaps.txt", sep="\t")
#metaTable("UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab", sep="\t")
#metaTable("UNEX_XJ_1000bp.lincRNA.overlaps.txt", templateFile="UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab.META", sep="\t")
#metaTable("UNEX_XJ.validCounts.intergenic.1000bpMerge.XieJennings.G15.min200bp20reads.contigTest.tab", templateFile="UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab.META", sep="\t")
# following gave problems.
#metaTable("UNEX_XJ.validCounts.intergenic.1000bpMerge.XieJennings.G15.tab", templateFile="UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab.META", sep="\t")

#metaTable("UNEX_XJ.validCounts.intergenic.1000bpMerge.XieJennings.G15.min200bp20reads.tab", templateFile="UNEX_XJ.validCounts.intergenic.1000bpMerge.G15.min200bp20reads.extra.tab.META", sep="\t")








