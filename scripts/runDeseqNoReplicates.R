#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################


########FUNCTIONS


###########################################

# plan to be able to call this script thus: Rscript script.R arg1 arg2 ..
# get args from call
Args <- commandArgs(TRUE)

if(length(Args) != 3)  { 
	stop("Usage: runDeseqNoReplicates.R [bedCounts1] [bedCounts2] [outputFile]\n
		Runs DESEQ on two bedCount files (chr	start	end	count) without headers.
		Outputs regions and p-values. 
	")
	
}

inFile.1 <- Args[1]
inFile.2 <- Args[2]
outputFile <- Args[3]
#sigFilter <- 0.1   # not actually used in this script


# source("http://www.bioconductor.org/biocLite.R")
# biocLite("DESeq")
library("DESeq")
old.o <- options(scipen=999) # required to stop genome co-ordinats being converted and written as scientific notation

#dir()

simpleBins.1 <- read.delim(inFile.1, sep="\t",header=F)
simpleBins.2 <- read.delim(inFile.2, sep="\t",header=F)
#should include identical() test here.
stopifnot(identical(simpleBins.1[1:3], simpleBins.2[1:3])
countsTable <- as.data.frame(cbind(simpleBins.1$V4,simpleBins.2$V4))
#head(countsTable)
conds <- c("Cond.1","Cond.2")
names(countsTable) <- conds
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
cds <- estimateVarianceFunctions( cds, pool=T )
res <- nbinomTest( cds, "Cond.1", "Cond.2" )
#resSig <- res[ res$padj < sigFilter, ]   # this currently does nothing, could be used to filter the results before output.
#head( res[ order(res$pval), ] )
islands <- simpleBins.1[,1:3]
names(islands) <- c("chr","start","end")
output <- cbind(islands, countsTable,res)
#head(output)
#head(output[,order(output$pval)])
#head(output[order(output$pval),])

write.table(output[order(output$pval),], file=outputFile, sep="\t",quote=F,row.names=F)
options(old.o)
savehistory("runDeseqNoReplicates.RHistory")
#quit(save="no")

