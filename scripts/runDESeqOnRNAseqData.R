#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################


##### Only run manually so far. Might work as an Rscript.

##### SET UP ENVIRONMENT AND PARAMETERS

old.o <- options(scipen=999) # required to stop genome co-ordinats being converted and written as scientific notation

# 
inputFile <- "C:/Users/dave/BioinfCore/saunders/Partek_Transcripts_UCSC_annot_0711.txt"
#getwd()
setwd("C:/Users/dave/BioinfCore/saunders/")
outputFile <- "DESeq_on_transcriptCounts_CvsE.tab"

############LOAD PACKAGES
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("DESeq")
library(DESeq)


##########LOAD DATA
count_data <- read.delim(inputFile ,header=T)
head(count_data)


###########INSPECT DATA
#hist(log(count_data$C1_tophat_dm3_chr_hetMU..Reads.), breaks=200)


################ CREATE DATA TABLE WITH COUNTS ONLY

countsTable <- subset(count_data, select=c( 
				"C1_tophat_dm3_chr_hetMU..Reads.",
				 "C2_tophat_dm3_chr_hetMU..Reads.",
				"E1_tophat_dm3_chr_hetMU..Reads.",
				"E2_tophat_dm3_chr_hetMU..Reads."))

row.names(countsTable) <- count_data$Transcript
head(countsTable)
countsTable <- round(countsTable)  # must be integer

conds <- factor( c( "C_type", "C_type", "E_type", "E_type" ) )


###########RUN DESEQ
cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )
#sizeFactors( cds )
#cds <- estimateDispersions( cds )
cds <- estimateVarianceFunctions( cds)
res <- nbinomTest( cds, "C_type", "E_type" )


############### CREATE OUTPUT TABLE AND WRITE TO FILE

names(res) <- paste("DESeq", names(res), sep=".")
output <- merge(count_data, res, by.x="Transcript", by.y="DESeq.id")

write.table(output[order(output$DESeq.padj),], file=outputFile, sep="\t",quote=F,row.names=F)


################ CLEAN UP AND LEAVE
historyFile <- paste(outputFile, "RHistory", sep=".")
options(old.o)
savehistory(file=historyFile)
quit(save="no")

###########################################

stopifnot(FALSE)
