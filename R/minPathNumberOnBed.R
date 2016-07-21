#!/usr/bin/Rscript

#in future use optparse library
Args <- commandArgs(TRUE)

if(length(Args) < 2)  {
        stop("Usage: minPathNumberOnBed.R [features bed file] [junctions bed file] [out file (optional)]\n
		For each feature in the feature file,  finds the overlapping set of junctions from the junctions file. 
		Then runs minPathNumber() over the subset of junction to calculate minimum number of separate transcripts
		required to explain the set of junctions.
		Output to STDOUT if no third argument given.  
		
		Requires that correct paths to two other R scripts are set within this script.
        ")

}
stdout <- FALSE

featuresFile <- Args[1]
junctionsFile <- Args[2]
if(is.na(Args[3]))  {
	stdout <- TRUE
} else {
	outputFile <- Args[3]
}


# TESTING # featuresFile <- "C:/Users/Dave/HanleyGroup/JenningsRNAseq/unexRegions/jenningsANDxie_UNEX_UNSTRANDED.merged.1000.intergenic.bed"
# TESTING # junctionsFile <- "C:/Temp/sampleJuncs.bed"

#outputFile 

source("C:/Users/Dave/utilsGerrardDT/minPathNumber.R")	# on PC
source("C:/Users/Dave/utilsGerrardDT/mapNearPeaksGenes.R")	# on PC

#source("/home/mqbssdgb/bin/minPathNumber.R")     # on KADMON
#source("/home/mqbssdgb/bin/mapNearPeaksGenes.R")     # on KADMON


### TODO allow strand specificity.

bedColNames <- c("chr", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")

features <- read.delim(featuresFile, header=F)
names(features) <- bedColNames[1:ncol(features)]
features <- features[order(features$chr, features$start), ]
#head(features)

# TESTING # features <- features[1:100,]
# TESTING # features <- data.frame(chr="chr1", start=c(310000, 320000, 665000, 710000, 877000), end=c(311000, 323000, 670000, 740000, 879000))

junctions <- read.delim(junctionsFile, header=F)
names(junctions) <- bedColNames[1:ncol(junctions )]
junctions <- junctions[order(junctions$chr, junctions$start), ]
#head(junctions)



global.results <- data.frame()

for(thisChrom in unique(features$chr))  {

	features.sub <- subset(features, chr == thisChrom)
	junctions.sub <- subset(junctions, chr == thisChrom) 

	chrom.results <- features.sub
	chrom.results$n.unique.juncs <- NA
	chrom.results$min.TC <- NA
	for(i in 1:nrow(features.sub))  {
		# for each region of interest, extract junctions from the junctions file 
		region.start <- features.sub$start[i]
		region.size <- features.sub$end[i] - features.sub$start[i]
		region.juncs <- junctions.sub[intersect( findPointFeaturesInWindow(region.start, junctions.sub$start, window_size=region.size, window_type="right"), findPointFeaturesInWindow(region.start, junctions.sub$end, window_size=region.size, window_type="right")), ]

		n.unique.juncs <-  nrow(unique(region.juncs[,c("chr", "start", "end")]))

		#! what about strand? 
		chrom.results$n.unique.juncs[i] <- n.unique.juncs 


		# calculate min number of transcripts (minPathLength)  'tc'
		if(n.unique.juncs  > 0 )  {
			chrom.results$min.TC[i] <- minPathNumber(region.juncs)
		}
		
		# out put feature (or just feature names) with tc number
		
	}
	global.results <- rbind(global.results, chrom.results)
}


if(stdout)  {
	write.table(global.results, file="", sep="\t", quote=F, row.names=F)

} else {
	write.table(global.results, file=outputFile, sep="\t", quote=F, row.names=F)

}



