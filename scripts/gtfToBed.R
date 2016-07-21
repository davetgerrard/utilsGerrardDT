
#gtf to bed converter  DAVE GERRARD, 2014.      
# simple one-to-one.
## FURTHER WORK: user selects which features to export
##          Can export bed12 of combined features. (e.g. group exons into transcripts).



require(getopt)

spec = matrix(c(
  'verbose', 'v', 2, "integer",
  'help'   , 'h', 0, "logical",
  'inFile'  , 'i', 1, "character",
  'outFile'  , 'o', 2, "character",
  'cLociFile' , 'l', 2, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);
#0=no argument, 1=required argument, 2=optional argument


# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$outFile    ) ) { 
  opt$outFile    <-  paste( dirname(opt$inFile), "/", sub("gtf$",  "bed", basename(opt$inFile), ignore.case=TRUE), sep="")    
}


gtf.file <- "C:/Users/Dave/HanleyGroup/BerryRNAseq_multiTissue/data/cufflinks/cuffComp_0.10/cuffC_all_0.10.class.u.gtf"
gtf.file <- opt$inFile

gtf.data <- read.delim(gtf.file, header=F)
#http://genome.ucsc.edu/FAQ/FAQformat.html#format3
names(gtf.data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "group")


want_list <- c("transcript_id",  "exon_number")

#resultsList <- list()

for(thisPattern in want_list) {
  pattern <- paste(".*",thisPattern," ([^;]+);.*",sep="")
  gtf.data[,thisPattern] <- sub(pattern, "\\1", gtf.data[grepl(pattern, gtf.data[,'group'], perl=T),'group'], ignore.case=T)
  
}

head(gtf.data )

gtf.data$name <- paste(gtf.data[,"transcript_id"], gtf.data[,"exon_number"], sep=":")

bed.data <- gtf.data[, c("seqname", "start", "end", "name", "score", "strand")]
bed.data$start <- bed.data$start - 1    # bed uses 0-based, gtf uses 1-based.


write.table(bed.data, opt$outFile, sep="\t", quote=F, row.names=F, col.names=F)

