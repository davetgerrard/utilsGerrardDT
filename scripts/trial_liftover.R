
setwd("C:/Temp")

### trial using liftover in rtracklayer package

require(rtracklayer)

chain.file.url <- "http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz"
chain.file.local.gz <- "C:/Users/Dave/data/UCSC/hg18ToHg19.over.chain.gz"
chain.file.local <- "C:/Users/Dave/data/UCSC/hg18ToHg19.over.chain"

## try separate download and open
#download.file(url=chain.file.url, destfile=basename(chain.file.url ), mode="wb")

chain=import.chain(chain.file.local)
# import bed, use as.bed() to add headers
bed.file.hg18 <- as.bed( read.delim("C:/Temp/Pasquali2014_Panc/FOXA2_HI_101_1e-10_peaks.bed", header=F))

bed.file.hg18.GR <- GRanges(bed.file.hg18$chr, IRanges(bed.file.hg18$start, bed.file.hg18$end), score=bed.file.hg18$score)

bed.file.hg19.GR <- liftOver(bed.file.hg18.GR, chain=chain)

bed.file.hg19.GR.single <- unlist(bed.file.hg19.GR ) 

stopifnot(FALSE)
## development#
# Strongly suspect that chain can be loaded directly from URL.



chain <- import.chain(con=gzcon(file(chain.file.local, "rb")))





chain <- import.chain(gzfile(chain.file.url, "rb"), format="chain")


con <- gzfile(chain.file.local.gz)  # works as a connection.
con <- gzcon(url(chain.file.url))

chain <- import.chain(chain.file.local, format="chain")
chain <- import.chain(txt, format="chain")
?import.chain





con <- gzfile(chain.file.local.gz)  # works as a connection.
txt <- readLines(con)  #


test <- import(text=txt, format="chain")
test <- import(chain.file.local, format="chain")
