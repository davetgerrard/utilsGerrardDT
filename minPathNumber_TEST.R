source("minPathNumber.R")

my.juncs.simple <- data.frame(chr="chr1", start=c(200,300,400,500), end=c(250,350,450,550))
my.juncs <- data.frame(chr="chr1", start=c(200,200,250,350), end=c(500,300,500,500))
my.juncs.complex <- data.frame(chr="chr1", start=c(200,200,250,350,350,550), end=c(500,300,500,500,600,700))




minPathNumber(my.juncs.simple)		# should be 1

minPathNumber(my.juncs)		# should be 3

minPathNumber(my.juncs.complex)		# should be 4 


