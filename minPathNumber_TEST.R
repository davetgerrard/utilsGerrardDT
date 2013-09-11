source("minPathNumber.R")

my.juncs.simple <- data.frame(chr="chr1", start=c(200,300,400,500), end=c(250,350,450,550))
my.juncs <- data.frame(chr="chr1", start=c(200,200,250,350), end=c(500,300,500,500))
my.juncs.complex <- data.frame(chr="chr1", start=c(200,200,250,350,350,550), end=c(500,300,500,500,600,700))
my.juncs.repeat <- data.frame(chr="chr1", start=c(200,300,400,500,500), end=c(250,350,450,550,550), extra=NA)


minPathNumber(my.juncs.simple)		# should be 1

minPathNumber(my.juncs.repeat)		# should be 1

minPathNumber(my.juncs)		# should be 3

minPathNumber(my.juncs.complex)		# should be 4 



### N.B. Works for small numbers of overlapping junctions (<10) but quickly becomes very slow to compute for anything more.

# Another test set with more junctions and overlaps.
testJuncs <- data.frame(chr="chr1", start=seq(1, 801, by=50), end= seq(121, 921, by=50))


# Don't do this it takes a VERY LONG TIME
# system.time(minPathNumber(testJuncs))


#instead, plot the times for 1:12 juncs
times <- numeric()

for(i in 1:12)  {
	times[i] <- system.time(minPathNumber(testJuncs[1:i,]))['elapsed']

}

plot(times, xlab="number of junctions", ylab="time")


# any kind of approximation for splcing complexity?  
# Number of overlaps not inlcuding self:- 
sum(getOverlapMatrix(testJuncs$start, testJuncs$end, half_matrix=TRUE), na.rm=T) - nrow(testJuncs)

# but can easily show that this number is not a good predictor of solvability with minPathNumber()
# e.g. compare the above testJuncs (which are very slow to compute) with a set where no junctions overlap EXCEPT two that overlap all others.
testJuncs.easy <- data.frame(chr="chr1", start=c(2,3, seq(1, 751, by=50)), end= c(820,821,seq(21, 771, by=50)))
sum(getOverlapMatrix(testJuncs.easy $start, testJuncs.easy $end, half_matrix=TRUE), na.rm=T) - nrow(testJuncs.easy )

#system.time(minPathNumber(testJuncs.easy))
# Errrr. not easy apparently!  This one also seems to take a very long time (i stopped it after 15 mins). 

##   GRRRR NP-COMPLETE!  http://en.wikipedia.org/wiki/Set_cover_problem
