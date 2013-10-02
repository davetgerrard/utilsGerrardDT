# geneId	a name for the gene (NOT USED!)
# bindingData	the full set of binding information for one or more binding factors. Data.frame with columns 'chr' 'start' 'end' 'name' 'score'. 'name' refers to the name of the factor and must be identical for all occurences of the same factor.
# focalChrom	the chromosome of the focal point. Must match format of bindingData (e.g. 'chr1')
# focalPoint	a single nucleotide position marking the centre of the region of interest (e.g. a transcription start site (tss))
# viewWindow	Collect data from factors bound up to this far from the focalPoint [Default = 10000]
# effectDistance  Binding scores are downweighted  using a sine-decay function between the focalPoint and this distance. Beyond this distance, scores are downweighted by a factor of 1/effectDistance. [Default == viewWindow]
# 
getFactorProfile <- function(geneId, bindingData, focalChrom, focalPoint, viewWindow=10000, effectDistance=viewWindow)  {
	factorProfile <- list()		#pre-declare so that empty list can be returned
	windowStart <- max(0,(focalPoint-viewWindow))
	windowEnd <- focalPoint+viewWindow

	validBindingIndex <- (allData$chr == focalChrom) & (allData$start > windowStart) & (allData$end < windowEnd )
	validBound <- allData[validBindingIndex ,]
	if(nrow(validBound) > 0)  {	

		validBound$midpoint <- ((validBound$end- validBound$start) /2 ) + validBound$start	
		validBound$distance <- abs(validBound$midpoint - focalPoint)

		## could be adapted to have plateau around the tss. 
		validBound$distWeight <- ifelse(validBound$distance >= effectDistance, 1/effectDistance,((sin(((pi) + (2*(validBound$distance/effectDistance) * pi))/2) ) + 1 ) / 2)
		validBound$weightedScore <- validBound$score * validBound$distWeight
		# collect aggregated results for each factor.
		factorProfile <- aggregate(validBound$weightedScore, list(factor=validBound$name), sum)
	}
	return(factorProfile)
}

#bindingData must have columns: chr, start, end, name, score.
getFactorPositionMatrixList <- function(bindingData, boundFactors, focalChroms, focalPoints, viewWindow=10000)  {
	stopifnot(length(focalChroms) == length(focalPoints))
	factorPositionMatrixList <- list()		#pre-declare so that empty list can be returned
	for(thisFactor in boundFactors)  {
		factorPositionMatrixList[[thisFactor]] <- data.frame()

	}
	
	for(i in 1:length(focalPoints))  {	
		cat(focalChroms[i],focalPoints[i],": ")
		windowStart <- max(0,(focalPoints[i]-viewWindow))
		windowEnd <- focalPoints[i]+viewWindow
	
		validBindingIndex <- (allData$chr == focalChroms[i]) & (allData$start > windowStart) & (allData$end < windowEnd )
		validBound <- allData[validBindingIndex ,]
		if(nrow(validBound) > 0)  {
			validBound$relPosition.start <- validBound$start - focalPoints[i]  
			validBound$relPosition.end  <- validBound$end - focalPoints[i] 
			validBound$relPosition.midpoint <- validBound$relPosition.start + floor(validBound$relPosition.end - validBound$relPosition.start)
			for(thisFactor in boundFactors) {
				cat(thisFactor, " ")
				thisFocalFactor.index <- validBound$name == thisFactor
				#cat(sum(thisFocalFactor.index), "\n")
				#cat(thisFocalFactor.index, "\n")
				thisRow <- data.frame(rel.start=validBound$relPosition.start[thisFocalFactor.index], rel.end=validBound$relPosition.end[thisFocalFactor.index], rel.midpoint=validBound$relPosition.midpoint[thisFocalFactor.index], score=validBound$score[thisFocalFactor.index])
				factorPositionMatrixList[[thisFactor]] <- rbind(factorPositionMatrixList[[thisFactor]], thisRow)
				cat(nrow(factorPositionMatrixList[[thisFactor]]), " ")
			}
		}
		cat("\n")
	}

	return(factorPositionMatrixList)
}

#test <- getFactorPositionMatrixList(allData, "CTCF", "chr13", 28494157)

# Stranded version
#bindingData must have columns: chr, start, end, name, score.
getFactorPositionMatrixListStranded <- function(bindingData, boundFactors, focalChroms, focalPoints, focalStrands='+', viewWindow=10000)  {
	stopifnot(length(focalChroms) == length(focalPoints))
	factorPositionMatrixList <- list()		#pre-declare so that empty list can be returned
	for(thisFactor in boundFactors)  {
		factorPositionMatrixList[[thisFactor]] <- data.frame()

	}
	
	for(i in 1:length(focalPoints))  {	
		cat(focalChroms[i],focalPoints[i],": ")
		windowStart <- max(0,(focalPoints[i]-viewWindow))
		windowEnd <- focalPoints[i]+viewWindow
	
		validBindingIndex <- (allData$chr == focalChroms[i]) & (allData$start > windowStart) & (allData$end < windowEnd )
		validBound <- allData[validBindingIndex ,]
		if(nrow(validBound) > 0)  {
			if(focalStrands[i] == '-')  {
				validBound$relPosition.start <- focalPoints[i] - validBound$start
				validBound$relPosition.end  <- focalPoints[i] - validBound$end
			} else  {
				validBound$relPosition.start <- validBound$start - focalPoints[i]  
				validBound$relPosition.end  <- validBound$end - focalPoints[i]  
			}
			validBound$relPosition.midpoint <- validBound$relPosition.start + floor(validBound$relPosition.end - validBound$relPosition.start)
			for(thisFactor in boundFactors) {
				cat(thisFactor, " ")
				thisFocalFactor.index <- validBound$name == thisFactor
				#cat(sum(thisFocalFactor.index), "\n")
				if(sum(thisFocalFactor.index ) > 0)  {
					#thisRow <- data.frame(rel.start=validBound$relPosition.start[thisFocalFactor.index], rel.end=validBound$relPosition.end[thisFocalFactor.index], rel.midpoint=validBound$relPosition.midpoint[thisFocalFactor.index], score=validBound$score[thisFocalFactor.index])
					#cat(thisFocalFactor.index)
					thisRow <- data.frame(focal.chr = as.character(focalChroms[i]), focal.point=as.integer(focalPoints[i]), rel.start=validBound$relPosition.start[thisFocalFactor.index], rel.end=validBound$relPosition.end[thisFocalFactor.index], rel.midpoint=validBound$relPosition.midpoint[thisFocalFactor.index], score=validBound$score[thisFocalFactor.index])
					factorPositionMatrixList[[thisFactor]] <- rbind(factorPositionMatrixList[[thisFactor]], thisRow)
				}
				cat(nrow(factorPositionMatrixList[[thisFactor]]), " ")
			}
		}
		cat("\n")
	}

	return(factorPositionMatrixList)
}


scoreTestOnFPML <-  function(fg.FPML, bg.FPML, test="wilcox", factors=names(fg.FPML))  {
	library(qvalue)
	testResults <- numeric()
	fg.means <- numeric()
	bg.means <- numeric()

	for(thisFactor in factors)  {
		#print(paste(thisFactor, ">"))
		fg.scores <- fg.FPML[[thisFactor]]	
		#print(fg.scores)
		bg.scores <- bg.FPML[[thisFactor]]
		#print(bg.scores)
		if( nrow(fg.scores) < 1  ||  nrow(bg.scores) < 1 )  {
			print(thisFactor)
			testResults[thisFactor] <- NA	
			fg.means[thisFactor] <- NA	
			bg.means[thisFactor] <- NA	
		} else  {
			testResults[thisFactor] <- wilcox.test(fg.scores[,'score'],bg.scores[,'score'], alternative="t")$p.value
			fg.means[thisFactor] <- mean(fg.scores[,'score'])
			bg.means[thisFactor] <- mean(bg.scores[,'score'])
		}
	}

	resultSummary <- data.frame(encodeFactor=names(testResults),wilcox.p.value=testResults, bg.mean=bg.means, fg.mean=fg.means)
	
	# cannot allow NA values to be submitted to qvalue(). Need to split the table and recombine after calculating qvalues.
	resultSummary$qvalue <- NA
	resultSummary.na <- resultSummary[is.na(resultSummary$wilcox.p.value), ]
	print(nrow(resultSummary.na))
	resultSummary.values <- resultSummary[!is.na(resultSummary$wilcox.p.value), ]	
	print(nrow(resultSummary.values))
	resultSummary.values$qvalue <- qvalue(resultSummary.values$wilcox.p.value)$qvalues
	resultSummary  <- rbind(resultSummary.values, resultSummary.na)
	resultSummary  <- resultSummary[order(resultSummary$qvalue),]
	
}

