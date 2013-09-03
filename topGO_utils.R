
#########################################
#					
#	Dave Gerrard 			
#	University of Manchester	
#			
#					
#########################################

##INFO: Utility functions to run multiple GO analyses over goData objects and store results in tables.

########## FUNCTIONS


listProtsInGoFromList <- function(goTerm,ontology,goDataList)  {	#INFO: list proteins for a GO term in a precalculated list of goData objects, accounting for ontology
	genesInTerm(goDataList[[ontology]],goTerm)[[1]]
}


topDiffGenes <- function(allScore) {	#INFO: the topGOdata objects require a method for gene/protein selection. For enrichment analysis, we keep all proteins.
	return(allScore )
}


resultSummaryAsText <-function (x)    #INFO: this was a hack to get results out of the GOdata object and used for text output. May be out-of-date.
{
    text <-     paste("\nDescription:", description(x), "\n",
    	"'", algorithm(x), "' algorithm with the '", testName(x), 
        "' test\n", 
    length(score(x)), "GO terms scored:", sum(score(x) <= 
        elimCutOff), "terms with p < ",elimCutOff ,"\n")
    xg <- geneData(x)
    if ("Annotated" %in% names(xg)) {
    	text <- paste(text,"    Annotated proteins:", xg["Annotated"], "\n")
    	}
    if ("Significant" %in% names(xg)) {
            text <- paste(text,"    Significant proteins:", xg["Significant"], "\n")
        }
    text <- paste(text,"    Min. no. of proteins annotated to a GO:", xg["NodeSize"],"\n")
    if ("SigTerms" %in% names(xg)) {
        text <- paste(text,"    Nontrivial nodes:", xg["SigTerms"], "\n")
        }
    return(text) ;
    #.printGeneData(geneData(x))
}


runDetectGoTests <- function(GOdataBase )  {	#INFO: runs a suite of GO detected/undetected tests over a topGOdata object and returns a result table
	require(gplots)
	thisGOgraph <- GOdataBase@ontology
	#GOdataBase <- goDataCollection.ubiq.binary[[thisGOgraph]]
	# not sure if updateGenes() required here. 
	test.stat <- new("classicCount", testStatistic = GOFisherTest, name="Fisher test")
	resultFisher <- getSigGroups(GOdataBase,test.stat)
	test.stat <- new("elimCount", testStatistic = GOFisherTest, name="elimFisher", cutOff = elimCutOff)
	resultElimFisher <- getSigGroups(GOdataBase,test.stat)
	resultElimFisherAdj <- resultElimFisher
	score(resultElimFisherAdj) <- qvalue(score(resultElimFisher))$qvalue
	allResFT <- GenTable(GOdataBase, classic=resultFisher,elim=resultElimFisher,qvalElim=resultElimFisherAdj,
					orderBy="qvalElim",topNodes=topTerms )
	#allResFT <- subset(allResFT, select=-"Rank in elim")
	#allResFT
	headText <- resultSummaryAsText(resultElimFisher)
	headText <- paste(headText, "\nTop ",topTerms,"GO terms shown")
	par(mar=c(1,1,1,1))
	layout(matrix(c(1,2), byrow=T),heights=c(1,3)) 
	textplot(headText,halign="left",mar=c(0,0,0,0))
	textplot(allResFT,mar=c(0,0,0,0))

	###  getMethods("GenTable")

	### create summmary table manually and add to single table for all results across ontologies

	manual.Fisher <- data.frame(Fisher=score(resultFisher),goTerm=names(score(resultFisher)))
	manual.elimFisher <- data.frame(elimFisher=score(resultElimFisher),goTerm=names(score(resultElimFisher)))
	manual.termStats <- termStat(GOdataBase,names(score(resultFisher)))
	manual.termStats$goTerm <- row.names(manual.termStats)

	thisGoResults <- NULL
	thisGoResults <- merge(manual.termStats,manual.Fisher,by="goTerm")
	thisGoResults <- merge(thisGoResults ,manual.elimFisher ,by="goTerm")

	thisGoResults$ontology <- thisGOgraph
	thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(as.character(thisGoResults$goTerm), ontology(GOdataBase),numChar=200))
	#goGroup <- as.character(thisGoResults$goTerm)
	#thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdataBase,x)[[1]])))

	#summaryDetectResults.ubiq <- rbind(summaryDetectResults.ubiq,thisGoResults)
	#rm(manual.Fisher, manual.elimFisher, manual.termStats,thisGoResults)
	#rm(resultElimFisherAdj, resultElimFisher,  resultFisher,allResFT,headText)
	return(thisGoResults)

}


## functions required by runScoreGoTests()   to find min rank of values from either extreme. 
#  	(first and last both rank 1, median value is ranked n/2)
minRankFromMedian <- function(vector)  {
	rankVector <- rank(vector)
	inv.rankVector <- (length(vector) +1 ) - rankVector
	minRankVector <- ifelse(rankVector < inv.rankVector, rankVector, inv.rankVector)	
	return(minRankVector)
}

runScoreGoTests <- function(GOdata,geneList,geneSelectionFun=topDiffGenes)  { #INFO: runs a suite of GO score tests over a topGOdata object and returns a result table
		thisGOgraph <- GOdata@ontology
	
		####INFO: define test statistics and apply over all GO terms 
		#INFO: resultWilcox = 2-sided wilcox test. Good for outliers in one direction
		test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
		resultWilcox <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcox = elim version of resultWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcox <- getSigGroups(GOdata, test.stat)

		test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
		resultWilcoxGreater <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcox = elim version of resultWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxGreater <- getSigGroups(GOdata, test.stat)

		test.stat <- new("classicScore", testStatistic = GOWilcoxTestLesser, name = "Wilcox tests") 
		resultWilcoxLesser <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcox = elim version of resultWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestLesser, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxLesser <- getSigGroups(GOdata, test.stat)

		#resultElimWilcoxAdj <- resultElimWilcox 
		#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue
		#INFO: resultKS = Kolmogorov smirnov test: may be good for general distribution changes.
		test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
		resultKS <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimKS = elim version of resultKS
		test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
		resultElimKS <- getSigGroups(GOdata, test.stat)

		#resultElimKSAdj <- resultElimKS
		#score(resultElimKSAdj ) <- qvalue(score(resultElimKS))$qvalue

		#INFO: change scores to absolute values for Greater than test 
		geneList2 <- abs(geneList)
		names(geneList2) <- names(geneList)
		GOdata <- updateGenes(GOdata,geneList2,topDiffGenes)
		#INFO: resultWilcoxAbs = test for outliers in both directions simultaneously. Less power than wilcox if true difference is unidireectional.
		test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
		resultWilcoxAbs <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcoxAbs = elim version of resultElimWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)


		#INFO: change scores to ranks from edge for Greater than test  (first and last both rank 1, median value is ranked n/2)
		geneList3 <- minRankFromMedian(geneList)
		names(geneList3) <- names(geneList)
		GOdata <- updateGenes(GOdata,geneList3,topDiffGenes)
		#INFO: resultWilcoxMedRank = test for outliers in both directions simultaneously. Less power than wilcox if true difference is unidireectional. Different to resultWilcoxAbs if distribution is skewed.
		test.stat <- new("classicScore", testStatistic = GOWilcoxTestLesser, name = "Wilcox tests") 
		resultWilcoxMedRank <- getSigGroups(GOdata, test.stat)
		#INFO: resultElimWilcoxAbs = elim version of resultElimWilcox
		test.stat <- new("elimScore", testStatistic = GOWilcoxTestLesser, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
		resultElimWilcoxMedRank <- getSigGroups(GOdata, test.stat)		



		#resultElimWilcoxAbsAdj <- resultElimWilcoxAbs 
		#score(resultElimWilcoxAbsAdj ) <- qvalue(score(resultElimWilcoxAbs ))$qvalue


		#INFO: collect and bind all the results together. 
		# This could be done with a function available in topGO but I wanted extra columns and control.
		Wilcox <- data.frame(Wilcox=score(resultWilcox ),goTerm=names(score(resultWilcox)))
		elimWilcox <- data.frame(elimWilcox=score(resultElimWilcox),goTerm=names(score(resultElimWilcox)))
		WilcoxGreater <- data.frame(WilcoxGreater=score(resultWilcoxGreater ),goTerm=names(score(resultWilcoxGreater)))
		elimWilcoxGreater <- data.frame(elimWilcoxGreater=score(resultElimWilcoxGreater),goTerm=names(score(resultElimWilcoxGreater)))
		WilcoxLesser <- data.frame(WilcoxLesser=score(resultWilcoxLesser ),goTerm=names(score(resultWilcoxLesser)))
		elimWilcoxLesser <- data.frame(elimWilcoxLesser=score(resultElimWilcoxLesser),goTerm=names(score(resultElimWilcoxLesser)))		
		KS <- data.frame(KS=score(resultKS),goTerm=names(score(resultKS)))
		elimKS <- data.frame(elimKS=score(resultElimKS),goTerm=names(score(resultElimKS)))
		absWilcox <- data.frame(absWilcox=score(resultWilcoxAbs),goTerm=names(score(resultWilcoxAbs)))
		elimAbsWilcox <- data.frame(elimAbsWilcox=score(resultElimWilcoxAbs),goTerm=names(score(resultElimWilcoxAbs)))
		medRankWilcox <- data.frame(medRankWilcox=score(resultWilcoxMedRank),goTerm=names(score(resultWilcoxMedRank)))
		elimMedRankWilcox <- data.frame(elimMedRankWilcox=score(resultElimWilcoxMedRank),goTerm=names(score(resultElimWilcoxMedRank)))


		thisGoResults <- NULL
		thisGoResults <- merge(Wilcox,elimWilcox,by="goTerm")
		thisGoResults <- merge(thisGoResults, WilcoxGreater,by="goTerm")
		thisGoResults <- merge(thisGoResults, elimWilcoxGreater,by="goTerm")
		thisGoResults <- merge(thisGoResults, WilcoxLesser,by="goTerm")
		thisGoResults <- merge(thisGoResults, elimWilcoxLesser,by="goTerm")
		thisGoResults <- merge(thisGoResults, KS ,by="goTerm")
		thisGoResults <- merge(thisGoResults, elimKS ,by="goTerm")
		thisGoResults <- merge(thisGoResults, absWilcox ,by="goTerm")
		thisGoResults <- merge(thisGoResults, elimAbsWilcox ,by="goTerm")
		thisGoResults <- merge(thisGoResults, medRankWilcox, by="goTerm")
		thisGoResults <- merge(thisGoResults, elimMedRankWilcox, by="goTerm")
		thisGoResults$ontology <- thisGOgraph
		thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(as.character(thisGoResults$goTerm), ontology(GOdata),numChar=200))

		goGroup <- as.character(thisGoResults$goTerm)
		thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdata,x)[[1]])))
		
		#INFO: bind the results from this GO ontology as rows to the table for the current PC in 'summaryResultsList'
		#summaryResultsList[[i]] <- rbind(summaryResultsList[[i]],thisGoResults)
		return(thisGoResults)
}


GOWilcoxTestGreater <- function (object,alternativeType="greater")	#INFO: Wilcox 1-sided (greater) test for use in topGO 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}


GOWilcoxTestLesser <- function (object,alternativeType="less")	#INFO: Wilcox 1-sided (less) test for use in topGO 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}

GOWilcoxTest2Sided <- function (object,alternativeType="two.sided")	#INFO: Wilcox 2-sided test for use in topGO  
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}


addQvalueToTable <- function(df, p.columns="p.value") {
	if(require(qvalue))  {
		for(thisCol in p.columns)  {
			colName <- paste("q",thisCol,sep=".")
			df[,colName] <- qvalue(df[,thisCol])$qvalues
		}
		return(df)
	}   else  {
		warning("Package qvalue not found!")
		return(df)
	}
}

#test <- data.frame(p.value = runif(20), other=runif(20))
#addQvalueToTable(test)
#addQvalueToTable(test, p.columns=c("p.value", "other"))


########### utils to cluster tables of GO results


####INFO: utility functions to combine a GO results table with functional clusters.


############GENERIC

listBestOverlappingCluster <- function(goProteins.valid, seedList, goContainProp)  { 	#INFO: give the index of the functional cluster which best represents this go term IF the cluster contains at least goContainProp of the proteins assigned to the go term.
	# which clusters contain greater than 'goContainProp' of the proteins annotated to this goTerm.
	maxValue <-  max(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))))
	clusHitList <- which(unlist(lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid))/length(goProteins.valid))) == maxValue )
	if((length(clusHitList) > 0) & (maxValue > goContainProp)) {
		clusHitList[1] 
	} else {NA}
}

listAllOverlappingClusters <- function(goProteins.valid, seedList, goContainProp)  { 	# give all clusters containing at least goContainProp of the proteins assigned to the go term as a singe text string.
	#length(goProteins.valid)
	#lapply(seedList,FUN = function(x) length(intersect(x,goProteins.valid)))
	# which clusters contain greater than 'goContainProp' of the proteins annotated to this goTerm.
	clusHitList <- which(lapply(seedList,FUN = function(x) (length(intersect(x,goProteins.valid))/length(goProteins.valid))) > goContainProp)
	if(length(clusHitList) > 0) {
		clusHitList <- paste(clusHitList,collapse=",")
	} else {NA}
}


#INFO: clusterTabel() filter a table based on a treshold and reorders by keeping members of the same cluster together. 
#INFO: the table is filtered, marked and re-ordered, grouping GO terms which represent the same functional cluster of proteins/genes
clusterTable <- function(goTable,orderBy,detectTableSigThreshold = 0.05, max.genes=500, ontology="all")  {
	# find a way to output the table (or a subsection) with (significant) cluster members grouped. 
	# Only works with single (best) cluster per GO term. Limit to terms below significance threshold.
	sigTable <- subset(goTable,goTable[,orderBy] < detectTableSigThreshold  & goTable[,"number"] <= max.genes)
	if(ontology != "all") {
		#sigTable <- subset(sigTable , ontology == ontology)
		#sigTable <- subset(sigTable , sigTable[,"ontology"] == ontology)
		sigTable <- sigTable[sigTable[,"ontology"] == ontology,]

	}
	# order the table to cluster lower terms with top term of each cluster.
	sigTable <- sigTable[order(sigTable[,orderBy]),]
	sigTable$outputRank <- sigTable$pureRank<- rank(sigTable[,orderBy])
	#levels(as.factor(as.character(sigTable$bestCluster)))
	# assign same rank to all significant members of same cluster.
	##sigTable$bestCluster <- as.factor(as.character(sigTable$bestCluster))
	##for(thisCluster in levels(sigTable$bestCluster))  {
	##	if(is.na(thisCluster)) {break}
	##	rankToSet <- min(na.omit(sigTable$outputRank[sigTable$bestCluster == thisCluster]))
	##	sigTable$outputRank[sigTable$bestCluster == thisCluster] <- rankToSet
	##}	
	sigTable$subClusterTerm <- ifelse((sigTable$outputRank == sigTable$pureRank),FALSE,TRUE)
	sigTable <- sigTable[order(sigTable[,orderBy]),]
}





############# additional utils to create vioplots.
## regquires gplots, lattice


addGroupScoresToGroupTable <- function(groupTable,groupIds,groupName,scoreTable,idColumn,scoreColum) {	#INFO: Gets data for a list of proteins and binds as rows onto a table. Resulting stack can be used as basis for large vioplot or boxplot of several groups.
	index <- na.omit(match(groupIds , scoreTable[,idColumn]))
	scores <- scoreTable[index ,scoreColum]
	groupTable <- rbind(groupTable,data.frame(score=scores,group=groupName))
	
}


bwPlotsAsPdf <- function(pc.scores.table,resultTableList,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, goDataList=goDataCollection.ubiq.scored, ontology="all", max.genes=500)  {	#INFO: Control plotPCsFromSummary() for output in a pdf
	pdfName <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,"terms",max.genes,"genes", ontology,"ontologies","pdf",sep=".")

	pdf(pdfName, paper="a4")
	frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results\n and ",max.genes, "genes per term.\nOntologies:",ontology,sep=" ")
	textplot(frontPageText)
	for(pc.i in 1:length(resultTableList)) {
		plotPCsFromSummary(pc.scores.table=pc.scores.table,resultTableList[[pc.i]],orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
	}
	dev.off()
}

bwPlotsAsFigs <- function(pc.scores.table,resultTableList,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, figType="tiff", goDataList=goDataCollection.ubiq.scored, ontology="all", max.genes=500)  { #INFO: Control plotPCsFromSummary() for output in a figure
	for(pc.i in 1:length(resultTableList)) {
		plotDir <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,"terms",max.genes,"genes", ontology,"ontologies",figType,sep="_")
		if(!file.exists(plotDir))  {dir.create(plotDir)}
		fileName <- paste("vioPlotsByPC",pc.i,orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep=".")
		fileName <- paste(plotDir,fileName,sep="/")
		tiff(fileName, compression="lzw",width=180, height=480,units="mm",res=300)	
		plotPCsFromSummary(pc.scores.table=pc.scores.table,resultTableList[[pc.i]],orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i, goDataList=goDataList)
		dev.off()	
	}
	#par(op)
}




## Currently very specific to this project and data objects. Could generalise by specifying scoreTable and score column as variables.
plotPCsFromSummary <- function(pc.scores.table, resultsTable,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, pc.i=1, goDataList=goDataCollection.ubiq.scored, ontology="all", max.genes=500)  {	#INFO: gets data and plots a nice set of violin plots
	#INFO: need to make a cluster table out of the resultsTable so that redundant terms are shown as a cluster.
	sigTable <- clusterTable(resultsTable, orderBy=orderBy.PC,detectTableSigThreshold = pcTableSigThreshold, ontology=ontology, max.genes= max.genes) 

	## need to count each cluster only once and each independent go term once. 
	# take the first instance of each outputRank
	plotTable <- sigTable[match(unique(sigTable$outputRank),sigTable$outputRank),]
	numbGraphResultsLocal <- min(nrow(plotTable),numbGraphResults)
	if(numbGraphResultsLocal < 1)  {	# are there going to be enough results to plot?
		warningMess <- paste("No sig results for PC", pc.i,"\n using", orderBy.PC, "below", pcTableSigThreshold)
		textplot(warningMess)
		return(NULL)
	}
	plotTable <- plotTable[1:numbGraphResultsLocal,]	# only going to try to plot the top terms

	# could try strwrp or substr
	compHead <- paste("Comp.",pc.i,sep="")
	#INFO: Want the scores for ALL proteins as first on plot for visual comparison.
	baseScore <- data.frame(score=pc.scores.table[,compHead],group="                            All Proteins")		# long name to handle padding of output
	groupTable <- baseScore
	plotList <- list()
	for(i in 1:nrow(plotTable))  {
		##if(is.na(plotTable$bestCluster[i]))  {
			#plotList[[i]] <- intersect(listProtsInGoFromList(goTerm=as.character(plotTable[i,"goTerm"]),ontology=as.character(plotTable[i,"ontology"]),goDataList),validProts)
			plotList[[i]] <- listProtsInGoFromList(goTerm=as.character(plotTable[i,"goTerm"]),ontology=as.character(plotTable[i,"ontology"]),goDataList)
			# some goDescriptions are very long. Currently wrapping the terms rather than truncating.
			goID <- paste(plotTable[i,"description"],plotTable[i,"goTerm"],sep=",")
			#goID <- paste(paste(strwrap(plotTable[i,"description"],width=40),collapse="\n"),plotTable[i,"goTerm"],sep="\n")
		##}	else {
			# this should be a cluster
		##	plotList[[i]] <- seedList[[plotTable$bestCluster[i]]]
		##	goID <- paste("CLUSTER", plotTable$bestCluster[i])
		#}
		# BEWARE: lots of dependencies here
		groupTable <- addGroupScoresToGroupTable(groupTable=groupTable,groupIds=plotList[[i]],
					groupName=goID,scoreTable=pc.scores.table,
					idColumn="gene_id",scoreColum=compHead)
	}


	#bwTitle <- compHead
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	bwTitle <- paste("Principal Component",pc.i,sep=" ")
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	plot(bwplot(bymedian ~ score, groupTable, main=bwTitle,  aspect="xy", cex.main=2, cex.lab=1.5,
		par.settings = list(layout.widths = list(axis.left = 0, ylab.axis.padding = 40)),
	       panel = function(..., box.ratio) {
		   panel.violin(..., col = "darkgrey",
				varwidth = FALSE, box.ratio = box.ratio)
		   panel.bwplot(..., fill = "lightgrey", box.ratio = .1)
	       } )
	)
}









