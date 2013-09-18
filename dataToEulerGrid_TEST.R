
#####EXAMPLES ########
source("dataToEulerGrid.R")


rawData <- data.frame(	sample.A.1=c(0.2, 0, 0, 0.4, 0.3, 0,0,0.8,0.6,0.1), 
			sample.A.2=c(0.1, 0, 0, 0.5, 0.3, 0.1,0.1,0.7,0.9,0.2),
			sample.B.1=c(0, 0.8, 0, 0.5, 0.4, 0,0,0.8,0.3,0), 
			sample.B.2=c(0, 0.9, 0, 0.6, 0.3, 0,0,0.9,0.8,0),
			sample.C.1=c(0.1, 0.2, 0, 0.2, 0.3, 0,0,0.7,0.2,0), 
			sample.C.2=c(0, 0, 0.1, 0.3, 0.3, 0,0,0.8,0.5,0)	)

rawDataLong <- rbind(rawData,rawData,rawData)




my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B"))

my.groupDesign <- data.frame(sampleName=c("sample.A.1", "sample.A.2", "sample.B.1", "sample.B.2", "sample.C.1", "sample.C.2"),
					groupName=c("group.A", "group.A", "group.B", "group.B", "group.C", "group.C"))

head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign))
head(reduceGroups(binarizeTable(rawDataLong), groupDesign =my.groupDesign,groupRule="all"))

scoreCardinalities(binarizeTable(rawDataLong))
rowCounts <-   scoreCardinalities(binarizeTable(rawDataLong))


plotEuler(rowCounts[,1:6], rowCounts$count, names(rowCounts[,1:6]))

# may want to ignore features present in ALL samples.  N.B. dropEmpty set does same for features absent from all samples but not relevant here.
plotEuler(rowCounts[,1:6], rowCounts$count, names(rowCounts[,1:6]), dropFullSet=T)



