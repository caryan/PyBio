# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 


library('XLConnect')
library('edgeR')

#Control
controlName <- 'MH'

#Drugs
drugNames <- c('Van','Cip', 'DaptoCaCl2', 'DaptoNO', 'Gent', 'Lin', 'Ox', 'Rif', 'ST', 'Cipin', 'Gentin', 'Linin', 'Oxin', 'Rifin', 'STin') 
# drugNames <- c('Cipin', 'Gentin', 'Linin', 'Oxin', 'Rifin', 'STin'  )
# drugNames = c('Van', 'Cip')
DEOutFile = 'DEins.xlsx'

load_data <- function(fileName){
	# #Load the appropriate columns from the Excel sheet
	# tmpData <- readWorksheetFromFile(fileName, sheet="MainSheet", rownames=2)[, c("Start", "Length", "ReadCount")]

	#Strip the troublesome intergenic line from the file
	system(paste("sed -n '/^intergenic/ !p' < ", fileName, " > tmpFile", sep=''))
	tmpData <- read.table("tmpFile", sep="\t", header=TRUE, row.names=2, na.strings="", quote="")[, c("Start", "Length", "ReadCount")]
	#Drop the intergenic
	# tmpData <- tmpData[c(rownames(tmpData) != "intergenic"),]
	#Replace empty values by zeros
	tmpData[is.na(tmpData)] <- 0
	#Sort by start position for convenience
	tmpData <- tmpData[with(tmpData, order(Start)),] 
}

#Load some data
controlData <- list()
drugData <- list()
controlData$A <- load_data(paste('NEWAGG/', controlName, "A.rdp", sep=''))
controlData$B <- load_data(paste('NEWAGG/', controlName, "B.rdp", sep=''))
for (drugName in drugNames){
	print(drugName)
	drugData[[drugName]] <- list()
	drugData[[drugName]]$A <- load_data(paste('NEWAGG/', drugName, "A.rdp", sep=''))
	drugData[[drugName]]$B <- load_data(paste('NEWAGG/', drugName, "B.rdp", sep=''))
}


#Use the read density as a function of position to figure out position normalization
add_norm_counts <- function(geneDataIn, plotName){
	#Remove the bottom/top 10% to reduce influence of outliers
	tmpData <- subset(geneDataIn, ReadCount > quantile(geneDataIn$ReadCount, 0.1) & ReadCount < quantile(geneDataIn$ReadCount, 0.9))
	tmpData$readDensity <- tmpData$ReadCount/tmpData$Length
	# X11()
	png(filename=paste(plotName , '_Bowl.png'))
	plot(tmpData$Start, tmpData$readDensity)
	#Make a loess fit
	loessFit <- loess(readDensity~position, data.frame(position=tmpData$Start, readDensity=tmpData$readDensity), 
				control = loess.control(statistics = c("approximate"),trace.hat = c("approximate"), iterations=3))
	tmpData$smoothed <- predict(loessFit, tmpData$Start)
	lines(tmpData$Start, tmpData$smoothed, col="red", lwd=2)
	dev.off()
	geneDataIn$normRatio <- predict(loessFit, geneDataIn$Start)
	#Flat-line extrapolate past the end points if necessary
	firstNonNAN <- which.min(is.na(geneDataIn$normRatio))
	geneDataIn$normRatio[1:firstNonNAN-1] <- geneDataIn$normRatio[firstNonNAN]
	lastNonNAN <- which.max(is.na(geneDataIn$normRatio))
	geneDataIn$normRatio[-(1:lastNonNAN-1)] <- geneDataIn$normRatio[lastNonNAN-1] 
	#Scale the normRatio so that the median is 1 and adjust the count data
	geneDataIn$normRatio <- geneDataIn$normRatio/median(geneDataIn$normRatio)
	geneDataIn$normCounts <- round(geneDataIn$ReadCount/geneDataIn$normRatio)

	return(geneDataIn)
}

controlData$A <- add_norm_counts(controlData$A, 'MH_A')
controlData$B <- add_norm_counts(controlData$B, 'MH_B')
for (drugName in drugNames){
	drugData[[drugName]]$A <- add_norm_counts(drugData[[drugName]]$A, drugName)
	drugData[[drugName]]$B <- add_norm_counts(drugData[[drugName]]$B, drugName)
}


#Now we can start to import the data into edgeR
edgeRData <- data.frame(MH_A=controlData$A$normCounts, MH_B=controlData$B$normCounts)
for (drugName in drugNames){
	tmpVar <- paste(drugName, '_A', sep='')
	edgeRData[[tmpVar]] <- drugData[[drugName]]$A$normCounts
	tmpVar <- paste(drugName, '_B', sep='')
	edgeRData[[tmpVar]] <- drugData[[drugName]]$B$normCounts
}


rownames(edgeRData) <- rownames(controlData$A)
group <- c(rep("MH",2))
for (drugName in drugNames){
	group <- append(group, rep(drugName, 2))
}
group <- factor(group)
cds <- DGEList(edgeRData, group=group)

#Remove low counts as per edgeR vignette
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds$samples$lib.size <- colSums(cds$counts)

#Now calculate normalization via TMM
cds <- calcNormFactors(cds)

plotMDS( cds )

cds <- estimateCommonDisp( cds, verbose=TRUE )

cds <- estimateTagwiseDisp( cds )

# # Open an Excel worksheet and write the data to it
# unlink(DEOutFile)
# wb <- loadWorkbook(DEOutFile, create=TRUE)

# for (drugName in drugNames){
# 	createSheet(wb, name=drugName)
# 	tmpET <- exactTest(cds, pair=c("MH", drugName))
# 	tmpData <- tmpET$table
# 	tmpData$PValue_adj <- p.adjust(tmpData$PValue, method='BH')
# 	tmpData$locus <- rownames(tmpData)
# 	writeWorksheet(wb, tmpData, sheet=drugName)
# }

# saveWorkbook(wb)
# rm(wb)