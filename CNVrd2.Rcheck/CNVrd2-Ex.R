pkgname <- "CNVrd2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('CNVrd2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CNVrd2-class")
### * CNVrd2-class

flush(stderr()); flush(stdout())

### Name: CNVrd2-class
### Title: Class '"CNVrd2"'
### Aliases: CNVrd2-class countReadInWindow,CNVrd2-method
###   plotCNVrd2,CNVrd2-method
### Keywords: classes

### ** Examples

showClass("CNVrd2")



cleanEx()
nameEx("calculateLDSNPandCNV")
### * calculateLDSNPandCNV

flush(stderr()); flush(stdout())

### Name: calculateLDSNPandCNV
### Title: calculateLDSNPandCNV
### Aliases: calculateLDSNPandCNV

### ** Examples

##Load data: fcgr3bMXL in CNVrd2 package############
data(fcgr3bMXL)
##Name a vcf file (vcfFile)
vcfFile <- system.file(package="CNVrd2", "extdata",
                      "chr1.161600000.161611000.vcf.gz")
##Make a data fame named sampleCNV including samples, CNs, population names

sampleCNV <- data.frame(copynumberGroups$allGroups[, c(1,2) ],rep("MXL", 58))

rownames(sampleCNV) <- substr(sampleCNV[, 1], 1, 7)
sampleCNV[, 1] <- rownames(sampleCNV)
##The first column must be the sample names 
tagSNPandINDELofMXL <- calculateLDSNPandCNV(sampleCNV = sampleCNV,
                                 vcfFile = vcfFile, cnvColumn = 2,
                                 population = "MXL", popColumn = 3,
                                 nChunkForVcf = 5, chr = "1",
                                 st = 161600000, en = 161611000,
                                 codeSNP= "Three", codeCNV = "ThreeGroup")
tagSNPandINDELofMXL[1:3,]







cleanEx()
nameEx("clusteringCNVs-class")
### * clusteringCNVs-class

flush(stderr()); flush(stdout())

### Name: clusteringCNVs-class
### Title: Class '"clusteringCNVs"'
### Aliases: clusteringCNVs-class objectCluster
###   emnormalCNV,clusteringCNVs-method groupCNVs,clusteringCNVs-method
###   searchGroupCNVs,clusteringCNVs-method
### Keywords: classes

### ** Examples

showClass("clusteringCNVs")



cleanEx()
nameEx("countReadInWindow-methods")
### * countReadInWindow-methods

flush(stderr()); flush(stdout())

### Name: countReadInWindow-methods
### Title: Method 'countReadInWindow'
### Aliases: countReadInWindow-methods
### Keywords: methods

### ** Examples

##data(fcgr3bMXL)
##readCountMatrix <- countReadInWindow(Object = objectCNVrd2, correctGC = TRUE)
##readCountMatrix[1:3, 1:3]




cleanEx()
nameEx("countReadInWindow")
### * countReadInWindow

flush(stderr()); flush(stdout())

### Name: countReadInWindow
### Title: Obtain read counts in constant windows.
### Aliases: countReadInWindow
### Keywords: methods

### ** Examples

## Not run: 
##D data(fcgr3bMXL)
##D bamFiles <- dir("Bam", pattern = ".bam$")
##D objectCNVrd2 <- new("CNVrd2", windows = 1000, chr = "chr1",
##D                    st = 161100001, en = 162100000,
##D                    dirBamFile = "Bam",
##D                    genes = c(161592986, 161601753),
##D                    geneNames = "3B")
##D 
##D readCountMatrix <- countReadInWindow(Object = objectCNVrd2, correctGC = TRUE)
##D readCountMatrix[1:3, 1:3]
## End(Not run)



cleanEx()
nameEx("emnormalCNV")
### * emnormalCNV

flush(stderr()); flush(stdout())

### Name: emnormalCNV
### Title: Implement the EM algorithm
### Aliases: emnormalCNV

### ** Examples

data(fcgr3bMXL)

sS <- resultSegment$segmentationScores
#########Histogram###########################
###View segmentation scores##################
hist(sS[, 1], 100)
############################################
##Number of components#######################
###Make an object of clusteringCNVs class######
objectCluster <- new("clusteringCNVs",
                     x = sS[, 1], k = 4, EV = TRUE)

set.seed(123)
copynumberGroups <- groupCNVs(Object = objectCluster)



cleanEx()
nameEx("groupBayesianCNVs")
### * groupBayesianCNVs

flush(stderr()); flush(stdout())

### Name: groupBayesianCNVs
### Title: groupBayesianCNVs
### Aliases: groupBayesianCNVs
### Keywords: methods

### ** Examples

## Not run: 
##D data(ccl3l1data)
##D 
##D xyEuro <- ccl3l1data[grep("CEU|TSI|IBS|GBR|FIN", ccl3l1data[, 2]), ]
##D 
##D names(yEuro) <- rownames(xyEuro)
##D 
##D ##Clustering European segmentation scores into group: 5 groups were chosen
##D 
##D objectClusterEuroCCL3L1 <- new("clusteringCNVs", x = yEuro, k = 5)
##D 
##D europeanCCL3L1Groups <- groupCNVs(Object = objectClusterEuroCCL3L1)
##D 
##D ##Obtain prior information
##D #Means
##D lambda0 <- as.numeric(europeanCCL3L1Groups$m)
##D #SD
##D sdEM <- as.numeric(europeanCCL3L1Groups$sigma)
##D #Proportions
##D pEM <- as.numeric(europeanCCL3L1Groups$p)
##D 
##D 
##D ###Calculate the distances between groups
##D for (ii in 2:5){print(lambda0[ii] - lambda0[ii-1])}
##D 
##D ###All segmentation scores
##D ccl3l1X <- ccl3l1data$SS
##D names(ccl3l1X) <- as.character(ccl3l1data$Name)
##D range(ccl3l1X)
##D 
##D 
##D  
##D ##Set prior information: 
##D #prior for the sd of the means of groups: 
##D #5 was set for the third group = 2 CN
##D sd <- c(1, 1, 5, 1, 1) 
##D ccl3l1X <- sort(ccl3l1X)
##D ###Data
##D xData <- ccl3l1X
##D ###Number of groups
##D nGroups <- 10 
##D ###prior for means of groups
##D lambda0 <- lambda0 
##D ###Prior for mixing proportions
##D alpha0 <-  c(3, 29, 44, 18, 7,  5, rep(2, nGroups -length(pEM) -1))
##D ##Prior for the distances between groups
##D distanceBetweenGroups = 0.485
##D 
##D sdEM = sdEM
##D 
##D 
##D ##Adjust standard deviation for the fifth group
##D sdEM[5] <- sdEM[4]
##D  
##D  set.seed(123)
##D  groupCCL3L1allPops <- groupBayesianCNVs(xData = xData, nGroups = nGroups,
##D                                          lambda0 = lambda0,
##D                                          sd0 = sdEM, alpha0 = alpha0,
##D                                          distanceBetweenGroups = distanceBetweenGroups,
##D                                          sdOftau = sd,
##D                                         rightLimit = 4)
##D 
##D 
##D 
## End(Not run)



cleanEx()
nameEx("groupCNVs")
### * groupCNVs

flush(stderr()); flush(stdout())

### Name: groupCNVs
### Title: Cluster segmentation scores into groups.
### Aliases: groupCNVs
### Keywords: methods

### ** Examples

data("fcgr3bMXL")
#resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)
objectCluster <- new("clusteringCNVs",
                     x = resultSegment$segmentationScores[, 1], k = 4, EV = TRUE)

#searchGroupCNVs(Object = objectCluster)
copynumberGroups <- groupCNVs(Object = objectCluster)




cleanEx()
nameEx("identifyPolymorphicRegion")
### * identifyPolymorphicRegion

flush(stderr()); flush(stdout())

### Name: identifyPolymorphicRegion
### Title: Identity polymorphic regions.
### Aliases: identifyPolymorphicRegion
### Keywords: methods

### ** Examples

## Not run: 
##D 
##D fcgr3PolymorphicRegion <- identifyPolymorphicRegion(Object = objectCNVrd2,
##D                                                     segmentObject = resultSegment, 
##D                                                     thresholdForPolymorphicRegions = c(0.75, 0.25),
##D                                                     plotLegend = FALSE)
##D 
##D 
## End(Not run)



cleanEx()
nameEx("plotCNVrd2")
### * plotCNVrd2

flush(stderr()); flush(stdout())

### Name: plotCNVrd2
### Title: Plot traces of samples.
### Aliases: plotCNVrd2
### Keywords: methods plots

### ** Examples

data(fcgr3bMXL)
##Obtain all information of CNVs
allGroups <- copynumberGroups$allGroups
###Obtain names of duplicate samples
duplicatedSamples <- rownames(allGroups[allGroups[, 2] > 2,])
###Plot the first duplicate samples
par(mfrow = c(3, 2))
for (ii in duplicatedSamples[1:6])  
plotCNVrd2(Object = objectCNVrd2,
           segmentObject = resultSegment,
           sampleName = ii)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plotPolymorphicRegion")
### * plotPolymorphicRegion

flush(stderr()); flush(stdout())

### Name: plotPolymorphicRegion
### Title: Plot polymorphic regions.
### Aliases: plotPolymorphicRegion
### Keywords: plot

### ** Examples

## Not run: 
##D 
##D plotPolymorphicRegion(Object = objectCNVrd2, polymorphicRegionObject = fcgr3PolymorphicRegion,
##D                       xlim = c(161300000, 161800000), drawThresholds = TRUE,
##D                       thresholdForPolymorphicRegions = c(0.75, 0.25))
##D 
##D 
##D 
##D ##Change thresholds
##D plotPolymorphicRegion(Object = objectCNVrd2, polymorphicRegionObject = fcgr3PolymorphicRegion,
##D                       xlim = c(161300000, 161800000), drawThresholds = TRUE,
##D                       thresholdForPolymorphicRegions = c(0.9, 0.1))
##D 
##D 
##D ##Plot standard deviation
##D 
##D plotPolymorphicRegion(Object = objectCNVrd2, polymorphicRegionObject = fcgr3PolymorphicRegion,
##D                       xlim = c(161300000, 161800000), typePlot = "SD",
##D                       thresholdForPolymorphicRegions = c(0.75, 0.25))
##D 
##D 
##D 
##D 
##D 
## End(Not run)



cleanEx()
nameEx("segmentSamples")
### * segmentSamples

flush(stderr()); flush(stdout())

### Name: segmentSamples
### Title: Implement the segmentation process
### Aliases: segmentSamples
### Keywords: methods

### ** Examples

data(fcgr3bMXL)
## Not run: resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
