
## ----setup, include=FALSE, cache=FALSE, eval = TRUE----------------------
library(knitr)
opts_chunk$set(fig.path='./figures/CNVrd2-', 
               fig.align='center', fig.show='asis', 
               eval = TRUE, fig.width = 6,
               fig.height = 6, size='small')
options(replace.assign=TRUE,width=90)


## ----options,echo=FALSE--------------------------------------------------
options(width=72)


## ----warning=FALSE, results="hide"---------------------------------------
library('CNVrd2')


## ----warning=FALSE, results="hide", eval=FALSE---------------------------
## objectCNVrd2 <- new("CNVrd2", windows = 1000, chr = "chr1",
##                     st = 161100001, en = 162100000,
##                     dirBamFile = "BamMXL",
##                     genes = c(161592986, 161601753),
##                     geneNames = "3B")
## 


## ----warnings=FALSE, results="hide", eval=FALSE--------------------------
## readCountMatrix <- countReadInWindow(Object = objectCNVrd2, correctGC = TRUE)
## 


## ----warning=FALSE, results="hide", fig.keep='none', eval=FALSE----------
## ##Obtain segmentation scores
## resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)


## ----'fcgr3bSS', fig.show='asis', results="hide", fig.cap='FCGR3B segmentation score.'----
##Load data into R
data(fcgr3bMXL)
##Reload readCountMatrix
readCountMatrix <- resultSegment$stdCntMatrix
##Take a quick look the data
readCountMatrix[1:2, 1:2]
##Make a CNVrd2 object 
objectCNVrd2 <- new("CNVrd2", windows = 1000, chr = "chr1",
                    st = 161100001, en = 162100000,
                    dirBamFile = "BamMXL",
                    genes = c(161592986, 161601753),
                    geneNames = "3B")
##Obtain segmentation scores
resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)
##View these segmentation results
sS <- resultSegment$segmentationScores
hist(sS[, 1], 100, xlab = 'Segmentation score', main = '')


## ----'fcgr3b4CNs', fig.show='asis', results="hide",warning=FALSE, fig.cap='FCGR3B CN groups.'----
objectCluster <- new("clusteringCNVs",
                     x = resultSegment$segmentationScores[, 1], k = 4, EV = TRUE)
#Cluster into 4 groups
copynumberGroups <- groupCNVs(Object = objectCluster)



## ----warning=FALSE-------------------------------------------------------
copynumberGroups$allGroups[1:3, ]


## ----'fcgr3b3CNs', fig.show='asis', results="hide",warning=FALSE, fig.cap = 'FCGR3B CN groups (rightLimit = 1.5).'----
#Set right limit = 1.5 to make values > 1.5 be into the largest group.
objectCluster <- new("clusteringCNVs",
                     x = resultSegment$segmentationScores[, 1], k = 3, EV = TRUE)
copynumberGroups <- groupCNVs(Object = objectCluster, rightLimit = 1.5)


## ----'fcgr3bdupSamples', fig.show='asis', warning=FALSE, fig.cap='MXL duplicated samples.', size='small'----
allGroups <- copynumberGroups$allGroups
###Obtain names of duplicate samples
duplicatedSamples <- rownames(allGroups[allGroups[, 2] > 2,])
###Plot 6 duplicate samples
par(mfrow = c(3, 2))
for (ii in duplicatedSamples[1:6])
    plotCNVrd2(Object = objectCNVrd2,
                         
               segmentObject = resultSegment,
               
                             
               sampleName = ii)



## ----warning=FALSE, results="hide"---------------------------------------
##Obtain VCF-file information in CNVrd2 package
vcfFile <- system.file(package="CNVrd2", "extdata",
                       "chr1.161600000.161611000.vcf.gz")
##Make a data frame named sampleCNV including samples, CNs, population names
sampleCNV <- data.frame(copynumberGroups$allGroups[, c(1,2) ],rep("MXL", dim(copynumberGroups$allGroups)[1]))
rownames(sampleCNV) <- substr(sampleCNV[, 1], 1, 7)
sampleCNV[, 1] <- rownames(sampleCNV)
##The first column must be the sample names and some samples should be in the vcf file
tagSNPandINDELofMXL <- calculateLDSNPandCNV(sampleCNV = sampleCNV,
                                            vcfFile = vcfFile, cnvColumn = 2,
                                            population = "MXL", popColumn = 3,
                                            nChunkForVcf = 5, chr = "1",
                                            st = 161600001, en = 161611000,
                                            codeSNP= "Three", codeCNV = "ThreeGroup")


## ----warning=FALSE-------------------------------------------------------
head(tagSNPandINDELofMXL)


## ----'ccl3l1Histogram', warning=FALSE, fig.cap = 'CCL3L1 segmentation score.'----
##Load data into R:
data(ccl3l1data)
head(ccl3l1data)
hist(ccl3l1data$SS, 100)


## ----'EUccl3l1Histogram', warning=FALSE, fig.cap = 'European-ancestry segmentation score.'----
xyEuro <- ccl3l1data[grep("CEU|TSI|IBS|GBR|FIN", ccl3l1data[, 2]), ]
yEuro <- xyEuro[, 3]
names(yEuro) <- rownames(xyEuro)
hist(yEuro, 100, xlab = '', main = '')



## ----'EUccl3l1results', warning=FALSE, fig.cap = 'Clustering results of European-ancestry sample sets.', results="hide"----
##Clustering European segmentation 
##scores into group: 5 groups were chosen

objectClusterEuroCCL3L1 <- new("clusteringCNVs", x = yEuro, k = 5)

europeanCCL3L1Groups <- groupCNVs(Object = objectClusterEuroCCL3L1)



## ----warning=FALSE, results="hide"---------------------------------------
#Means
lambda0 <- as.numeric(europeanCCL3L1Groups$m)
#SD
sdEM <- as.numeric(europeanCCL3L1Groups$sigma)
#Proportions
pEM <- as.numeric(europeanCCL3L1Groups$p)


## ----warning=FALSE-------------------------------------------------------
lambda0
sdEM
pEM
###Calculate the distances between groups
for (ii in 2:5){print(lambda0[ii] - lambda0[ii-1])}

###All segmentation scores
ccl3l1X <- ccl3l1data$SS
names(ccl3l1X) <- as.character(ccl3l1data$Name)
range(ccl3l1X)



## ------------------------------------------------------------------------
##Set prior information:
#prior for the sd of the means of groups: 
#5 was set for the third group = 2 CN
sd <- c(1, 1, 5, 1, 1) 
ccl3l1X <- sort(ccl3l1X)
###Data
xData <- ccl3l1X
###Number of groups
nGroups <- 10 
###prior for means of groups
lambda0 <- lambda0 
###Prior for mixing proportions
alpha0 <-  c(3, 29, 44, 18, 7,  5, rep(2, nGroups -length(pEM) -1))
##Prior for the distances between groups
distanceBetweenGroups = 0.485

sdEM = sdEM



## ------------------------------------------------------------------------
##Adjust standard deviation for the fifth group
sdEM[5] <- sdEM[4]
 


## ----eval = FALSE, echo = TRUE-------------------------------------------
## set.seed(123)
## groupCCL3L1allPops <- groupBayesianCNVs(xData = xData, nGroups = nGroups,
##                                         lambda0 = lambda0,
##                                         sd0 = sdEM, alpha0 = alpha0,
##                                         distanceBetweenGroups = distanceBetweenGroups,
##                                         sdOftau = sd,
##                                         rightLimit = 4)


## ----warning=FALSE, results="hide"---------------------------------------
rownames(ccl3l1data) <- ccl3l1data[, 1]



## ----warning=FALSE, results="hide"---------------------------------------
##Obtain vcf-file information in CNVrd2
vcfFileCCL3L1 <- system.file(package="CNVrd2", "extdata",
                       "chr17.34800000.34830000.vcf.gz")
##Set populations we would like to identify tagSNPs
allPops <- c("TSI", "CEU", "GBR", "FIN", "IBS")


## ----warning = FALSE, results="hide"-------------------------------------
##Identify tag SNPs/INDELs
tagSNPandINDELofCCL3L1 <- calculateLDSNPandCNV(sampleCNV = ccl3l1data,
                                            vcfFile = vcfFileCCL3L1, cnvColumn = 4,
                                            population = allPops, popColumn = 2,
                                            nChunkForVcf = 5, chr = "17",
                                            st = 34800000, en = 34830000 )



## ----warnings=FALSE, size='tiny'-----------------------------------------
lapply(tagSNPandINDELofCCL3L1, head)


## ----warnings=FALSE------------------------------------------------------
fcgr3PolymorphicRegion <- identifyPolymorphicRegion(Object = objectCNVrd2,
                                                    segmentObject = resultSegment, 
                                                    plotLegend = FALSE)



## ----'FCGR3polymorphc1', warning=FALSE, fig.cap = 'CN polymorphic region at FCGR3 locus, represented by quantiles of the distribution of segmentation scores across samples.'----
plotPolymorphicRegion(Object = objectCNVrd2, polymorphicRegionObject = fcgr3PolymorphicRegion,
                      xlim = c(161300000, 161800000), drawThresholds = TRUE,
                      typePlot = "SD")



## ----'FCGR3polymorphc2', warning=FALSE, fig.cap = 'CN polymorphic region at FCGR3 locus, represented by quantiles of the distribution of segmentation scores across samples.'----
plotPolymorphicRegion(Object = objectCNVrd2, polymorphicRegionObject = fcgr3PolymorphicRegion,
                      xlim = c(161300000, 161800000), sdThreshold = 0.05, drawThresholds = TRUE,
                      typePlot = "SD")


## ------------------------------------------------------------------------
sessionInfo()


