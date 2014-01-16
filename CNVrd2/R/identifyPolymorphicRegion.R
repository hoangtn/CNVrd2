setGeneric("identifyPolymorphicRegion",
           function(Object, ...){standardGeneric("identifyPolymorphicRegion")})

setMethod("identifyPolymorphicRegion", "CNVrd2",
          function(Object, segmentObject = NULL,
                                      xlim = NULL,
                                      quantileValue = c(0.1, 0.5, 0.9),
                                      quantileColor = NULL,
                                      thresholdForPolymorphicRegions = 
                                      c(0.975, 0.025),
                                      drawThresholds = FALSE,
                   VstTest = FALSE,
                   popName = NULL,
                   thresholdVST = NULL,

                                      geneColor = 'lightpink',
                                      cex = 1.5, lwd = 1.5,
                                      verticalAdjustment = 0.3,
                                      plotLegend = TRUE,
                   plotPolymorphicRegion = TRUE){
    
##Check parameters
  st = Object@st
  en = Object@en
  if (!is.null(xlim)){
    if ((xlim[1] < st) | (xlim[2] > en))
      stop("Please choose coordinates in the region")

    outputST <- xlim[1]
    outputEND <- xlim[2]
    }
  else {
  outputST = st
  outputEND = en}
  
   windows = Object@windows
  chr = Object@chr
  genes <- Object@genes
  geneNames = Object@geneNames
  quantileValue <- sort(quantileValue)

    nQuantile <- length(quantileValue)


  thresholdForPolymorphicRegions <- sort(thresholdForPolymorphicRegions)

  if (is.null(quantileColor))
      quantileColor <- seq(length(quantileValue)) + 1
#####################################
  genes <- matrix(genes, nrow = 2)

  segmentResults <- segmentObject$segmentResults
  sampleid <- sapply(segmentResults, function(x) x[1, 1])
  cna.out <- do.call(rbind, segmentResults)

  cna.out[, 4] <- cna.out[, 4] + windows

  ##Obtain sub-regions#########
  subRegionData <- sort(unique(c(cna.out[, 3], cna.out[, 4])))
  subRegionData <- data.frame(subRegionData[-length(subRegionData)], subRegionData[-1])

###Segmentation scores for sub-regions
  subRegionMatrix <- matrix(0, nrow = length(sampleid),
                          ncol = dim(subRegionData[1]))

  for (ii in 1:length(sampleid)){
      subCNA <- cna.out[grep(paste("^", sampleid[ii], "$", sep = ""), cna.out[, 1]), c(3, 4, 6)]

    for (jj in 1:dim(subRegionData)[1]){
        subRegionMatrix[ii, jj] <-  subCNA[(subCNA[, 1] <= subRegionData[jj, 1]) &(subCNA[, 2] >= subRegionData[jj, 2]), 3]

        }}


#########Standard deviation for sub-regions
  sdQ <- apply(subRegionMatrix, 2, sd)
  mSD <- data.frame(subRegionData, sdQ)
  colnames(mSD) <- c("st", "en", "sd")


  rownames(subRegionMatrix) <- sampleid
  nCol <- dim(subRegionMatrix)[2]
  ###########Make a matrix of quantiles
  mQuantile <- matrix(0, ncol = nCol, nrow = length(quantileValue))
  mQuantile <- apply(subRegionMatrix, 2, function(x) quantile(x, quantileValue))

#########Identify max, min of all quantile values
  minQuantile <- min(apply(mQuantile, 1, min))
  maxQuantile <- max(apply(mQuantile, 1, max))


  ###Make a list to store quantile values
  listQ <- list()
  for (ii in 1:length(quantileValue)){
      mQ <- data.frame(subRegionData, mQuantile[ii, ])
      listQ[[ii]] <- mQ
      mQ <- mQ[(mQ[, 1] >= outputST) & (mQ[, 2] <=outputEND ), ]
  }

  thresholdsTogetPolymorphicRegions <- c(quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]),
                                         quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]))

  
  ######################Identify polymorphic regions##################
  ##Extract the first and the last rows of quantile matrixes
  lowBoundary   <- listQ[[1]][listQ[[1]][, 3] <= quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]), ]
  highBoundary <- listQ[[nQuantile]][listQ[[nQuantile]][, 3] >=
                                     quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]), ]

  lowBoundary <- reduce(IRanges(lowBoundary[, 1], lowBoundary[, 2]))
  highBoundary <- reduce(IRanges(highBoundary[, 1], highBoundary[, 2]))

  unionBoundary <- reduce(union(lowBoundary, highBoundary))


  ##########Calculate Vst #################################################
  valueVST <- NULL
  if (VstTest){
      message("Calculate Vst")
      popAll <- as.character(popName)

      popNames <- names(table(popAll))
      nPops <- length(popNames)

      valueVST <- apply(subRegionMatrix, 2, function(x){
          bTemp <- data.frame(data = x, pop = popAll)
          nameValueVST <- c()
          vResults <- c()
          indexVST = 1

          for (iST in 1:(nPops-1)){
              for (jST in (iST+1):nPops){
                  bTempIJ <- bTemp[(bTemp[, 2] == popNames[iST]) |
                             (bTemp[, 2] == popNames[jST]), ]
                  bVT <- sd(bTempIJ[, 1])

                  bVSTtemp <- sapply(split(bTempIJ, bTempIJ$pop), function(x) sd(x[, 1]))

                  bVSTtemp <- data.frame(bVSTtemp, table(popAll))
                  bVSTtemp <- bVSTtemp[!is.na(bVSTtemp[, 1]), ]


                  bVS <- sum(bVSTtemp[, 1]*bVSTtemp[, 3])/sum(bVSTtemp[, 3])
                  vResults[indexVST] <- (bVT - bVS)/bVT

                  nameValueVST[indexVST] <- paste(popNames[iST], "-", popNames[jST], sep = "")
                  indexVST <- indexVST + 1
              }}
          names(vResults) <- nameValueVST

          return(vResults)
          })

      ########################
      maxVST <- valueVST
      if (is.matrix(maxVST))
	      maxVST <- apply(valueVST, 2, function(x) max(x))
      	
      if (is.null(thresholdVST))
          thresholdVST <- quantile(maxVST, 0.95)
      regionVST <- data.frame(subRegionData, maxVST)

      shortregionVST <- regionVST[regionVST[, 3] >= thresholdVST, ]
      
      shortregionVST <- reduce(IRanges(shortregionVST[, 1], shortregionVST[, 2]))
      #########Plot######################################

      ###############################################################################
      ##############################################################################

      unionBoundary <- reduce(intersect(shortregionVST, unionBoundary))

  }
  unionBoundary <- unionBoundary[width(unionBoundary) > 1]
  #########Obtain segmentation scores of putative regions#####################################
  message("Calculate segmentation scores for polymorphic regions")
  SSofunionBoundary <- as.data.frame(unionBoundary)[, c(1, 2)]
  SSofunionBoundary <- data.frame(SSofunionBoundary, matrix(0, nrow = dim(SSofunionBoundary)[1],
                                                          ncol = dim(subRegionMatrix)[1]))
  b1 <- as.data.frame(subRegionData)
  for (ii in 1:dim(SSofunionBoundary)[1]){
    b11 <- b1[b1[, 1] >= SSofunionBoundary[ii, 1], ]
    b11 <- b11[b11[, 2] <= SSofunionBoundary[ii, 2], ]
    mA11 <- as.matrix(subRegionMatrix[, as.numeric(rownames(b11))])
    mA11s <- apply(mA11, 1, mean)
    SSofunionBoundary[ii, -c(1, 2)] <- mA11s
    }

  colnames(SSofunionBoundary)[-c(1, 2)] <- rownames(subRegionMatrix)
  SSofunionBoundary[, 2] <- SSofunionBoundary[, 2] -1
  end(unionBoundary) <- end(unionBoundary) - 1
  #########################################################################
  if (plotPolymorphicRegion){

      polymorphicRegionObject <- list(putativeBoundary = unionBoundary, SSofPolymorphicRegions = SSofunionBoundary,
                                      subRegionMatrix = subRegionMatrix,
                                      subRegion = subRegionData,
                                      mQuantile = mQuantile,
                                      Vst = valueVST,
                                      mSD = mSD)
      plotPolymorphicRegion(Object = Object,
                            polymorphicRegionObject = polymorphicRegionObject,
                            xlim = c(outputST, outputEND),
                            drawThresholds = drawThresholds,
                            thresholdForPolymorphicRegions = thresholdForPolymorphicRegions ,
                            quantileValue = quantileValue,
                            quantileColor = quantileColor,
                            geneColor = geneColor,
                            cex = cex, lwd = lwd,
                            verticalAdjustment = verticalAdjustment,
                            plotLegend = plotLegend)

      }
      

  #########################################################################

  return(list(putativeBoundary = unionBoundary, SSofPolymorphicRegions = SSofunionBoundary,
              subRegionMatrix = subRegionMatrix,
              subRegion = subRegionData,
              mQuantile = mQuantile,
              mSD = mSD, 
              Vst = valueVST))
  
  
})
      

