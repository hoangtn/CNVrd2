setGeneric("plotPolymorphicRegion",
           function(Object, ...){standardGeneric("plotPolymorphicRegion")})

setMethod("plotPolymorphicRegion", "CNVrd2",
          function(Object, polymorphicRegionObject = NULL,
                                      xlim = NULL,
                                      quantileValue = c(0.1, 0.5, 0.9),
                                      quantileColor = NULL,
                                      thresholdForPolymorphicRegions =
                                      c(0.975, 0.025),
                                      drawThresholds = FALSE,
                                      geneColor = 'lightpink',
                                      cex = 1.5, lwd = 1.5,
                                      verticalAdjustment = 0.3,
                                      plotLegend = TRUE){
    
##Checking parameters
  st = Object@st
  en = Object@en
  if (!is.null(xlim)){
    if ((xlim[1] < st) | (xlim[2] > en))
      stop("Please choose coordinates in the region")
     if ((xlim[1] > en ) | (xlim[2] < st ))
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

  mxInterAll <- polymorphicRegionObject$subRegionMatrix
  subRegionData <- polymorphicRegionObject$subRegion
  nCol <- dim(subRegionData)[1]
  mQuantile <- polymorphicRegionObject$mQuantile


  ###################Plot####################################
#########Identify max, min of all quantile values
  minQuantile <- min(apply(mQuantile, 1, min))
  maxQuantile <- max(apply(mQuantile, 1, max))

  #########Plot the region#############################
  plot(subRegionData[, 1], mQuantile[1, ], col = 'white',
       xlab = paste(chr, ":", outputST, "-", outputEND, sep = ""),
       ylab = "Segmentation score",
       xlim = c(outputST, outputEND),
       ylim = c(minQuantile - verticalAdjustment, maxQuantile + verticalAdjustment))

  ##############Plot genes
   for(k in 1:ncol(genes)){
       rect(genes[1, k], minQuantile - verticalAdjustment - 0.3,
            genes[2, k], maxQuantile + verticalAdjustment + 0.3,
            col= geneColor)
       
            if (is.null(geneNames))
                text(genes[2, k],
                     maxQuantile + verticalAdjustment,
                     paste("Gene", k, sep = ""), 
                     cex = cex, lwd = lwd)
            else
              text(genes[2, k],
                   maxQuantile + verticalAdjustment,
                   geneNames[k], 
                   cex = 0.4*cex, lwd = 0.4*lwd
                   )
            }



  ###Make a list to store quantile values
  listQ <- list()
  for (ii in 1:length(quantileValue)){
      mQ <- data.frame(subRegionData, mQuantile[ii, ])
      listQ[[ii]] <- mQ
      
      index1 <- as.integer(rownames(mQ[(mQ[, 1] <= outputST) & (mQ[, 2]>=outputST), ]))[1]
          
      index2 <- as.integer(rownames(mQ[(mQ[, 1] <= outputEND) & (mQ[, 2]>=outputEND), ]))[1]
      

      mQ <- mQ[index1:index2, ]

      apply(mQ, 1, function(x)
            lines(c(x[1], x[2]), c(x[3], x[3]), col = quantileColor[ii],
                  lwd = lwd, cex = cex))
  }

  thresholdsTogetPolymorphicRegions <- c(quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]),
                                         quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]))
  if (drawThresholds)
      abline(h = thresholdsTogetPolymorphicRegions, col = 'grey', lty = 3)
  
  if (plotLegend)
      legend("topright", paste("Quantile: ", quantileValue, sep = ""), lty = 1,
         col = quantileColor)

  
  ######################Identity polymorphic regions##################
  ##Extract the first and the last rows of quantile matrixes
  lowBoundary   <- listQ[[1]][listQ[[1]][, 3] <= quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]), ]
  highBoundary <- listQ[[nQuantile]][listQ[[nQuantile]][, 3] >=
                                     quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]), ]

  lowBoundary <- reduce(IRanges(lowBoundary[, 1], lowBoundary[, 2]))
  highBoundary <- reduce(IRanges(highBoundary[, 1], highBoundary[, 2]))
  unionBoundary <- reduce(union(lowBoundary, lowBoundary))

  end(unionBoundary) <- end(unionBoundary) - 1


   return(putativeBoundary = unionBoundary)



})
