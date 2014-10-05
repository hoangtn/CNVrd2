setGeneric("plotPolymorphicRegion",
           function(Object, ...){standardGeneric("plotPolymorphicRegion")})

setMethod("plotPolymorphicRegion", "CNVrd2",
          function(Object, polymorphicRegionObject = NULL,
                                      xlim = NULL, typePlot = c("All", "Quantile", "SD", "negativeSum", "positiveSum"),
                   xlab = NULL, ylabQuantile = NULL,
                   ylabSD = NULL, ylabNSum = NULL, ylabPSum = NULL,
                                      quantileValue = c(0.1, 0.5, 0.9),
                                      quantileColor = NULL, textSize = 0.7,
                   polymorphicCriteria = c("SD", "Quantile", "negativeSum", "positiveSum"),
		sdThreshold = 0.1,
                                      thresholdForPolymorphicRegions =
                                      c(0.975, 0.025),
                                      drawThresholds = FALSE,
                                      geneColor = 'lightpink',
                                      cex = 1.2, lwd = 1.5,
                                      verticalAdjustment = 0.3, yGeneNameSD = NULL,
                   yGeneNameQuantile = NULL, yGeneNameNSum = NULL, yGeneNamePSum = NULL,
                                      plotLegend = TRUE){
    
##Checking parameters
  st = Object@st
  en = Object@en
  if (!is.null(xlim)){
          if (xlim[1] < st)
              xlim[1] <- st
          if (xlim[2] > en)
              xlim[2] <- en
          outputST <- xlim[1]
          outputEND <- xlim[2]


    }
  else {
  outputST = st
  outputEND = en}
  
   windows = Object@windows
  chr = Object@chr
  genes <- Object@genes
  geneNames <- Object@geneNames
  
  genes <- matrix(genes, nrow = 2)

  
  geneNames = Object@geneNames
  quantileValue <- sort(quantileValue)

    nQuantile <- length(quantileValue)

  typePlot <- match.arg(typePlot)
  polymorphicCriteria <- match.arg(polymorphicCriteria)


  if (is.null(ylabQuantile))
      ylabQuantile <- "Quantile"
  if (is.null(ylabSD))
      ylabSD <- "SD"
  if (is.null(ylabNSum))
      ylabNSum <- "negative Sum"
  if (is.null(ylabPSum))
      ylabPSum <- "positive Sum"
  
  if (is.null(xlab))
      xlab <- paste(chr, ":", outputST, "-", outputEND, sep = "")
      
  thresholdForPolymorphicRegions <- sort(thresholdForPolymorphicRegions)

  if (is.null(quantileColor))
      quantileColor <- seq(length(quantileValue)) + 1
#####################################
  genes <- matrix(genes, nrow = 2)

  mxInterAll <- polymorphicRegionObject$subRegionMatrix
  subRegionData <- polymorphicRegionObject$subRegion
  mSD <- polymorphicRegionObject$mSD
  nCol <- dim(subRegionData)[1]
  
  mQuantile <- apply(mxInterAll, 2, function(x)
                     quantile(x, quantileValue))

  negativeSum <- apply(mxInterAll, 2, function(x) sum(x[x < 0]))
  positiveSum <- apply(mxInterAll, 2, function(x) sum(x[x > 0]))

  
  ###################Plot####################################
#########Identify max, min of all quantile values
  minQuantile <- min(apply(mQuantile, 1, min))
  maxQuantile <- max(apply(mQuantile, 1, max))

  ###Function to transform data############################
  transtoDataFrame <- function(x, y, colour = 1){

    tempF <- data.frame(x1= sort(c(x[, 1], x[, 2])),
               x2 = as.numeric(apply(matrix(y, ncol = 1), 1,
               function(x) rep(x, 2))))
    tempF <- data.frame(tempF, Quantile = rep(colour, dim(tempF)[1]))


    return(tempF)

    }
  ########################################


  dfQuantile <- NULL
    ###Make a list to store quantile values
  listQ <- list()
  for (ii in 1:length(quantileValue)){
      mQ <- data.frame(subRegionData, mQuantile[ii, ])
      listQ[[ii]] <- mQ
      index1 <- as.integer(rownames(mQ[(mQ[, 1] <= outputST) & (mQ[, 2]>=outputST), ]))[1]
          
      index2 <- as.integer(rownames(mQ[(mQ[, 1] <= outputEND) & (mQ[, 2]>=outputEND), ]))[1]

      mQ <- mQ[index1:index2, ]
      dfQuantile <- rbind(dfQuantile, transtoDataFrame(x = mQ[, c(1, 2)],
                                                       y = mQ[, 3],
                                                       colour = as.character(rownames(mQuantile)[ii])))
  }
  

  #########Plot the region#############################
  thresholdsTogetPolymorphicRegions <- c(quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]),
                                         quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]))
      p1 <- ggplot() + geom_line(data= dfQuantile, aes(x=x1, y=x2, group = Quantile, colour=Quantile), size=cex) +
          coord_cartesian(xlim = c(outputST, outputEND)) + theme(legend.position="top") + xlab(xlab) + ylab(ylabQuantile)
##############Plot genes
      if (!is.null(genes)){
          genesDataFrame <- data.frame(xmin = genes[1, ], xmax = genes[2, ],
                                 ymin = rep(min(dfQuantile[, 2]) - 0.1, length(genes[1, ])),
                                 ymax = rep(max(dfQuantile[, 2]) + 0.1 , length(genes[1, 1])))

          p1 <- p1 + geom_rect(data = genesDataFrame,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),alpha = 0.4, col = 'pink', fill = 'pink')
           }
      if (!is.null(geneNames)){
          if (is.null(yGeneNameQuantile))
              yGeneNameQuantile <- min(dfQuantile[, 2]) + 0.05
          geneNamesDF <- data.frame(x = genes[1, ],
                                y = rep(yGeneNameQuantile , length(genes[1, 1])),
                                label = geneNames)
      p1 <- p1 + geom_text(data = geneNamesDF, aes(x = x, y = y, label = label), angle = 90, size =textSize)
          }
      if (drawThresholds == TRUE){
          dfQuantileThreshold1 <- data.frame(x = c(mQ[1, 1],  mQ[dim(mQ)[1], 2]),
                                             y = c(rep(thresholdsTogetPolymorphicRegions[1], 2)))
          dfQuantileThreshold2 <- data.frame(x = c(mQ[dim(mQ)[1], 2],mQ[1, 1]),
                                             y = rep(thresholdsTogetPolymorphicRegions[2], 2))

          p1 <- p1 + geom_line(data = dfQuantileThreshold1, aes(x = x, y = y), alpha = 0.3, fill = 'green', col = 'green', size = 0.8*cex)  +
              geom_line(data = dfQuantileThreshold2, aes(x = x, y = y), alpha = 0.3, fill = 'green', col = 'green', size = 0.8*cex)
          }


  #############Plot regions' sds####################################
      mSD <- mSD[index1:index2,]
      dfSD <- data.frame(transtoDataFrame(x = mSD[, c(1, 2)],
                                          y = mSD[, 3]))

      p2 <- ggplot() + geom_line(data= dfSD, aes(x=x1, y=x2), size=cex, col = 'blue') +
          coord_cartesian(xlim = c(outputST, outputEND)) + xlab(xlab) + ylab(ylabSD)
      if (!is.null(genes)){
          genesDataFrame <- data.frame(xmin = genes[1, ], xmax = genes[2, ],
                                       ymin = rep(min(dfSD[, 2]) - 0.1, length(genes[1, ])),
                                       ymax = rep(max(dfSD[, 2]) + 0.1 , length(genes[1, 1])))
          p2 <- p2 + geom_rect(data = genesDataFrame,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),alpha = 0.4, col = 'pink', fill = 'pink')
           }
      if (!is.null(geneNames)){
          if (is.null(yGeneNameSD))
              yGeneNameSD <- min(dfSD[, 2]) +  0.05
          geneNamesDF <- data.frame(x = genes[1, ],
                                    y = rep(yGeneNameSD , length(genes[1, 1])),
                                    label = geneNames)

          p2 <- p2 + geom_text(data = geneNamesDF, aes(x = x, y = y, label = label), angle = 90, size = textSize)

           if (drawThresholds == TRUE){
                dfSDtemp <- data.frame(x = c(mSD[1, 1],  mSD[dim(mSD)[1], 2]),
                                             y = c(rep(quantile(mSD[, 3], 1 - sdThreshold), 2)))
                 p2 <- p2 + geom_line(data = dfSDtemp, aes(x = x, y = y), alpha = 0.3, fill = 'green', col = 'green', size = 0.8*cex)
                
                

           }
          }
  ######################################################################
  #################Plot negative Sum####################################

  mNSum <- data.frame(subRegionData, negativeSum)
  mNSum <- mNSum[index1:index2,]
      dfNSum <- data.frame(transtoDataFrame(x = mNSum[, c(1, 2)],
                                          y = mNSum[, 3]))

      pNSum <- ggplot() + geom_line(data= dfNSum, aes(x=x1, y=x2), size=cex, col = 'blue') +
          coord_cartesian(xlim = c(outputST, outputEND)) + xlab(xlab) + ylab(ylabNSum)
      if (!is.null(genes)){
          genesDataFrame <- data.frame(xmin = genes[1, ], xmax = genes[2, ],
                                       ymin = rep(min(dfNSum[, 2]) - 0.1, length(genes[1, ])),
                                       ymax = rep(max(dfNSum[, 2]) + 0.1 , length(genes[1, 1])))
          pNSum <- pNSum + geom_rect(data = genesDataFrame,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),alpha = 0.4, col = 'pink', fill = 'pink')
           }
      if (!is.null(geneNames)){
          if (is.null(yGeneNameNSum))
              yGeneNameNSum <- min(dfNSum[, 2]) +  0.05
          geneNamesDF <- data.frame(x = genes[1, ],
                                    y = rep(yGeneNameNSum , length(genes[1, 1])),
                                    label = geneNames)

          pNSum <- pNSum + geom_text(data = geneNamesDF, aes(x = x, y = y, label = label), angle = 90, size = textSize)

          }
  

   #################Plot positive Sum####################################

  mPSum <- data.frame(subRegionData, positiveSum)
  mPSum <- mPSum[index1:index2,]
      dfPSum <- data.frame(transtoDataFrame(x = mPSum[, c(1, 2)],
                                          y = mPSum[, 3]))

      pPSum <- ggplot() + geom_line(data= dfPSum, aes(x=x1, y=x2), size=cex, col = 'blue') +
          coord_cartesian(xlim = c(outputST, outputEND)) + xlab(xlab) + ylab(ylabPSum)
      if (!is.null(genes)){
          genesDataFrame <- data.frame(xmin = genes[1, ], xmax = genes[2, ],
                                       ymin = rep(min(dfPSum[, 2]) - 0.1, length(genes[1, ])),
                                       ymax = rep(max(dfPSum[, 2]) + 0.1 , length(genes[1, 1])))
          pPSum <- pPSum + geom_rect(data = genesDataFrame,
                               aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),alpha = 0.4, col = 'pink', fill = 'pink')
           }
      if (!is.null(geneNames)){
          if (is.null(yGeneNamePSum))
              yGeneNamePSum <- min(dfPSum[, 2]) +  0.05
          geneNamesDF <- data.frame(x = genes[1, ],
                                    y = rep(yGeneNamePSum , length(genes[1, 1])),
                                    label = geneNames)

          pPSum <- pPSum + geom_text(data = geneNamesDF, aes(x = x, y = y, label = label), angle = 90, size = textSize)

          }
 

  ########################################################################
    if (typePlot == "Quantile")
        print(p1)
  else {if (typePlot == "SD")
            
      print(p2)
  else { if (typePlot == "negativeSum")
             print(pNSum)
  else { if (typePlot == "positiveSum")
             print(pPSum)
     
  else {
      library('gridExtra')
      sidebysideplot <- grid.arrange(p2, p1, pPSum, pNSum, ncol=1, nrow = 4)
  }}}}

        


            
  ######################Identity polymorphic regions##################
  ##Extract the first and the last rows of quantile matrixes
  if (polymorphicCriteria == "Quantile"){
      lowBoundary   <- listQ[[1]][listQ[[1]][, 3] <= quantile(listQ[[1]][, 3],
                                            thresholdForPolymorphicRegions[1]), ]
      highBoundary <- listQ[[nQuantile]][listQ[[nQuantile]][, 3] >=
                                     quantile(listQ[[nQuantile]][, 3],
                                              thresholdForPolymorphicRegions[2]), ]

      lowBoundary <- reduce(IRanges(lowBoundary[, 1], lowBoundary[, 2]))
      highBoundary <- reduce(IRanges(highBoundary[, 1], highBoundary[, 2]))
      unionBoundary <- reduce(union(lowBoundary, lowBoundary))
  }
  else {
      tempSD <- as.data.frame(mSD)
      tempSD <- tempSD[tempSD[, 3] >= quantile(tempSD[, 3], 1 - sdThreshold),]
      unionBoundary <- reduce(IRanges(tempSD[, 1], tempSD[, 2]))


  }

  end(unionBoundary) <- end(unionBoundary) - 1

  
   return(putativeBoundary = unionBoundary )



})
