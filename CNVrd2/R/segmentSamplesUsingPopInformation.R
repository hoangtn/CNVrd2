#Using population information to segment multiple samples
#pops: a data.frame including 2 columns: First column is
#population information and Second column is samples' names
#which are the same as the names in the matrix of read counts.

#Read-count matrix must be raw read-count matrix (not transformed
#or standardized

#ss: segmentation score

setGeneric("segmentSamplesUsingPopInformation",
           function(Object, ...){standardGeneric("segmentSamplesUsingPopInformation")})


setMethod("segmentSamplesUsingPopInformation", "CNVrd2",
          function(Object, rawReadCountMatrix, pops = NULL,
                   entireGene = FALSE, inputBamFile = FALSE,
                    bThresholds = NULL,
		alpha = 0.01, nperm = 10000, p.method = "hybrid", 
		   min.width = 2, kmax = 25, nmin = 200, eta = 0.05, 
		   trim = 0.025, undo.splits = "none", 
	           undo.prune = 0.05, undo.SD = 3, verbose = 1){

              
                pops <- pops
                entireGene <- entireGene
                inputBamFile <- inputBamFile
                alpha <- alpha
                min.width <- min.width
                bThresholds <- bThresholds
                ObjectAllPops <- Object
                
                rawReadCountMatrix <- rawReadCountMatrix

                

                if (!is.null(pops)){

                    pops <- pops[order(pops[, 1]),]
                    rawReadCountMatrix <- rawReadCountMatrix[pmatch(pops[, 2], rownames(rawReadCountMatrix)),]
                    

                }

                
                stdCntMatrix <- apply(rawReadCountMatrix, 1, function(x)
                                      x/median(x))
                stdCntMatrix <- apply(stdCntMatrix, 1, function(x){
                        if (sd(x)==0)
                            ab = x
                        else
                            ab = (x - median(x))/sd(x)
                        return(ab)
                        })

                
 ################Obtain segmentation scores for all populations                    
                allNotAdjustedResults <-
                    segmentSamples(Object = ObjectAllPops , stdCntMatrix = stdCntMatrix, entireGene = entireGene,
                               inputBamFile = inputBamFile, bThresholds = bThresholds,
                               alpha = alpha, min.width = min.width, 
                                nperm = nperm, p.method = p.method,
			kmax = kmax, nmin = nmin, eta = eta,
			trim = trim, undo.splits = undo.splits, 
                   undo.prune = undo.prune, undo.SD = undo.SD, verbose = verbose)

                ssAllPopResults <- allNotAdjustedResults$segmentationScores
                ssAllPopResultsAdjusted <- ssAllPopResults
#################Segment for single population################
                
                if (!is.null(pops)){
                    ssSinglePopResults <- list()
                    popNames <- names(table(pops[, 1]))
                    for (singlePop in popNames){
                        samplePopNames <- as.character(pops[grep(singlePop, pops[, 1]), 2])

                        popRawMatrixReadCount <- rawReadCountMatrix[pmatch(samplePopNames, rownames(rawReadCountMatrix)), ]

                        popStdMatrixReadCount <- apply(popRawMatrixReadCount, 1, function(x)
                                                      x/median(x))
                        popStdMatrixReadCount <- apply(popStdMatrixReadCount, 1, function(x){
                            if (sd(x)==0)
                                ab = x
                            else
                                ab = (x - median(x))/sd(x)
                            return(ab)
                            })

                        ssSinglePopResults[[singlePop]] <-
                            segmentSamples(Object = ObjectAllPops, stdCntMatrix = popStdMatrixReadCount,
                                           entireGene = entireGene, inputBamFile = FALSE, bThresholds = bThresholds,
					alpha = alpha, min.width = min.width, nperm = nperm, p.method = p.method,
					kmax = kmax, nmin = nmin, eta = eta,
					trim = trim, undo.splits = undo.splits, 
			                   undo.prune = undo.prune, undo.SD = undo.SD, verbose = verbose)$segmentationScores
                        
		ssSinglePopFromAllPops <- ssAllPopResults[pmatch(rownames(ssSinglePopResults[[singlePop]]),
                                                                         rownames(ssAllPopResults)), ]

                        if (is.null(dim(ssSinglePopFromAllPops))){
                            ssSinglePopAdjusted <- fitted(lm(ssSinglePopFromAllPops ~ ssSinglePopResults[[singlePop]][, 1]))
                            names(ssSinglePopAdjusted) <- names(ssSinglePopFromAllPops)

                            ssAllPopResultsAdjusted[pmatch(names(ssSinglePopAdjusted), rownames(ssAllPopResultsAdjusted)), 1] <- ssSinglePopAdjusted
                            

                            
                        }
                        else
                            for (colGene in 1:dim(ssSinglePopFromAllPops)[2]){
                                ssSinglePopAdjusted <- fitted(lm(ssSinglePopFromAllPops[, colGene] ~ ssSinglePopResults[[singlePop]][, colGene]))
                                names(ssSinglePopAdjusted) <- rownames(ssSinglePopFromAllPops)
                                ssAllPopResultsAdjusted[pmatch(names(ssSinglePopAdjusted), rownames(ssAllPopResultsAdjusted)), colGene] <- ssSinglePopAdjusted
                                                      
                            

                        }
                                              }
                }
                    
################Using a linear regression to adjust segmentation scores for
###multiple populations
                allResults <- list(segmentationScores = ssAllPopResultsAdjusted,
                                   segmentationScoresFromSinglePops = ssSinglePopResults,
                                   unAdjustedsegmentationScores = allNotAdjustedResults$segmentationScores,
                                   segmentResults = allNotAdjustedResults$segmentResults,
                                   observedReadCountRatios = allNotAdjustedResults$observedReadCountRatios,
                                   stdCntMatrix = stdCntMatrix)



                return(allResults)


          }
)
