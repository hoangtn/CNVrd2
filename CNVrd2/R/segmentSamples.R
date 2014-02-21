setGeneric("segmentSamples",
           function(Object, ...){standardGeneric("segmentSamples")})

setMethod("segmentSamples", "CNVrd2",
          function(Object, stdCntMatrix, entireGene = FALSE, inputBamFile = FALSE,
                   testThreshold2Merge = 0.25,
                   bThresholds = NULL, 
		   alpha = 0.01, nperm = 10000, p.method = "hybrid", 
		   min.width = 2, kmax = 25, nmin = 200, eta = 0.05, 
		   trim = 0.025, undo.splits = "none", 
	           undo.prune = 0.05, undo.SD = 3, verbose = 1){
##testThreshold2Merge: to merge two consecutive segments in the gene/locus being considered

##Checking parameters
  st = Object@st
  en = Object@en
  
  dirCoordinate = Object@dirCoordinate
  windows = Object@windows
  chr = Object@chr
  genes <- Object@genes
  if (is.null(bThresholds))
      bThresholds = 0.5*windows
  geneNames = Object@geneNames

   if (is.na(dirCoordinate)){
       dirCoordinate <- "TempAll"}
                  
       
                  
         if (substr(dirCoordinate, length(dirCoordinate), 1) != "/")
          dirCoordinate <- paste(dirCoordinate, "/", sep = "")

#####################################
  genes <- matrix(genes, nrow = 2)
  numberofWindows <- ceiling((en - st +1)/windows)

          
  stdCntMatrix <- stdCntMatrix
  if (is.null(testThreshold2Merge)){
    testThreshold2Merge <- 0}
  
  testThreshold2Merge <- abs(testThreshold2Merge)


  bamFile <- rownames(stdCntMatrix)
  nnn <- dim(stdCntMatrix)[1]

  #####Score of Genes########################
  sampleGenes <- matrix(0, nr = nnn, ncol = dim(genes)[2])
  rownames(sampleGenes) <- rownames(stdCntMatrix)

  ####################################################
  ######Testing read-count###########################
  testReadCount <- function(st1, en1, st, en, sampleTest){
    bamFile <- rownames(stdCntMatrix)
    nnn <- length(bamFile)
    fileTests <- c()
    for (ii in 1:nnn){
      bb1 <- scan(paste(dirCoordinate,  bamFile[ii], ".coordinate.txt", sep = ""))
      bbD <- bb1[(bb1 >= st1) & (bb1 <= en1)] ##Observed read counts
      lambda <- (en1 - st1 + 1)*length(bb1)/(en - st + 1)##Expected read counts
      fileTests[ii] <- length(bbD)/lambda ##Observed/Expected ratios
    }
    names(fileTests) <- bamFile
  ################################################################
  
fileTests1 <- (fileTests - median(fileTests))/sd(fileTests)
return(fileTests1[sampleTest])
}
#########################################################################
  ########Calculate observed read-count ratios for genes##################
  observedReadCountRatio <- matrix(0, ncol = ncol(genes), nrow = nnn)
  rownames(observedReadCountRatio) <- bamFile
  if (inputBamFile == TRUE){
    for (cc in 1:ncol(genes)){
      for (ii in 1:nnn){
        st1 <- genes[1, cc]
        en1 <- genes[2, cc]
        bb1 <- scan(paste(dirCoordinate,  bamFile[ii], ".coordinate.txt", sep = ""))
        bbD <- bb1[(bb1 >= st1) & (bb1 <= en1)] ##Observed read counts
        lambda <- (en1 - st1)*length(bb1)/(en - st)##Expected read counts
        observedReadCountRatio[ii, cc] <- length(bbD)/lambda ##Observed/Expected ratios
        }
        }}
    else
      observedReadCountRatio <- NA
      
 ########################Get Segmentation scores##########
  getSegScores <- function(x, geneLength){
    if (missing(geneLength))
      geneLength <- max(x[, 2]) - min(x[, 1]) + 1
    return(sum(x[, 3]*(x[, 2] - x[, 1] +1))/geneLength)}
  #######################TestSign################################3
  signTest <- function(x){
    signTest <- ifelse(x[1] > 0, 1, -1)
    for (i in 1:length(x)){
      if (x[i]*x[1] < 0){
        signTest <- 0
        break     }}
    return(signTest)
    }
   ##################################################
  

  ######################################################################
  #######Obtaining segmentation scores##################################
  #################Main function#######################################
    
  runFunction <- function(kk){
    # create 'CNA' Object for use with DNAcopy
    cna.obj0 <-DNAcopy::CNA(as(stdCntMatrix[kk,], "numeric"),
                 chrom=as(rep(1, numberofWindows), "numeric"),
                 seq(st, en, by = windows)[1:numberofWindows], data.type="logratio",
                 sampleid= paste(rownames(stdCntMatrix)[kk]))
    cna.obj <- DNAcopy::smooth.CNA(cna.obj0)
    
    # perform segmentation analysis & extract output
    cna.seg<-DNAcopy::segment(cna.obj, alpha = alpha, min.width = min.width,
                              nperm = nperm, p.method = p.method,
                              kmax = kmax, nmin = nmin, eta = eta,
                              trim = trim, undo.splits = undo.splits,
                              undo.prune = undo.prune, undo.SD = undo.SD, verbose = verbose)
    
    cna.out<-cna.seg$output
    cna.out
}

cnaListOut <- list()
for (kk in 1:nnn){

  cnaListOut[[kk]] <- runFunction(kk)

  if (entireGene == TRUE){
    cna.out <- cnaListOut[[kk]]
    for (jj in 1:ncol(genes)){

            
      cnaL <- dim(cna.out)[1]
      list1 <- data.frame(list(Start=rep(genes[1, jj], cnaL), End=rep(genes[2, jj], cnaL)))
      list2 <- data.frame(list(Start=cna.out[, 3], End=cna.out[, 4] + windows))
      cnaInter <- data.frame(Start = pmax(list1$Start, list2$Start),
                             End = pmin(list1$End, list2$End))
      cnaInter <- data.frame(cnaInter, seg.mean =cna.out[, 6])
      cnaInter[cnaInter$Start > cnaInter$End, ] <- NA
      cnaInter <- cnaInter[!is.na(cnaInter$Start),]
      if (dim(cnaInter)[1] > 1)

        segScores <- 0
      else
        segScores <- cnaInter[, 3]
      
       sampleGenes[kk, jj] <- segScores

    }
    }

  else{

     ############################
    ###Get scores of genes
    for (jj in 1:ncol(genes)){
      geneLength <- genes[2, jj] - genes[1, jj]

      #######################################################################################
      
      cnaL <- dim(cnaListOut[[kk]])[1]
      
      list1 <- data.frame(list(Start=rep(genes[1, jj], cnaL), End=rep(genes[2, jj], cnaL)))
      list2 <- data.frame(list(Start=cnaListOut[[kk]][, 3], End=cnaListOut[[kk]][, 4] + windows))
      
      cnaInter <- data.frame(Start = pmax(list1$Start, list2$Start),
                             End = pmin(list1$End, list2$End))
      cnaInter <- data.frame(cnaInter, seg.mean =cnaListOut[[kk]][, 6])
      cnaInter[cnaInter$Start > cnaInter$End, ] <- NA
      cnaInter <- cnaInter[!is.na(cnaInter$Start),]
      testSegment <- TRUE
      segScores <- 0
      dimRow <- dim(cnaInter)[1] ##Number of rows of results
######Sign of segmentation scores
      signOfSS <- ifelse(as.numeric(cnaInter[, 3]) > testThreshold2Merge, 1,
                         ifelse(as.numeric(cnaInter[, 3]) <  -testThreshold2Merge, -1, 0))
##Checking the same sign
      if (signTest(signOfSS) == 0)
          testSegment <- FALSE

      if (testSegment == TRUE)
        segScores <- getSegScores(cnaInter)
      else {
        testSegment <- TRUE
          if (dimRow == 2){
           if (((cnaInter[1, 2] - cnaInter[1, 1]) <= bThresholds) &
                ((cnaInter[2, 2] - cnaInter[2, 1]) > bThresholds)){
              testSup <- testReadCount(st1 = cnaInter[1, 1], en1 = cnaInter[1, 2],
                                       st = st, en = en, sampleTest = kk
                                       )
              
              segScores <- ifelse(testSup*signOfSS[2] > 0,
                                  getSegScores(cnaInter, geneLength = geneLength), 0)
            }
            else if (((cnaInter[1, 2] - cnaInter[1, 1]) > bThresholds) &
                ((cnaInter[2, 2] - cnaInter[2, 1]) <= bThresholds)){
              testSup <- testReadCount(st1 = cnaInter[2, 1], en1 = cnaInter[2, 2],
                                       st = st, en = en, sampleTest = kk)
              
              segScores <- ifelse(testSup*signOfSS[1] > 0,
                                  getSegScores(cnaInter, geneLength = geneLength), 0)
            }

          }
        #########################################################
        ######################## > 2 rows########################
        else  {
            if ((cnaInter[1, 2] - cnaInter[1, 1]) <= bThresholds)
            testSup1 <- testReadCount(st1 = cnaInter[1, 1], en1 = cnaInter[1, 2],
                                      st = st, en = en, sampleTest = kk)
            
            if ((cnaInter[dimRow, 2] - cnaInter[dimRow, 1]) <= bThresholds)
            testSup2 <- testReadCount(st1 = cnaInter[dimRow, 1], en1 = cnaInter[dimRow, 2],
                                      st = st, en = en, sampleTest = kk)

            
            if (signTest(signOfSS[1])*signTest(signOfSS[-1]) < 0){
              if ((cnaInter[1, 2] - cnaInter[1, 1]) <= bThresholds){
                segScores <- ifelse(testSup1*signTest(signOfSS[-1]) > 0,
                                    getSegScores(cnaInter[-1, ], geneLength = geneLength), 0)}
              
            }
            else if (signTest(signOfSS[length(signOfSS)])*signTest(signOfSS[-length(signOfSS)]) < 0){
              if ((cnaInter[dimRow, 2] - cnaInter[dimRow, 1]) <= bThresholds){
              segScores <- ifelse(testSup2*cnaInter[dimRow, 3] > 0,
                                  getSegScores(cnaInter[-dimRow, ], geneLength = geneLength), 0)}
              
            }
            else if (signTest(signOfSS[length(signOfSS)])*signTest(signOfSS[1]) > 0){
              if (signTest(signOfSS[-c(1, length(signOfSS))])*signTest(signOfSS[1]) == 0)
                segScores <- 0
              else if (((cnaInter[1, 2] - cnaInter[1, 1]) > bThresholds) |
                       ((cnaInter[dimRow, 2] - cnaInter[dimRow, 1]) > bThresholds))
                segScores <- 0
              else {

              segScores <- ifelse((testSup1*signTest(signOfSS[2:(length(signOfSS)-1)]) > 0) &
                                  (testSup2*signTest(signOfSS[2:(length(signOfSS)-1)]) > 0),
                                   getSegScores(cnaInter[2:(dimRow -1),], geneLength = geneLength), 0)
            }
            }
          }
             
       }
      sampleGenes[kk, jj] <- segScores
}    
}}
return(list(segmentationScores = sampleGenes, segmentResults = cnaListOut, observedReadCountRatios = observedReadCountRatio,
            stdCntMatrix = stdCntMatrix))
})
