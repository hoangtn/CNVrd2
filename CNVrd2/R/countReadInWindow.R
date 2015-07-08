#This method counts number of reads in constant windows
setMethod("countReadInWindow", "CNVrd2",
          function(Object, correctGC = FALSE, standardizingAllSamples = TRUE,
                   qualityThreshold = 0,
                   rawReadCount = FALSE, byGCcontent = 1, useRSamtoolsToCount = FALSE,
			writeCoordinate = TRUE,
                   referenceGenome = "BSgenome.Hsapiens.UCSC.hg19",reference_fasta=NULL){

              if (correctGC){
                  if(is.null(reference_fasta)){
                  library(referenceGenome, character.only=TRUE)
                    variable_name <- strsplit(referenceGenome,'.',fixed=T)[[1]][2]
                    do.call("<-",list(referenceGenome,get(variable_name)))
                  } else {
                    referenceGenome = readDNAStringSet(reference_fasta,format="fasta")
                }
             }
              
              windows = Object@windows
              chr = Object@chr
              st = Object@st
              en = Object@en
              dirBamFile = Object@dirBamFile
              dirCoordinate <- Object@dirCoordinate

              if (is.na(dirCoordinate)){
                  dirCoordinate <- "TempAll"
                  dir.create(dirCoordinate)
              }
              if (is.na(dirBamFile))
                  dirBamFile <- "./"
              if (substr(dirBamFile, length(dirBamFile), 1) != "/")
                  dirBamFile <- paste(dirBamFile, "/", sep = "")
              if (substr(dirCoordinate, length(dirCoordinate), 1) != "/")
                  dirCoordinate <- paste(dirCoordinate, "/", sep = "")

              bamFile <- dir(path = dirBamFile, pattern = ".bam$")

              ##########################################################
              #######Function to divide reads into windows##############
              getWindows <- function(data, windows, st){
                  data <- data - st + 1
                  data <- data[data > 0]


                  return(table(ceiling(data/windows)))
                  }
              what <- c("pos", "mapq")
              param <- Rsamtools::ScanBamParam( what = what)
              numberofWindows <- ceiling((en - st + 1)/windows)
              seqStart <- seq(st, en, by = windows)[1:numberofWindows]

              ###Function to read Bam files and write out coordinates###############
              countReadForBamFile <- function(x){
                    bam <- Rsamtools::scanBam(paste(dirBamFile, bamFile[x], sep = ""),  param=param)[[1]]
                    allPos <- cbind(bam[[1]], bam[[2]])
                    allPos <- allPos[as.numeric(allPos[, 2]) >= qualityThreshold,]
                    bam <- allPos[, 1]
                    
                    bam <- bam[!is.na(bam)]
                    
                    bam <- bam[(bam >= st) & (bam <= en)]
                    if (writeCoordinate)
                        write.table(bam, paste(dirCoordinate, bamFile[x], ".coordinate.txt", sep = ""),
                                col.names = FALSE, quote = FALSE, row.names = FALSE)
                    aa <- getWindows(data = bam, windows = windows, st = st)
                    if (length(aa) > numberofWindows)
                        aa <- aa[1:numberofWindows]
                    names(aa) <- as.integer(names(aa))
                    
                    tempRow <- rep(0, numberofWindows)
                    names(tempRow) <- as.integer(c(1:numberofWindows))
                    tempRow[names(tempRow) %in% names(aa)] <- aa
                    message("Reading file: ", bamFile[x])
                    return(tempRow)
                    }
              ####Read all Bam files#########################################
              if (useRSamtoolsToCount == TRUE){

                  fileName <- paste(dirBamFile, bamFile, sep = "")
                  what <- c("pos")
                  which <- IRanges::RangesList('2' = IRanges(seq(objectCNVrd2@st, objectCNVrd2@en, by = objectCNVrd2@windows),
                            seq(objectCNVrd2@st, objectCNVrd2@en, by = objectCNVrd2@windows) + objectCNVrd2@windows))

                  names(which) <- as.character(as.name(gsub("chr", "", objectCNVrd2@chr)))
                  param <- ScanBamParam( what = what, which = which)
                  aa1 <- lapply(fileName, function(x) {
                      message("Reading: ", x)
                      return(countBam(x, param = param)$records)
                         })
                  readCountMatrix <-   do.call(rbind, aa1)

              } else {
                  readCountMatrix <- do.call(rbind, lapply(1:length(bamFile), countReadForBamFile))
              }
              rownames(readCountMatrix) <- bamFile
              message("=============================================")
              message(dim(readCountMatrix)[1], " bam files were read")
              message("=============================================")
########################################################################################
########Correct GC content###############################################################

              if (correctGC){
                  gcContent <- function(){
                      message("Correcting the GC content")
                      chr <- as.character(chr)
                      if(is.null(reference_fasta)){
                      tempG <- unmasked(Hsapiens[[chr]])[(st):en]} else{
                          names(referenceGenome) <- toupper(names(referenceGenome))

                          message("names(referenceGenome) ")
                          print(names(referenceGenome))
                          
                          positionChr <- grep(toupper(chr), names(referenceGenome))[1]

                          message("position: ")
                          print(positionChr)
                          tempG <- referenceGenome[positionChr]
                          names(tempG) <- toupper(chr)
                          tempG <- tempG[[toupper(chr)]][st:en]
    }
                      sT <- seq(1, length(tempG), by = windows)

                      endPos <- sT + windows - 1
                      ####Change endPos if endPos is over end
                      if (endPos[length(endPos)] > length(tempG))
                          endPos[length(endPos)] <- length(tempG)
                      
                      gcSegment <- Views(tempG, sT, endPos)
                      gcContentInSegment <- apply(letterFrequency(gcSegment, letters = c("G", "C")), 1, sum)

                      gcContentInSegment <- ifelse(is.na(gcContentInSegment), 0, gcContentInSegment)

                      gcContentInSegment <- gcContentInSegment/windows

                      return(gcContentInSegment)
                      }
################################################################################
                       gcn <- gcContent()
                  
                  gcn <- gcn[1:numberofWindows]


                  gcList <- list()


                       readCountMatrix <- readCountMatrix
                       nnn <- dim(readCountMatrix)[1]
                       readCountMatrix <- as.matrix(readCountMatrix)
                       cnt1 <- readCountMatrix
###Normalization GC content by median
                  gcContentInSegment <- gcn*100

                  windowByGCcontent <- seq(0, 100, by = byGCcontent)

                  if (windowByGCcontent[length(windowByGCcontent)] <= 100)
                      windowByGCcontent[length(windowByGCcontent)] <- 101
                  dataframeToCorrect <- data.frame(windowByGCcontent[-length(windowByGCcontent)],
                                 windowByGCcontent[-1])

                  gcList <- apply(dataframeToCorrect, 1, function(x){
                      tempMap <- which((gcContentInSegment >= x[1]) & (gcContentInSegment < x[2]))
                      return(tempMap)
                  })

                  gcList <- unique(gcList)

                  gcList <- gcList[lapply(gcList,length)>0]
                  
                  lengthGC <- length(gcList)
                  
                  correctGCforRow <- function(xRow){
                      medianAll <- median(xRow)
                      for (jj in 1:lengthGC){
                          x = gcList[[jj]]
                          if (!is.null(x)){
                              x1 = xRow[x]
                              medianRegion <- median(x1)
                              if (medianRegion != 0){
                                  xRow[x] <- x1*medianAll/medianRegion
                              }
                              else
                                  xRow[x] <- xRow[x]
                          }    }
                      return(xRow)
                  }

                  readCountMatrix <- t(apply(cnt1, 1, correctGCforRow))
              }
 ################################################################################
 ###################Transfer to the same coverage
              if (rawReadCount == FALSE){
                       readCountMatrix <- t(apply(readCountMatrix, 1, function(x) x <- x/median(x)))
 ###############################################################################
#####################Standardize across samples
                       if (standardizingAllSamples == TRUE){
                         readCountMatrix <- apply(readCountMatrix, 2, function(x) ifelse(is.na(x), 0, x))
                         
                         readCountMatrix <- apply(readCountMatrix, 2, function(x) ifelse(is.infinite(x), 1, x))
                         readCountMatrix <- apply(readCountMatrix, 2, function(x) (x -median(x))/sd(x))
                                                  readCountMatrix <- apply(readCountMatrix, 2, function(x) ifelse(is.nan(x), 0, x))
                       }
                 }
return(readCountMatrix)
}     )
