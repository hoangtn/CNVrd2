#This method counts number of reads in constant windows
setMethod("countReadInWindow", "CNVrd2",
          function(Object, correctGC = FALSE, standardizingAllSamples = TRUE,
                   rawReadCount = FALSE, byGCcontent = 5,
                   referenceGenome = "BSgenome.Hsapiens.UCSC.hg19"){

              if (correctGC){
                  library(referenceGenome, character.only=TRUE)
                  }
              
              windows = Object@windows
              chr = Object@chr
              st = Object@st
              en = Object@en
              dirBamFile = Object@dirBamFile
              dirCoordinate <- Object@dirCoordinate
              
              if (is.na(dirCoordinate)){
                  dir.create("TempAll")
                  dirCoordinate <- "TempAll"}
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
              what <- c("pos")
              param <- ScanBamParam( what = what)
              numberofWindows <- ceiling((en - st + 1)/windows)
              seqStart <- seq(st, en, by = windows)[1:numberofWindows]

              ###Function to read Bam files and write out coordinates###############
              countReadForBamFile <- function(x){
                    bam <- scanBam(paste(dirBamFile, bamFile[x], sep = ""),  param=param)[[1]]$pos
                    bam <- bam[!is.na(bam)]

                    
                    bam <- bam[(bam >= st) & (bam <= en)]
                    write.table(bam, paste(dirCoordinate, bamFile[x], ".coordinate.txt", sep = ""),
                                col.names = FALSE, quote = FALSE, row.names = FALSE)
                    aa <- getWindows(data = bam, windows = windows, st = st)
                    if (length(aa) > numberofWindows)
                        aa <- aa[1:numberofWindows]
                    names(aa) <- as.integer(names(aa))
                    
                    tempRow <- rep(0, numberofWindows)
                    names(tempRow) <- as.integer(c(1:numberofWindows))
                    tempRow[names(tempRow) %in% names(aa)] <- aa
                    cat("Reading file: ", bamFile[x], "\n")
                    return(tempRow)
                    }
              ####Read all Bam files#########################################
              readCountMatrix <- do.call(rbind, lapply(1:length(bamFile), countReadForBamFile))
              rownames(readCountMatrix) <- bamFile
              cat("=============================================\n")
              cat(dim(readCountMatrix)[1], " bam files were read",  "\n")
              cat("=============================================\n")
########################################################################################
########Correct GC content###############################################################

              if (correctGC){
                  gcContent <- function(){
                      cat("Correcting the GC content\n")
                      chr <- as.character(chr)
                      tempG <- unmasked(Hsapiens[[chr]])[(st):en]
                      gc <- c()
                      temp <- seq(1, length(tempG), by = windows)
                      for (ii in 1:length(temp)){
                          if (temp[ii] < (length(tempG) - windows))
                                 gc[ii] <- sum(alphabetFrequency(tempG[temp[ii]:(temp[ii+1] - 1)], baseOnly= TRUE)[2:3])/windows
                              else
                              gc[ii] <- sum(alphabetFrequency(tempG[temp[ii]:length(tempG)], baseOnly= TRUE)[2:3])/windows
                          }
                      gc <- ifelse(is.na(gc), 0, gc)
                      return(gc)
                      }
################################################################################
                       gcn <- 100*gcContent()
                  gcn <- gcn[1:numberofWindows]
                                          

                       readCountMatrix <- readCountMatrix
                       nnn <- dim(readCountMatrix)[1]
                       readCountMatrix <- as.matrix(readCountMatrix)
                       cnt1 <- readCountMatrix
###Normalization GC content by median
                  gcSeq <- seq(0, 100, by = byGCcontent)
                  gcList <- list()
                  for (ii in 2:length(gcSeq)){
                      tempGC <- gcn[(gcn >= gcSeq[ii - 1]) & (gcn < gcSeq[ii])]
                      if (length(tempGC) > 0){
                          gcList[[ii]] <- pmatch(tempGC, gcn)
                          }   }
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
