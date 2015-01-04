setMethod("plotCNVrd2", "CNVrd2",
          function(Object, segmentObject = NULL,
                   sampleName = NULL, xlim = NULL, ylim = NULL, 
			data1000Genomes = TRUE, geneColor = 'lightpink',
			xlab = NULL, ylab = NULL, main = NULL, lwd = 0.5,
                   cex = 1, segmentCEX = 3, segmentLWD = 3, segmentCOLOUR = 'green'){
##Checking parameters
  st = Object@st
  en = Object@en
  if (!is.null(xlim)){
    if (xlim[1] < st){
        xlim[1] <- st
        message("Left coordinate is less than the start position")
    }
    if (xlim[2] > en){
        xlim[2] <- en
        message("Right coordinate is larger than the end position")
    }
      
  outputST <- xlim[1]
  outputEND <- xlim[2]
}
else {
  outputST = st
  outputEND = en}

  if (is.null(xlab))
      xlab <- "Coordinate"
  if (is.null(ylab))
      ylab <- "Standardized read count"



   windows = Object@windows
  chr = Object@chr
  genes <- Object@genes
  geneNames = Object@geneNames
  
#####################################
  genes <- matrix(genes, nrow = 2)

  
  stdCntMatrix <- segmentObject$stdCntMatrix
  
  kk <- which(rownames(stdCntMatrix) == as.character(sampleName))
  xSample <- as.numeric(stdCntMatrix[kk, ])
  xSample <- ifelse(is.nan(xSample), 0, xSample)
  minGene <- as.numeric(min(xSample) + 0.3)
  maxGene <- as.numeric(max(xSample) - 0.3)
  
  if (is.null(ylim))
	ylim <- c(minGene, maxGene)
  nnn <- dim(stdCntMatrix)[1]
  cna.out <- segmentObject$segmentResults[[kk]]
  
  yCC <- xSample
  xCC <- seq(st, en, by = windows)[1:length(yCC)]

  if (data1000Genomes)
    sampleName <- substr(sampleName, 1, 7)
  

  if (is.null(main))
	main <- paste(sampleName, "\nwindow = ", windows, sep = "")
#####################################################3
###Plot#######################3
    plot(xCC, yCC, type='l',
         col = 'white',
         xlab= xlab,
         ylab= ylab,
         main = main,
         ylim = ylim,
         xlim = c(outputST, outputEND)
         )

          for(k in 1:ncol(genes)){
            rect(genes[1, k], ylim[1], genes[2, k], ylim[2] ,col= geneColor, border = NA)
            if (is.null(geneNames))
              text(genes[2, k], maxGene, paste("Gene", k, sep = ""), col = 'blue', cex = 0.7, lwd = 1.1, srt= 90, pos = 2)
            else
              text(genes[2, k], maxGene, geneNames[k], col = 'blue', cex = 0.7, lwd = 1.1, srt= 90, pos = 2)
            }

          lines(xCC, yCC, lwd = lwd, cex = cex)

          abline(h = median(yCC))
          if (nrow(cna.out) > 0) {
            for (ii in 1:nrow(cna.out))
              lines(cna.out[ii, 3:4], rep(cna.out[ii, 6], 2), col = segmentCOLOUR, cex = segmentCEX, lwd = segmentLWD)
            }

        }

)





