setGeneric("plotCNVrd2",
           function(Object, ...){standardGeneric("plotCNVrd2")})

setMethod("plotCNVrd2", "CNVrd2",
          function(Object, segmentObject = NULL,
                   sampleName = NULL, xlim = NULL, ylim = NULL, 
			data1000Genomes = TRUE, geneColor = 'lightpink'){
##Checking parameters
  st = Object@st
  en = Object@en
  if (!is.null(xlim)){
    if ((xlim[1] < st) | (xlim[2] > en))
      stop("Please choose coordinates in the region")
    if ((xlim[1] > en) | (xlim[2] < st))
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
  

#####################################################3
###Plot#######################3
    plot(xCC, yCC, type='l',
         col = 'white',
         xlab= 'Coordinate',
         ylab= 'Standardized read count',
         main = paste(sampleName, "\nwindow = ", windows, sep = ""),
         ylim = ylim,
         xlim = c(outputST, outputEND)
         )

          for(k in 1:ncol(genes)){
            rect(genes[1, k], minGene, genes[2,k], maxGene ,col= geneColor)
            if (is.null(geneNames))
              text(genes[2, k], maxGene, paste("Gene", k, sep = ""), col = 'blue', cex = 0.7, lwd = 1.1, srt= 90, pos = 2)
            else
              text(genes[2, k], maxGene, geneNames[k], col = 'blue', cex = 0.7, lwd = 1.1, srt= 90, pos = 2)
            }

          lines(xCC, yCC, lwd =0.5, cex = 0.5)

          abline(h = median(yCC))
          if (nrow(cna.out) > 0) {
            for (ii in 1:nrow(cna.out))
              lines(cna.out[ii, 3:4], rep(cna.out[ii, 6], 2), col = 'green', cex = 2.5, lwd = 2.5)
            }

        }

)





