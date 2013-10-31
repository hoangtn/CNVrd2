#Use vectorORfactor to correct errors that can happen while reading a vcf file
#Link: https://stat.ethz.ch/pipermail/bioconductor/2013-January/050620.html
calculateLDSNPandCNV <- function(sampleCNV = NULL,
                                vcfFile = NULL, matrixGenotype = NULL,
                                cnvColumn = NULL,  popColumn = NULL,
                                population = NULL, chr = NULL, hg = "hg19", 
				st = NULL, en = NULL, nChunkForVcf = 10,
                                codeSNP = c("Two", "Three"),
                                codeCNV = c("CN", "ThreeGroup"),
                                typeTest = c("All", "Dup", "Del"),
                                parallel = FALSE){
###################################################################
##Test arguments###################################################
    if (is.null(sampleCNV))
      stop("Sample file is not empty")
    if (is.null(popColumn))
      stop("Population column  is not empty")
    if (is.null(cnvColumn))
        stop("Must input copy-number column")
    if (is.null(chr) | is.null(st) | is.null(en))
      stop("Chromosome, coordinates are not empty")
    codeSNP <- match.arg(codeSNP)
    codeCNV <- match.arg(codeCNV)
    typeTest <- match.arg(typeTest)

    sampleCNV <- sampleCNV
                                         
####################################################################
####################################################################
####################Main program####################################
    if (is.null(population))
      population <- as.character(sampleCNV[, popColumn][1])

    
#######################################################################
                        
########Read chunks of VCF files#########################################
    if (!is.null(matrixGenotype))
        snp.matrix <- matrixGenotype
    else {
        chunKs <- seq.int(st, en, length.out = nChunkForVcf)
        xMatrix <- cbind(round(chunKs[-length(chunKs)], 0), round(chunKs[-1], 0))
        xMatrix[-1, 1] <- xMatrix[-1, 1] + 1
        xMatrix <- split(xMatrix, row(xMatrix))
        tabixFile <- TabixFile(vcfFile)
        ###Function to read chunk
        readChunkVCF <- function(x){
            rangeVCF <- GRanges(seqnames = chr, ranges = IRanges(
                                                start = x[1],
                                                end = x[2]))
            message("VCF file: ", x[1], " to ", x[2])
            chunkVCF <- readVcf(tabixFile, hg, rangeVCF)
            return(geno(chunkVCF)$GT)
            }
        message("Reading the VCF file ",chr, ":", st, "-", en, " with ", nChunkForVcf, " blocks each")
        if (parallel == FALSE)
            snp.matrix <- do.call(rbind, lapply(xMatrix, readChunkVCF))
        else
            snp.matrix <- do.call(rbind, mclapply(xMatrix, readChunkVCF))
        }
    #############Remove duplicate rows################################
    snp.matrix <- snp.matrix[!duplicated(snp.matrix),]
    ############Match samples measured CNVs with samples in the VCF file

    calcPandR2 <- function(kk) {
      samples <- sampleCNV[sampleCNV[, popColumn] == kk,]
      if (codeCNV == "ThreeGroup"){
        medianCNV <- median(samples[, cnvColumn])
        samples[, cnvColumn] <- ifelse(samples[, cnvColumn] > medianCNV, 3,
                                   ifelse(samples[, cnvColumn] < medianCNV, 1, 2))
        if (typeTest == "Dup")
          samples <- samples[samples[, cnvColumn] >= 2,]
        if (typeTest == "Del")
          samples <- samples[samples[, cnvColumn] <= 2,]
        }
      if (dim(samples)[1] < 1)
        stop(paste("No samples of the population ", kk, "\n Please check the column cnvColumn", sep = ""))
      message("Calculating p and r2 values for ", kk, " population.")
      samples <- samples[order(samples[, cnvColumn]), ]
      colMatch <- pmatch(samples[, 1], colnames(snp.matrix))
      colMatch <- colMatch[!is.na(colMatch)]
    ############Keep only one population##############################
      snpMatrixPop <- snp.matrix[, colMatch]
      samples <- samples[pmatch(colnames(snpMatrixPop), samples[, 1]), ]
    #######################################################################
    ######################################################################
    ####Calculate p and r2 values########################################
      pv0 <- c()
      cov <- c()
    #####################################################################
      nameSamples <- colnames(snpMatrixPop)
      cnv.number <- samples[, cnvColumn]
      names(cnv.number) <- nameSamples
      snpMatrixPopT <- apply(snpMatrixPop, 2, function(x)
                             ifelse((x == "0|0") | (x == "0/0"), 0,
                                    ifelse((x == "1|1") | (x == "1/1"), 2, 1)))
      if (codeSNP == "Two")
        snpMatrixPopT <- apply(snpMatrixPopT, 2, function(x)
                               ifelse(x == 0, 0, 1))
      LDandP <- function(x){
        x1 <- ifelse(x == 0, 0, 1)
        tableCode <- table(x1, cnv.number)
        valueReturn <- c(NA, NA)
        if (nrow(tableCode) > 1){
          if (ncol(tableCode) > 1){
            result <- try(pv0 <- fisher.test(tableCode, workspace=2e+7,hybrid=TRUE)$p.value)
            if(class(result) == "try-error") {
              pv0 <- chisq.test(tableCode)$p.value
                    }

            cov0 <- cor(x, cnv.number, method = "spearman")
                }
          valueReturn <- c(pv0, cov0)
          }
        return(valueReturn)
        }
      resultPandR2 <- t(apply(snpMatrixPopT, 1, LDandP))
      rownames(resultPandR2) <- rownames(snpMatrixPopT)
    #########################################################
      resultPandR2 <- resultPandR2[!is.na(resultPandR2[, 1]),]
      pv0 <- resultPandR2[, 1]
      pv1 <- p.adjust(pv0, method = "BH")
      pv <- data.frame(resultPandR2, pv1, resultPandR2[, 2]^2)
      rownames(pv) <- rownames(resultPandR2)
      pv[, c(2, 4)] <- round(pv[, c(2, 4)], 2)
      pv[, c(1, 3)] <- format(signif(pv[, c(1, 3)], 2), scientific = TRUE)

      snpMatrixPopT <- snpMatrixPopT[pmatch(rownames(pv), rownames(snpMatrixPopT)), ]
      snpMatrixPopT <- snpMatrixPopT[rowSums(snpMatrixPopT) > 0 , ]

      cnvTable <- apply(snpMatrixPopT, 1, function(x){
        bTemp <- table(x, cnv.number)
        cnvFrequency <- bTemp[rownames(bTemp) != "0", ]/colSums(bTemp)
        if (is.matrix(cnvFrequency))
            cnvFrequency <- colSums(cnvFrequency)
        
        return(cnvFrequency)
        })

      
      cnvTable <- apply(cnvTable, 1, function(x) round(100*x, 2))
      pv <- cbind(cnvTable[pmatch(rownames(pv), rownames(cnvTable)),], pv, rep(kk, dim(pv)[1]))
      colnames(pv) <- c(paste(names(table(cnv.number)), "CN_(n=", table(cnv.number), ")", sep = ""),
                        "p.values", "r", "p.valuesAdjusted", "r2", "POP")
      pv <- pv[order(pv[, dim(pv)[2]-1] , decreasing = TRUE),]
    }
    if (length(population) == 1)
      PandR2list <- calcPandR2(population)
    else
      PandR2list <- lapply(population, calcPandR2)

    return(r2Andpvalues = PandR2list)
    }

