\name{calculateLDSNPandCNV}
\alias{calculateLDSNPandCNV}

\title{
calculateLDSNPandCNV
}
\description{
Identifying SNPs/INDELs being in linkage disequilibrium with CNV.
}
\usage{
calculateLDSNPandCNV(sampleCNV = NULL, vcfFile = NULL,
	matrixGenotype = NULL, cnvColumn = NULL, popColumn = NULL,
	population = NULL, chr = NULL, hg = "hg19", 
	st = NULL, en = NULL,	nChunkForVcf = 10, 
	codeSNP = c("Two", "Three"), codeCNV = c("CN", "ThreeGroup"),
	typeTest = c("All", "Dup", "Del"), parallel = FALSE)

}

\arguments{
  \item{sampleCNV}{
A data frame with no missing data; The first column is samples, the
other columns are copy number, population names. This object can be obtained from the clustering step (allGroups).

}
\item{vcfFile}{
  Name of zipped vcf file including SNPs/INDELs.

}

  \item{matrixGenotype}{
A matrix of SNPs/INDELs coded: 0 for 0|0 or 0/0; 1 for 0/1, 1/0, 0|1,
  1|0; 2 for 1/1, 1|1. Rows are SNPs/INDELs and columns are samples. If
  users use this argument then the argument vcfFile is not necessary.
}
  \item{cnvColumn}{
A number indicating the column of CNV.
}
\item{popColumn}{
A number indicating the column of a population being calculated LD values.

}
  \item{population}{
A character() vector indicating names of populations.  
}
  \item{chr}{
A character string indicating a chromosome of genes/SNPs/INDELs.
}
  \item{hg}{
A character string indicating the version of a reference genome (default: hg19).

}
  \item{st}{
A number indicating a starting coordinate to read the vcf file.
}
\item{en}{
A number indicating an ending coordinate to read the vcf file.

}
\item{nChunkForVcf}{
A number indicating how many chunks users would like to divide the vcf file to read into R. It depends on users' computers. We usually use 10 to 50 for this argument. 
}
\item{codeSNP}{

A character string indicating a way to code unphased/phased SNPs/INDELs into numeric values (Two: 0, 1 or Three: 0, 1, 2).

}
\item{codeCNV}{
A character string partial matching to one of: \emph{\dQuote{All}} will test copy-number counts and SNPs/INDELs, \emph{\dQuote{ThreeGroup}}
  will divide samples into three groups: deletion, normality
  and duplication.
}
  \item{typeTest}{
A character string partial matching to one of: \dQuote{All} will test all copy-number status, \dQuote{Dup/Del} will test for
only duplicated/deleted versus normal.
}
\item{parallel}{
A logical value indicating whether multicores are used or NOT (default: FALSE).

}
}

\value{
  \code{r2Andpvalues}: 
  A data frame (or a list) with each row for a SNP/INDEL: all information includes p-values adjusted by the method of Benjamini and Hocheberg (1995) and r2 values
  between the SNP/INDEL and copy-number status. 

}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\seealso{
\code{\link{fcgr3bMXL}}
}
\note{
\emph{st} and \emph{en} must not be outside the coordinates of the VCF file.
}
\references{
Benjamini, Y., and Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society. Series B (Methodological), 57, 289-300.

}
\examples{
##Load data: fcgr3bMXL in CNVrd2 package############
data(fcgr3bMXL)
##Name a vcf file (vcfFile)
vcfFile <- system.file(package="CNVrd2", "extdata",
                      "chr1.161600000.161611000.vcf.gz")
##Make a data fame named sampleCNV including samples, CNs, population names

sampleCNV <- data.frame(copynumberGroups$allGroups[, c(1,2) ],rep("MXL", 58))

rownames(sampleCNV) <- substr(sampleCNV[, 1], 1, 7)
sampleCNV[, 1] <- rownames(sampleCNV)
##The first column must be the sample names 
tagSNPandINDELofMXL <- calculateLDSNPandCNV(sampleCNV = sampleCNV,
                                 vcfFile = vcfFile, cnvColumn = 2,
                                 population = "MXL", popColumn = 3,
                                 nChunkForVcf = 5, chr = "1",
                                 st = 161600000, en = 161611000,
                                 codeSNP= "Three", codeCNV = "ThreeGroup")
tagSNPandINDELofMXL[1:3,]




}
