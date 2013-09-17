\name{countReadInWindow}
\alias{countReadInWindow}

\title{
Obtain read counts in constant windows.
}
\description{
Counting, transfering and standardizing read counts for all windows of
samples. If correctGC = TRUE then all read-count windows will be
corrected by the method of Yoon et al. (2009). 

}
\usage{
countReadInWindow(Object, ...)
}

\arguments{
  \item{Object}{
An object of class CNVrd2.
}
  \item{\dots}{
Further aguments.
}
}
\value{
\code{readCountMatrix}: a matrix of read counts for all samples (rows).

}
\references{
Yoon, S., Xuan, Z., Makarov, V., Ye, K., Sebat, J., 2009. Sensitive and accurate
detection of copy number variants using read depth of coverage. Genome research 19 (9), 1586-1592.
}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}

}

\examples{
data(fcgr3bMXL)
#bamFiles <- dir("Bam", pattern = ".bam$")
objectCNVrd2 <- new("CNVrd2", windows = 1000, chr = "chr1",
                   st = 161100001, en = 162100000,
                   dirBamFile = "Bam",
                   genes = c(161592986, 161601753),
                   geneNames = "3B")

##readCountMatrix <- countReadInWindow(Object = objectCNVrd2, correctGC = TRUE)
##readCountMatrix[1:3, 1:3]
}
\keyword{methods}

