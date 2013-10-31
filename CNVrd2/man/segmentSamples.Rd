\name{segmentSamples}
\alias{segmentSamples}

\title{
Implement the segmentation process
}
\description{
Segment read-count windows into region having similar signal values by
using the DNAcopy package (Venkatraman and Olshen, 2007) and refine this process 
to obtain segmentation scores at genes.


}
\usage{
segmentSamples(Object, ...)
}

\arguments{
  \item{Object}{
An object of class CNVrd2.
}
  \item{\dots}{
Optional arguments.
}
}

\value{

  \item{segmentResults}{All results of the segmentation process.}
  \item{segmentationScores}{Segmentation scores of the gene(s) being measured.}
  \item{observedReadCountRatios}{Observed read-count ratios of genes. This value is a matrix of observed read-count ratios at genes if (only inputBamFile = TRUE).}
  \item{stdCntMatrix}{Matrix of read counts (standardized).}

}
\references{
Venkatraman, E., Olshen, A. B., 2007. A faster circular binary segmentation algorithm for the analysis of array chg data. Bioinformatics 23 (6), 657-663.}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\seealso{
\code{\link{countReadInWindow}}, \code{DNAcopy}
}
\examples{
data(fcgr3bMXL)
\dontrun{resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)}
}
\keyword{methods}
