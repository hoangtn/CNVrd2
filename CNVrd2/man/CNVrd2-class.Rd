\name{CNVrd2-class}
\Rdversion{1.1}
\docType{class}
\alias{CNVrd2-class}
\alias{countReadInWindow,CNVrd2-method}
\alias{plotCNVrd2,CNVrd2-method}
\title{Class \code{"CNVrd2"}}
\description{
A class of reading BAM files into R and grouping read-count windows into similar segments.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNVrd2", ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{windows}:}{Object of class \code{"numeric"} indicating
      a window size. }
    \item{\code{chr}:}{Object of class \code{"character"} indicating the
    chromosome of the region.}
    \item{\code{st}:}{Object of class \code{"numeric"} indicating the
      starting coordinate of the region. }
    \item{\code{en}:}{Object of class \code{"numeric"} indicating the
      ending coordinate of the region. }
    \item{\code{dirBamFile}:}{Object of class \code{"character"}
      indicating a directory of BAM files. }
    \item{\code{dirCoordinate}:}{Object of class \code{"character"}
      indicating a directory where all the positions of mapped reads
      will be written out to prepare for the segmentation process. }
    \item{\code{genes}:}{Object of class \code{"numeric"} indicating
      gene coordinates. }
    \item{\code{geneNames}:}{Object of class \code{"character"}
      indicating names of genes. }
  }
}
\section{Methods}{
  \describe{
    \item{countReadInWindow}{\code{signature(Object = "CNVrd2")}: Count
      reads in windows. }
    \item{plotCNVrd2}{\code{signature(Object = "CNVrd2")}: Plot traces
      of samples. }
    \item{segmentSamples}{\code{signature(Object = "CNVrd2")}: Cluster
      windows of read counts into regions having similar signal values. }
	 }
}

\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\examples{
showClass("CNVrd2")
}
\keyword{classes}
