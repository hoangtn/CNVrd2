\name{identifyPolymorphicRegion}
\alias{identifyPolymorphicRegion}
\title{
Identity polymorphic regions.
}
\description{
Using quantile values to identify polymorphic regions.
}
\usage{
identifyPolymorphicRegion(Object, ...)
}
\arguments{
  \item{Object}{An object of class CNVrd2}
  \item{\dots}{
    \code{polymorphicRegionObject}{An object obtained from the process of
  identifying poplymorphic regions.}

  Optional arguments.
 
}

}

\value{
  \item{putativeBoundary}{resultant boundaries based on quantile
    thresholds.}
  \item{subRegionMatrix}{segmentation-score matrix of sub-regions derived from the
    segmentation process.}
  \item{subRegion}{sub-regions derived from the segmentation process.}
  \item{mQuantile}{matrix of quantile values.}
  \item{mSD}{data frame including subregions and their standard deviations.}
  \item{Vst}{data frame of Vst between populations if VstTest = TRUE.}
  \item{SSofPolymorphicRegions}{segmentation scores of polymorphic regions.}

}

\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\note{
Users can choose various quantile values and adjust different thresholds
to obtain polymorphic regions.
}

\seealso{
\code{\link{plotPolymorphicRegion}}
}
\keyword{methods}
