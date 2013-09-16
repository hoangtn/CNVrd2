\name{plotPolymorphicRegion}
\alias{plotPolymorphicRegion}
\title{
Plot polymorphic regions.
}
\description{
Plot polymorphic regions based on coordinates set by users.
}
\usage{
plotPolymorphicRegion(Object, ...)
}

\arguments{
  \item{Object}{
    An object of class CNVrd2
}
\item{\dots}{
  \code{polymorphicRegionObject}{An object obtained from the process of
  identifying poplymorphic regions.}

  Optional arguments.
}
}

\value{
  \item{putativeBoundary}{Putative boundaries of polymorphic regions
    based on quantile values}

}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\note{
Users can choose various quantile values and adjust different thresholds
to obtain polymorphic regions.
}


\seealso{
 \code{\link{identifyPolymorphicRegion}}
}
\keyword{plot}
