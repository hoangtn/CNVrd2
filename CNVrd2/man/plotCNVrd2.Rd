\name{plotCNVrd2}
\alias{plotCNVrd2}

\title{
Plot traces of samples.
}
\description{
Plot traces of samples.
}
\usage{
plotCNVrd2(Object, ...)
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
Plot
}

\author{

  Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}

}
\note{
Users can plot multiple samples simultaneously.
}

\examples{
data(fcgr3bMXL)
##Obtain all information of CNVs
allGroups <- copynumberGroups$allGroups
###Obtain names of duplicate samples
duplicatedSamples <- rownames(allGroups[allGroups[, 2] > 2,])
###Plot the first duplicate samples
par(mfrow = c(3, 2))
for (ii in duplicatedSamples[1:6])  
plotCNVrd2(Object = objectCNVrd2,
           segmentObject = resultSegment,
           sampleName = ii)

}
\keyword{methods}
\keyword{plots}
