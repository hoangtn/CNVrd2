\name{groupCNVs}
\alias{groupCNVs}

\title{
Group segmentation scores into groups.
}
\description{
  Use the EM algorithm (Dempster et al., 1977) to group segmentation
  scores into various groups.

}
\usage{
groupCNVs(Object, ...)
}

\arguments{
  \item{Object}{
An object of class clusteringCNVs.
}
 
  \item{\dots}{
Further arguments.
}
}
\details{
Users can set limits of segmentation scores: values being smaller than the
  left limit will be assigned to the smallest group and values being larger
  than righ limit will be assigned to the largest group.
}
\value{

  \item{allGroups}{Samples and their corresponding groups}
  \item{means}{Means of groups.}
  \item{sigma}{Variances of groups.}
  \item{p}{Proportions of groups in all data set.}
  \item{loglk}{Value of loglikehood function.}
}
\references{
Dempster, A. P., Laird, N. M., Rubin, D. B., 1977. Maximum likelihood from incomplete data via the em algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 1-38.
}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}

\seealso{
\code{\link{emnormalCNV}}, \code{\link{searchGroupCNVs}}
}
\examples{
data("fcgr3bMXL")
#resultSegment <- segmentSamples(Object = objectCNVrd2, stdCntMatrix = readCountMatrix)
objectCluster <- new("clusteringCNVs",
                     x = resultSegment$segmentationScores[, 1], k = 4, EV = TRUE)

#searchGroupCNVs(Object = objectCluster)
copynumberGroups <- groupCNVs(Object = objectCluster)

}
\keyword{methods}
