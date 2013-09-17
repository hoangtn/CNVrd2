\name{emnormalCNV}
\alias{emnormalCNV}
\title{
Implement the EM algorithm
}
\description{
This function is used to obtain the maximization likelihood estimation
of normal mixture model by using the EM algorithm (Demster et al., 1977). 
}
\usage{
emnormalCNV(Object, ...)
}

\arguments{
  \item{Object}{
An object of class clusteringCNVs.
}
  \item{\dots}{
Optional arguments
}
}
\note{
In the package, the distance between two initial means of the two nearest neighbor groups was set \code{groupDistance} = 0.25 as 
a default value to obtain initial values (using the \code{kmeans} function in R).
}

\value{
  \item{loglk}{Value of the likelihood function.}
  \item{p}{Proportions of groups.}
  \item{m}{Means of groups.}
  \item{sigma}{Standard deviations of groups.}
  \item{count}{A number of iteration to obtain convergence stage.}
  \item{bic}{See \code{\link{searchGroupCNVs}}.}
  \item{z}{Data frame of proportions of data in mixture components.}

}
\references{
Dempster, A. P., Laird, N. M., Rubin, D. B., 1977. Maximum likelihood from incomplete data via the em algorithm. Journal of the Royal Statistical Society. Series B (Methodological), 1-38.
}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\seealso{
\code{\link{searchGroupCNVs}}, \code{\link{groupCNVs}}

}
\examples{
data(fcgr3bMXL)

sS <- resultSegment$segmentationScores
#########Histogram###########################
###View segmentation scores##################
hist(sS[, 1], 100)
############################################
##Number of components#######################
###Make an object of clusteringCNVs class######
objectCluster <- new("clusteringCNVs",
                     x = sS[, 1], k = 4, EV = TRUE)

set.seed(123)
copynumberGroups <- groupCNVs(Object = objectCluster)
}
