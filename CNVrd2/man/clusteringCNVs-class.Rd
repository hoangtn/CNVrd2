\name{clusteringCNVs-class}
\Rdversion{1.1}
\docType{class}
\alias{clusteringCNVs-class}
\alias{objectCluster}
\alias{emnormalCNV,clusteringCNVs-method}
\alias{groupCNVs,clusteringCNVs-method}
\alias{searchGroupCNVs,clusteringCNVs-method}
\title{Class \code{"clusteringCNVs"}}
\description{
This class is used to cluster segmentation scores into copy-number groups.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("clusteringCNVs", ...)}.
 
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}.}
    \item{\code{k}:}{Object of class \code{"numeric"} indicating a number of groups. }
    \item{\code{p}:}{Object of class \code{"numericOrNULL"} indicating
      groups' proportions.}
    \item{\code{m}:}{Object of class \code{"numericOrNULL"} indicating groups' means. }
    \item{\code{sigma}:}{Object of class \code{"numericOrNULL"} indicating groups' standard deviations.}
    \item{\code{small}:}{Object of class \code{"numeric"} indicating the
    value to stop the iteration process of the EM algorithm.}
    \item{\code{nMax}:}{Object of class \code{"numeric"} indicating a
       maximum number of iterations. }
    \item{\code{EV}:}{Object of class \code{"logical"} indicating
      whether all groups having equal variances or not (default). }
    \item{\code{eee}:}{Object of class \code{"numeric"} indicating a pseudo value of 0.}
    \item{\code{nmaxInit}:}{Object of class \code{"numeric"} indicating
      a number of iterations to obtain initial values. }
    \item{\code{nChangeVariance}:}{Object of class \code{"numeric"}
      indicating a number of times to change from unequal variances to
      equal variances (\emph{\dQuote{this option is used to avoid the EM algorithm being broken down
      if there is one (or a few) sample in a group}}).}
    \item{\code{verbose}:}{Object of class \code{"logical"} indicating
      whether printing out all loops. }
    \item{\code{groupDistance}:}{Object of class \code{"numericOrNULL"} indicating the distance between groups.}
  }
}
\section{Methods}{
  \describe{
    \item{emnormalCNV}{\code{signature(Object = "clusteringCNVs")}: run
      the EM algorithm. }
    \item{groupCNVs}{\code{signature(Object = "clusteringCNVs")}:
      cluster segmentation scores into groups. }
    \item{searchGroupCNVs}{\code{signature(Object = "clusteringCNVs")}:
      identify a number of groups.}
	 }
}

\author{
  Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyen@gmail.com}

}
\examples{
showClass("clusteringCNVs")
}
\keyword{classes}
