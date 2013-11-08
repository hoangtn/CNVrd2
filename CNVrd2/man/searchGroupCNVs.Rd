\name{searchGroupCNVs}
\alias{searchGroupCNVs}

\title{
Choose a number of CN groups 
}
\description{
Choose a number of CN groups by using Bayesian Information Criterion (Schwarz, 1978).

}
\usage{
searchGroupCNVs(Object, ...)
}
\arguments{
  \item{Object}{
An object of class \emph{clusteringCNVs}.
}
  \item{\dots}{
Optional arguments.
  }
}

\value{
\item{Groups}{General information.}
\item{nComponents}{A suitable number of components.}

}
\references{
Schwarz, G. , 1978. Estimating the dimension of a model. The Annals of Statistics 6(2), 461-464.
}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}



\keyword{methods}

