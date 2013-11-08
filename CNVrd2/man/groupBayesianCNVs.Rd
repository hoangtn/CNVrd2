\name{groupBayesianCNVs}
\alias{groupBayesianCNVs}

\title{
groupBayesianCNVs
}
\description{
  Cluster segmentation scores into different groups by using prior
  information from one population.

}
\usage{
groupBayesianCNVs(xData, nGroups, lambda0, sd0, alpha0, distanceBetweenGroups, inits = NULL, 
		precisionOfGroupMeans = 3000, sdOftau = NULL, n.adapt = 100, 
		nUpdate = 1000, n.iter = 20000, thin = 5, n.chains = 1, 
		heidel.diag = FALSE, leftLimit = NULL, rightLimit = NULL)
}

\arguments{
  \item{xData}{
    a numeric vector of observations (segmentation scores).
}
\item{nGroups}{
  an integer indicating a number of groups.

}
\item{lambda0}{
  Prior means of groups. 

}
\item{sd0}{
  Prior standard deviations of groups.

}
\item{alpha0}{
  Prior parameters for mixing proportions.

}
\item{distanceBetweenGroups}{
  Prior value for the distance between groups.

}
\item{inits}{
  A list of initial values of parameters.
}
\item{precisionOfGroupMeans}{
  Prior parameter of group means (default = 3000).

}
\item{sdOftau}{
  Prior parameter of the standard deviations of group precisions.

}
\item{n.adapt}{
  the number of iterations for adaptation (\code{rjags}'s parameter)

}
\item{nUpdate}{
  the number of iterations for burn-in process.
%%     ~~Describe \code{nUpdate} here~~
}
\item{n.iter}{
  the number of iterations for sampling (\code{rjags}'s parameter)

}
\item{thin}{
  thinning interval for monitors (\code{rjags}'s parameter)

}
\item{n.chains}{
  the number of parallel chains for the model (\code{rjags}'s parameter,
  default = 1)

}
\item{heidel.diag}{
If heidel.diag = TRUE then Heidelberger and Welch's convergence diagnostic is used.
}
\item{leftLimit}{
Values which are less than this value will be allocated to the smallest group.
}
\item{rightLimit}{
Values which are larger than this value will be allocated to the largest group.
}

}
\details{
  This function assumes that users already know the information of
  groups' means, standard deviations; the distances between groups.

  
}
\value{
  \item{mcmcChains}{A list of mcarray objects for means, standard
    deviations, proportions}
  \item{m1}{Means of groups}
  \item{s1}{Standard deviations of groups}
  \item{p1}{Proportions of groups}
  \item{allGroups}{A data.frame includes samples and their corresponding
    groups}
  \item{hTest}{Results of Heidelberger and Welch's convergence diagnostic}
}
\references{
  Martyn Plummer (2013). rjags: Bayesian graphical models using MCMC. R
  package version 3-10. \url{http://CRAN.R-project.org/package=rjags}.

  Lunn, David J., et al. WinBUGS-a Bayesian modelling framework:
  concepts, structure, and extensibility. Statistics and Computing 10.4
  (2000): 325-337.

}
\author{
Hoang Tan Nguyen, Tony R Merriman and MA Black. \email{hoangtannguyenvn@gmail.com}
}
\note{
#####
}


\seealso{
  \link{groupCNVs}

}
\examples{
\dontrun{
data(ccl3l1data)

xyEuro <- ccl3l1data[grep("CEU|TSI|IBS|GBR|FIN", ccl3l1data[, 2]), ]

names(yEuro) <- rownames(xyEuro)

##Clustering European segmentation scores into group: 5 groups were chosen

objectClusterEuroCCL3L1 <- new("clusteringCNVs", x = yEuro, k = 5)

europeanCCL3L1Groups <- groupCNVs(Object = objectClusterEuroCCL3L1)

##Obtain prior information
#Means
lambda0 <- as.numeric(europeanCCL3L1Groups$m)
#SD
sdEM <- as.numeric(europeanCCL3L1Groups$sigma)
#Proportions
pEM <- as.numeric(europeanCCL3L1Groups$p)


###Calculate the distances between groups
for (ii in 2:5){print(lambda0[ii] - lambda0[ii-1])}

###All segmentation scores
ccl3l1X <- ccl3l1data$SS
names(ccl3l1X) <- as.character(ccl3l1data$Name)
range(ccl3l1X)


 
##Set prior information: 
#prior for the sd of the means of groups: 
#5 was set for the third group = 2 CN
sd <- c(1, 1, 5, 1, 1) 
ccl3l1X <- sort(ccl3l1X)
###Data
xData <- ccl3l1X
###Number of groups
nGroups <- 10 
###prior for means of groups
lambda0 <- lambda0 
###Prior for mixing proportions
alpha0 <-  c(3, 29, 44, 18, 7,  5, rep(2, nGroups -length(pEM) -1))
##Prior for the distances between groups
distanceBetweenGroups = 0.485

sdEM = sdEM


##Adjust standard deviation for the fifth group
sdEM[5] <- sdEM[4]
 
 set.seed(123)
 groupCCL3L1allPops <- groupBayesianCNVs(xData = xData, nGroups = nGroups,
                                         lambda0 = lambda0,
                                         sd0 = sdEM, alpha0 = alpha0,
                                         distanceBetweenGroups = distanceBetweenGroups,
                                         sdOftau = sd,
                                        rightLimit = 4)



}
}
\keyword{methods}

