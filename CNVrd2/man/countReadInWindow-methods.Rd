\name{countReadInWindow-methods}
\docType{methods}
\alias{countReadInWindow-methods}
\title{ Method \code{countReadInWindow} }

\description{
Method to count reads in windows  
}
\seealso{\code{\link{countReadInWindow}}}
\section{Methods}{
\describe{

\item{\code{signature(Object = "CNVrd2")}}{
Count, stransfer and standardize read count in windows of samples
}
}}
\examples{
##data(fcgr3bMXL)
##readCountMatrix <- countReadInWindow(Object = objectCNVrd2, correctGC = TRUE)
##readCountMatrix[1:3, 1:3]

}
\keyword{methods}
