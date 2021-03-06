% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs2.R
\docType{data}
\name{bearings}
\alias{bearings}
\title{Fatigue lifetime data}
\format{
A vector containing 10 observations.
}
\usage{
data(bearings)
}
\description{
Several data sets related to the BS distribution are available in the RBS package, which have been taken from the literature on this topic.
}
\details{
McCool (1974) gives data related to the fatigue life in hours of 10 bearings of a certain type.
}
\examples{
## Load data sets
data(bearings)
## Histogram for bearings
hist(bearings)

}
\references{
McCool, J. I. (1974). Inferential techniques for Weibull populations. Aerospace Research Laboratories Report ARL TR74-0180, Wright-Patterson Air Force Base, Dayton, OH.
}
\keyword{datasets}
