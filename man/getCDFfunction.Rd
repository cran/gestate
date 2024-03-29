% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_architecture.R
\name{getCDFfunction}
\alias{getCDFfunction}
\title{Method for returning the CDF function for a Curve object}
\usage{
getCDFfunction(theObject, ...)
}
\arguments{
\item{theObject}{The name of the Curve Object}

\item{...}{Pass-through arguments}
}
\description{
This retrieves the full CDF function of the Curve object as a string
}
\examples{
getCDFfunction(Weibull(100,1))
}
