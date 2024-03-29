% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_architecture.R
\name{getParam}
\alias{getParam}
\title{Method for returning a single parameter from a Curve object}
\usage{
getParam(theObject, ...)
}
\arguments{
\item{theObject}{The name of the Curve Object}

\item{...}{Pass-through arguments}
}
\description{
This retrieves a single parameter from a Curve object
}
\examples{
getParam(Weibull(100,1),1)
}
