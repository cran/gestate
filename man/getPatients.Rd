% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_architecture.R
\name{getPatients}
\alias{getPatients}
\title{Method for calculating expected number of recruited patients at a given time from an RCurve}
\usage{
getPatients(theObject, ...)
}
\arguments{
\item{theObject}{The name of the RCurve Object}

\item{...}{Pass-through arguments}
}
\description{
This calculates the expected number of recruited patients at a given time based upon an RCurve
}
\examples{
getPatients(LinearR(12,100,100),7)
}
