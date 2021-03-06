% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s4_architecture.R
\name{plotSF}
\alias{plotSF}
\title{Method for plotting the Survival Function of a Curve object}
\usage{
plotSF(theObject, ...)
}
\arguments{
\item{theObject}{The name of the RCurve Object}

\item{...}{Pass-through arguments}
}
\description{
This plots a Curve Survival Function
}
\examples{
plotSF(Weibull(100,1))
plotSF(Weibull(100,1),xlab="Test x label",maxT=60)
plotSF(Weibull(80,0.8),overlay=TRUE,col=2,lty=2)
}
