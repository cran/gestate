% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curve_constructors.R
\name{PieceExponential}
\alias{PieceExponential}
\title{Piecewise Exponential Curve constructor function}
\usage{
PieceExponential(start, lambda)
}
\arguments{
\item{start}{Vector of start times for each period. First element must be 0. Must be same length as lambda vector.}

\item{lambda}{Vector of rate parameters from the corresponding respective time defined in start vector until the start of the next period. Must be same length as start vector.}
}
\description{
This creates a Curve object for a Piecewise Exponential distribution.\cr
Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
Parameterisation follows that used by pexp etc. This function requires a vector of start times (beginning with 0) and a corresponding vector of rates. See Details for more information on parameterisation.
}
\details{
The piecewise exponential distribution with rates lambda_1 to lambda_n and start times t_1 to t_n has parameterisation:\cr
Product(x=1:length(lambda)) of (e^(-lambda[x].t[x])) where t[x] = min(start[x+1],max(0,t-start[x])). start[x+1] is defined as Inf if otherwise undefined.
}
\examples{
PieceExponential(start=c(0,6,24),lambda=c(0.05,0.01,0.001))
}
\author{
James Bell
}
