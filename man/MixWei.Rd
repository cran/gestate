% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/curve_constructors.R
\name{MixWei}
\alias{MixWei}
\title{Mixture Weibull Curve constructor function}
\usage{
MixWei(props, alphas, betas = rep(1, length(props)))
}
\arguments{
\item{props}{Vector of length x for the probabilities of the two subpopulations. Must sum to 1.}

\item{alphas}{Vector of length x for the scale parameters for the corresponding subpopulations define by props.}

\item{betas}{Vector of length x for the shape parameters for the corresponding subpopulations define by props. Default is rep(1,length(props)), i.e. all exponential distributions.}
}
\description{
This creates a Curve object for a Mixture Weibull distribution.\cr
Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
Parameterisation follows that used by pweibull etc. See Details for more information on parameterisation.
}
\details{
The mixture distribution with scales alpha1 and alpha2 etc, shapes beta1 and beta2 etc, and prevalences p1 and p2 etc has parameterisation:\cr
f(x) = p1 (beta1/alpha1) (x/alpha1)^(beta1-1) exp(- (x/alpha1)^beta1) + p2 (beta2/alpha2) (x/alpha2)^(beta2-1) exp(- (x/alpha2)^beta2)+...\cr
F(x) = p1 (1 - exp(- (x/alpha1)^beta1) + p2 (1 - exp(- (x/alpha2)^beta2)+...
}
\examples{
MixWei(props=c(0.8,0.2),alphas=c(100,10),betas=c(1.1,0.9))
}
\author{
James Bell
}
