
##################################################################################################################################
# Curve constructors - modify to add new curve shapes/types
#
##################################################################################################################################

#################################################################################################
# Curve constructors
#  These are the methods used to create particular types of curve
#  If you have a new curve type to add, create a new constructor in this section.
#  'type' contains the name of the type of curve being specified
#  The PDF entry must be the name of a function containing the PDF of the curve
#  Likewise, the CDF entry has the name of the CDF function
#  The parameter names go in the vector 'pnames'
#  The parameters go into the list 'pvalues'. Note that the order should correspond to that in pnames.
#  'paramno' should be the length of the pvalue list.
#
#################################################################################################


# Note that the Weibull curve object uses the parameterisation found in the default R functions rweibull/qweibull/pweibull
# use '?rweibull' for further information.
#' Weibull Curve constructor function
#'
#' This creates a Curve object for a Weibull distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by pweibull etc. See Details for more information on parameterisation.
#' @param alpha Scale parameter for Weibull distribution.
#' @param beta Shape parameter for Weibull distribution. Default is 1; an exponential distribution.
#' @details The Weibull distribution with shape parameter beta and scale parameter alpha has parameterisation:\cr
#' f(x) = (beta/alpha) (x/alpha)^(beta-1) exp(- (x/alpha)^beta)\cr
#' F(x) = 1 - exp(- (x/alpha)^beta)
#' @author James Bell
#' @examples
#' Weibull(alpha=100,beta=0.8)
#' @export
Weibull <- function(alpha,beta=1){
  if(length(beta)!=1){stop("beta parameter must be a single value")}
  if(length(alpha)!=1){stop("alpha parameter must be a single value")}
  new("Curve",type="Weibull",PDF="dweibull",CDF="pweibull",RF="rweibull",inverse="qweibull",paramno=2,pnames=c("scale","shape"),pvalue=list(alpha,beta))
}

#' Log-normal Curve constructor function
#'
#' This creates a Curve object for a Log-normal distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by plnorm etc. See Details for more information on parameterisation.
#' @param mu Mean (on log scale) parameter for Log-normal distribution.
#' @param sigma Standard deviation (on log scale) parameter for Log-normal distribution.
#' @details The log normal distribution has parameterisation:\cr
#' f(x) = 1/(sqrt(2*pi) sigma  x) e^-((log x - mu)^2 / (2 sigma^2))\cr
#' F(x) = 0.5(1 + erf((log(x)-mu)/(sigma sqrt(2))))\cr
#' where erf is the error function.
#' @author James Bell
#' @examples
#' Lognormal(mu=5,sigma=1.2)
#' Lognormal(6)
#' @export
Lognormal <- function(mu,sigma=1){
  if(length(mu)!=1){stop("mu parameter must be a single value")}
  if(length(sigma)!=1){stop("sigma parameter must be a single value")}
  new("Curve",type="Lognormal",PDF="dlnorm",CDF="plnorm",RF="rlnorm",inverse="qlnorm",paramno=2,pnames=c("meanlog","sdlog"),pvalue=list(mu,sigma))
}

#' Exponential Curve constructor function
#'
#' This creates a Curve object for a Exponential distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by pexp etc. Note that lambda = 1/alpha from the Weibull constructor. See Details for more information on parameterisation.
#' @param lambda Rate parameter for Exponential distribution.
#' @details The exponential distribution has parameterisation:\cr
#' f(x) = lambda e^(- lambda x)\cr
#' F(x) = 1 - e^(- lambda x)
#' @author James Bell
#' @examples
#' Exponential(0.01)
#' @export
Exponential <- function(lambda){
  if(length(lambda)!=1){stop("lambda parameter must be a single value")}
  new("Curve",type="Exponential",PDF="dexp",CDF="pexp",RF="rexp",inverse="qexp",paramno=1,pnames="rate",pvalue=list(lambda))
}

#' Blank Curve constructor function
#'
#' This creates a Curve object for a 'Blank' pseudo-distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' This distribution is 0 by definition for all times. It is not therefore a true probability distribution.
#' @details The blank pseudo-distribution is used for impossible events, notably where censoring is not possible/allowed.\cr
#' f(x) = 0\cr
#' F(x) = 0
#' @author James Bell
#' @examples
#' Blank()
#' @export
Blank <- function(){
  new("Curve",type="Blank",PDF="pmin",CDF="pmin",RF="INF",inverse="INF",paramno=1,pnames="Zero",pvalue=list(0))
}

#' Piecewise Exponential Curve constructor function
#'
#' This creates a Curve object for a Piecewise Exponential distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by pexp etc. This function requires a vector of start times (beginning with 0) and a corresponding vector of rates. See Details for more information on parameterisation.
#' @param start Vector of start times for each period. First element must be 0. Must be same length as lambda vector.
#' @param lambda Vector of rate parameters from the corresponding respective time defined in start vector until the start of the next period. Must be same length as start vector.
#' @details The piecewise exponential distribution with rates lambda_1 to lambda_n and start times t_1 to t_n has parameterisation:\cr
#' Product(x=1:length(lambda)) of (e^(-lambda[x].t[x])) where t[x] = min(start[x+1],max(0,t-start[x])). start[x+1] is defined as Inf if otherwise undefined.
#' @author James Bell
#' @examples PieceExponential(start=c(0,6,24),lambda=c(0.05,0.01,0.001))
#' @export
PieceExponential <- function(start,lambda){
  if(length(start)!=length(lambda)){stop("Piecewise exponential curve has mismatched length 'start' and 'lambda' vectors")}
  if(start[1]!=0){stop("First element of piecewise exponential curve must start at 0")}
  if(is.unsorted(start,strictly=TRUE)){stop("Start times must be in ascending order with no duplicates")}
  new("Curve",type="PieceExponential",PDF="dpieceexp",CDF="ppieceexp",RF="rpieceexp",inverse="qpieceexp",paramno=2,pnames=c("start","rate"),pvalue=list(start,lambda))
}

#' Mixture Exponential Curve constructor function
#'
#' This creates a Curve object for a Mixture Exponential distribution, commonly used for modelling distributions with subpopulations.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by pexp etc. See Details for more information on parameterisation.
#' @param props Vector of length x for the probabilities of the subpopulations. Must sum to 1.
#' @param lambdas Vector of length x for the rate parameters for the corresponding subpopulations define by props.
#' @details The mixture distribution with rates lambda1 to lambda2 etc and prevalence p1 and p2 etc has parameterisation:\cr
#' f(x) = p1 lambda1 e^(- lambda1 x) + p2 lambda2 e^(- lambda2 x)+...\cr
#' F(x) = p1 (1 - e^(- lambda1 x)) + p2 (1 - e^(- lambda2 x))+...
#' @author James Bell
#' @examples MixExp(props=c(0.8,0.2),lambdas=c(0.01,0.1))
#' @export
MixExp <- function(props,lambdas){
  if(length(props)!=length(lambdas)){stop("Mixture exponential curve has mismatched length 'props' and 'lambdas' vectors")}
  if(sum(props)!=1){stop("Proportions must sum to 1!")}
  new("Curve",type="MixExp",PDF="dmixexp",CDF="pmixexp",RF="rmixexp",inverse="qmixexp",paramno=2,pnames=c("props","lambdas"),pvalue=list(props,lambdas))
}

#' Mixture Weibull Curve constructor function
#'
#' This creates a Curve object for a Mixture Weibull distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' Parameterisation follows that used by pweibull etc. See Details for more information on parameterisation.
#' @param props Vector of length x for the probabilities of the two subpopulations. Must sum to 1.
#' @param alphas Vector of length x for the scale parameters for the corresponding subpopulations define by props.
#' @param betas Vector of length x for the shape parameters for the corresponding subpopulations define by props. Default is rep(1,length(props)), i.e. all exponential distributions.
#' @details The mixture distribution with scales alpha1 and alpha2 etc, shapes beta1 and beta2 etc, and prevalences p1 and p2 etc has parameterisation:\cr
#' f(x) = p1 (beta1/alpha1) (x/alpha1)^(beta1-1) exp(- (x/alpha1)^beta1) + p2 (beta2/alpha2) (x/alpha2)^(beta2-1) exp(- (x/alpha2)^beta2)+...\cr
#' F(x) = p1 (1 - exp(- (x/alpha1)^beta1) + p2 (1 - exp(- (x/alpha2)^beta2)+...
#' @author James Bell
#' @examples MixWei(props=c(0.8,0.2),alphas=c(100,10),betas=c(1.1,0.9))
#' @export
MixWei <- function(props,alphas,betas=rep(1,length(props))){
  if(length(props)!=length(alphas)){stop("Mixture weibull curve has mismatched length 'props' and 'alphas' vectors")}
  if(length(props)!=length(betas)){stop("Mixture weibull curve has mismatched length 'props' and 'betas' vectors")}
  if(sum(props)!=1){stop("Proportions must sum to 1!")}
  new("Curve",type="MixWei",PDF="dmixwei",CDF="pmixwei",RF="rmixwei",inverse="qmixwei",paramno=3,pnames=c("props","betas","alphas"),pvalue=list(props,betas,alphas))
}

#' Log-logistic Curve constructor function
#'
#' This creates a Curve object for a Log-logistic distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' See Details for information on parameterisation.
#' @param theta Scale parameter for Log-logistic distribution.
#' @param eta Shape parameter for Log-logistic distribution.
#' @details The log-logistic distribution has parameterisation:\cr
#' f(x) = eta (theta^beta) x^(eta-1) (theta^eta + x^eta)^-2\cr
#' F(x) = (x^eta) /(theta^eta+x^eta)
#' @author Jasmin Ruehl
#' @examples LogLogistic(theta=20,eta=2)
#' @export
LogLogistic <- function(theta, eta){
  new('Curve', type='LogLogistic', PDF='dloglog', CDF='ploglog', RF='rloglog',inverse='qloglog',
      paramno=2, pnames=c('scale', 'shape'), pvalue=list(theta, eta))
}

#' Gompertz Curve constructor function
#'
#' This creates a Curve object for a Gompertz distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' See Details for information on parameterisation.
#' @param theta Scale parameter for Log-logistic distribution.
#' @param eta Shape parameter for Log-logistic distribution.
#' @details The Gompertz distribution has parameterisation:\cr
#' f(x) = theta eta e^(eta + theta x - eta e^(theta x))\cr
#' F(x) = 1 - exp(eta - eta e^(theta x))
#' @author Jasmin Ruehl
#' @examples Gompertz(theta=0.02,eta=2)
#' @export
Gompertz <- function(theta, eta){
  new('Curve', type='Gompertz', PDF='dgompertz', CDF='pgompertz', RF='rgompertz',inverse='qgompertz',
      paramno=2, pnames=c('scale', 'shape'), pvalue=list(theta, eta))
}

#' Generalised Gamma Curve constructor function
#'
#' This creates a Curve object for a Generalised Gamma distribution.\cr
#' Curve objects contain all necessary information to describe a distribution, including functions and parameters describing it.\cr
#' See Details for information on parameterisation.
#' @param theta Scale parameter for Generalised Gamma distribution.
#' @param eta Shape parameter for Generalised Gamma distribution.
#' @param rho Family parameter for Generalised Gamma distribution.
#' @details The Generalised Gamma distribution has parameterisation:\cr
#' f(x) = (rho x^((rho eta)-1) e^(-(x/theta)^rho) theta^(-rho eta) )/Gamma(eta)\cr
#' F(x) = LPGamma(eta,(x/theta)^rho)/Gamma(eta)\cr
#' where Gamma is the gamma function, and LPGamma is the lower partial gamma function.\cr
#' As of v1.4.0, all values of eta are now fully supported.\cr
#' @author Jasmin Ruehl
#' @references Tadikamalla PR, Random Sampling from the Generalized Gamma Distribution. Computing, 1979, 23(2), 199-203.
#' @examples GGamma(theta=20,eta=2,rho=0.7)
#' @export
GGamma <- function(theta, eta, rho){
  new('Curve', type='GGamma', PDF='dggamma', CDF='pggamma', RF='rggamma',inverse="qggamma",
      paramno=3, pnames=c('scale', 'shape', 'family'), pvalue=list(theta, eta, rho))
}

###################################################################################################
# Recruitment Curve constructors
#  These are the methods used to create particular types of recruitment curve
#  If you have a new recruitment curve type to add, create a new constructor in this section.
#  See the text for the Curve constructors.
#  There are additional slots for N, Nactive, Ncontrol and Ratio, referring to the overall recruitment numbers/ratios
#  There are several possible ways to specify these parameters; it would be possible to write different constructors that create the same
#     recruitment curve types but using different inputs. Currently only one parameterisation method is provided per curve type.
#
####################################################################################################
#' LinearR RCurve constructor function
#'
#' This creates a RCurve object for a linear recruitment distribution.\cr
#' RCurve objects contain all necessary information to describe a recruitment distribution. They are a particular type of Curve object containing additional recruitment-related information, including patient numbers and the randomisation ratio.\cr
#' @param rlength Length of recruitment.
#' @param Nactive Number of patients recruited in the active arm.
#' @param Ncontrol Number of patients recruited in the control arm.
#' @details This RCurve is used when it is expected that patients enter a trial at a constant rate until the required number is achieved.
#' @author James Bell
#' @examples LinearR(rlength=12,Nactive=100,Ncontrol=100)
#' @export
LinearR <- function(rlength,Nactive,Ncontrol){
  new("RCurve",type="LinearR",PDF="linear_recruitPDF",CDF="linear_recruit",RF="linear_sim",inverse="NULL",paramno=1,pnames="rlength",pvalue=list(rlength),N=Nactive+Ncontrol,Nactive=Nactive,Ncontrol=Ncontrol,Ratio=Nactive/Ncontrol,Length=rlength,maxF=Inf)
}

#' InstantR RCurve constructor function
#'
#' This creates a RCurve object for an instant recruitment distribution.\cr
#' RCurve objects contain all necessary information to describe a recruitment distribution. They are a particular type of Curve object containing additional recruitment-related information, including patient numbers and the randomisation ratio.\cr
#' @param Nactive Number of patients recruited in the active arm.
#' @param Ncontrol Number of patients recruited in the control arm.
#' @details This RCurve is used when either all patients enter at the same time, or a fixed-length follow-up design is used. Note that a PDF function is not provided for this RCurve type, but is not required for standard calculations.
#' @author James Bell
#' @examples InstantR(Nactive=100,Ncontrol=100)
#' @export
InstantR <- function(Nactive,Ncontrol){
  new("RCurve",type="InstantR",PDF="NULL",CDF="instant_recruit",RF="instant_sim",inverse="NULL",paramno=1,pnames="Dummy",pvalue=list(0),N=Nactive+Ncontrol,Nactive=Nactive,Ncontrol=Ncontrol,Ratio=Nactive/Ncontrol,Length=0,maxF=Inf)
}

#' PieceR RCurve constructor function
#'
#' This creates a RCurve object for a piecewise-linear recruitment distribution.\cr
#' RCurve objects contain all necessary information to describe a recruitment distribution. They are a particular type of Curve object containing additional recruitment-related information, including patient numbers and the randomisation ratio.\cr
#' @param recruitment 2-column matrix with recruitment parameters. First column gives the lengths of each period of recruitment. Second column gives the corresponding rates of recruitment for each period.
#' @param ratio Randomisation ratio; active arm divided by control arm.
#' @details This RCurve is used when it is expected that patients enter a trial at a rate that varies over time.
#' @author James Bell
#' @examples
#' rmatrix <- matrix(c(rep(4,3),5,10,15),ncol=2)
#' rmatrix
#' PieceR(rmatrix,1)
#' @export
PieceR <- function(recruitment,ratio){
  lengths <- recruitment[,1]
  rates <- recruitment[,2]
  N <- sum(rates*lengths)
  Nactive <- N*(ratio/(ratio+1))
  Ncontrol <- N-Nactive
  new("RCurve",type="PieceR",PDF="piece_recruitPDF",CDF="piece_recruit",RF="piece_sim",inverse="NULL",paramno=2,pnames=c("lengths","rates"),pvalue=list(lengths,rates),N=N,Ratio=ratio,Nactive=Nactive,Ncontrol=Ncontrol,Length=sum(lengths),maxF=Inf)
}

#' PieceR RCurve constructor function
#'
#' This creates a RCurve object for a piecewise-linear recruitment distribution with a fixed (maximum) per-patient follow-up time.\cr
#' RCurve objects contain all necessary information to describe a recruitment distribution. They are a particular type of Curve object containing additional recruitment-related information, including patient numbers and the randomisation ratio.\cr
#' @param recruitment 2-column matrix with recruitment parameters. First column gives the lengths of each period of recruitment. Second column gives the corresponding rates of recruitment for each period.
#' @param ratio Randomisation ratio; active arm divided by control arm.
#' @param maxF Fixed follow-up time per patient, i.e. maximum time a patient will be at risk independent of length of study. 
#' @details This RCurve is used when it is expected that patients enter a trial at a rate that varies over time and there is a fixed maximum follow-up time per patient.
#' @author James Bell
#' @examples
#' rmatrix <- matrix(c(rep(4,3),5,10,15),ncol=2)
#' rmatrix
#' PieceRMaxF(recruitment=rmatrix,ratio=1,maxF=12)
#' @export
PieceRMaxF <- function(recruitment,ratio,maxF){
  lengths <- recruitment[,1]
  rates <- recruitment[,2]
  N <- sum(rates*lengths)
  Nactive <- N*(ratio/(ratio+1))
  Ncontrol <- N-Nactive
  new("RCurve",type="PieceR",PDF="piece_recruitPDFMaxF",CDF="piece_recruitMaxF",RF="piece_simMaxF",inverse="NULL",paramno=3,pnames=c("lengths","rates","maxF"),pvalue=list(lengths,rates,maxF),N=N,Ratio=ratio,Nactive=Nactive,Ncontrol=Ncontrol,Length=sum(lengths),maxF=maxF)
}

