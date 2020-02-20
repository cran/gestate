##################################################################################################################################
# Functions for defining PDFs, CDFs, and RFs
# 	These are required for setting up non-trivial Curve/RCurve objects
##################################################################################################################################

########################################################################################################
# Function for the censoring CDF of a linear recruitment
# Note that this is in trial time
# Takes as input:
#  q: 	time
#  assess:	time of assessment
#  rlength: length of recruitment period
########################################################################################################
linear_recruit <- function(q,assess,rlength){
  punif(q=q,min=assess-rlength,max=assess)
}

########################################################################################################
# Function for the censoring PDF of a linear recruitment
# Note that this is in trial time
# Takes as input:
#  x: 	time
#  assess:	time of assessment
#  rlength: length of recruitment period
########################################################################################################
linear_recruitPDF <- function(x,assess,rlength){
  dunif(x=x,min=assess-rlength,max=assess)
}

########################################################################################################
# Function for the RF of a linear recruitment
# Note that this is in patient-time
# Takes as input:
#  n: 	how many values to create
#  rlength: length of recruitment period
########################################################################################################
linear_sim <- function(n,rlength){
  runif(n=n,min=0,max=rlength)
}

########################################################################################################
# Function for the RF of a Blank function (used for censoring)
# Takes as input:
#  n: 	how many values to create
#  Zero:    dummy parameter
########################################################################################################
INF <- function(n,Zero){
  rep(Inf,n)
}

########################################################################################################
# Function for the censoring CDF of a piecewise linear recruitment
#  Note that the vectors for lengths and rates should be the same length
#  This can be assured by using an RCurve object as the source, loaded using a recruitment matrix
#
# Takes as input:
#  q: 	time
#  assess:	time of assessment
#  lengths: vector of lengths of recruitment period
#  rates:	vector of recruitment rates of recruitment periods
########################################################################################################
piece_recruit <- function(q,assess,lengths,rates){
  rows <- length(lengths)
  rates <- c(0,rates[rows:1],0)
  lengths <- c(assess-sum(lengths),c(lengths[rows:1]),max(q))
  cumlengths <- cumsum(lengths)
  cumpats <- c(0,cumsum(rates*lengths))
  location <- findInterval(q,cumlengths)+1
  cumlengths <- c(0,cumlengths)
  out <- (cumpats[location]+rates[location]*(q-cumlengths[location]))/sum(rates*lengths)
  return(out)
}

########################################################################################################
# Function for the censoring PDF of a piecewise linear recruitment
#  Note that the vectors for lengths and rates should be the same length
#  This can be assured by using an RCurve object as the source, loaded using a recruitment matrix
#
# Takes as input:
#  x: 	time
#  assess:	time of assessment
#  lengths: vector of lengths of recruitment period
#  rates:	vector of recruitment rates of recruitment periods
########################################################################################################
piece_recruitPDF <- function(x,assess,lengths,rates){
  rows <- length(lengths)
  rates <- c(0,rates[rows:1],0)
  lengths <- c(assess-sum(lengths),c(lengths[rows:1]),max(x))
  cumlengths <- cumsum(lengths)
  location <- findInterval(x,cumlengths)+1
  out <- rates[location]/sum(rates*lengths)
  return(out)
}

########################################################################################################
# Function for the RF of a piecewise linear recruitment
# Takes as input:
#  n: 	number of patients to simulate
#  lengths: vector of lengths of recruitment rate periods
#  rates:   vector of recruitment rates per period
########################################################################################################
piece_sim <- function(n,lengths,rates){
  lengths <- c(0,lengths)
  rates <- c(0,rates)
  numbers <- lengths*rates
  total <- sum(numbers)
  rand <- runif(n=n,min=0,max=total)
  cumtime <- cumsum(lengths)
  cumnumbers <- cumsum(numbers)
  locations <- findInterval(rand,cumnumbers)
  output <- cumtime[locations]+(rand-cumnumbers[locations])/rates[locations+1]
  return(output)
}

########################################################################################################
# Function for the censoring CDF of a piecewise linear recruitment with maximum follow-up time
#  Note that the vectors for lengths and rates should be the same length
#  This can be assured by using an RCurve object as the source, loaded using a recruitment matrix
#
# Takes as input:
#  q: 	time
#  assess:	time of assessment
#  lengths: vector of lengths of recruitment period
#  rates:	vector of recruitment rates of recruitment periods
########################################################################################################
piece_recruitMaxF <- function(q,assess,lengths,rates,maxF){
  rows <- length(lengths)
  rates <- c(0,rates[rows:1],0)
  lengths <- c(assess-sum(lengths),c(lengths[rows:1]),max(q))
  cumlengths <- cumsum(lengths)
  cumpats <- c(0,cumsum(rates*lengths))
  location <- findInterval(q,cumlengths)+1
  cumlengths <- c(0,cumlengths)
  out <- (cumpats[location]+rates[location]*(q-cumlengths[location]))/sum(rates*lengths)
  out[q >= maxF] <- 1
  return(out)
}

########################################################################################################
# Function for the censoring PDF of a piecewise linear recruitment with maximum follow-up time
#  Note that the vectors for lengths and rates should be the same length
#  This can be assured by using an RCurve object as the source, loaded using a recruitment matrix
#
# Takes as input:
#  x: 	time
#  assess:	time of assessment
#  lengths: vector of lengths of recruitment period
#  rates:	vector of recruitment rates of recruitment periods
########################################################################################################
piece_recruitPDFMaxF <- function(x,assess,lengths,rates,maxF){
  rows <- length(lengths)
  rates <- c(0,rates[rows:1],0)
  lengths <- c(assess-sum(lengths),c(lengths[rows:1]),max(x))
  cumlengths <- cumsum(lengths)
  location <- findInterval(x,cumlengths)+1
  out <- rates[location]/sum(rates*lengths)
  maximum <- findInterval(maxF-0.1,cumlengths)+1
  cumpats <- c(0,cumsum(rates*lengths))
  patsatmax <- cumpats[maximum] + rates[maximum]*((maxF-0.1)-cumlengths[maximum-1])
  out[x >= (maxF-0.1) & x <= maxF] <- 10*(sum(rates*lengths)-patsatmax)/sum(rates*lengths)
  out[x > maxF] <- 0
  return(out)
}

########################################################################################################
# Function for the RF of piecewise linear recruitment with maximum follow-up time
# Takes as input:
#  n: 	number of patients to simulate
#  lengths: vector of lengths of recruitment rate periods
#  rates:   vector of recruitment rates per period
########################################################################################################
piece_simMaxF <- function(n,lengths,rates,maxF){
  lengths <- c(0,lengths)
  rates <- c(0,rates)
  numbers <- lengths*rates
  total <- sum(numbers)
  rand <- runif(n=n,min=0,max=total)
  cumtime <- cumsum(lengths)
  cumnumbers <- cumsum(numbers)
  locations <- findInterval(rand,cumnumbers)
  output <- cumtime[locations]+(rand-cumnumbers[locations])/rates[locations+1]
  output <- pmin(output,maxF)
  return(output)
}

########################################################################################################
# Function for the RF of instant recruitment
# Takes as input:
#  n: 	number of patients to simulate
#  ...:     placeholder for dummy arguments
########################################################################################################
instant_sim <- function(n, ...){
  rep(0,n)
}

########################################################################################################
# Function for the censoring CDF of instant recruitment
# Takes as input:
#  q: 	time
#  assess:	time of assessment
#  ...:     placeholder for dummy arguments
########################################################################################################
instant_recruit <- function(q,assess,...){
  rep(0,length(q))
}

########################################################################################################
# Function for the censoring PDF of instant recruitment
# Note that due to the difficulties with numerically-integrating impulse functions,
# it is assumed that the PDF is that of a 0.1 month linear recruitment period.
# This will lead to a small amount of error being introduced into calculations using this function
# Takes as input:
#  x: 	time
#  assess:	time of assessment
#  ...:     placeholder for dummy arguments
########################################################################################################
instant_recruitPDF <- function(x,assess,...){
  linear_recruitPDF(x,assess,0.1)
}

########################################################################################################
# Function for the CDF of a piecewise exponential distribution
# Takes as input:
#  q: 	time
#  start:   vector of lengths of pieces (first element must be 0, must be strictly ascending)
#  rate:	vector of rates of pieces (length must match length of start)
########################################################################################################
ppieceexp <- function(q,start,rate){
  lengths <- length(start)
  transitions <- cumprod(c(1,exp(-rate[1:(lengths-1)]*diff(start))))
  places <- findInterval(q,start)
  output <- 1-transitions[places]*exp(-rate[places]*(q-start[places]))
  return(output)
}

########################################################################################################
# Function for the PDF of a piecewise exponential distribution
# Takes as input:
#  x: 	time
#  start:   vector of lengths of pieces (first element must be 0, must be strictly ascending)
#  rate:	vector of rates of pieces (length must match length of start)
########################################################################################################
dpieceexp <- function(x,start,rate){
  lengths <- length(start)
  transitions <- cumprod(c(1,exp(-rate[1:(lengths-1)]*diff(start))))
  places <- findInterval(x,start)
  output <- rate[places]*transitions[places]*exp(-rate[places]*(x-start[places]))
  return(output)
}

########################################################################################################
# Function for the RF of a piecewise exponential distribution
# Takes as input:
#  n: 	number of random draws to make
#  start:   vector of lengths of pieces (first element must be 0, must be strictly ascending)
#  rate:	vector of rates of pieces (length must match length of start)
########################################################################################################
rpieceexp <- function(n,start,rate){
  return(qpieceexp(runif(n),start,rate))
}

########################################################################################################
# Function for the inverse-CDF of a piecewise exponential distribution
# Takes as input:
#  p: 	probability
#  start:   vector of lengths of pieces (first element must be 0, must be strictly ascending)
#  rate:	vector of rates of pieces (length must match length of start)
########################################################################################################
qpieceexp <- function(p,start,rate){
  lengths <- length(start)
  transitions <- cumprod(c(1,exp(-rate[1:(lengths-1)]*diff(start))))
  places <- (lengths+1) - findInterval(1-p,c(0,rev(transitions)))
  output <- start[places]+(log(transitions[places]) - log(1-p))/rate[places]
  return(output)
}

########################################################################################################
# Function for the CDF of a mixed exponential distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  q: 	time
#  props:   vector of proportions (Must sum to 1. NB: these are not checked!)
#  lambdas:	vector of rates (length must match length props)
########################################################################################################
pmixexp <- function(q,props,lambdas){
  width <- length(props)
  leng <- length(q)
  qmat <- matrix(rep(q,width),ncol=width)
  propsmat <- matrix(rep(props,each=leng),ncol=width)
  lambdasmat <- matrix(rep(lambdas,each=leng),ncol=width)
  rowSums(propsmat*pexp(q=qmat,rate=lambdasmat))
}

########################################################################################################
# Function for the PDF of a mixed exponential distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  x: 	time
#  props:   vector of proportions (Must sum to 1. NB: these are not checked!)
#  lambdas:	vector of rates (length must match length props)
########################################################################################################
dmixexp <- function(x,props,lambdas){
  width <- length(props)
  leng <- length(x)
  xmat <- matrix(rep(x,width),ncol=width)
  propsmat <- matrix(rep(props,each=leng),ncol=width)
  lambdasmat <- matrix(rep(lambdas,each=leng),ncol=width)
  rowSums(propsmat*dexp(x=xmat,rate=lambdasmat))
}

########################################################################################################
# Function for the RF of a mixed exponential distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  n: 	number of random draws to make
#  props:   vector of proportions (must sum to 1. NB: these are not checked!)
#  lambdas:	vector of rates (length must match length props)
########################################################################################################
rmixexp <- function(n,props,lambdas){
  assignment <- sample.int(n=length(props),size=n,replace=TRUE,prob=props)
  output <- rexp(n=n,rate=lambdas[assignment])
  return(output)
}

########################################################################################################
# Function for the inverse-CDF of a mixed exponential distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  p: 	probability
#  props:   vector of proportions (must sum to 1. NB: these are not checked!)
#  lambdas:	vector of rates (length must match length props)
########################################################################################################
qmixexp <- function(p,props,lambdas){
  q <- rep(NA,length(p))
  rootfunction <- function(t,p,props,lambdas){
    sum(props*exp(-lambdas*t))+p-1
  }
  for (i in 1:length(p)){
    limits <- qexp(p[i],lambdas)
    q[i] <- uniroot(rootfunction,lower=min(limits),upper=max(limits),p=p[i],props=props,lambdas=lambdas)$root
  }
  return(q)
}

########################################################################################################
# Function for the CDF of a mixed weibull distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  q: 	time
#  props:   vector of proportions (Must sum to 1. NB: these are not checked!)
#  betas:	vector of shape parameters (length must match length props)
#  alphas:	vector of scale parameters (length must match length props)
########################################################################################################
pmixwei <- function(q,props,betas,alphas){
  width <- length(props)
  leng <- length(q)
  qmat <- matrix(rep(q,width),ncol=width)
  propsmat <- matrix(rep(props,each=leng),ncol=width)
  alphasmat <- matrix(rep(alphas,each=leng),ncol=width)
  betasmat <- matrix(rep(betas,each=leng),ncol=width)
  rowSums(propsmat*pweibull(q=qmat,shape=betasmat,scale=alphasmat))
}

########################################################################################################
# Function for the PDF of a mixed weibull distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  x: 	time
#  props:   vector of proportions (Must sum to 1. NB: these are not checked!)
#  betas:	vector of shape parameters (length must match length props)
#  alphas:	vector of scale parameters (length must match length props)
########################################################################################################
dmixwei <- function(x,props,betas,alphas){
  width <- length(props)
  leng <- length(x)
  xmat <- matrix(rep(x,width),ncol=width)
  propsmat <- matrix(rep(props,each=leng),ncol=width)
  alphasmat <- matrix(rep(alphas,each=leng),ncol=width)
  betasmat <- matrix(rep(betas,each=leng),ncol=width)
  rowSums(propsmat*dweibull(x=xmat,shape=betasmat,scale=alphasmat))
}

########################################################################################################
# Function for the RF of a mixed weibull distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  n: 	number of random draws to make
#  props:   vector of proportions (Must sum to 1. NB: these are not checked!)
#  betas:	vector of shape parameters (length must match length props)
#  alphas:	vector of scale parameters (length must match length props)
########################################################################################################
rmixwei <- function(n,props,betas,alphas){
  assignment <- sample.int(n=length(props),size=n,replace=TRUE,prob=props)
  output <- rweibull(n=n,shape=betas[assignment],scale=alphas[assignment])
  return(output)
}

########################################################################################################
# Function for the inverse-CDF of a mixed weibull distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  p: 	probability
#  props:   vector of proportions (must sum to 1. NB: these are not checked!)
#  betas:	vector of shape parameters (length must match length props)
#  alphas:	vector of scale parameters (length must match length props)
########################################################################################################
qmixwei <- function(p,props,betas,alphas){
  q <- rep(NA,length(p))
  rootfunction <- function(t,p,props,betas,alphas){
    sum(props*exp(-(t/alphas)^betas))+p-1
  }
  for (i in 1:length(p)){
    limits <- qweibull(p[i],betas,alphas)
    q[i] <- uniroot(rootfunction,lower=min(limits),upper=max(limits),p=p[i],props=props,betas=betas,alphas=alphas)$root
  }
  return(q)
}

########################################################################################################
# Function for the CDF of a Log-Logistic distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  q: 	time
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
ploglog <- function(q, scale, shape){
  q ^ shape / (scale ^ shape + q ^ shape)
}

########################################################################################################
# Function for the PDF of a Log-Logistic distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  x: 	time
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
dloglog <- function(x, scale, shape){
  shape * scale ^ shape * x ^ (shape - 1) / (scale ^ shape + x ^ shape) ^ 2
}

########################################################################################################
# Function for the inverse-CDF of a Log-Logistic distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  p: 	probability
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
qloglog <- function(p, scale, shape){
  (p * scale ^ shape / (1 - p)) ^ (1/shape)
}

########################################################################################################
# Function for the RF of a Log-Logistic distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  n: 	number of random draws to make
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
rloglog <- function(n, scale, shape){
  qloglog(runif(n),scale,shape)
}

########################################################################################################
# Function for the PDF of a Gompertz distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  x: 	time
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
dgompertz <- function(x, scale, shape){
  scale * shape * exp(shape + scale * x - shape * exp(scale * x))
}

########################################################################################################
# Function for the CDF of a Gompertz distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  q: 	time
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
pgompertz <- function(q, scale, shape){
  1 - exp(shape - shape * exp(scale * q))
}
########################################################################################################
# Function for the inverse-CDF of a Gompertz distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  p: 	probability
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
qgompertz <- function(p, scale, shape){
  log(1 - log(1 - p) / shape) / scale
}

########################################################################################################
# Function for the RF of a Gompertz distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  n: 	number of random draws to make
#  scale:   scale parameter
#  shape:	shape parameter
########################################################################################################
rgompertz <- function(n, scale, shape){
  qgompertz(runif(n),scale,shape)
}

########################################################################################################
# Function for the PDF of a Generalised Gamma distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  x: 	time
#  scale:   scale parameter
#  shape:	shape parameter
#  family:  family parameter
########################################################################################################
dggamma <- function(x, scale, shape, family){
  family * x ^ (family * shape - 1) * exp(-(x / scale) ^ family) / (scale ^ (family * shape) * gamma(shape))
}

########################################################################################################
# Function for the CDF of a Generalised Gamma distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  q: 	time
#  scale:   scale parameter
#  shape:	shape parameter
#  family:  family parameter
########################################################################################################
pggamma <- function(q, scale, shape, family){
  pgamma((q/scale) ^ family, shape)
}

########################################################################################################
# Function for the RF of a Generalised Gamma distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# To maximise speed, the function uses a fast rejection sampling approach where shape > 1, 
#   and a slower qgamma-based approach where shape <= 1. 
# Takes as input:
#  n: 	number of random draws to make
#  scale:   scale parameter
#  shape:	shape parameter
#  family:  family parameter
########################################################################################################
rggamma <- function(n, scale, shape, family){
# This first approach can apply to any valid parameter combination, but is slower than rejection sampling
# Rejection sampling only works when shape >1 however, so we use this slower method to cover the missing parameter space. 

  ##qgamma-based approach
  if(shape <= 1){
    draws <- runif(n)
    return(qggamma(draws,scale,shape,family))
  }
  ## Rejection sampling approach
  ## N: number of samples required
  K = 4 * shape ^ shape * exp(-shape) / (sqrt(2 * shape - 1) * gamma(shape)) # expected number of samples required for 1 success
  N = n * K + 10 * sqrt(n * (K - 1) / K ^ 2) # total number of samples required (add 10*SE in order to guarantee a sufficient number)
  ## generate samples
  u1 <- runif(N)
  u2 <- runif(N)
  V = log(u1 / (1 - u1)) / (family * sqrt(2 * shape - 1))
  ## select the first n successful samples
  index <- which(log(u1 * u2 * (1 - u1)) <= shape + log(1/4) + shape * family * V - shape * exp(family * V))
  return(scale * shape ^ (1 / family) * exp(V[index[1:n]]))
}


########################################################################################################
# Function for the inverse-CDF of a Generalised Gamma distribution
# Meant for use solely within GESTATE's architecture: input checking done at object level
# Takes as input:
#  p: 	probability
#  scale:   scale parameter (theta)
#  shape:	shape parameter (eta)
#  family:  family parameter (rho)
########################################################################################################
qggamma <- function(p,scale, shape, family){
  q <- scale*qgamma(p,shape,1)^(1/family)
}
