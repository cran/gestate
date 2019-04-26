
##################################################################################################################################
# Functions for converting power to sample size and vice-versa
#
##################################################################################################################################

#'Calculate Schoenfeld or Event Proportion power from number of events
#'
#' @param events Number of events.
#' @param HR Hazard Ratio. Values below 1 indicate a benefit to the active arm vs control.
#' @param ratio Randomisation ratio (Schoenfeld), or event ratio (event proportion method). Default=1.
#' @param alpha1 1-sided alpha. Default=0.025.
#' @return Power as a decimal
#' @keywords internal
events2power <- function(events,HR,ratio=1,alpha1=0.025){
  power <- pnorm(sqrt(events*ratio)*-log(HR)/(ratio+1)- qnorm(1-alpha1))
  return(power)
}

#'Calculate Schoenfeld or Event Proportion based event numbers from power
#'
#' @param power Desired power (decimal).
#' @param HR Hazard Ratio. Values below 1 indicate a benefit to the active arm vs control.
#' @param ratio Randomisation ratio (Schoenfeld), or event ratio (event proportion method). Default=1.
#' @param alpha1 1-sided alpha. Default=0.025.
#' @return Required event number
#' @keywords internal
power2events <- function(power,HR,ratio=1,alpha1=0.025){
  events <- (((ratio+1)^2/ratio)*(qnorm(1-alpha1)+qnorm(power))^2)/(log(HR)^2)
  events[HR > 1] <- Inf
  return(events)
}

########################################################################################################
#'Calculate Normal distribution test power based on Effect size and Variance
#'
#' @param Z effect size.
#' @param V variance. Default=1.
#' @param alpha1 1-sided alpha. Default=0.025.
#' @return Power as a decimal
#' @keywords internal
ZV2power <- function(Z,V=1,alpha1=0.025){
  pnorm(abs(Z)/sqrt(V) - qnorm(1-alpha1))
}

########################################################################################################
#'Calculate Frontier power from number of events
#'
#' @param events Number of events.
#' @param HR Hazard Ratio. Values below 1 indicate a benefit to the active arm vs control.
#' @param Eratio Event ratio.
#' @param Rratio Randomisation ratio. Default=1.
#' @param startpower Initial estimate of power. Default=0.5.
#' @param alpha1 1-sided alpha. Default=0.025.
#' @param iter Number of iterations to perform. Default=10.
#' @return Power as a decimal
#' @author James Bell
#' @references Bell J, Power Calculations for Time-to-Event Trials Using Predicted Event Proportions, 2019, paper under review.
#' @examples frontierpower(events=300,HR=0.7,Eratio=1.2,Rratio=1.5,alpha1=0.025)
#' @export
frontierpower <- function(events,HR,Eratio,Rratio=1,startpower=0.5,alpha1=0.025,iter=10){
  HR <- pmin(pmax(1,events),pmax(pmin(1,1/events),HR))
  P <- Rratio/(Rratio+1)
  pi <- Eratio/(Eratio+1)
  power <- startpower
  cap <- pmin(0.1,P,1/P)
  for(i in 1:iter){
    pi2 <- pmin(1-cap,pmax(cap,P+(pi-P)*(qnorm(1-alpha1)/(qnorm(1-alpha1)+qnorm(power)))))
    ratio <- pi2/(1-pi2)
    power <- events2power(events=events,HR=HR,ratio=ratio,alpha1=alpha1)
  }
  return(power)
}

########################################################################################################
#'Calculate Freedman power from number of events
#'
#' @param events Number of events.
#' @param HR Hazard Ratio.
#' @param ratio Randomisation ratio. Default=1.
#' @param alpha1 1-sided alpha. Default=0.025.
#' @return Power as a decimal
#' @keywords internal
freedmanpower <- function(events,HR,ratio=1,alpha1=0.025){
  power <- pnorm((sqrt(events*ratio)*(1-HR)/(1+HR*ratio))- qnorm(1-alpha1))
  return(power)
}
