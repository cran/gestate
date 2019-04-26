##################################################################################################################################
# Composite function definitions, e.g. defining expected events
#
##################################################################################################################################

########################################################################################################
# Function for constructing the expected event function ('E')
#  Note that the D and N functions are complex, containing recruitment, events and censoring components
#
# Takes as input:
#  tim: 	time
#  assess:	time of assessment
#  Nfunction_control:	Function for the number at risk in the control arm
#  Dfunction_control:	Function for the number of events in the control arm
#  Nfunction_active:	Function for the number at risk in the active arm
#  Dfunction_active:	Function for the number of events in the active arm
########################################################################################################
Efunction_active <- function(tim,assess=assess,Nfunction_control,Dfunction_control,Nfunction_active,Dfunction_active){
  Nco <- Nfunction_control(tim=tim,assess=assess)
  Nac <- Nfunction_active(tim=tim,assess=assess)
  N <- pmax(Nco+Nac,0.0000000000001)
  D <- Dfunction_control(tim=tim,assess=assess)+Dfunction_active(tim=tim,assess=assess)
  Nac*D/N
}

########################################################################################################
# Function for constructing the expected variance function ('V')
#  Note that the D and N functions are complex, containing recruitment, events and censoring components
#
# Takes as input:
#  tim: 	time
#  assess:	time of assessment
#  Nfunction_control:	Function for the number at risk in the control arm
#  Dfunction_control:	Function for the number of events in the control arm
#  Nfunction_active:	Function for the number at risk in the active arm
#  Dfunction_active:	Function for the number of events in the active arm
########################################################################################################
Vfunction <- function(tim,assess=assess,Nfunction_control,Dfunction_control,Nfunction_active,Dfunction_active){
  Nco <- Nfunction_control(tim=tim,assess=assess)
  Nac <- Nfunction_active(tim=tim,assess=assess)
  N <- pmax(Nco+Nac,0.0000000000001)
  D <- Dfunction_control(tim=tim,assess=assess)+Dfunction_active(tim=tim,assess=assess)
  Nco*Nac*D/(N^2)
}

########################################################################################################
# Function for constructing the Greenwood function
#  Note that the D and N functions are complex, containing recruitment, events and censoring components
#
# Takes as input:
#  tim: 	time
#  assess:	time of assessment
#  Nfunction:	Function for the number at risk
#  Dfunction:	Function for the number of events
########################################################################################################
Greenwoodfunction <- function(tim,assess=assess,Nfunction,Dfunction){
  N <- pmax(Nfunction(tim=tim,assess=assess),0.0000000000001)
  D <- Dfunction(tim=tim,assess=assess)
  D/(N^2)
}

########################################################################################################
# Function for constructing the expected RMST for an arm
#
# Takes as input:
#  Sfunction:	Survival Function of the events
#  restriction:	Time of restriction
########################################################################################################
RMST <- function(Sfunction,restriction){
  integrate(Sfunction,lower=0,upper=restriction)$value
}

########################################################################################################
# Function for constructing the expected RMST variance function
#  Note that the D and N functions are complex, containing recruitment, events and censoring components
#
# Takes as input:
#  n: 	number of patients in the arm
#  Dfunction:	Function for the number of events
#  Sfunction:	Survival Function of the events
#  restriction:	Time of restriction
#  assess:	time of assessment
#
#  Output is a list, with:
#  [[1]] = underlying variance
#  [[2]] = censoring adjusted n
#  [[3]] = sample variance
#  [[4]] = number of patients required to produce number of events observed if no censoring
########################################################################################################
RMST_Var <- function(n,Dfunction,Sfunction,restriction,assess){
  #Formula for variance of area under the curve (V), see Royston 2013
  v1 <- function(r){
    r*Sfunction(r)
  }
  V1 <- integrate(v1,lower=0,upper=restriction)$value
  V2 <- integrate(Sfunction,lower=0,upper=restriction)$value
  V <- 2*V1-V2^2
  
  # Set of nested helper functions used for the n* calculation requiring double integration
  # Overall formula is: 
  # n*Integrate(0-->restrict){2*Integrate(0-->r){obsPDF}*CDF*(r-Integrate(0-->r){CDF})/(1-CDF)}/V
  RMST_nstar3 <- function(r,Dfunction,Sfunction,assess){
    2*integrate(Dfunction_wrapper,lower=0,upper=r,Dfunction=Dfunction,assess=assess)$value*Sfunction(r)*(r-integrate(Sfunction_wrapper,lower=0,upper=r,Sfunction=Sfunction)$value)/(1-Sfunction(r))
  }
  RMST_nstar2 <- function(r,Dfunction,Sfunction,assess){
    sapply(r,RMST_nstar3,Dfunction=Dfunction,Sfunction=Sfunction,assess=assess)
  }
  RMST_nstar <- function(n,V,Dfunction,Sfunction,restriction,assess){
    integrate(RMST_nstar2,lower=0,upper=restriction,Dfunction=Dfunction,Sfunction=Sfunction,assess=assess)$value/V
  }

  #Calculation of n*, the censoring-adjusted 'n'
  nstar <- RMST_nstar(n=n,V=V,Dfunction=Dfunction,Sfunction=Sfunction,restriction=restriction,assess=assess)
  #Calculation of number of patients required to reproduce number of events if no censoring were present
  neff <- integrate(Dfunction,lower=0,upper=restriction,assess=assess)$value/(1-Sfunction(restriction))
  
  out <- vector("list",4)
  out[[1]] <- V
  out[[2]] <- nstar
  out[[3]] <- V/nstar
  out[[4]] <- neff
  
  return(out)
}

########################################################################################################
# Dfunction_wrapper, Sfunction_wrapper are two simple functions that wrap Sfunction and Dfunction
# within sapply to allow for double integration - used for calculation of n*
#
########################################################################################################
Dfunction_wrapper <- function(r,Dfunction,assess){
  sapply(r,Dfunction,assess=assess)
}

Sfunction_wrapper <- function(r,Sfunction){
  sapply(r,Sfunction)
}
