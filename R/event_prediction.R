if(getRversion() >= "2.15.1")  utils::globalVariables(c("Ialphafunction", "Ibetafunction","Dfunction","Cfunction","IalphaNARfunction","IbetaNARfunction"))

###############################################################################################################################
#'Fit Weibull and Log-Normal survival curves to Kaplan Meier estimates
#'
#' This is a function to fit Weibull and log-normal curves to Survival data in life-table form using non-linear regression.\cr
#' By default it fits both, then picks the best fit based on the lowest (un)weighted residual sum of squares.\cr
#' Alternatively, just one shape may be fitted, by changing the 'type' argument to either "Weibull" or "Lognormal".
#' Weighted or unweighted fitting are possible. In general, weighted fitting using the number at risk as the weights seems to work best.\cr
#' This function is primarily used by event_prediction function, but also useful for general KM curve fitting.\cr
#' One useful aspect of this is for fitting the 'inverse KM', where drop-outs are events, while events and 'time-outs' are censored.
#'    This allows for finding a suitable parameterisation for the censoring curve.\cr
#' Primary advantage over likelihood-based methods is ability to use aggregated, rather than patient-level data.
#' Primary disadvantage is that the covariance matrix is heavily underestimated (essentially unusable) due to strong correlation between the input data points going into the regression.
#' @param KMcurve The dataframe object containing the survival data in lifetable form
#' @param Survival The column name for the survival function (i.e. the probabilities). Default is "Survival"
#' @param Time The column name for the times. Default is "Time"
#' @param type Type of event curve to fit.Default is "Automatic", fitting both Weibull and Log-normal curves.
#'     Alternatively accepts "Weibull" or "Lognormal" to force the type.
#' @param weighting Boolean for whether to use weighting. Default=TRUE as it greatly improves curve fitting.
#' @param Weights Name of Weights column. Default="Weights". Optional if weighting=FALSE. Recommended to use number at risk or remaining.
#' @param Weight_power Power to raise the weights to. Useful in large trials to give added weight to later points where numbers may still be high. Default=1 (Use weights as specified)
#' @param startbeta Starting value for the Weibull beta (shape) parameter to be used in the non-linear regression. Default=1 (exponential).
#' @param startsigma Starting value for the Lognormal sigma (sd) parameter to be used in the non-linear regression. Default=1.
#' @return Returns a 3-item list providing information needed to define a Curve object:
#' \itemize{
#'  \item{"Item 1"}{The type of Curve object fitted.}
#'  \item{"Item 2"}{A list of fitted parameters for the curve type.}
#'  \item{"Item 3"}{A vector containing the covariance-matrix parameters for the curve type. Note that these are heavily underestimated as the data going into the regression is not independent.}
#'  \item{"Item 4"}{A data frame containing the goodness of fit metrics for each curve type.}
#' }
#' @author James Bell
#' @examples recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' example_data_short <- simulate_trials(active_ecurve=Weibull(50,0.8),
#' control_ecurve=Weibull(50,0.8), rcurve=recruit,assess=10,iterations=1,seed=12345,
#' detailed_output=TRUE)
#'
#' library(survival)
#'
#' temp1  <- summary(survfit(Surv(example_data_short[,"Time"],1-example_data_short[,"Censored"])~ 1,
#' error="greenwood"))
#' out1 <- cbind(temp1$time,temp1$n.risk,temp1$surv,temp1$std.err)
#' out1 <- rbind(c(0,out1[1,2],1,0),out1)
#' colnames(out1) <- c("Time","NAR","Survival","Std.Err")
#' x1 <- ceiling(max(out1[,"Time"]))
#' example_lifetable <- out1[findInterval(0:x1,out1[,"Time"]),]
#' example_lifetable[,"Time"] <- 0:x1
#'
#' fit_KM(KMcurve=example_lifetable,Survival="Survival",Time="Time",Weights="NAR",type="automatic")
#' @export
fit_KM <- function(KMcurve,Survival="Survival",Time="Time",type=c("automatic","Weibull","Lognormal"),weighting=TRUE,Weights="Weights",Weight_power=1,startbeta=1,startsigma=1){
  #Multiple choice
  type <- match.arg(type)

  if(missing(KMcurve))stop("Please specify the life-table using the 'KMcurve' argument.")

  if(!(Time %in% colnames(KMcurve)))stop("Error: Specified name of Time column does not exist.")
  if(!(Survival %in% colnames(KMcurve)))stop("Error: Specified name of Survival column does not exist.")

  if(!is.logical(weighting))stop("Error: weighting argument must be boolean: default=TRUE (weighting is to be used).")
  if(weighting){
    if(!(Weights %in% colnames(KMcurve)))stop("Error: Weighting specified, but name of Weights column does not exist.")
    if(!is.numeric(Weight_power))stop("Error: Weight_power argument must be numerical. default=1.")
  }
  if(!is.numeric(startbeta))stop("Error: startbeta argument must be numerical. default=1.")
  if(!is.numeric(startsigma))stop("Error: startbeta argument must be numerical. default=1.")

  # Y is the Survival function
  Y <- KMcurve[,Survival]
  # T is the time variable
  T <- KMcurve[,Time]

  # If Weibull or automatic fitting specified, do a Weibull fit
  if(type=="automatic"|type=="Weibull"){
  # Automatic start function based upon startbeta and the final time point
    if(any(Y==0.5)){
      median <- min(T[Y==0.5])
      startalpha <- median /((log(2)^(1/startbeta)))
    } else if(min(Y)<0.5){
      Y2 <- max(Y[Y < 0.5])
      Y1 <- min(Y[Y > 0.5])
      T2 <- min(T[Y < 0.5])
      T1 <- max(T[Y > 0.5])
      median <- T2*(0.5-Y2)/((Y1-0.5)+(0.5-Y2))+T1*(Y1-0.5)/((Y1-0.5)+(0.5-Y2))
      startalpha <- median /((log(2)^(1/startbeta)))
    } else{
      startalpha= max(T)/(-log(min(Y)))^(1/startbeta)
    }
  #Perform a weighted or unweighted Weibull curve fit. If no weighting, you should only supply point with sufficient patients  at risk
  # If weighting, recommended to start with the number of patients at risk as weighting
    fitW <- Inf
    if(weighting){
      #weights are labelled W
      W <- (KMcurve[,Weights])^Weight_power

      tryCatch(fitW <- nls(Y ~ 1-pweibull(q=T,scale=alpha,shape=beta),weights=W, start=list(alpha=startalpha,beta=startbeta))
        ,error=function(e){warning("Non-linear regression unable to fit Weibull curve\n")})

    } else{
      tryCatch(fitW <- nls(Y ~ 1-pweibull(q=T,scale=alpha,shape=beta),start=list(alpha=startalpha,beta=startbeta))
        ,error=function(e){warning("Non-linear regression unable to fit Weibull curve\n")})
    }

    WSS <- tryCatch(deviance(fitW),error=function(e){return(Inf)})
  }
  # If Lognormal or automatic fitting specified, do a Lognormal fit
  if(type=="automatic"|type=="Lognormal"){
  # Automatic start function based upon startbeta and the final time point
    if(any(Y==0.5)){
      startmu <- log(min(T[Y==0.5]))
    } else if(min(Y)<0.5){
      Y2 <- max(Y[Y < 0.5])
      Y1 <- min(Y[Y > 0.5])
      T2 <- min(T[Y < 0.5])
      T1 <- max(T[Y > 0.5])
      startmu<- log(T2*(0.5-Y2)/((Y1-0.5)+(0.5-Y2))+T1*(Y1-0.5)/((Y1-0.5)+(0.5-Y2)))
    } else{
      startmu= log(-0.5*max(T)/(min(Y)-1))
    }

  #Perform a weighted or unweighted Lognormal curve fit. If no weighting, you should only supply point with sufficient patients  at risk
  # If weighting, recommended to start with the number of patients at risk as weighting
    fitLN <- Inf
    if(weighting){
      #weights are labelled W
      W <- (KMcurve[,Weights])^Weight_power
      tryCatch(fitLN <- nls(Y ~ 1-plnorm(T,mu,sigma),weights=W, start=list(mu=startmu,sigma=startsigma))
        ,error=function(e){warning("Non-linear regression unable to fit Lognormal curve\n")})
    } else{
      tryCatch(fitLN <- nls(Y ~ 1-plnorm(T,mu,sigma), start=list(mu=startmu,sigma=startsigma))
        ,error=function(e){warning("Non-linear regression unable to fit Lognormal curve\n")})
    }
    LNSS <- tryCatch(deviance(fitLN),error=function(e){return(Inf)})
  }

  metric <- ifelse(weighting,"Weighted Sum of Squares", "Sum of Squares")
  if(type=="Weibull"){
    curve <- "Weibull"
    stats <- WSS
  } else if (type=="Lognormal"){
    curve <- "Lognormal"
    stats <- LNSS
  } else if(type=="automatic"){
    curve <- c("Weibull","Lognormal")
    metric <- rep(metric,2)
    stats <- c(WSS,LNSS)
    if(WSS == LNSS && LNSS == Inf)stop("Unable to fit KM curve using any method. More data points probably needed.")
    if(WSS < LNSS){
      type="Weibull"
    } else{
      type="Lognormal"
    }
  }
  
  Fit <- cbind(curve,metric,stats)
  colnames(Fit) <- c("Distribution","Metric","Value")

  if(type=="Weibull"){
    if(WSS == Inf)stop("Unable to fit KM curve using specified Weibull method. Try Lognormal; otherwise more data points probably needed.")
  # Extract the Weibull coefficients
    coefs <- coef(summary(fitW))[,1]
    VCov <- c(vcov(fitW)[1,1],vcov(fitW)[2,2],vcov(fitW)[1,2])
    names(coefs) <- c("Alpha","Beta")
    names(VCov) <- c("Alpha_Var","Beta_Var","Covariance")
  }else if(type=="Lognormal"){
    if(LNSS == Inf)stop("Unable to fit KM curve using specified Lognormal method. Try Weibull; otherwise more data points probably needed.")
  # Extract the Lognormal coefficients
    coefs <- coef(summary(fitLN))[,1]
    VCov <- c(vcov(fitLN)[1,1],vcov(fitLN)[2,2],vcov(fitLN)[1,2])
    names(coefs) <- c("Mu","Sigma")
    names(VCov) <- c("Mu_Var","Sigma_Var","Covariance")
  }
  output <- list(Curvetype=type,Parameters=coefs,VCov=VCov,Fit=Fit)
  return(output)
}

###############################################################################################################################
#'Fit Weibull and Log-Normal survival curves to patient-level time-to-event data
#'
#' This is a function to fit Weibull and log-normal curves to patient-level Survival data using maximum likelihood estimation.\cr
#' By default it fits both, then picks the best fit based on the log-likelihood (and implicitly the AIC).\cr
#' Alternatively, just one shape may be fitted, by changing the 'type' argument to either "Weibull" or "Lognormal".
#' This function is primarily used by event_prediction_data function, but also useful for general Survival function curve fitting.\cr
#' One useful aspect of this is for fitting the 'inverse KM', where drop-outs are events, while events and 'time-outs' are censored.
#'    This allows for finding a suitable parameterisation for the censoring curve.\cr
#' Where patient-level data is available, this function will typically perform substantially better than fit_KM, with lower variability of point estimates (and more accurate quantification of it).
#' @param data The dataframe object containing the patient-level survival data
#' @param Time The column name for the times. Default is "Time"
#' @param Event The column name for the events column (i.e. the binary variable denoting events vs censorings). Default is "Event"
#' @param censoringOne Specify whether censoring is denoted in the Event column by a one (TRUE) or zero (FALSE). Default=FALSE (censorings denoted by 0, events by 1)
#' @param type Type of event curve to fit. Default is "Automatic", fitting both Weibull and Log-normal curves.
#'     Alternatively accepts "Weibull" or "Lognormal" to force the type.
#' @param init Vector of starting values for parameter values; useful if survreg experiences convergence issues. Default=NULL (no values specified)
#' @return Returns a 3-item list providing information needed to define a Curve object:
#' \itemize{
#'  \item{"Item 1"}{The type of Curve object fitted.}
#'  \item{"Item 2"}{A list of fitted parameters for the curve type.}
#'  \item{"Item 3"}{A vector containing the covariance-matrix parameters for the curve type.}
#'  \item{"Item 4"}{A data frame containing the goodness of fit metrics for each curve type.}
#' }
#' @author James Bell
#' @examples recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' example_data <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),
#' rcurve=recruit, assess=10,iterations=1,seed=12345,detailed_output=TRUE)
#'
#' fit_tte_data(data=example_data,Time="Time",Event="Censored",censoringOne=TRUE,type="automatic")
#' @export
fit_tte_data <- function(data,Time="Time",Event="Event",censoringOne=FALSE,type=c("automatic","Weibull","Lognormal"),init=NULL){
   #Multiple choice
  type <- match.arg(type)

  #data
  if(missing(data))stop("Error: Please specify the data using the 'data' argument.")

  #Strings / column names
  if(!is.character(Event))stop("'Event' argument must be a string")
  if(!is.character(Time))stop("'Time' argument must be a string")

  if(!(Time %in% colnames(data)))stop("Error: Specified name of Time column does not exist.")
  if(!(Event %in% colnames(data)))stop("Error: Specified name of Event column does not exist.")
  if(!is.logical(censoringOne))stop("Error: censoringOne argument must be boolean: default=FALSE (0 is censored, 1 is event).")

  #Booleans
  if(!is.logical(censoringOne))stop("Error: censoringOne argument must be boolean: default=FALSE (0 is censored, 1 is event).")

  #Numbers
  if(!is.null(init) && (!is.numeric(init) || min(init < 0)))stop("Error: 'init' argument must be NULL or a vector of two positive numbers corresponding to the initiation values for the MLE. This should only be specified if there are numerical convergence issues. default=NULL (automatic starting values used).")

  # T is the time variable
  T <- data[,Time]
  E <- data[,Event]
  if(censoringOne){E <- 1-E}

  if(type=="automatic"|type=="Lognormal"){
    fitLN <- survreg(Surv(T,E)~1,dist="lognormal",init=init)
    LNLL <- fitLN$loglik[1]
  }
  if(type=="automatic"|type=="Weibull"){
    fitW <- survreg(Surv(T,E)~1,dist="weibull",init=init)
    WLL <- fitW$loglik[1]
  }

  metric <- "Log-Likelihood"
  if(type=="Weibull"){
    curve <- "Weibull"
    stats <- WLL
  } else if (type=="Lognormal"){
    curve <- "Lognormal"
    stats <- LNLL
  } else if(type=="automatic"){
    curve <- c("Weibull","Lognormal")
    metric <- rep(metric,2)
    stats <- c(WLL,LNLL)
    if(WLL >= LNLL){
      type="Weibull"
    } else{
      type="Lognormal"
    }
  }

  Fit <- cbind(curve,metric,stats)
  colnames(Fit) <- c("Distribution","Metric","Value")

  if(type=="Weibull"){
  # Extract the Weibull coefficients
    rawcoefs <- fitW$icoef
    coefs <-c(exp(rawcoefs[1]),1/exp(rawcoefs[2]))
    VCov <- c(fitW$var[1,1]*coefs[1]^2,fitW$var[2,2]*coefs[2]^2,fitW$var[1,2]*-coefs[1]*coefs[2])
    names(coefs) <- c("Alpha","Beta")
    names(VCov) <- c("Alpha_Var","Beta_Var","Covariance")
  }else if(type=="Lognormal"){
  # Extract the Lognormal coefficients
    rawcoefs <- fitLN$icoef
    coefs <- c(rawcoefs[1],exp(rawcoefs[2]))
    VCov <- c(fitLN$var[1,1],fitLN$var[2,2]*coefs[2]^2,fitLN$var[1,2]*coefs[2])
    names(coefs) <- c("Mu","Sigma")
    names(VCov) <- c("Mu_Var","Sigma_Var","Covariance")
  }
  output <- list(Curvetype=type,Parameters=coefs,VCov=VCov,Fit=Fit)
  return(output)
}

#####################################################################################################################
#'Event prediction using patient-level survival data and a recruitment RCurve
#'
#' This is a function to perform event prediction\cr
#' It uses the fit_KM_tte_data function to perform MLE regression of Weibull and log-normal curves to the provided survival data.\cr
#' It creates an event Curve object from this, and combines it with a recruitment RCurve and an optional dropout(censoring) Curve.\cr
#' Using the same numerical integration approach as nph_curve_trajectories it performs an unconditional event prediction.\cr
#' If a conditioning time, event number (and preferably number at risk) are provided, a conditional event prediction is also calculated.\cr
#' Analytic standard errors for conditional and unconditional event numbers are provided for the whole trajectory.\cr
#' SEs calculated by propagating parameter estimate errors through the integrals by the delta method and then invoking a beta-binomial distribution.\cr\cr
#' For event prediction, conditional predictions with the Conditional SE of Prediction are most accurate and appropriate.\cr
#' Unconditional predictions should be close to conditional ones but technically relate to predictions if the trial were rerun, rather than this specific instance.
#' Point estimates are usually very close to the unconditional ones, but the CIs are typically much wider than necessary.
#' The conditional and unconditional SEs of fitting relate to the accuracy of the estimated mean event number at a given time, rather than the spread of future observations.
#' The conditional and unconditional SEs of prediction relate to the accuracy of prediction of future observations, and should therefore be used for event prediction.
#' Note that the Prediction SEs are wider than the Fitting SEs as they also take into account the binomial uncertainty of events occurring (beta-binomial model).
#' @param data The dataframe object containing the patient-level survival data
#' @param Time The column name for the times. Default is "Time"
#' @param Event The column name for the events column (i.e. the binary variable denoting events vs censorings). Default is "Event"
#' @param censoringOne Specify whether censoring is denoted in the Event column by a one (TRUE) or zero (FALSE). Default=FALSE (censorings denoted by 0, events by 1)
#' @param type Type of event curve to fit. Default is "Automatic", fitting both Weibull and Log-normal curves.
#'     Alternatively accepts "Weibull" or "Lognormal" to force the type.
#' @param rcurve Observed and/or expected recruitment distribution as an RCurve object. This should typically be of PieceR type (piecewise linear recruitment).
#' @param max_time Maximum time to predict events up to.
#' @param dcurve Dropout/censoring distribution as a Curve object. This is Blank() by default, i.e. no dropout.
#' @param CI Number between 0 and 1 for the size of CI to calculate. Default is 0.95 (95 percent CI).
#' @param condition Boolean whether to also do a conditional event prediction. Default=FALSE
#'     Note that If all conditioning options are left as defaults, conditioned calculation will equal the unconditional one.
#' @param cond_Events Number of events to condition on. Default=0. Optional unless condition=TRUE.
#' @param cond_NatRisk Number of patients at risk to condition on. Default=-1.
#'     By default, the program will estimate the number at risk assuming no censoring. It is highly recommended to specify this if conditioning. Optional unless condition=TRUE.
#' @param cond_Time Time, in months, to condition on. Default=0. Optional unless condition=TRUE.
#' @param units Units that the KM-curve is specified in. Accepts "Days", "Months". Default="Days".
#' @param init Vector of starting values for parameter values; useful if survreg experiences convergence issues. Default=NULL (no values specified)
#' @param discountHR Hazard ratio for discounting events e.g. used to predict adjudicated events from unadjudicated data where patients remain 'at risk' after an event is adjudicated not to have occurred.
#'     Values below 1 indicate fewer events will occur than predicted by the curve-fitting.
#'     Note that changing this argument is only allowed if type="Weibull" since log-normal curves are not compatible with proportional hazards.
#'     Default=1 (No discounting)
#' @return Returns a list object with the fitted ecurve, the dcurve, the rcurve and a summary table with one row per month up to max_time containing the following columns:
#' \itemize{
#'  \item{"Assessment_Time"}{Time of assessment.}
#'  \item{"Patients"}{Number of patients recruited by the assessment time.}
#'  \item{"Predicted_Events"}{Number of events unconditionally predicted at the assessment time.}
#'  \item{"SE_Fitting"}{SE of the estimate of the fitted mean. Note that this corresponds to the accuracy of the estimate of the underlying parameter, not future observed event numbers.}
#'  \item{"SE_Prediction"}{SE of event prediction.}
#'  \item{"Prediction_Lower"}{Lower bound of X percent CI of unconditional event prediction, where X is the 'CI' argument. This CI is based on the quantiles of the beta-binomial distribution and hence is asymmetric about the predicted event number.}
#'  \item{"Prediction_Upper"}{Upper bound of X percent CI of unconditional event prediction, where X is the 'CI' argument. This CI is based on the quantiles of the beta-binomial distribution and hence is asymmetric about the predicted event number.}
#'  \item{"Conditioned_Events"}{Number of events unconditionally predicted at the assessment time (Column present only if conditioning specified).}
#'  \item{"Cond_SE_Fitting"}{SE of the estimate of the fitted conditional mean. Note that this corresponds to the accuracy of the estimate of the underlying parameter, not future observed event numbers (Column present only if conditioning specified).}
#'  \item{"Cond_SE_Prediction"}{SE of the conditional event prediction (Column present only if conditioning specified).}
#'  \item{"Cond_Prediction_Lower"}{Lower bound of X percent CI of conditional event prediction, where X is the 'CI' argument. This CI is based on the quantiles of the beta-binomial distribution and hence is asymmetric about the predicted event number.}
#'  \item{"Cond_Prediction_Upper"}{Upper bound of X percent CI of conditional event prediction, where X is the 'CI' argument. This CI is based on the quantiles of the beta-binomial distribution and hence is asymmetric about the predicted event number.}
#' }
#' @author James Bell
#' @references Bell J, unpublished work.
#' @examples recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' trial_short <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),
#' rcurve=recruit, assess=10,iterations=1,seed=12345,detailed_output=TRUE)
#'
#' predictions <- event_prediction(data=trial_short, Event="Censored", censoringOne=TRUE, 
#' type="Weibull", rcurve=recruit, max_time=60, condition=TRUE, cond_Events=49, cond_NatRisk=451, 
#' cond_Time=10, units="Months")
#' @export
event_prediction <- function(data, Time="Time", Event="Event",censoringOne=FALSE, type=c("automatic","Weibull","Lognormal"), rcurve,max_time=100, dcurve=Blank(), CI=0.95,condition=FALSE,cond_Events=0,cond_NatRisk=NULL,cond_Time=0,units=c("Days","Months"),init=NULL,discountHR=1){
# Input checks to ensure valid inputs
# Firstly check for missing variables
  if(missing(data))stop("Please supply the name of the event dataset using the 'KMcurve' argument")
  if(missing(rcurve))stop("Please specify the recruitment distribution using the 'rcurve' argument. Please note that this should be an RCurve object; Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")

# Secondly, check arguments are all of the correct type
  #Multiple choices
  type <- match.arg(type)
  units <- match.arg(units)

  #Strings / column names
  if(!is.character(Event))stop("'Event' argument must be a string")
  if(!is.character(Time))stop("'Time' argument must be a string")
  if(!(Time %in% colnames(data)))stop("Error: Specified name of Time column does not exist.")
  if(!(Event %in% colnames(data)))stop("Error: Specified name of Event column does not exist.")

  #Curves
  if(class(rcurve)[1]!= "RCurve") stop("Argument 'rcurve' must be an RCurve object in order to define a recruitment distribution. Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")
  if(class(dcurve)[1]!= "Curve") stop("Argument 'dcurve' must be a Curve object in order to define the censoring curve. By default, no censoring is specified.")

  #Booleans
  if(!is.logical(censoringOne))stop("Error: censoringOne argument must be boolean: default=FALSE (0 is censored, 1 is event).")

  #Numbers
  if(!is.numeric(max_time) || max_time%%1!=0 || max_time < 1 ) stop("Error: Please specify a positive integer for the maximum time to predict events up to, using the 'max_time' argument")
  if(!is.numeric(discountHR) || discountHR < 0 ) stop("Error: Please specify a positive number for the 'discountHR'. Default = 1 (no discounting)")
  if(!is.numeric(CI) || !(CI > 0) || !(CI < 1)) stop("Error: Please specify a value between 0 and 1 (exclusive) for the 'CI'. Default = 0.95 (95% CI)")

  if(!is.null(init) && (!is.numeric(init) || min(init < 0)))stop("Error: 'init' argument must be NULL or a vector of two positive numbers corresponding to the initiation values for the MLE. This should only be specified if there are numerical convergence issues. default=NULL (automatic starting values used).")
  if(!(discountHR ==1) && !(type=="Weibull"))stop("Error: 'DiscountHR' may only be set when 'type='Weibull'':Discounting of events is only possible with Weibull curves as it requires a proportional hazards assumption.")

  #Conditioning
  if(!is.logical(condition))stop("Error: 'condition' argument must be boolean: default=FALSE (no conditioning).")
  if(condition){
   if(!is.numeric(cond_Events) || cond_Events%%1!=0 || cond_Events < 0 || cond_Events > getN(rcurve))stop("Error: 'cond_Events' argument must be a positive integer no greater than the number of patients. default=0 (0 events).")
    if(!is.numeric(cond_Time) || cond_Time%%1!=0 || cond_Time < 0 || cond_Time >= max_time)stop("Error: 'cond_Time' argument must be a positive integer smaller than max_time. default=0 (conditioning from start time).")
    if(!is.null(cond_NatRisk) && (!is.numeric(cond_NatRisk) || cond_NatRisk%%1!=0 || cond_NatRisk < 0 || cond_NatRisk > getN(rcurve)))stop("Error: 'cond_NatRisk' argument must be NULL or a positive integer no greater than the number of patients. default=NULL (conditioning on expected number at risk).")
  }

### Warnings and notes

  if(condition){
    if(cond_Events == 0){
      warning("Conditional event prediction has been requested but no observed events have been specified - results may not be reasonable.\n")
    }
    if(is.null(cond_NatRisk)){
      message("NOTE: No number at risk at the time of conditioning has been specified. The program by default estimates this as the predicted number minus the difference between the observed and predicted events at the time of conditioning. This will be accurate only if the predicted censoring distribution is also accurate. It is strongly recommended to supply the observed number at risk at the time of conditioning.\n")
    }
    if(cond_Time == 0){
      warning("No conditioning time has been specified. The program by default will condition from the start of the trial (time 0).\n")
    }
  }

####
# Required functions

qbb2 <- function(p,mu,sigma2,n){
  alpha <- max(0,((1-mu)/sigma2 - 1/mu)*mu^2)
  beta <- alpha*(1/mu - 1)
  range <- 0:n
  x <- cumsum(exp(lbeta(range+alpha, n-range+beta) - lbeta(alpha, beta) + lchoose(n, range)))
  a <- which.max(x > p)
  a <- ifelse(length(a)==0,NA,a)
  return(a-1)
}

# Function for finding quantiles of beta-binomial distribution
qbetabin <- Vectorize(qbb2)

# Function for calculating variance of a beta-binomial distribution
betabinomialVar <- function(mu,sigma2,n){
  alpha <- ((1-mu)/sigma2 - 1/mu)*mu^2
  beta <- alpha*(1/mu - 1)
  var <- sigma2*n*(alpha+beta+n)
}

# Set of functions for Weibull event prediction CIs
Weibull_d_alpha <- function(x,alpha,beta){
  y <- (x/alpha)^beta
  d_alpha <- exp(-y)*y*(y-1)*(beta^2)/(x*alpha)
}

Weibull_d_beta <- function(x,alpha,beta){
  y <- (x/alpha)^beta
  d_beta <- (exp(-y)*y/x)*(1+log(y)-y*log(y))
}

Weibull_ST_d_alpha <- function(x,alpha,beta){
  y <- (x/alpha)^beta
  CDFd_alpha <- exp(-y)*y*beta/alpha
}

Weibull_ST_d_beta <- function(x,alpha,beta){
  y <- (x/alpha)^beta
  CDFd_beta <- -exp(-y)*log(y)*y/beta
}

# Set of functions for Lognormal event prediction CIs
LN_d_mu <- function(x,mu,sigma){
  y <- (log(x)-mu)/(sigma*sqrt(2))
  d_mu <- y*exp(-(y^2))/(x*sqrt(pi)*sigma^2)
}

LN_d_sigma <- function(x,mu,sigma){
  y <- (log(x)-mu)/(sigma*sqrt(2))
  d_sigma <- (2*y^2-1)*exp(-(y^2))/(x*sqrt(2*pi)*sigma^2)
}

LN_ST_d_mu <- function(x,mu,sigma){
  y <- (log(x)-mu)/(sigma*sqrt(2))
  CDFd_mu <- -sqrt(2)*exp(-(y^2))/(sqrt(pi)*sigma)
}

LN_ST_d_sigma <- function(x,mu,sigma){
  y <- (log(x)-mu)/(sigma*sqrt(2))
  CDFd_sigma <- -2*y*exp(-(y^2))/(sqrt(pi)*sigma)
}
#####

  #Fit chosen distribution(s) to patient-level.

  fitted <- fit_tte_data(data=data,Time=Time,Event=Event,censoringOne=censoringOne,type=type,init=init)

  if(units=="Days"){
    conversion_factor <- 365.25/12
  } else if(units=="Months"){
    conversion_factor <- 1
  } else {
    warning("Unrecognised unit specified for units:",units, "\nPlease specify either \"Days\" or \"Months\".Assuming units of days by default.\n")
    conversion_factor <- 365.25/12
  }
  if(fitted[[1]]=="Weibull"){
    discount <- discountHR^(1/fitted[[2]][2])
    fitted[[2]][1] <- fitted[[2]][1]/(conversion_factor*discount)
    fitted[[3]][1] <- fitted[[3]][1]/((conversion_factor*discount)^2)
    fitted[[3]][3] <- fitted[[3]][3]/(conversion_factor*discount)
    ecurve <- Weibull(alpha=fitted[[2]][1],beta=fitted[[2]][2])
  } else if(fitted[[1]]=="Lognormal") {
    if(discountHR !=1){warning("A discount HR for event occurrence has been specified but a log-normal curve has been fitted - no discounting will be applied due to incompatibility.\n")}
    fitted[[2]][1] <- fitted[[2]][1]-log(conversion_factor)
    ecurve <- Lognormal(mu=fitted[[2]][1],sigma=fitted[[2]][2])
  }

  N <- getN(rcurve)

# Self-modifying R code! These 2 statements create the two required functions tailored specifically to the combinations of events, censoring and recruitment
# found within each arm to calculate the events and number at risk. Relatively fast as only needs to do this once per run, and allows code to adapt to all sorts of input.

  eval(parse(text=paste(
    "Dfunction <- function(tim,assess){",
      N,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*",getPDFfunction(ecurve,x="tim"),
    "}"
  )))

# Write function for numbers at risk in the control arm
  eval(parse(text=paste(
    "Cfunction <- function(tim,assess){",
      N,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(ecurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),
    "}"
  )))

# Set up variance integrals

# More self-modifying R code! These 2 statements create the two required functions tailored specifically to the combinations of events, censoring and recruitment
# found within each arm to calculate the events and number at risk. Relatively fast as only needs to do this once per run, and allows code to adapt to all sorts of input.

  if(fitted[[1]]=="Weibull"){
    eval(parse(text=paste(
      "Ialphafunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*Weibull_d_alpha(x=tim,alpha=fitted[[2]][1],beta=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "Ibetafunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*Weibull_d_beta(x=tim,alpha=fitted[[2]][1],beta=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "IalphaNARfunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),"*Weibull_ST_d_alpha(x=tim,alpha=fitted[[2]][1],beta=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "IbetaNARfunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),"*Weibull_ST_d_beta(x=tim,alpha=fitted[[2]][1],beta=fitted[[2]][2])}"
    )))

  } else if(fitted[[1]]=="Lognormal") {
    eval(parse(text=paste(
      "Ialphafunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*LN_d_mu(x=tim,mu=fitted[[2]][1],sigma=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "Ibetafunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*LN_d_sigma(x=tim,mu=fitted[[2]][1],sigma=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "IalphaNARfunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),"*LN_ST_d_mu(x=tim,mu=fitted[[2]][1],sigma=fitted[[2]][2])}"
    )))

    eval(parse(text=paste(
      "IbetaNARfunction <- function(tim,assess){",
        "(1-",getAssessCDFfunction(rcurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),"*LN_ST_d_sigma(x=tim,mu=fitted[[2]][1],sigma=fitted[[2]][2])}"
    )))
  }

# Set up list to store output, and a 'single line' temporary storage vector
  output_storage <- vector("list", max_time)
  temp_storage <- data.frame(t(rep(NA,7)))
  colnames(temp_storage) <- c(
    "Time",
    "Patients",
    "Predicted_Events",
    "SE_Fitting",
    "SE_Prediction",
    "Prediction_Lower",
    "Prediction_Upper"
   )

  Ialpha <- rep(NA,max_time)
  Ibeta <- rep(NA,max_time)

  #SE_fitting is the variance associated with the curve fitting only
  for(i in 1:max_time){
    Ialpha[i] <- tryCatch(integrate(Ialphafunction,lower=0,upper=i,assess=i)$value,error=function(e){return(Inf)})
    Ibeta[i]  <- tryCatch(integrate(Ibetafunction ,lower=0,upper=i,assess=i,abs.tol=1e-30)$value,error=function(e){return(Inf)})
  }


  isBetaInf <- is.infinite(Ibeta)
  if(length(isBetaInf)>0){
    for(i in c(1:max_time)[isBetaInf]){
      below <- c(1:max_time) < i
      above <- c(1:max_time) > i
      lower <- Ibeta[!isBetaInf & below]
      upper <- Ibeta[!isBetaInf & above]
      if(length(upper)>=1 & length(lower)>=1){
        Ibeta[i] <- max(lower[length(lower)],upper[1])
      }
    }
  }

  for(i in 1:max_time){
    temp_storage$Time <- i
    temp_storage$Patients <- getPatients(rcurve,i)
    temp_storage$Predicted_Events <- integrate(Dfunction,lower=0,upper=i,assess=i)$value
    temp_storage$SE_Fitting <- tryCatch(N*pmin(sqrt(Ialpha[i]^2*fitted[[3]][1]+Ibeta[i]^2*fitted[[3]][2]+2*Ialpha[i]*Ibeta[i]*fitted[[3]][3]),1/sqrt(12)),error=function(e){return(N/sqrt(12))})
    output_storage[[i]] <- temp_storage
  }

  output <- do.call("rbind",output_storage)
  # SE_prediction is the variance of the actually observed event numbers: it includes both the curve fitting error, and the binomial error associated with the random 'sampling'
  # It is based on a beta-binomial model with size corresponding to the number of patients recruited to date
  output[,"SE_Prediction"] <- sqrt(betabinomialVar(mu=output[,"Predicted_Events"]/output[,"Patients"],sigma2=(output[,"SE_Fitting"]/output[,"Patients"])^2,n=output[,"Patients"]))
  output[,"Prediction_Lower"] <- qbetabin(p=(1-CI)/2,mu=output[,"Predicted_Events"]/output[,"Patients"],sigma2=(output[,"SE_Fitting"]/output[,"Patients"])^2,n=output[,"Patients"])
  output[,"Prediction_Upper"] <- pmin(qbetabin(p=1-((1-CI)/2),mu=output[,"Predicted_Events"]/output[,"Patients"],sigma2=(output[,"SE_Fitting"]/output[,"Patients"])^2,n=output[,"Patients"]),output[,"Patients"])

# Perform conditioned event prediction if conditioning required
  if(condition){
    if(cond_Time <= 0){
      NatRisk <- N
      est_events_at_condition <- 0
    } else {
      NatRisk <- (N-(integrate(Dfunction,lower=0,upper=cond_Time,assess=cond_Time)$value+integrate(Cfunction,lower=0,upper=cond_Time,assess=cond_Time)$value))
      est_events_at_condition <- integrate(Dfunction,lower=0,upper=cond_Time,assess=cond_Time)$value
    }
    cond_Events <- max(cond_Events,0)
    if(is.null(cond_NatRisk)){
      cond_NatRisk <- NatRisk+est_events_at_condition-cond_Events
    }
    atRiskRatio <- cond_NatRisk/NatRisk

# Create conditioning output
    Conditioned_Events <- Cond_SE_Fitting <- Cond_SE_Prediction <- Cond_Prediction_Lower <- Cond_Prediction_Upper <- rep(NA,nrow(output))
    if(cond_Time < 1){
      Conditioned_Events <- output[,"Predicted_Events"]
      Cond_SE_Fitting <- output[,"SE_Fitting"]
      Cond_SE_Prediction <- output[,"SE_Prediction"]
      Cond_Prediction_Lower <- output[,"Prediction_Lower"]
      Cond_Prediction_Upper <- output[,"Prediction_Upper"]
    }else{
    # Events exactly known at time of conditioning
      Cond_Prediction_Lower[cond_Time] <- Cond_Prediction_Upper[cond_Time] <- Conditioned_Events[cond_Time] <- cond_Events
      Cond_SE_Fitting[cond_Time] <- Cond_SE_Prediction[cond_Time] <-  0

    # Calculate values at time of conditioning
      Ialphacond <- tryCatch(integrate(Ialphafunction,lower=0,upper=cond_Time,assess=cond_Time)$value,error=function(e){return(Inf)})
      Ibetacond <- tryCatch(integrate(Ibetafunction,lower=0,upper=cond_Time,assess=cond_Time)$value,error=function(e){return(Inf)})
      Econd <- est_events_at_condition/N
      NARcond <- NatRisk/N
      IalphaNARcond <- -Ialphacond-tryCatch(integrate(IalphaNARfunction,lower=0,upper=cond_Time,assess=cond_Time)$value,error=function(e){return(Inf)})
      IbetaNARcond <- -Ibetacond-tryCatch(integrate(IbetaNARfunction,lower=0,upper=cond_Time,assess=cond_Time)$value,error=function(e){return(-Inf)})

    # Conditioned number of events is observed number + differences in predictions at the time and the conditioning time, multiplied by the at-risk-ratio to reflect the relative number of patients known to be at risk
    # mu is then the proportion of the patients known to be at risk when conditioned that went on to have an event by that time
      Conditioned_Events <- cond_Events + (output[,"Predicted_Events"]-est_events_at_condition)*atRiskRatio
      if(cond_Time>1){Conditioned_Events[1:(cond_Time-1)] <- NA}
      after_cond <- (cond_Time+1):length(Cond_SE_Fitting)
      ALPHA <- ((Ialpha-Ialphacond)*NARcond-IalphaNARcond*(output[,"Predicted_Events"]/N-Econd))/(NARcond)^2
      BETA <- ((Ibeta-Ibetacond)*NARcond-IbetaNARcond*(output[,"Predicted_Events"]/N-Econd))/(NARcond)^2
      Cond_SE_Fitting[after_cond] <- cond_NatRisk*sqrt(ALPHA^2*fitted[[3]][1]+BETA^2*fitted[[3]][2]+2*ALPHA*BETA*fitted[[3]][3])[after_cond]
      mu <- (Conditioned_Events-cond_Events)/cond_NatRisk
      sigma2 <- (Cond_SE_Fitting/cond_NatRisk)^2
      Cond_SE_Prediction[after_cond] <- sqrt(betabinomialVar(mu=mu,sigma2=sigma2,n=cond_NatRisk))[after_cond]
      Cond_Prediction_Lower[after_cond] <- cond_Events + qbetabin(p=(1-CI)/2,mu=mu,sigma2=sigma2,n=cond_NatRisk)[after_cond]
      Cond_Prediction_Upper[after_cond] <- cond_Events + pmin(cond_NatRisk,qbetabin(p=1-((1-CI)/2),mu=mu,sigma2=sigma2,n=cond_NatRisk)[after_cond])
    }
    output <- cbind(output,Conditioned_Events,Cond_SE_Fitting,Cond_SE_Prediction,Cond_Prediction_Lower,Cond_Prediction_Upper)
    output[,"Conditioned_Events"] <- round(output[,"Conditioned_Events"],3)
    output[,"Cond_SE_Fitting"] <- round(output[,"Cond_SE_Fitting"],4)
    output[,"Cond_SE_Prediction"] <- round(output[,"Cond_SE_Prediction"],4)
  }

  output[,"Predicted_Events"] <- round(output[,"Predicted_Events"],3)
  output[,"SE_Fitting"] <- round(output[,"SE_Fitting"],4)
  output[,"SE_Prediction"] <- round(output[,"SE_Prediction"],4)

  outputlist <- list(ecurve=ecurve,dcurve=dcurve,rcurve=rcurve,CI=CI,Fitted=fitted,Summary=output)
  # Return condensed trajectory output
  return(outputlist)
}

#####################################################################################################################
#'Event prediction using a lifetable and a recruitment RCurve
#'
#' This is a function to perform event prediction\cr
#' It uses the fit_KM function to perform non-linear regression of Weibull and log-normal curves to the provided survival data.\cr
#' It creates an event Curve object from this, and combines it with a recruitment RCurve and an optional dropout(censoring) Curve.\cr
#' Using the same numerical integration approach as nph_curve_trajectories it performs an unconditional event prediction.\cr
#' If a conditioning time, event number (and preferably number at risk) are provided, a conditional event prediction is also calculated.\cr
#' @param KMcurve The dataframe object containing the survival data
#' @param Survival The column name for the survival function (i.e. the probabilities). Default is "Survival"
#' @param Time The column name for the times. Default is "Time"
#'     Alternatively accepts "Weibull" or "Lognormal" to force the type.
#' @param weighting Boolean for whether to use weighting. Default=TRUE as it greatly improves curve fitting.
#' @param Weights Name of Weights column. Default="Weights". Optional if weighting=FALSE. Recommended to use number at risk or remaining.
#' @param Weight_power Power to raise the weights to. Useful in large trials to give added weight to later points where numbers may still be high. Default=1 (Use weights as specified).
#' @param rcurve Observed and/or expected recruitment distribution as an RCurve object. This should typically be of PieceR type (piecewise linear recruitment).
#' @param max_time Maximum time to predict events up to.
#' @param dcurve Dropout/censoring distribution as a Curve object. This is Blank() by default, i.e. no dropout.
#' @param type Type of event curve to fit.Default is "Automatic", fitting both Weibull and Log-normal curves.
#' @param startbeta Starting value for the Weibull beta (shape) parameter to be used in the non-linear regression. Default=1 (exponential).
#' @param startsigma Starting value for the Lognormal sigma (sd) parameter to be used in the non-linear regression. Default=1.
#' @param condition Boolean whether to also do a conditional event prediction. Default=FALSE
#'     Note that If all conditioning options are left as defaults, conditioned calculation will equal the unconditional one.
#' @param cond_Events Number of events to condition on. Default=0. Optional unless condition=TRUE.
#' @param cond_NatRisk Number of patients at risk to condition on. By default, the program will estimate the number at risk assuming no censoring.
#'        It is highly recommended to specify this if conditioning. Default=NULL(takes value of N - cond_Events). Optional unless condition=TRUE.
#' @param cond_Time Time, in months, to condition on. Default=0. Optional unless condition=TRUE.
#' @param units Units that the KM-curve is specified in. Accepts "Days", "Months". Default="Days".
#' @param discountHR Hazard ratio for discounting events e.g. used to predict adjudicated events from unadjudicated data where patients remain 'at risk' after an event is adjudicated not to have occurred.
#'     Values below 1 indicate fewer events will occur than predicted by the curve-fitting.
#'     Note that changing this argument is only allowed if type="Weibull" since log-normal curves are not compatible with proportional hazards.
#'     Default=1 (No discounting)
#' @return Returns a list object with the fitted ecurve, the dcurve, the rcurve, the fitting details, and a summary table with one row per month up to max_time containing the following columns:
#' \itemize{
#'  \item{"Time"}{Time of assessment.}
#'  \item{"Patients"}{Number of patients recruited by the assessment time.}
#'  \item{"Predicted_Events"}{Number of events unconditionally predicted at the assessment time.}
#'  \item{"Conditioned_Events"}{Number of events unconditionally predicted at the assessment time (Column present only if conditioning specified).}
#' }
#' @author James Bell
#' @examples recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' example_data_short <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),
#' rcurve=recruit, assess=10,iterations=1,seed=12345,detailed_output=TRUE)
#'
#' library(survival)
#'
#' temp1 <- summary(survfit(Surv(example_data_short[,"Time"],1-example_data_short[,"Censored"])~ 1,
#' error="greenwood"))
#' out1 <- cbind(temp1$time,temp1$n.risk,temp1$surv,temp1$std.err)
#' out1 <- rbind(c(0,out1[1,2],1,0),out1)
#' colnames(out1) <- c("Time","NAR","Survival","Std.Err")
#' x1 <- ceiling(max(out1[,"Time"]))
#' example_lifetable <- out1[findInterval(0:x1,out1[,"Time"]),]
#' example_lifetable[,"Time"] <- 0:x1
#'
#' event_prediction_KM(KMcurve=example_lifetable, weighting=TRUE, Weights="NAR", rcurve=recruit,
#' max_time=60, type="automatic", condition=TRUE, cond_Events=49, cond_NatRisk=451, cond_Time=10, 
#' units="Months")
#' @export
event_prediction_KM <- function(KMcurve, Survival="Survival", Time="Time", weighting=FALSE, Weights="Weights",Weight_power=1,rcurve,max_time=100, dcurve=Blank(), type=c("automatic","Weibull","Lognormal"), startbeta=1, startsigma=1,
                                condition=FALSE,cond_Events=0,cond_NatRisk=NULL,cond_Time=0,units=c("Days","Months"),discountHR=1){
# Input checks to ensure valid inputs
# Firstly check for missing variables
  if(missing(KMcurve))stop("Please supply the name of the dataset with the Kaplan Meier curve to be modelled using the 'KMcurve' argument")
  if(missing(rcurve))stop("Please specify the recruitment distribution using the 'rcurve' argument. Please note that this should be an RCurve object; Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")

  #Multiple choices
  type <- match.arg(type)
  units <- match.arg(units)

# Secondly, check arguments are all of the correct type
  #Characters/columns
  if(!is.character(Survival))stop("'Survival' argument must be a string")
  if(!is.character(Time))stop("'Time' argument must be a string")
  if(!(Time %in% colnames(KMcurve)))stop("Error: Specified name of Time column does not exist.")
  if(!(Survival %in% colnames(KMcurve)))stop("Error: Specified name of Survival column does not exist.")

  #Weighting
  if(!is.logical(weighting))stop("Error: 'weighting' argument must be boolean: default=TRUE (weighting is to be used).")
  if(weighting){
    if(!(Weights %in% colnames(KMcurve)))stop("Error: Weighting specified, but name of Weights column does not exist.")
    if(!is.numeric(Weight_power))stop("Error: Weight_power argument must be numerical. default=1.")
  }

  #Numeric
  if(!is.numeric(startbeta) || startbeta < 0 )stop("Error: 'startbeta' argument must be a positive number")
  if(!is.numeric(startsigma) || startsigma < 0)stop("Error: 'startsigma' argument must be a positive number")
  if(!is.numeric(max_time) || max_time%%1!=0 || max_time < 1 ) stop("Please specify a positive integer for the maximum time to predict events up to, using the 'max_time' argument")
  if(!is.numeric(discountHR) || discountHR < 0 ) stop("Error: Please specify a positive number for the 'discountHR'. Default = 1 (no discounting)")
  if(!(discountHR ==1) && !(type=="Weibull"))stop("Error: 'DiscountHR' may only be set when 'type='Weibull'':Discounting of events is only possible with Weibull curves as it requires a proportional hazards assumption.")

  #Curve
  if(class(rcurve)[1]!= "RCurve") stop("Argument 'rcurve' must be an RCurve object in order to define the combined observed/predicted recruitment distribution. For event prediction, this will typically be a piecewise distribution (create using PieceR constructor function)")
  if(class(dcurve)[1]!= "Curve") stop("Argument 'dcurve' must be a Curve object in order to define the censoring distribution. By default, no censoring is specified.")

  #Conditioning
  if(!is.logical(condition))stop("Error: 'condition' argument must be boolean: default=FALSE (no conditioning).")
  if(condition){
    if(!is.numeric(cond_Events) || cond_Events%%1!=0 || cond_Events < 0 || cond_Events > getN(rcurve))stop("Error: 'cond_Events' argument must be a positive integer no greater than the number of patients. default=0 (0 events).")
    if(!is.numeric(cond_Time) || cond_Time%%1!=0 || cond_Time < 0 || cond_Time >= max_time)stop("Error: 'cond_Time' argument must be a positive integer smaller than max_time. default=0 (conditioning from start time).")
    if(!is.null(cond_NatRisk) && (!is.numeric(cond_NatRisk) || cond_NatRisk%%1!=0 || cond_NatRisk < 0 || cond_NatRisk > getN(rcurve)))stop("Error: 'cond_NatRisk' argument must be NULL or a positive integer no greater than the number of patients. default=NULL (conditioning on expected number at risk).")
  }

### Warnings and notes

  if(condition){
    if(cond_Events == 0){
      warning("Conditional event prediction has been requested but no observed events have been specified - results may not be reasonable.\n")
    }
    if(is.null(cond_NatRisk)){
      message("NOTE: No number at risk at the time of conditioning has been specified. The program by default estimates this as the predicted number minus the difference between the observed and predicted events at the time of conditioning. This will be accurate only if the predicted censoring distribution is also accurate. It is strongly recommended to supply the observed number at risk at the time of conditioning.\n")
    }
    if(cond_Time == 0){
      warning("No conditioning time has been specified. The program by default will condition from the start of the trial (time 0).\n")
    }
  }

  #Fit selected distribution(s) to KM curve.

  fitted <- fit_KM(KMcurve=KMcurve,Survival=Survival,Time=Time,weighting=weighting,Weights=Weights,Weight_power=Weight_power,startbeta,startsigma,type=type)

  if(units=="Days"){
    conversion_factor <- 365.25/12
  }else if(units=="Months"){
    conversion_factor <- 1
  } else {
    warning("Unrecognised units type specified:",units, "\nPlease specify either \"Days\" or \"Months\". Assuming units of days by default.")
    conversion_factor <- 365.25/12
  }

  if(fitted[[1]]=="Weibull"){
    discount <- discountHR^(1/fitted[[2]][2])
    fitted[[2]][1] <- fitted[[2]][1]/(conversion_factor*discount)
    fitted[[3]][1] <- fitted[[3]][1]/((conversion_factor*discount)^2)
    fitted[[3]][3] <- fitted[[3]][3]/(conversion_factor*discount)
    ecurve <- Weibull(alpha=fitted[[2]][1],beta=fitted[[2]][2])
  } else if(fitted[[1]]=="Lognormal") {
    if(discountHR !=1){warning("A discount HR for event occurrence has been specified but a log-normal curve has been fitted - no discounting will be applied due to incompatibility.\n")}
    fitted[[2]][1] <- fitted[[2]][1]-log(conversion_factor)
    ecurve <- Lognormal(mu=fitted[[2]][1],sigma=fitted[[2]][2])
  }

  N <- getN(rcurve)

# Self-modifying R code! These 2 statements create the two required functions tailored specifically to the combinations of events, censoring and recruitment
# found within each arm to calculate the events and number at risk. Relatively fast as only needs to do this once per run, and allows code to adapt to all sorts of input.

    eval(parse(text=paste(
      "Dfunction <- function(tim,assess){",
        N,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(dcurve,q="tim"),")*",getPDFfunction(ecurve,x="tim"),
      "}"
    )))

# Write function for numbers at risk in the control arm
    eval(parse(text=paste(
      "Cfunction <- function(tim,assess){",
        N,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(ecurve,q="tim"),")*",getPDFfunction(dcurve,x="tim"),
      "}"
    )))

# Set up list to store output, and a 'single line' temporary storage vector
  output_storage <- vector("list", max_time)
  temp_storage <- data.frame(t(rep(NA,4)))
  colnames(temp_storage) <- c(
    "Time",
    "Patients",
    "Predicted_Events",
    "Conditioned_Events"
   )
# Perform conditioned event prediction if conditioning required
  if(condition){
    if(cond_Time == 0){
      atRiskRatio <- 1
      est_events_at_condition <- 0
    } else {
      est_events_at_condition <- integrate(Dfunction,lower=0,upper=cond_Time,assess=cond_Time)$value
      NatRisk <- (N-(integrate(Dfunction,lower=0,upper=cond_Time,assess=cond_Time)$value+integrate(Cfunction,lower=0,upper=cond_Time,assess=cond_Time)$value))
      if(is.null(cond_NatRisk)){
        cond_NatRisk <- NatRisk+est_events_at_condition-cond_Events
      }
      atRiskRatio <- cond_NatRisk/NatRisk
    }

# Create output if conditioning required
    for(i in 1:max_time){
      temp_storage$Time <- i
      temp_storage$Patients <- getPatients(rcurve,i)
      temp_storage$Predicted_Events <- integrate(Dfunction,lower=0,upper=i,assess=i)$value
      if(i >= cond_Time){
        temp_storage$Conditioned_Events <- cond_Events + (temp_storage$Predicted_Events-est_events_at_condition)*atRiskRatio
      } else {
        temp_storage$Conditioned_Events <- NA
      }
      output_storage[[i]] <- temp_storage
    }
# Create output if conditioning not required
  } else{
    for(i in 1:max_time){
      temp_storage$Time <- i
      temp_storage$Patients <- getPatients(rcurve,i)
      temp_storage$Predicted_Events <- integrate(Dfunction,lower=0,upper=i,assess=i)$value
      temp_storage$Conditioned_Events <- NA
      output_storage[[i]] <- temp_storage
    }
  }
  output <- do.call("rbind",output_storage)

  output[,"Predicted_Events"] <- round(output[,"Predicted_Events"],3)
  output[,"Conditioned_Events"] <- round(output[,"Conditioned_Events"],3)

  outputlist <- list(ecurve=ecurve,dcurve=dcurve,rcurve=rcurve,Fitted=fitted,Summary=output)
  # Return condensed trajectory output
  return(outputlist)
}


