if(getRversion() >= "2.15.1")  utils::globalVariables(c("Dfunction","Cfunction"))

#####################################################################################################################
#'Event prediction using a lifetable and a recruitment RCurve
#'
#' This is a function to perform event prediction using summary-level data. As of 1.4.0 this function is deprecated as event_prediction gives improved accuracy.\cr
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
#' @param type Type of event curve to fit. Default is "Automatic", fitting both Weibull and Log-normal curves.
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
event_prediction_KM <- function(KMcurve, Survival="Survival", Time="Time", weighting=FALSE, Weights="Weights",Weight_power=1,rcurve,max_time=100, dcurve=Blank(), type=c("automatic","Weibull","Lognormal","Exponential"), startbeta=1, startsigma=1,
                                condition=FALSE,cond_Events=0,cond_NatRisk=NULL,cond_Time=0,units=c("Days","Months"),discountHR=1){
.Deprecated(new="event_prediction", package="gestate",msg="event_prediction_KM is now deprecated: The event_prediction function performs the same role more accurately and allows creation of prediction intervals. Patient-level data should be available for event prediction, so no need is foreseen for event prediction using only summary-level information.")
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
  } else if(fitted[[1]]=="Exponential") {
    fitted[[2]][1] <- fitted[[2]][1]*conversion_factor*discount
    ecurve <- Exponential(lambda=fitted[[2]][1])
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









