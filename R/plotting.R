#'Plot output from nph_traj
#'
#' This function plots the output from nph_traj().\cr
#' By default, it produces 6 plots:
#' \itemize{
#'  \item{"KM plot"}{ Kaplan Meier plot for events. This is in patient time.}
#'  \item{"Censoring plot"}{ Plot of CDFs for censoring functions. This is in patient time.}
#'  \item{"Recruitment plot"}{ Number of patients expected to have been recruited over time. This is in trial time.}
#'  \item{"Event plot"}{ Total number of events expected to occur over time. This is in trial time.}
#'  \item{"log(HR) plot"}{ Expected log(HR), with expected confidence interval, over time. This is in trial time.}
#'  \item{"Power plot"}{ Expected power over time for various methods. This is in trial time.}
#' }
#' Plots may be omitted via arguments.\cr
#' All calculated powers automatically plotted unless specified otherwise.\cr
#' @param data Full output list from nph_traj()
#' @param KM Boolean to include KM plot (Default = TRUE)
#' @param censor Boolean to include censoring plot (Default = TRUE)
#' @param recruitment Boolean to include recruitment plot (Default = TRUE)
#' @param events Boolean to include events plot (Default = TRUE)
#' @param logHR Boolean to include log(HR) plot (Default = TRUE)
#' @param power Boolean to include power plot (Default = TRUE)
#' @param include_frontier Boolean to include frontier power curve in power plot (Default = TRUE)
#' @param include_RMST Boolean to include RMST power curve in power plot if available (Default = TRUE)
#' @param include_landmark Boolean to include landmark power curve in power plot if available (Default = TRUE)
#' @param alpha1 One-sided alpha to use for estimation of log(HR) confidence intervals (Default = 0.025)
#' @param legend_position String with any of "top_left","top_right" or "bottom_right", corresponding to legend position in power plot. Default is "top_left".
#' @return Returns NULL
#' @author James Bell
#' @examples 
#' trial <- nph_traj(Weibull(100,1),Weibull(70,1),rcurve=LinearR(12,100,100),RMST=20,
#'   landmark=20,max_assessment=30)
#'
#' plot_npht(trial)
#' plot_npht(data=trial,KM=FALSE,censor=FALSE,recruitment=FALSE)
#' plot_npht(data=trial,KM=FALSE,censor=FALSE,recruitment=FALSE,events=FALSE,logHR=FALSE,
#' include_frontier=FALSE, include_RMST=FALSE,include_landmark=FALSE,legend_position="top_right")
#' 
#' @import graphics
#' @export
plot_npht <- function(data,KM=TRUE,censor=TRUE,recruitment=TRUE,events=TRUE,logHR=TRUE,power=TRUE,include_frontier=TRUE,include_RMST=TRUE,include_landmark=TRUE,alpha1=0.025,legend_position=c("top_left","top_right","bottom_right")){
  if(missing(data))stop("Error: Please use the 'data' argument to specify the name of a table created by nph_curve_trajectories.")
  if(!("Summary" %in% names(data)))stop("Error: Cannot find Summary table")
  #Booleans
  if(!is.logical(KM)) stop("Error:'KM' must be a boolean argument.")
  if(!is.logical(censor)) stop("Error: 'censor' must be a boolean argument.")
  if(!is.logical(recruitment)) stop("Error: 'recruitment' must be a boolean argument.")
  if(!is.logical(events)) stop("Error:'events' must be a boolean argument.")
  if(!is.logical(logHR)) stop("Error:'logHR' must be a boolean argument.")
  if(!is.logical(power)) stop("Error:'Power' must be a boolean argument.")
  if(!is.logical(include_frontier)) stop("Error: 'include_frontier' must be a boolean argument.")
  if(!is.logical(include_RMST)) stop("Error: 'include_RMST' must be a boolean argument.")
  if(!is.logical(include_landmark)) stop("Error: 'include_landmark' must be a boolean argument.")

  #Multiple choice
  legend_position <- match.arg(legend_position)
  #Numeric
  if(alpha1 <= 0 | alpha1 > 0.5 ) stop("Please specify a positive one-sided alpha less than 0.5 using the 'alpha1' argument.")

  # Set graphical display at 2 rows of 3 graphs
  plotnum <- KM+censor+recruitment+events+logHR+power
  if(plotnum <= 1){par(mfrow=c(1,1))
  } else if(plotnum == 2){par(mfrow=c(1,2))
  } else if(plotnum == 3){par(mfrow=c(1,3))
  } else if(plotnum == 4){par(mfrow=c(2,2))
  } else {par(mfrow=c(2,3))}

  dat <- data$Summary
  maxtime <- max(dat[,"Time"])
  if(KM){
    plotSF(data$control_ecurve,overlay=FALSE,maxT=maxtime,xlab="Patient Time (months)",main="KM Curves for \nActive (red) and Control (black) Arms")
    plotSF(data$active_ecurve,overlay=TRUE,maxT=maxtime,col=2)
  }
  if(censor){
    plotCDF(data$control_dcurve,overlay=FALSE,maxT=maxtime,xlab="Patient Time (months)",main="Censoring Curves for \nActive (red) and Control (black) Arms")
    plotCDF(data$active_dcurve,overlay=TRUE,maxT=maxtime,col=2)
  }
  if(recruitment){
    plot(dat[,"Time"],dat[,"Patients"],type="l",xlab="Trial Time (months)",ylab="Patients",main="Expected Recruited Patients over Time")
  }
  if(events){
  plot(dat[,"Time"],dat[,"Events_Total"],type="l",xlab="Trial Time (months)",ylab="Observed Events",main="Expected Events over Time")
  }
  if(logHR){
    plot(dat[,"Time"],dat[,"LogHR"],type="l",xlab="Trial Time (months)",ylab="log(HR)",main="Expected log(HR) over Time",ylim=c(min(-1,dat[,"LogHR"]),max(1,dat[,"LogHR"])))
    lines(dat[,"Time"],dat[,"LogHR"]+qnorm(1-alpha1)*dat[,"LogHR_SE"],type="l",col=2,lty=3)
    lines(dat[,"Time"],dat[,"LogHR"]-qnorm(1-alpha1)*dat[,"LogHR_SE"],type="l",col=2,lty=3)
    lines(dat[,"Time"],dat[,"Time"]*0,type="l",col=3)
  }
  if(power){
    #Set legend locations
    if(legend_position=="top_left"){
      xpos <- 0
    } else {xpos <- 0.7}
    if(legend_position=="top_right" || legend_position=="top_left" ){
      ypos <- 1
    } else {ypos <- 0.2}

    plot(dat[,"Time"],dat[,"Schoenfeld_Power"],type="l",xlab="Trial Time (months)",ylab="Power",main="Expected Power over Time",ylim=c(0,1))
    legend(x=xpos*maxtime,y=ypos,legend=c("Schoenfeld"),cex=0.7,lwd=c(1),lty=c(1),col=c(1),bty = "n")
    ypos <- ypos-0.05

    if(include_frontier){
      lines(dat[,"Time"],dat[,"Frontier_Power"],type="l",col=2)
      legend(x=xpos*maxtime,y=ypos,legend=c("Frontier"),cex=0.7,lwd=c(1),lty=c(1),col=c(2),bty = "n")
      ypos <- ypos-0.05
    }
    if(include_RMST){
      if("RMST_Power" %in% names(dat)){
        lines(dat[,"Time"],dat[,"RMST_Power"],type="l",col=3)
        legend(x=xpos*maxtime,y=ypos,legend=c(paste("RMST",dat[1,"RMST_Restrict"])),cex=0.7,lwd=c(1),lty=c(1),col=c(3),bty = "n")
        ypos <- ypos-0.05
      }
    }
    if(include_landmark){
      if("LM_Power" %in% names(dat)){
        lines(dat[,"Time"],dat[,"LM_Power"],type="l",col=4)
        legend(x=xpos*maxtime,y=ypos,legend=c(paste("Landmark",dat[1,"LM_Time"])),cex=0.7,lwd=c(1),lty=c(1),col=c(4),bty = "n")
        ypos <- ypos-0.05
      }
    }
  }
}

#'Plot event prediction output
#'
#' This function plots the output from event_prediction().\cr
#' By default, produces a plot of predicted events over time, with confidence intervals if available.\cr
#' Alternatively, produces a plot with the fitting CI over time with same percentage as prediction interval.\cr
#' By default, both conditional and unconditional trajectories are plotted (if conditional event prediction is available).\cr
#' Options are available to customise inclusion.\cr
#' @param data Full output list from event_prediction().
#' @param trajectory String, choice of "both","conditional","unconditional", for which trajectory to plot. (Default="both")
#' @param which_PI String, choice of "both","conditional","unconditional","none", for which PIs to plot. (Default="both")
#' @param prediction_fitting String, choice of "prediction" or "fitting" (Default = "prediction"). Determines the nature of the intervals; PIs are relevant for prediction of the observation of future trajectories (sample level), whereas fitting CIs concern the interval for the mean trajectory itself (population level).
#' @param observed Optional trajectory of observed event numbers. If vector specified, will plot values at integer times starting from 1. If 2-column matrix specified, will take x-values from column 1 and y-values from column 2. (Default=NULL; not plotted).
#' @param target Optional target number of events to plot. (Default=NULL; not plotted)
#' @param max_time Optional maximum time to plot up to if you do not want to plot full trajectory. (Default=NULL; maximum time determined by input data)
#' @param max_E Optional maximum number of events to plot up to. (Default=NULL; maximum event number is the number of patients)
#' @param legend_position String with any of "top_left", or "bottom_right", corresponding to legend position in power plot. (Default="top_left").
#' @param no_legend Boolean to turn off legend. Default is FALSE; legend shown.
#' @param ... Additional graphical parameters.
#' @return Returns NULL
#' @author James Bell
#' @examples 
#' recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' trial_long <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),
#' rcurve=recruit,fix_events=200, iterations=1,seed=12345,detailed_output=TRUE)
#' trial_short <- set_assess_time(data=trial_long,time=10,detailed_output = FALSE)
#'
#' maxtime <- max(ceiling(trial_long[,"Assess"]))
#' events <- rep(NA,maxtime)
#' for (i in 1:maxtime){events[i] <- sum(1-set_assess_time(trial_long,i)[,"Censored"])}
#'
#' predictions <- event_prediction(data=trial_short, Event="Censored", censoringOne=TRUE, 
#' type="Weibull", rcurve=recruit, max_time=60, cond_Events=49, cond_NatRisk=451, 
#' cond_Time=10, units="Months")
#'
#' plot_ep(predictions,trajectory="conditional",which_PI="conditional",max_time=40,observed=events,
#' target=200,max_E=200)
#'
#' plot_ep(predictions,trajectory="unconditional",which_PI="unconditional",max_time=40,
#' observed=events,target=200,max_E=200)
#'
#' plot_ep(predictions,trajectory="conditional",which_PI="none",observed=events[1:10],max_time=20,
#' max_E=150)
#' 
#' @export
plot_ep <- function(data,trajectory=c("both","conditional","unconditional"),which_PI=c("both","conditional","unconditional","none"),prediction_fitting=c("prediction","fitting"),observed=NULL,target=NULL,max_time=NULL,max_E=NULL,legend_position=c("top_left","bottom_right"),no_legend=FALSE,...){
  if(missing(data))stop("Error: Please use the 'data' argument to specify the name of a table created by nph_traj().")
  if(!("Summary" %in% names(data)))stop("Error: Cannot find Summary table")
  if(!("rcurve" %in% names(data)))stop("Error: Cannot find the Rcurve")
  if(!("Predicted_Events" %in% names(data$Summary)))stop("Error: Cannot find Predicted_Events column - incorrect Summary table format")

  #Multiple choice
  trajectory <- match.arg(trajectory)
  legend_position <- match.arg(legend_position)
  which_PI <- match.arg(which_PI)
  prediction_fitting <- match.arg(prediction_fitting)

 # observed
  if(!is.null(max_time) && (!is.numeric(max_time) || max_time < 1) ) stop("A maximum time to plot has been specified that is not a positive number of at least 1.")
  if(!is.null(max_E) && (!is.numeric(max_E) || max_E < 1) ) stop("A maximum event number to plot has been specified that is not a positive number of at least 1.")
  if(!is.null(target) && (target%%1!=0 | target < 1 )) stop("A target event number has been specified, but it is not a positive integer.")
  if(!is.null(observed) && is.vector(observed)) {observed <- cbind(1:length(observed),observed)}
  if(!is.null(observed) && (ncol(observed) != 2 || !is.numeric(observed) ))stop("observed must be entirely numeric and either a vector or a 2-column matrix/data frame")

  dat <- data$Summary
  maxtime <- ceiling(min(max(dat[,"Time"]),max_time))

  plotU <- FALSE
  plotC <- FALSE
  plotUPI <- FALSE
  plotCPI <- FALSE

  plotU <- trajectory != "conditional"
  plotC <- (trajectory != "unconditional" && "Conditioned_Events" %in% names(dat))
  plotUPI <- (plotU && which_PI != "conditional" && which_PI != "none" && "SE_Prediction" %in% names(dat))
  plotCPI <- (plotC && which_PI != "unconditional" && which_PI != "none" && "Cond_SE_Prediction" %in% names(dat))

  maxE <- ceiling(min(getN(data$rcurve),max_E))

  #Set legend locations
  if(legend_position=="top_left"){
    xpos <- 0
    ypos <- 1
  } else {
    xpos <- 0.7
    ypos <- 0.3
  }
  percentage <- 100*data$PI
  alpha1 <- (1-data$PI)/2
  times <- dat[,"Time"]

  plot(NULL, xlim=c(0,maxtime), ylim=c(0,maxE), xlab="Trial Time (months)", ylab="Events",...)
  if(plotU){
    lines(times,dat[,"Predicted_Events"],type="l",col=4)
    if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=c("Unconditional Prediction"),cex=0.7,lwd=1,lty=1,col=4,bty = "n")}
    ypos <- ypos-0.05
  }
  if(plotUPI){
    if(prediction_fitting == "prediction"){
      lines(times,dat[,"Prediction_Upper"],type="l",col=4,lty=2)
      lines(times,dat[,"Prediction_Lower"],type="l",col=4,lty=2)
      if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=paste("Unconditional ",percentage,"% PI",sep=""),cex=0.7,lwd=1,lty=2,col=4,bty = "n")}
    } else{
      lines(times,dat[,"Predicted_Events"]+qnorm(1-alpha1)*dat[,"SE_Fitting"],type="l",col=4,lty=2)
      lines(times,dat[,"Predicted_Events"]-qnorm(1-alpha1)*dat[,"SE_Fitting"],type="l",col=4,lty=2)
      if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=paste("Unconditional ",percentage,"% Fitting CI",sep=""),cex=0.7,lwd=1,lty=2,col=4,bty = "n")}
    }
    ypos <- ypos-0.05
  }
  if(plotC){
    lines(times,dat[,"Conditioned_Events"],type="l",col=2)
    if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=c("Conditional Prediction"),cex=0.7,lwd=1,lty=1,col=2,bty = "n")}
    ypos <- ypos-0.05
  }
  if(plotCPI){
    if(prediction_fitting == "prediction"){
      lines(times,dat[,"Cond_Prediction_Upper"],type="l",col=2,lty=2)
      lines(times,dat[,"Cond_Prediction_Lower"],type="l",col=2,lty=2)
      if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=paste("Conditional ",percentage,"% PI",sep=""),cex=0.7,lwd=1,lty=2,col=2,bty = "n")}
    } else{
      lines(times,dat[,"Conditioned_Events"]+qnorm(1-alpha1)*dat[,"Cond_SE_Fitting"],type="l",col=2,lty=2)
      lines(times,dat[,"Conditioned_Events"]-qnorm(1-alpha1)*dat[,"Cond_SE_Fitting"],type="l",col=2,lty=2)
      if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=paste("Conditional ",percentage,"% Fitting CI",sep=""),cex=0.7,lwd=1,lty=2,col=2,bty = "n")}
    }
    ypos <- ypos-0.05
  }
  if(!is.null(observed)){
    lines(observed[,1],observed[,2],lwd=2,lty=1,col=1)
    if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend="Observed Events",cex=0.7,lwd=2,lty=1,col=1,bty = "n")}
    ypos <- ypos-0.05
  }
  if(!is.null(target)){
    lines(c(-10,times),rep(target,length(times)+1),lwd=0.5,lty=3,col=1)
    if(!no_legend){legend(x=xpos*maxtime,y=ypos*maxE,legend=paste("Target: ",target," Events",sep=""),cex=0.7,lwd=0.5,lty=3,col=1,bty = "n")}
    ypos <- ypos-0.05
  }
}

#'Kaplan Meier Plot of Curve-Fit
#'
#' This function creates a Kaplan Meier plot with the fitted curve from the output of event_prediction(), fit_tte_data() or fit_KM().\cr
#' Where available, it will include fitting confidence intervals based upon the variance derived by the delta method.\cr
#' Options are available to customise inclusion.\cr
#' @param fit Full output list from event_prediction(), fit_tte_data() or fit_KM().
#' @param data Name of patient-level data set, used to generate the KM plot.
#' @param Time The column name for the times. Default is "Time"
#' @param Event The column name for the events column (i.e. the binary variable denoting events vs censorings). Default is "Event"
#' @param censoringOne Specify whether censoring is denoted in the Event column by a one (TRUE) or zero (FALSE). Default=FALSE (censorings denoted by 0, events by 1)
#' @param CI Number between 0 and 1 for the size of Kaplan Meier confidence interval to calculate. Default is 0.95 (95 percent confidence interval).
#' @param colour_CI Boolean for whether to colour the fitting confidence interval area. Default=TRUE (colour area)
#' @param maxT Maximum time to calculate point estimate and CIs up to. Default=NULL (Calculate up to last time in patient data)
#' @param xlim Standard graphical parameter for x-axis limits. Default=NULL (Plots from 0 to maximum patient time)
#' @param ylim Graphical parameter for y-axis limits. Default=c(0,1) (Plots survival function from 0 to 1)
#' @param main String for plot title. Default="Kaplan Meier Curve Fit Plot"
#' @param fit_col Colour for fitting curve Default=2 (red)
#' @param km_col Colour for km curve Default=1 (black)
#' @param area_col Colour for CI area Default="skyblue" (sky blue)
#' @param CI_col Colour for CI Default=4 (blue)
#' @param CI_lty Line type for CI Default=2 (dashed)
#' @param no_legend Boolean to turn off legend. Default is FALSE; legend shown.
#' @param legend_position String with "top_right", or "bottom_left", corresponding to legend position in power plot. (Default="bottom_left").
#' @param overlay Boolean whether to overlay on existing plot (vs start a new one). Default=FALSE
#' @param ... Additional graphical parameters.
#' @return Returns NULL
#' @author James Bell
#' @examples 
#' recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)
#' trial_long <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),
#' rcurve=recruit,fix_events=200, iterations=1,seed=12345,detailed_output=TRUE)
#' trial_short <- set_assess_time(data=trial_long,time=10,detailed_output = FALSE)
#'
#' maxtime <- max(ceiling(trial_long[,"Assess"]))
#' events <- rep(NA,maxtime)
#' for (i in 1:maxtime){events[i] <- sum(1-set_assess_time(trial_long,i)[,"Censored"])}
#'
#' predictions <- event_prediction(data=trial_short, Event="Censored", censoringOne=TRUE, 
#' type="Weibull", rcurve=recruit, max_time=60, cond_Events=49, cond_NatRisk=451, 
#' cond_Time=10, units="Months")
#'
#' plot_km_fit(fit=predictions,data=trial_short,Event="Censored",censoringOne=TRUE,maxT=70)
#' 
#' @export
plot_km_fit <- function(fit,data,Time="Time",Event="Event",censoringOne=FALSE,CI=0.95,colour_CI=TRUE,maxT=NULL,xlim=NULL,ylim=c(0,1),main="Kaplan Meier Curve Fit Plot",fit_col=2,km_col=1,area_col="skyblue",CI_col=4,CI_lty=2,no_legend=FALSE,legend_position=c("bottom_left","top_right"),overlay=FALSE,...){

# Input checks to ensure valid inputs
# Firstly check for missing variables
  if(missing(data))stop("Please supply the name of the patient-level data set using the 'data' argument")
  if(missing(fit))stop("Please supply the name of the event_prediction(), fit_tte_data() or fit_KM() output containing the curve fitting")

  #Strings / column names
  if(!is.character(Event))stop("'Event' argument must be a string")
  if(!is.character(Time))stop("'Time' argument must be a string")
  if(!(Time %in% colnames(data)))stop("Error: Specified name of Time column does not exist in 'data'.")
  if(!(Event %in% colnames(data)))stop("Error: Specified name of Event column does not exist in 'data'.")
  #Booleans
  if(!is.logical(censoringOne))stop("Error: censoringOne argument must be boolean: default=FALSE (0 is censored, 1 is event).")
  if(!is.logical(colour_CI))stop("Error: colour_CI argument must be boolean: default=TRUE (0 is censored, 1 is event).")
  #Numbers
  if(!is.numeric(CI) || !(CI > 0) || !(CI < 1)) stop("Error: Please specify a value between 0 and 1 (exclusive) for the confidence interval. Default = 0.95 (95% CI)")
  if(!is.numeric(ylim) || length(ylim) != 2 || ylim[1] >= ylim[2] || !(ylim[1] >= 0) || !(ylim[2] <= 1)) stop("Error: Please specify a two-element vector with values between 0 and 1 for ylim; the minimum and maximum values for the y axis plot. Default = c(0,1) (plot survival function from 0 to 1)")
  if(!is.null(xlim)){
    if(!is.numeric(xlim) || length(xlim) != 2 || xlim[1] >= ylim[2] || !(xlim[1] >= 0)) stop("Error: Please specify a two-element vector with values greater than 0 for the x-axis limits, or leave blank for automatic selection. Default = NULL (automatically plot full range of observed times)")
  }
  if(!is.null(maxT) && (!is.numeric(maxT) || !(maxT > 0))) stop ("Error: maxT must be a positive number or left null.")
  #Multiple choice
  legend_position <- match.arg(legend_position)

#Check that input fit data is recognised, and converge formats
  if("Fitted" %in% names(fit)){
    sourcefit <- fit$Fitted
  } else if("Curvetype" %in% names(fit)){
    sourcefit <- fit
  } else stop("Unrecognised input for argument 'fit'. 'fit' must be output from event_prediction, fit_tte_data or fit_KM")

  NAvar <- any(is.na(sourcefit$VCov))

### Functions for calculating Survival fitting variance
  VWeibull <- function(x,alpha,beta,Va,Vb,Cov){
    y <- (x/alpha)^beta
    a <- beta*y*exp(-y)/alpha
    b <- -exp(-y)*y/beta*log(y)
    V <- a^2*Va+b^2*Vb+2*a*b*Cov
  }
  VExp <- function(x,lambda,Var){
    V <- (-x*exp(-lambda*x))^2*Var
  }
  VLognormal <- function(x,mu,sigma,Vm,Vs,Cov){
    m <- 1/sqrt(2*pi)*exp(-((log(x)-mu)^2)/(2*sigma^2))/sigma
    s <- (log(x)-mu)/(sigma)
    V <- m^2*(Vm+s^2*Vs+2*s*Cov)
  }

  estBetaCIs <- function(mu, var, CI) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    lower <- suppressWarnings(qbeta((1-CI)/2,alpha,beta))
    upper <- suppressWarnings(qbeta(1-((1-CI)/2),alpha,beta))
    lower[is.nan(lower)] <- 1
    upper[is.nan(upper)] <- 1
    return(list(upper=upper,lower=lower))
  }

# Do the KM plot
  times <- data[,Time]
  events <- data[,Event]
  if(censoringOne){events = 1-events}
  maxT <- ifelse(is.null(maxT),max(times),maxT)
  xlim <- if(is.null(xlim)){
            c(0,maxT)
          } else {xlim}
  temp  <- survfit(Surv(times,events)~ 1,error="greenwood")
  x <- seq(from=0,to=maxT,by=maxT/10000)

#Create ecurve from fitting data
  if(sourcefit$Curvetype=="Exponential"){
    ecurve <- Exponential(sourcefit$Parameters[1])
    variance <- VExp(x,sourcefit$Parameters,sourcefit$VCov)
  } else if(sourcefit$Curvetype=="Weibull"){
    ecurve <- Weibull(alpha=sourcefit$Parameters[1],beta=sourcefit$Parameters[2])
    variance <- VWeibull(x=x,alpha=sourcefit$Parameters[1],beta=sourcefit$Parameters[2],Va=sourcefit$VCov[1],Vb=sourcefit$VCov[2],Cov=sourcefit$VCov[3])
  } else if(sourcefit$Curvetype=="Lognormal"){
    ecurve <- Lognormal(mu=sourcefit$Parameters[1],sigma=sourcefit$Parameters[2])
    variance <- VLognormal(x=x,mu=sourcefit$Parameters[1],sigma=sourcefit$Parameters[2],Vm=sourcefit$VCov[1],Vs=sourcefit$VCov[2],Cov=sourcefit$VCov[3])
  } else stop("Error: Curve type from curve fitting not supported by plotting function.")


  central <- 1-evaluateCDFfunction(ecurve,x)
#Estimate CIs using beta distribution with known moments
  limits <- estBetaCIs(central,variance,CI)

  if(overlay){
    lines(temp,conf.int=FALSE,col=km_col,...)
  } else {
    plot(temp,xlab="Patient Time",ylab="S(t)",xlim=xlim,ylim=ylim,conf.int=FALSE,main=main,col=km_col,...)
  }

  if(legend_position=="top_right"){
    xpos <- 0.6
    ypos <- 1
  } else {
    xpos <- 0
    ypos <- 0.2
  }
    if(!no_legend){legend(x=xpos*xlim[2],y=ypos*ylim[2],legend=c("KM curve"),cex=1,lwd=1,lty=1,col=km_col,bty = "n")}
    ypos <- ypos-0.05
    if(!no_legend){legend(x=xpos*xlim[2],y=ypos*ylim[2],legend=paste("Fitted ",getType(ecurve)," curve",sep=""),cex=1,lwd=1,lty=1,col=fit_col,bty = "n")}
    ypos <- ypos-0.05
    if(!no_legend && !NAvar){legend(x=xpos*xlim[2],y=ypos*ylim[2],legend=paste("Fitted curve ",CI*100,"% CI",sep=""),cex=1,lwd=1,lty=CI_lty,col=CI_col,bty = "n")}
    ypos <- ypos-0.05

  if(colour_CI){
    polygon(c(x,rev(x)),c(limits[[1]],rev(limits[[2]])),col=area_col,border=CI_col,lty=CI_lty)
  } else {
    lines(x,limits[[1]],col=CI_col,lty=CI_lty,...)
    lines(x,limits[[2]],col=CI_col,lty=CI_lty,...)
  }
  lines(temp,conf.int=FALSE,col=km_col,...)
  lines(x,central,col=fit_col,...)
}




