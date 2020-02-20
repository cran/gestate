## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gestate)
library(survival)
# Simulate a single data set based upon a Weibull curve with shape 0.8 and scale 100:
# Distributions
example_distribution <- Weibull(alpha=100,beta=0.8)
example_recruitment  <- LinearR(rlength=24,Nactive=600,Ncontrol=0)
# Create the full trajectory
maxtime <- 100
example_data_long <- simulate_trials(active_ecurve=example_distribution,control_ecurve=example_distribution,rcurve=example_recruitment,assess=maxtime,iterations=1,seed=1234567,detailed_output=TRUE)
# Shrink the trajectory to an early time point at which event prediction will occur
example_data_short <- set_assess_time(example_data_long,18)

######################
# Create life-table with 7-day 'sampling'
# Note that this is to create an example lifetable and would not normally be done as it is typically better to use patient-level data if available
temp1  <- summary(survfit(Surv(example_data_short[,"Time"],1-example_data_short[,"Censored"])~ 1,error="greenwood"))
out1 <- cbind(temp1$time,temp1$n.risk,temp1$surv,temp1$std.err)
out1 <- rbind(c(0,out1[1,2],1,0),out1)
colnames(out1) <- c("Time","NAR","Survival","Std.Err")
interval <- 1
x1 <- ceiling(max(out1[,"Time"]/interval))*interval
times1 <- seq(from = 0, to = x1,by=interval)
keys1 <- findInterval(times1,out1[,"Time"])
example_lifetable <- out1[keys1,]
example_lifetable[,"Time"] <- times1

  

## -----------------------------------------------------------------------------
# Curve-fit the example lifetable in automatic mode
fit1 <- fit_KM(KMcurve=example_lifetable,Survival="Survival",Time="Time",Weights="NAR",type="automatic")
fit1

## -----------------------------------------------------------------------------
# Curve-fit the example lifetable in automatic mode
fit2 <- fit_tte_data(data=example_data_short,Time="Time",Event="Censored",censoringOne=TRUE,type="automatic")
fit2

## -----------------------------------------------------------------------------
# Create an RCurve incorporating both observed and predicted recruitment with arbitrary randomisation ratio of 1
recruit <- PieceR(matrix(c(rep(1,12),10,15,25,30,45,60,55,50,65,60,55,30),ncol=2),1)

## -----------------------------------------------------------------------------
# Create the full length trial showing what 'actually' happens
trial_long <- simulate_trials(active_ecurve=Weibull(50,0.8),control_ecurve=Weibull(50,0.8),rcurve=recruit,fix_events=200,iterations=1,seed=12345,detailed_output=TRUE)
# Create a shortened version with data for use in the event prediction
trial_short <- set_assess_time(data=trial_long,time=10,detailed_output = FALSE)
# Create the trajectory of events
maxtime <- max(ceiling(trial_long[,"Assess"]))
events <- rep(NA,maxtime)
for (i in 1:maxtime){
  events[i] <- sum(1-set_assess_time(trial_long,i)[,"Censored"])
}


## -----------------------------------------------------------------------------
predictions <- event_prediction(data=trial_short, Event="Censored", censoringOne=TRUE, type="Weibull", rcurve=recruit,max_time=60, cond_Events=49, cond_NatRisk=451, cond_Time=10, units="Months")

## ---- fig.show='hold',fig.height = 5, fig.width = 7, fig.align = "center"-----
# Plot observed events and conditional predictions
plot_ep(predictions,trajectory="conditional",which_PI="conditional",max_time=40,observed=events,target=200,max_E=200)

## -----------------------------------------------------------------------------
trial_less_short <- set_assess_time(data=trial_long,time=14,detailed_output = FALSE)
predictions2 <- event_prediction(data=trial_less_short, Event="Censored",censoringOne=TRUE, type="Weibull", rcurve=recruit,max_time=60, cond_Events=101,cond_NatRisk=399,cond_Time=14, units="Months")

trial_mature <- set_assess_time(data=trial_long,time=18,detailed_output = FALSE)
predictions3 <- event_prediction(data=trial_mature, Event="Censored",censoringOne=TRUE, type="Weibull", rcurve=recruit,max_time=60, cond_Events=148,cond_NatRisk=352,cond_Time=18, units="Months")

## ---- echo=FALSE, fig.show='hold',fig.height = 5, fig.width = 7, fig.align = "center"----
plot_ep(predictions2,trajectory="conditional",which_PI="conditional",max_time=40,observed=events,target=200,max_E=200)
plot_ep(predictions3,trajectory="conditional",which_PI="conditional",max_time=40,observed=events,target=200,max_E=200)


## ---- echo=FALSE, fig.show='hold',fig.height = 5, fig.width = 7, fig.align = "center"----
plot_km_fit(fit=predictions,data=trial_short, Event="Censored",censoringOne=TRUE,main="Kaplan Meier Curve Fit Plot: 10 Months",maxT=45,xlim=c(0,40))
plot_km_fit(fit=predictions2,data=trial_less_short, Event="Censored",censoringOne=TRUE,main="Kaplan Meier Curve Fit Plot: 14 Months",maxT=45,xlim=c(0,40))
plot_km_fit(fit=predictions3,data=trial_mature, Event="Censored",censoringOne=TRUE,main="Kaplan Meier Curve Fit Plot: 18 Months",maxT=45,xlim=c(0,40))


