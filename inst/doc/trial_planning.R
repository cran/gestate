## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gestate)
# Define a Weibull curve with shape 0.8 and scale 200:
curve1 <- Weibull(alpha=200,beta=0.8)
# Show the specified curve
curve1

## -----------------------------------------------------------------------------
# Define a linear recruitment over 12 months with 100 patients in each arm:
rcurve1 <- LinearR(rlength=12,Nactive=100,Ncontrol=100)
rcurve1

# Define a more complex recruitment with increasing rates
pieces <- matrix(c(1,2,3,4,5,10,15,32.5),ncol=2)
#Matrix to create a piecewise recruitment distribution. First column should be duration, second column should be rate.
pieces
rcurve2 <- PieceR(recruitment=pieces,ratio=1)

# Show the specified piecewise curve
rcurve2

## -----------------------------------------------------------------------------
# Define exponential distributions for active and control arms, with HR of 0.7
controlCurve <- Exponential(lambda=0.1)
activeCurve <- Exponential(lambda=0.1*0.7)

# Define a linear recruitment of 400 patients
rcurve3 <- LinearR(rlength=12,Nactive=200,Ncontrol=200)

# Run nph_traj with default settings to calculate expected properties up to 10 months
output <- nph_traj(active_ecurve=activeCurve,control_ecurve=controlCurve,rcurve=rcurve3,max_assessment=10)


## -----------------------------------------------------------------------------
# Show output from previous section
output

## -----------------------------------------------------------------------------
# Define exponential distributions for active and control arms, with HR of 0.7
censorCurveA <- Exponential(lambda=0.001)
censorCurveC <- Exponential(lambda=0.002)

# Run nph_traj with default settings to calculate expected properties up to 10 months
output1 <- nph_traj(active_ecurve=activeCurve,control_ecurve=controlCurve,active_dcurve=censorCurveA,control_dcurve=censorCurveC,rcurve=rcurve3,max_assessment=10,detailed_output = TRUE,required_power = 0.9)
# Display output
output1


## -----------------------------------------------------------------------------
# Define exponential distributions for active and control arms, with HR of 0.7
activeCurveNPH <- Weibull(alpha=100,beta=0.8)
controlCurveNPH <-  Weibull(alpha=50,beta=1)

# Run nph_traj with default settings to calculate expected properties up to 30 months
output2 <- nph_traj(active_ecurve=activeCurveNPH,control_ecurve=controlCurveNPH,rcurve=rcurve3,max_assessment=30,RMST=20,landmark=20)
# Display output
output2


## ---- fig.show='hold',fig.height = 5, fig.width = 7, fig.align = "center"-----
# Plot output from nph_traj
plot_npht(output2)


## -----------------------------------------------------------------------------
# Simple simulation corresponding to first example of nph_traj
# Number of iterations kept unrealistically low to reduce processing time
simulation1 <- simulate_trials(active_ecurve=activeCurve,control_ecurve=controlCurve,rcurve=rcurve3,assess=10,iterations=100,seed=1234)


## -----------------------------------------------------------------------------
# Simple simulation corresponding to first example of nph_traj
# Number of iterations kept unrealistically low to reduce processing time
simulation2 <- simulate_trials(active_ecurve=activeCurveNPH,control_ecurve=controlCurveNPH,active_dcurve=censorCurveA,control_dcurve=censorCurveC,rcurve=rcurve3,fix_events=100,iterations=100,seed=12345,detailed_output = TRUE,output_type = "list")


## -----------------------------------------------------------------------------
# Define a Curve list for the active arm but keep just a single value for the control arm, to represent a predictive covariate.
activeCurveStrata <- c(activeCurve, Weibull(alpha=100,beta=0.8))

simulation3 <- simulate_trials_strata(stratum_probs=c(0.5,0.5),active_ecurve=activeCurveStrata,control_ecurve=controlCurve,rcurve=rcurve3, assess=10,iterations=100,seed=1234)


## -----------------------------------------------------------------------------
# Increase the number of events for simulation 2
simulation2a <- set_event_number(data=simulation2,events=120)

#Rerun simulation 1 with detailed output specified
simulation1a <- simulate_trials(active_ecurve=activeCurve,control_ecurve=controlCurve,rcurve=rcurve3,assess=10,iterations=100,seed=1234,detailed_output=TRUE)
# Retain the assessment time for simulation 1 at 10, but convert to a list format with simple output
simulation1b <- set_assess_time(data=simulation1a,time=10,output_type="list",detailed_output = FALSE)


## -----------------------------------------------------------------------------
# Perform log-rank test and Cox regression on simulated data
analysis1 <- analyse_sim(simulation1)
head(analysis1)


## -----------------------------------------------------------------------------
# Perform stratified log-rank test and covariate-adjusted Cox regression on simulated data
analysis3 <- analyse_sim(data=simulation3, stratum="Stratum")
head(analysis3)


## -----------------------------------------------------------------------------
# Perform log-rank test, Cox regression, landmark analysis and RMST analysis on simulated data using parallel processing.
analysis2 <- analyse_sim(data=simulation2, RMST=10, landmark=10, parallel_cores=2)
head(analysis2)


## -----------------------------------------------------------------------------
# Summarise analysis from analysis2 output
summarise_analysis(analysis2)


## -----------------------------------------------------------------------------
# Run analytic planning and select a suitable time
outputX <- nph_traj(active_ecurve=activeCurveNPH, control_ecurve=controlCurveNPH, rcurve=rcurve3, max_assessment=50, RMST=40, landmark=40)
outputX$Summary[41:50,]

## -----------------------------------------------------------------------------
# Run simulation, analysis, summary at chosen time
simulationX <- simulate_trials(active_ecurve=activeCurveNPH, control_ecurve=controlCurveNPH, rcurve=rcurve3, assess=47, seed=123456, iterations = 2000, output_type="list")
analysisX   <- analyse_sim(simulationX,RMST=40,landmark=40)
summaryX    <- summarise_analysis(analysisX)


## -----------------------------------------------------------------------------
# Run simulation, analysis, summary at chosen time
A <- outputX$Summary[47,]
B <- summaryX
Events_Active <- c(A$Events_Active,B$Events_Active)
Events_Control <- c(A$Events_Control,B$Events_Control)
Events_Total <- c(A$Events_Total,B$Events_Total)
HR <- c(A$HR,B$HR)
LogHR <- c(A$LogHR,B$LogHR)
LogHR_SE <- c(A$LogHR_SE,B$LogHR_SE)
LR_Power <- c(A$Schoenfeld_Power,B$Power_LR)
RMST <-c(A$RMST_Delta,B$RMST_Delta)
RMST_SE <-c(A$RMST_SE,B$RMST_D_SE)
RMST_Power <- c(A$RMST_Power,B$RMST_Power)
Landmark <- c(A$LM_Delta,B$LM_Delta)
Landmark_SE <- c(A$LM_D_SE,B$LM_D_SE)
Landmark_Power <- c(A$LM_Power,B$LM_Power)
collated <- rbind(Events_Active,Events_Control,Events_Total,HR,LogHR,LogHR_SE,LR_Power,RMST, RMST_SE, RMST_Power, Landmark,Landmark_SE,Landmark_Power)
colnames(collated) <- c("Analytic","Simulated")

# Display comparison table
collated

