if(getRversion() >= "2.15.1")  utils::globalVariables(c("Dfunction_control","Dfunction_active","Nfunction_control","Sfunction_control","Sfunction_active","Nfunction_active","N_Control","N_Active","E_Events_Active","E_Events_Control","HR_CI_Upper","HR_CI_Lower","Peto_LogHR","Expected_Z","Expected_P","Log_Rank_Stat","Variance","V_Pike_Peto","Event_Ratio","Event_Prop_Power","Z_Power"))

#'Calculate analytic time-to-event trial properties under non-proportional hazards or other complex assumptions
#'
#' This function calculates the expected parameters/outputs/properties for a 2-arm Time-To-Event trial under complex assumptions.\cr
#' It is designed to work with non-proportional hazards and ought to be able to accommodate any distributional assumptions for
#' events, censoring and recruitment, so long as they are correctly detailed in Curve or RCurve objects.\cr
#' It performs power calculation and hence sample size planning. However, identifying the optimum assessment time is also key.\cr
#' It uses numerical integration across event, censoring and recruitment functions to calculate expected observed and expected expected event numbers.\cr
#' From these, it estimates an expected HR using the Pike method, with the same interpretation as that found using Cox regression.\cr
#' The expected event numbers and HR can then be used in calculating power by one of several methods, including the Schoenfeld and Frontier methods.\cr
#' A separate, direct, power calculation can also be performed using the log-rank test formula and its Z-distribution.\cr
#' To assist sample size finding, it will also estimate the required sample size to reach a given power keeping all variables other than recruitment.\cr
#' @param active_ecurve Event distribution for the active arm, specified as a Curve object
#' @param control_ecurve Event distribution for the control arm, specified as a Curve object
#' @param active_dcurve Dropout/censoring distribution for the active arm, specified as a Curve object. By default, a Blank() object, i.e. no dropout.
#' @param control_dcurve Dropout/censoring distribution for the control arm, specified as a Curve object. By default, a Blank() object, i.e. no dropout.
#' @param rcurve Recruitment distribution, specified as an RCurve object
#' @param max_assessment Maximum assessment time to calculate properties up to
#' @param landmark (Optional) Time in months of landmark analysis, if required. Otherwise NULL (Not calculated; default).
#' @param RMST (Optional) Restriction time for RMST analysis in months, if required. Otherwise NULL (Not calculated; default).
#' @param alpha1 One-sided alpha required, as a decimal. 0.025 by default. Requires 0 < alpha1 <= 0.5.
#' @param required_power (Optional) Power required for estimated sample sizes. Otherwise NULL (not calculated; default).
#' @param detailed_output Boolean to require a more detailed output table, including Peto LogHR, expectations of various quantities and alternative power calculations. Default = FALSE (detailed outputs omitted).
#' @import stats
#' @return Returns a table with one row per assessment time. Table contains both all input parameters as well as the following expected quantities:
#' \itemize{
#'  \item{"Time"}{ Time at which the assessment is made}
#'  \item{"Patients"}{ Number of patients expected to be recruited to date}
#'  \item{"Events_Active"}{ Expected number of observed events in active arm}
#'  \item{"Events_Control"}{ Expected number of observed events in control arm}
#'  \item{"Events_Total"}{ Expected number of events across both arms}
#'  \item{"HR"}{ Expected Hazard Ratio (using the Pike method)}
#'  \item{"LogHR"}{ Log of the expected Hazard Ratio}
#'  \item{"LogHR_SE"}{ SE of the log of the expected Hazard Ratio}
#'  \item{"Schoenfeld_Power"}{ Estimated power based on Schoenfeld formula}
#'  \item{"Frontier_Power"}{ Estimated power based on Frontier method, using estimated event ratio at 0.5 power}
#' }
#' In addition, if the detailed_output argument is set to TRUE, the following additional columns are provided:
#' \itemize{
#'  \item{"E_Events_Active"}{ Expected number of expected events in active arm}
#'  \item{"E_Events_Control"}{ Expected number of expected events in control arm}
#'  \item{"HR_CI_Upper"}{ Estimated Upper Bound of the CI for the Hazard Ratio}
#'  \item{"HR_CI_Lower"}{ Estimated Lower Bound of the CI for the Hazard Ratio}
#'  \item{"Peto_LogHR"}{ Expected Log Hazard Ratio using the Peto method}
#'  \item{"Expected_Z"}{ Estimated Z-score based on expected quantities for O, E and V, and log-rank test formula}
#'  \item{"Expected_P"}{ Estimated p-value based on estimated Z-score}
#'  \item{"Log_Rank_Stat"}{ Expected log-rank statistic}
#'  \item{"Variance"}{ Expected variance of LR-statistic by integration of V_function}
#'  \item{"V_Pike_Peto"}{ Expected variance based upon Pike and Peto approximations}
#'  \item{"Event_Ratio"}{ Expected ratio of events between arms; active divided by control}
#'  \item{"Event_Prop_Power"}{ Estimated power based on event proportion method, using event ratio rather than randomisation ratio}
#'  \item{"Z_Power"}{ Estimated power based on expected value of Z}
#' }
#' If RMST calculations are requested, the following columns are included:
#' \itemize{
#'  \item{"RMST_Restrict"}{ Specified RMST restriction time}
#'  \item{"RMST_Active"}{ Expected RMST for active arm}
#'  \item{"RMST_Control"}{ Expected RMST for control arm}
#'  \item{"RMST_Delta"}{ Absolute difference in expected RMSTs between arms (active minus control)}
#'  \item{"RMST_SE"}{ Estimated SE of the RMST delta}
#'  \item{"RMST_Z"}{ Estimated RMST Z score}
#'  \item{"RMST_Failure"}{ Estimated probability of RMST difference being uncomputable for the specified restriction time}
#'  \item{"RMST_Power"}{ Estimated RMST Power}
#' }
#' If landmark calculations are requested, the following columns are included:
#' \itemize{
#'  \item{"LM_Time"}{ Time of landmark analysis}
#'  \item{"LM_Active"}{ Expected Kaplan Meier estimate of active arm at landmark time }
#'  \item{"LM_Control"}{ Expected Kaplan Meier estimate of control arm at landmark time }
#'  \item{"LM_Delta"}{ Expected absolute difference in Kaplan Meiers estimates at landmark time (active-control) }
#'  \item{"LM_A_SE"}{ Estimated Greenwood SE for active arm at landmark time }
#'  \item{"LM_C_SE"}{ Estimated Greenwood SE for control arm at landmark time }
#'  \item{"LM_D_SE"}{ Estimated Greenwood SE for delta at landmark time }
#'  \item{"LM_Z"}{ Estimated landmark analysis Z-score based on Greenwood SE }
#'  \item{"LM_Power"}{ Estimated landmark analysis power based on Greenwood SE }
#' }
#' If a required power is requested, the following column is included:
#' \itemize{
#'  \item{"Estimated_SS"}{ Estimated sample size required, keeping constant all parameters other than rate of recruitment and total sample size}
#' }
#' @author James Bell
#' @references Bell J, Accurate Sample Size Calculations in Trials with Non-Proportional Hazards, 2018, presentation at PSI Conference.
#' \url{https://www.psiweb.org/docs/default-source/default-document-library/james-bell-slides.pdf?sfvrsn=3324dedb_0}
#' Bell J, Power Calculations for Time-to-Event Trials Using Predicted Event Proportions, 2019, paper under review.
#' Ruehl J, Sample Size Calculation in Time-To-Event Trials with Non-Proportional Hazards Using GESTATE, 2018, BSc thesis at University of Ulm.
#' Pike MC, Contribution to discussion in Asymptotically efficient rank invariant test procedures by Peto R and Peto J,
#' Journal of the Royal Statistical Society Series A, 135(2), 201-203.
#' @examples nph_traj(max_assessment=100,rcurve=LinearR(12,100,100),control_ecurve=Weibull(100,1),
#' active_ecurve=Weibull(250,0.8))
#' @export
nph_traj <- function(active_ecurve,control_ecurve,active_dcurve=Blank(),control_dcurve=Blank(),rcurve,max_assessment=100,landmark=NULL,RMST=NULL,alpha1=0.025,required_power=NULL,detailed_output=FALSE){
  # Perform checks of inputs
  # Firstly, check for missing arguments
  if(missing(active_ecurve))stop("Please specify the active event curve using the 'active_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(control_ecurve))stop("Please specify the control event curve using the 'control_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(rcurve))stop("Please specify the recruitment distribution using the 'rcurve' argument. Please note that this should be an RCurve object; Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")

  # Secondly, check for improper arguments, in particular that inputs are Curve objects
  #Curve and RCurve objects
  if(class(active_ecurve)[1]!= "Curve") stop("Argument 'active_ecurve' must be a Curve object in order to define the active event curve. Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(control_ecurve)[1]!= "Curve") stop("Argument 'control_ecurve' must be a Curve object in order to define the control event curve. Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(active_dcurve)[1]!= "Curve") stop("Argument 'active_dcurve' must be a Curve object in order to define the active censoring curve. Create one using a Curve constructor function, e.g. active_dcurve <- Weibull(beta=1,lambda=1)")
  if(class(control_dcurve)[1]!= "Curve") stop("Argument 'control_dcurve' must be a Curve object in order to define the control censoring curve. Create one using a Curve constructor function, e.g. control_dcurve <- Weibull(beta=1,lambda=1)")
  if(class(rcurve)[1]!= "RCurve") stop("Argument 'rcurve' must be an RCurve object in order to define a recruitment distribution. Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")
  #Positive integers
  if(max_assessment%%1!=0 | max_assessment < 1 ) stop("Please specify a positive integer for the maximum assessment time using the 'max_assessment' argument")
  if(!is.null(landmark) && (landmark%%1!=0 | landmark < 1 )) stop("Landmark calculations specified, but landmark time is not a positive integer.")
  if(!is.null(RMST) && (RMST%%1!=0 | RMST < 1 )) stop("RMST calculations specified, but restriction time is not a positive integer.")
  #Positive continuous
  if(!is.numeric(alpha1) || alpha1 <= 0 || alpha1 > 0.5 ) stop("Please specify a positive one-sided alpha less than 0.5 using the 'alpha1' argument.")
  if(!is.null(required_power) && (!(required_power > 0) | !(required_power < 1) )) stop("Sample size for a required power is specified, but power is not a value between 0 and 1 (exclusive).")
  #Boolean
  if(!is.logical(detailed_output))stop("Error: detailed_output argument must be boolean: default=FALSE (simplified output).")

  #As the expected HR method is only accurate, not exact, where simple PH can be detected, GESTATE will
  #preferentially use the easily-calculated exact theoretical value. Exp and Weibull curve detection supported
  #This typically affects power by less than 1%, but is nonetheless an improvement.

  #Check for theoretical PH in Exponential/Weibull case
  PH <- FALSE
  PHa <- -1234567
  PHb <- -9876543

  if(getType(active_ecurve) == "Exponential"){
    PHa <- 1
    lambdaA <- getParams(active_ecurve)[[1]]
  } else if(getType(active_ecurve) == "Weibull"){
    PHa <- getParams(active_ecurve)[[2]]
    lambdaA <- getParams(active_ecurve)[[1]]^-getParams(active_ecurve)[[2]]
  }
  if(getType(control_ecurve) == "Exponential"){
    PHb <- 1
    lambdaC <- getParams(control_ecurve)[[1]]
  } else if(getType(control_ecurve) == "Weibull"){
    PHb <- getParams(control_ecurve)[[2]]
    lambdaC <- getParams(control_ecurve)[[1]]^-getParams(control_ecurve)[[2]]
  }
  if(PHa == PHb){
    PH <- TRUE
    HRtemp <- lambdaA/lambdaC
  }
  N_active <- getNactive(rcurve)
  N_control <- getNcontrol(rcurve)
  Ratio <- getRatio(rcurve)

  # Set up the output data frame list
  output_storage <- vector("list", max_assessment)
  temp_storage <- data.frame(t(rep(NA,23)))

  colnames(temp_storage) <- c(
    "Time",
    "Patients",
    "Events_Active",
    "Events_Control",
    "E_Events_Active",
    "E_Events_Control",
    "Events_Total",
    "HR",
    "LogHR",
    "LogHR_SE",
    "HR_CI_Upper",
    "HR_CI_Lower",
    "Peto_LogHR",
    "Expected_Z",
    "Expected_P",
    "Log_Rank_Stat",
    "Variance",
    "V_Pike_Peto",
    "Event_Ratio",
    "Schoenfeld_Power",
    "Event_Prop_Power",
    "Z_Power",
    "Frontier_Power"
  )
  # If landmark enabled, set up additional vector to contain Greenwood results
  if(!is.null(landmark)){
    green_active <- green_control <- rep(NA,max_assessment)
  }

  # If RMST enabled, set up additional infrastructure to contain RMST results
  if(!is.null(RMST)){
    RMST_storage <- vector("list", max_assessment)
    RMST_temp <- data.frame(t(c(RMST,rep(NA,7))))
    colnames(RMST_temp) <- c(
      "RMST_Restrict",
      "RMST_Active",
      "RMST_Control",
      "RMST_Delta",
      "RMST_SE",
      "RMST_Z",
      "RMST_Power",
      "RMST_Failure"
    )
  }

  # Self-modifying R code! This statement creates the expected observed event function tailored specifically to the combinations of events, censoring and recruitment
  # found within each arm to calculate the events and number at risk. Relatively fast as only needs to do this once per run, and allows code to adapt to all sorts of input.
  # Note that is requires S4 curve objects for the event and censoring distributions, as well as an S4 RCurve object for the recruitment distribution
  # The three 'functions' (getAssessCDFfunction, getCDFfunction and getPDFfunction) are all methods of the Curve and/or RCurve classes that return the appropriate function call.
  # eval(parse(text=paste(XXXXXXXXXX))) is used to 'run' the text of the function that is written.

  # Write function for numbers at risk in the control arm
  eval(parse(text=paste(
    "Nfunction_control <- function(tim,assess){",
    N_control,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(control_dcurve,q="tim"),")*(1-",getCDFfunction(control_ecurve,q="tim"),")
}"
    )))

  # Write function for observed events in the control arm
  eval(parse(text=paste(
    "Dfunction_control <- function(tim,assess){",
    N_control,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(control_dcurve,q="tim"),")*",getPDFfunction(control_ecurve,x="tim"),
    "}"
  )))

  # Write function for numbers at risk in the active arm
  eval(parse(text=paste(
    "Nfunction_active <- function(tim,assess){",
    N_active,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(active_dcurve,q="tim"),")*(1-",getCDFfunction(active_ecurve,q="tim"),")
}"
    )))

  # Write function for observed events in the active arm
  eval(parse(text=paste(
    "Dfunction_active <- function(tim,assess){",
    N_active,"*(1-",getAssessCDFfunction(rcurve,q="tim"),")*(1-",getCDFfunction(active_dcurve,q="tim"),")*",getPDFfunction(active_ecurve,x="tim"),
    "}"
  )))

  # Write CDF for events in the control arm
  eval(parse(text=paste(
    "Sfunction_control <- function(tim){
    (1-",
    getCDFfunction(control_ecurve,q="tim"),
    ")
}"
      )))

  # Write CDF for events in the active arm
  eval(parse(text=paste(
    "Sfunction_active <- function(tim){
    (1-",
    getCDFfunction(active_ecurve,q="tim"),
    ")
}"
      )))

  # Loop through all the possible (integer) assessment times
  # For each, calculate the expected observed number of events in each arm, and the expected expected number of events in the active arm

  for(i in 1:max_assessment){
    O1 <- integrate(Dfunction_control,lower=0,upper=i,assess=i)$value
    O2 <- integrate(Dfunction_active,lower=0,upper=i,assess=i)$value
    E2 <- integrate(Efunction_active,lower=0,upper=i,assess=i,Nfunction_control=Nfunction_control,Dfunction_control=Dfunction_control,Nfunction_active=Nfunction_active,Dfunction_active=Dfunction_active)$value
    V  <- integrate(Vfunction,lower=0,upper=i,assess=i,Nfunction_control=Nfunction_control,Dfunction_control=Dfunction_control,Nfunction_active=Nfunction_active,Dfunction_active=Dfunction_active)$value

    # Calculations for RMST
    if(!is.null(RMST)){
      if(i > RMST){
        restrict <- min(RMST,i)
        RMST_temp$RMST_Control <- RMST(Sfunction_control,restriction=restrict)
        RMST_temp$RMST_Active <- RMST(Sfunction_active,restriction=restrict)
        RMST_temp$RMST_Failure <- (1-(1-((N_control-Nfunction_control(tim=restrict,assess=i))/N_control)^N_control)*(1-((N_active-Nfunction_active(tim=restrict,assess=i))/N_active)^N_active))

        RMSTV1 <- tryCatch(
          {RMST_Var(n=N_Control,Dfunction=Dfunction_control,Sfunction=Sfunction_control,restriction=restrict,assess=i)[[3]]
          },error=function(e){NA}
        )
        RMSTV2 <- tryCatch(
          {RMST_Var(n=N_Active,Dfunction=Dfunction_active,Sfunction=Sfunction_active,restriction=restrict,assess=i)[[3]]
          },error=function(e){NA}
        )
        RMST_temp$RMST_SE <- sqrt(RMSTV1+RMSTV2)
        RMST_storage[[i]] <- RMST_temp
      } else{
        RMST_storage[[i]] <- c(NA,NA,NA,1)
      }
    }

    # Store the observed and expected events, as well as the assessment time and the expected number of recruited patients
    # Then copy into the temporary storage list, output_storage
    temp_storage$Time <- i
    temp_storage$Events_Control <- O1
    temp_storage$Events_Active <- O2
    temp_storage$E_Events_Active <- E2
    temp_storage$Patients <- getPatients(rcurve,i)
    temp_storage$Variance <- V

    output_storage[[i]] <- temp_storage

    # If landmark SE calculations required, integrate the continuised Greenwood formula. Note: will return NA if none at risk at landmark
    if(!is.null(landmark)){
      if(landmark <= i){
        green_active[i] <- tryCatch(
          {integrate(Greenwoodfunction,lower=0,upper=landmark,assess=i,Nfunction=Nfunction_active,Dfunction=Dfunction_active)$value
          },error=function(e){NA}
        )
        green_control[i] <- tryCatch(
          {integrate(Greenwoodfunction,lower=0,upper=landmark,assess=i,Nfunction=Nfunction_control,Dfunction=Dfunction_control)$value
          },error=function(e){NA}
        )
      } else{
        green_active[i] <- NA
        green_control[i] <- NA
      }
    }
  }

  # Bind together the list of single-row outputs into the output data frame
  output <- do.call("rbind",output_storage)

  # Across all times, calculate derived quantities
  # First, define quantities to simplify future code
  O1 <- output$Events_Control
  O2 <- output$Events_Active
  E2 <- output$E_Events_Active

  # Round numbers
  output$Events_Control <- round(output$Events_Control,3)
  output$Events_Active <- round(output$Events_Active,3)
  output$E_Events_Active <- round(output$E_Events_Active,3)

  # Next, calculate some common derived quantities
  # Note that HR will use exact value for detectable PH (see earlier), or estimated value for NPH
  # Careful that rounding does not occur before values are used for a later calculation
  E1 <- O1+O2-E2
  if(PH) {HR <- HRtemp} else { HR <- (O2*E1)/(O1*E2) }
  Z <- (O2-E2)/sqrt(output$Variance)
  output$Patients <- round(output$Patients,2)
  output$E_Events_Control <- round(E1,3)
  output$HR <- round(HR,4)
  output$LogHR <- round(log(HR),4)
  output$LogHR_SE <- round(abs(log(HR)/Z),5)
  output$HR_CI_Upper <- round(exp(log(HR)+qnorm(1-alpha1)*abs(log(HR)/Z)),4)
  output$HR_CI_Lower <- round(exp(log(HR)-qnorm(1-alpha1)*abs(log(HR)/Z)),4)
  output$Events_Total <- round(O1 + O2,3)
  output$Expected_Z <- round(Z,4)
  output$Expected_P <- round(pnorm(Z),4)
  output$Event_Ratio <- round(O2/O1,4)
  output$Log_Rank_Stat <- round(O2-E2,4)
  output$Peto_LogHR <- round((O2-E2)/output$Variance,4)
  output$V_Pike_Peto <- round(1/(1/E1+1/E2),3)
  output$Variance <- round(output$Variance,3)
  output$Schoenfeld_Power <- round(events2power(events=O1+O2,HR=HR,ratio=N_active/N_control,alpha1=alpha1),4)
  output$Event_Prop_Power <- round(events2power(events=O1+O2,HR=HR,ratio=O2/O1,alpha1=alpha1),4)
  output$Z_Power <- round(ZV2power(Z=Z,V=1,alpha1=alpha1),4)
  output$Frontier_Power <- round(frontierpower(events=O1+O2,HR=HR,Eratio=O2/O1,Rratio=N_active/N_control,startpower=output$Event_Prop_Power,alpha1=alpha1,iter=10),4)

  if(!is.null(required_power) && required_power>0 && required_power <1){
    Estimated_SS <- ceiling((N_active+N_control)*power2events(power=required_power,HR=HR,ratio=N_active/N_control,alpha1=alpha1)/(O1+O2))
    output <- cbind(output,Estimated_SS)
  }

  if(!is.null(RMST)){
    # Bind together the list of single-row outputs into the output data frame
    RMSToutput <- do.call("rbind",RMST_storage)
    RMSToutput$RMST_Restrict <- RMST
    RMSToutput$RMST_Delta <- RMSToutput$RMST_Active-RMSToutput$RMST_Control
    RMSToutput$RMST_Z <- RMSToutput$RMST_Delta/RMSToutput$RMST_SE
    RMSToutput$RMST_Power <- round(ZV2power(Z=RMSToutput$RMST_Z,V=1,alpha1=alpha1),4)
    RMSToutput$RMST_Active <- round(RMSToutput$RMST_Active,4)
    RMSToutput$RMST_Control <- round(RMSToutput$RMST_Control,4)
    RMSToutput$RMST_Delta <- round(RMSToutput$RMST_Delta,4)
    RMSToutput$RMST_SE <- round(RMSToutput$RMST_SE,4)
    RMSToutput$RMST_Z <- round(RMSToutput$RMST_Z,4)
    RMSToutput$RMST_Failure <- round(RMSToutput$RMST_Failure,4)
    output <- cbind(output,RMSToutput)
  }

  # Multiply greenwood integral by the expected KM curve for each arm to give the
  # greenwood SE for each arm
  if(!is.null(landmark)){
    LM_Time <- landmark
    LM_Control <- Sfunction_control(landmark)
    LM_Active <- Sfunction_active(landmark)
    LM_Delta <- LM_Active-LM_Control

    LM_A_SE <- sqrt(green_active)*LM_Active
    LM_C_SE <- sqrt(green_control)*LM_Control
    LM_D_SE <- sqrt(LM_A_SE^2+LM_C_SE^2)
    LM_Z <- LM_Delta/LM_D_SE
    LM_Power <- round(ZV2power(Z=LM_Z,V=1,alpha1=alpha1),4)

    LM_A_SE <- round(LM_A_SE,4)
    LM_C_SE <- round(LM_C_SE,4)
    LM_D_SE <- round(LM_D_SE,4)
    LM_Z <- round(LM_Z,4)
    LM_Control <- round(LM_Control,4)
    LM_Active <- round(LM_Active,4)
    LM_Delta <- round(LM_Delta,4)

    output <- cbind(output, LM_Time, LM_Active, LM_Control, LM_Delta, LM_A_SE, LM_C_SE, LM_D_SE, LM_Z, LM_Power)
    output[is.na(LM_C_SE),"LM_Control"] <- NA
    output[is.na(LM_A_SE),"LM_Active"] <- NA
    output[is.na(LM_D_SE),"LM_Delta"] <- NA
  }
  if(!detailed_output){
    output <- subset(output, select = -c(E_Events_Active, E_Events_Control,HR_CI_Upper,HR_CI_Lower,Peto_LogHR,Expected_Z,Expected_P,Log_Rank_Stat,Variance,V_Pike_Peto,Event_Ratio,Event_Prop_Power,Z_Power))
  }
  outputlist <- list(active_ecurve=active_ecurve,control_ecurve=control_ecurve,active_dcurve=active_dcurve,control_dcurve=control_dcurve,rcurve=rcurve,Summary=output)
  return(outputlist)
}
