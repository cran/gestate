if(getRversion() >= "2.15.1")  utils::globalVariables(c("i","ETime","Rec_Time","Assess","RCTime","CTime","Iter"))

#######################################################################################################
#'Perform simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for simulating generalised two-arm time-to-event trial data for NPH trials with arbitrary event, censoring and recruitment distributions.\cr
#' Event and censoring distributions specified via Curve objects. Recruitment specified via an RCurve object.
#' As it uses same architecture and similar syntax to nph_curve_trajectories, results ought to be directly comparable.\cr
#' Can be used to validate outputs from nph_curve_trajectories.\cr
#' Data sets from this are set up to be automatically analysed with the analyse_sim function.\cr
#' @param active_ecurve Event distribution for the active arm, specified as a Curve object
#' @param control_ecurve Event distribution for the control arm, specified as a Curve object
#' @param active_dcurve Dropout/censoring distribution for the active arm, specified as a Curve object. By default, a Blank() object, i.e. no dropout.
#' @param control_dcurve Dropout/censoring distribution for the control arm, specified as a Curve object. By default, a Blank() object, i.e. no dropout.
#' @param rcurve Recruitment distribution, specified as an RCurve object
#' @param assess Positive number for the assessment time at which administrative censoring will be performed.
#' @param fix_events Positive integer for the number of events to fix (if required), letting the assessment time vary. Alternatively, NULL for fixed time assessment with variable event numbers. Notes: Fixing event numbers overrides any specified assessment time and slows simulation considerably. Default = NULL (fixed analysis time)
#' @param iterations Number of simulations to perform. Depending on trial size, 10,000-20,000 is typically OK to analyse on 8GB RAM.
#' @param seed Seed number to use. Numerical, although if "Rand" is specified, a system-time-derived number will be used.
#' @param detailed_output Boolean to require full details of timings of competing processes. If FALSE, the simplified data only includes the *'ed output columns - this approximately halves RAM requirements. Default=FALSE (simplified).
#' @param output_type "matrix" or "list" specifying the type of output required. "matrix" requests a single matrix with a column "iter" to denote the simulation, while "list" creates a list with one entry per simulation. Default="matrix".
#' @import survival
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom survRM2 rmst2
#' @return Returns a table with one row per patient per simulation. Table contains the following columns:
#' \itemize{
#'  \item{"Time"}{ Simulated actually observed (patient) time of event or censoring: This is the main column of interest for analysis*}
#'  \item{"Censored"}{ Simulated censoring indicator: 1 denotes censoring (administrative or dropout), 0 denotes an event*}
#'  \item{"Trt"}{ Treatment group number - 1 is active, 2 is control*}
#'  \item{"Iter"}{ Simulation number*}
#'  \item{"ETime"}{ Simulated actual event (patient) time (may or may not be observed)}
#'  \item{"CTime"}{ Simulated actual censoring/dropout (patient) time (may or may not be observed)}
#'  \item{"Rec_Time"}{ Simulated (trial) time of recruitment}
#'  \item{"Assess"}{ Prespecified (trial) time of assessment}
#'  \item{"RCTime"}{ Simulated actual administrative censoring (patient) time (may or may not be observed)}
#' }
#' @author James Bell
#' @examples
#' example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=20, iterations=5,seed=12345)
#' @export
simulate_trials <- function(active_ecurve,control_ecurve,active_dcurve=Blank(),control_dcurve=Blank(),rcurve,assess=NULL,fix_events=NULL,iterations,seed,detailed_output=FALSE,output_type=c("matrix","list")){

  # Perform input checks
  # Firstly, check for missing arguments
  if(missing(active_ecurve))stop("Please specify the active event curve using the 'active_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(control_ecurve))stop("Please specify the control event curve using the 'control_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(rcurve))stop("Please specify the recruitment distribution using the 'rcurve' argument. Please note that this should be an RCurve object; Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")
  if(missing(seed)) stop("Please specify a value for the random seed using the 'seed' argument. For a system-time-derived random seed number, use 0")
  if(missing(iterations)) stop("Please specify the required number of simulations using the 'iterations' argument")
  if(is.null(assess) & is.null(fix_events)) stop("Please specify either an assessment time (assess) or a fixed number of events (fix_events)")

  # Secondly, check argument forms
  #Curves
  if(class(control_ecurve)[1]!= "Curve") stop("Argument 'control_ecurve' must be a Curve object in order to define the control event curve. Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(active_ecurve)[1]!= "Curve") stop("Argument 'active_ecurve' must be a Curve object in order to define the active event curve. Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(control_dcurve)[1]!= "Curve") stop("Argument 'control_dcurve' must be a Curve object in order to define the control censoring curve. Create one using a Curve constructor function, e.g. control_dcurve <- Weibull(beta=1,lambda=1)")
  if(class(active_dcurve)[1]!= "Curve") stop("Argument 'active_dcurve' must be a Curve object in order to define the active censoring curve. Create one using a Curve constructor function, e.g. active_dcurve <- Weibull(beta=1,lambda=1)")
  if(class(rcurve)[1]!= "RCurve") stop("Argument 'rcurve' must be an RCurve object in order to define a recruitment distribution. Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")

  #Positive integers
  if(iterations%%1!=0 || iterations < 1 ) stop("Please specify a positive integer for the number of simulations using the 'iterations' argument")
  if(is.null(assess) ||(!is.null(fix_events) & assess < getLength(rcurve))){assess=ceiling(getLength(rcurve)+1)}
  if(!is.numeric(assess) || assess < 1 ) stop("Please specify a positive number for the assessment time using the 'assess' argument")
  if(seed%%1!=0 || seed < 0 ) stop("Please specify a positive integer for the seed (or 0 for a random one) using the 'seed' argument")
  if(!is.null(fix_events) && (fix_events%%1!=0 || fix_events < 1 )) stop("Fixed event assessment specified, but fix_events is not a positive integer.")

  #Multiple choice
  output_type <- match.arg(output_type)

  #Boolean
  if(!is.logical(detailed_output))stop("Error: detailed_output argument must be boolean: default=FALSE (simplified output).")

  if(is.null(fix_events) && assess<getLength(rcurve)){
    message("Note: assessment time is shorter than the length of the recruitment period:",getLength(rcurve),". Any attempt to increase the assessment time at a later date will result in missing patients!\n")
  }

  # Set the random seed based on input argument
  if(seed==0){
    seed <- as.numeric(Sys.time())
  }
  set.seed(seed)

  # Write function for simulating events in the control arm
  Econtrol_sim <- createRFfunction(control_ecurve)

  # Write function for simulating events in the active arm
  Eactive_sim <- createRFfunction(active_ecurve)

  # Write function for simulating censorings in the control arm
  Dcontrol_sim <- createRFfunction(control_dcurve)

  # Write function for simulating censorings in the active arm
  Dactive_sim <- createRFfunction(active_dcurve)

  # Write function for simulating recruitment
  R_sim <- createRFfunction(rcurve)

  columns <- 9
  nactive <- ceiling(getNactive(rcurve))
  ncontrol <- floor(getNcontrol(rcurve))
  n <- nactive+ncontrol

  if(nactive > 0){
    activematrix <- matrix(rep(0,nactive*iterations*columns),nrow=nactive*iterations)
    # Create events columns
    activematrix[,5] <- Eactive_sim(nactive*iterations)
    # Create censoring columns
    activematrix[,6] <- Dactive_sim(nactive*iterations)
    # Create treatment columns
    activematrix[,3] <- 2
    # Creates iterations columns
    activematrix[,4] <- rep(1:iterations, each=nactive)
  }

  if(ncontrol > 0){
    controlmatrix <- matrix(rep(0,ncontrol*iterations*columns),nrow=ncontrol*iterations)
    # Create events columns
    controlmatrix[,5] <- Econtrol_sim(ncontrol*iterations)
    # Create censoring columns
    controlmatrix[,6] <- Dcontrol_sim(ncontrol*iterations)
    # Create treatment columns
    controlmatrix[,3] <- 1
    # Creates iterations columns
    controlmatrix[,4] <- rep(1:iterations, each=ncontrol)
  } else{output <- activematrix}

  #Combine treatments now that the treatment-specific bits are done
  if(ncontrol > 0 & nactive > 0){
    output <- rbind(activematrix,controlmatrix)
  } else if(nactive <= 0) {output <- controlmatrix}

  # Create recruitment time
  colnames(output) <- c("Time","Censored","Trt","Iter","ETime", "CTime", "Rec_Time", "Assess", "RCTime")
  output[,"Rec_Time"] <- R_sim(n*iterations)
  output <- update_times_assess(output,assess)
  output <- output[output[,"RCTime"]>=0,]
  # This section handles event fixing and choice of output formatting
  # Note that output formatting is also outsourced to the event-fixing function if event-fixing requested
  if(!is.null(fix_events)){
    output <- set_event_number(data=output,events=fix_events,output_type=output_type,detailed_output=detailed_output)
  } else {
    if(!detailed_output){
      output <- output[,c("Time","Censored","Trt","Iter")]
    }
    if(output_type=="list"){
      output <- split.data.frame(as.matrix(subset(output,select=-Iter)),output[,"Iter"])
    }
  }
  return(output)
}


#######################################################################################################
#'Perform multi-strata simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for simulating generalised two-arm multi-strata time-to-event trial data for NPH trials with arbitrary event, censoring and recruitment distributions.\cr
#' Acts as a wrapper for simulate_trials.\cr
#' Vector of strata proportions supplies number of strata. Event and censoring distributions specified via lists of Curve objects. If only one Curve supplied then assumed to be common to all strata. Recruitment specified via a single RCurve object.\cr
#' As it uses same architecture and similar syntax to nph_curve_trajectories, results ought to be directly comparable to e.g. use of MixExp or MixWei distributions.\cr
#' Can be used to validate outputs from nph_curve_trajectories.\cr
#' Data sets from this are set up to be automatically analysed with the analyse_sim function (including stratified analysis if you provide it the name of stratum column).\cr
#' @param stratum_probs Vector of probabilities that patients belong to each stratum. Must sum to 1. Its length determines the number of strata.
#' @param active_ecurve List of event distributions for the active arm, specified as a list of Curve objects. If single Curve is specified, will be used for all strata.
#' @param control_ecurve List of event distributions for the control arm, specified as a list of Curve objects. If single Curve is specified, will be used for all strata.
#' @param active_dcurve List of dropout/censoring distribution for the active arm, specified as a Curve object. If single Curve is specified, will be used for all strata. By default, a Blank() object, i.e. no dropout in any stratum.
#' @param control_dcurve List of dropout/censoring distribution for the control arm, specified as a Curve object. If single Curve is specified, will be used for all strata. By default, a Blank() object, i.e. no dropout in any stratum.
#' @param rcurve Recruitment distribution, specified as a single RCurve object.
#' @param assess Positive number for the assessment time at which administrative censoring will be performed.
#' @param fix_events Positive integer for the number of events to fix (if required), letting the assessment time vary. Alternatively, NULL for fixed time assessment with variable event numbers. Notes: Fixing event numbers overrides any specified assessment time and slows simulation considerably. Default = NULL (fixed analysis time)
#' @param stratum_name Name of the column defining the stratum. Default="Stratum".
#' @param iterations Number of simulations to perform. Depending on trial size, 10,000-20,000 is typically OK to analyse on a laptop. 100,000 typically requires a system with more RAM.
#' @param seed Seed number to use. Numerical, although if "Rand" is specified, a system-time-derived number will be used.
#' @param detailed_output Boolean to require full details of timings of competing processes. If FALSE, the simplified data only includes the *'ed output columns - this approximately halves RAM requirements. Default=FALSE (simplified).
#' @param output_type "matrix" or "list" specifying the type of output required. "matrix" requests a single matrix with a column "iter" to denote the simulation, while "list" creates a list with one entry per simulation. Default="matrix".
#' @import survival
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom survRM2 rmst2
#' @return Returns a table with one row per patient per simulation. Table contains the following columns:
#' \itemize{
#'  \item{"Time"}{ Simulated actually observed (patient) time of event or censoring: This is the main column of interest for analysis*}
#'  \item{"Censored"}{ Simulated censoring indicator: 1 denotes censoring (administrative or dropout), 0 denotes an event*}
#'  \item{"Trt"}{ Treatment group number - 1 is active, 2 is control*}
#'  \item{"Iter"}{ Simulation number*}
#'  \item{"ETime"}{ Simulated actual event (patient) time (may or may not be observed)}
#'  \item{"CTime"}{ Simulated actual censoring/dropout (patient) time (may or may not be observed)}
#'  \item{"Rec_Time"}{ Simulated (trial) time of recruitment}
#'  \item{"Assess"}{ Prespecified (trial) time of assessment}
#'  \item{"RCTime"}{ Simulated actual administrative censoring (patient) time (may or may not be observed)}
#'  \item{"Stratum"}{ Stratum number. Column name will be the value of the stratum_name argument.)}
#' }
#' @author James Bell
#' @examples
#' example_strat_sim <- simulate_trials_strata(stratum_probs=c(0.5,0.5),
#' active_ecurve=c(Weibull(250,0.8),Weibull(100,1)), control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100),assess=20,iterations=5,seed=12345)
#' @export
simulate_trials_strata <- function(stratum_probs,active_ecurve,control_ecurve,active_dcurve=Blank(),control_dcurve=Blank(),rcurve,assess=NULL,fix_events=NULL,stratum_name="Stratum",iterations,seed,detailed_output=FALSE,output_type=c("matrix","list")){
# Input checks and define number of strata
  output_type <- match.arg(output_type)

  # Perform input checks
  # Firstly, check for missing arguments
  if(missing(stratum_probs))stop("Please specify a vector of probabilities for strata membership (summing to 1) using the 'stratum_probs' argument.")
  if(missing(active_ecurve))stop("Please specify the active event curve using the 'active_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(control_ecurve))stop("Please specify the control event curve using the 'control_ecurve' argument. Please note that this should be a Curve object; Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(missing(rcurve))stop("Please specify the recruitment distribution using the 'rcurve' argument. Please note that this should be an RCurve object; Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")
  if(missing(seed)) stop("Please specify a value for the random seed using the 'seed' argument. For a system-time-derived random seed number, use \"random\"")
  if(missing(iterations)) stop("Please specify the required number of simulations using the 'iterations' argument")
  if(is.null(assess) & is.null(fix_events)) stop("Please specify either an assessment time (assess) or a fixed number of events (fix_events)")

  # Secondly, check argument forms
  if(class(rcurve)[1]!= "RCurve") stop("Argument 'rcurve' must be an RCurve object in order to define a recruitment distribution. Create one using an RCurve constructor function, e.g. rcurve <- LinearR(rlength=10,Nactive=100,Ncontrol=100)")
  #Positive integers
  if(is.null(assess) ||(!is.null(fix_events) & assess < getLength(rcurve))){assess=ceiling(getLength(rcurve)+1)}
  if(iterations%%1!=0 | iterations < 1 ) stop("Please specify a positive integer for the number of simulations using the 'iterations' argument")
  if(!is.numeric(assess) || assess < 1 ) stop("Please specify a positive number for the assessment time using the 'assess' argument")
  if(seed%%1!=0 | seed < 1 ) stop("Please specify a positive integer for the seed using the 'seed' argument")
  if(!is.null(fix_events) && (fix_events%%1!=0 | fix_events < 1 )) stop("Fixed event assessment specified, but fix_events is not a positive integer.")

  #Boolean
  if(!is.logical(detailed_output))stop("Error: detailed_output argument must be boolean: default=FALSE (simplified output).")


  #String
  if(!is.character(stratum_name))stop("Error: stratum_name argument must be a string: default='Stratum'.")

  # Check stratum-related things
  nstrata <- length(stratum_probs)
  if(sum(stratum_probs)!=1)stop("stratum_probs must sum to 1")
  if(length(control_ecurve)!=nstrata & length(control_ecurve)!=1)stop("Number of control_ecurves specified is inconsistent with number of strata")
  if(length(active_ecurve)!=nstrata & length(active_ecurve)!=1)stop("Number of active_ecurves specified is inconsistent with number of strata")
  if(length(control_dcurve)!=nstrata & length(control_dcurve)!=1)stop("Number of control_dcurves specified is inconsistent with number of strata")
  if(length(active_dcurve)!=nstrata & length(active_dcurve)!=1)stop("Number of active_dcurves specified is inconsistent with number of strata")
# Catch all single-Curve inputs, apply them to all strata
  if(length(control_ecurve)==1){
    control_ecurve <- list(control_ecurve)[rep(1,nstrata)]
  }
  if(length(active_ecurve)==1){
    active_ecurve <- list(active_ecurve)[rep(1,nstrata)]
  }
  if(length(control_dcurve)==1){
    control_dcurve <- list(control_dcurve)[rep(1,nstrata)]
  }
  if(length(active_dcurve)==1){
    active_dcurve <- list(active_dcurve)[rep(1,nstrata)]
  }
  rcurve <- list(rcurve)[rep(1,nstrata)]

  #Curves
  if(class(control_ecurve[[1]])[1]!= "Curve") stop("Argument 'control_ecurve' must be a Curve object (or list thereof) in order to define the control event curve. Create one using a Curve constructor function, e.g. control_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(active_ecurve[[1]])[1]!= "Curve") stop("Argument 'active_ecurve' must be a Curve object (or list thereof) in order to define the active event curve. Create one using a Curve constructor function, e.g. active_ecurve <- Weibull(beta=1,lambda=1)")
  if(class(control_dcurve[[1]])[1]!= "Curve") stop("Argument 'control_dcurve' must be a Curve object (or list thereof) in order to define the control censoring curve. Create one using a Curve constructor function, e.g. control_dcurve <- Weibull(beta=1,lambda=1)")
  if(class(active_dcurve[[1]])[1]!= "Curve") stop("Argument 'active_dcurve' must be a Curve object (or list thereof) in order to define the active censoring curve. Create one using a Curve constructor function, e.g. active_dcurve <- Weibull(beta=1,lambda=1)")

# Go through all strata in turn, simulating them. Add a column for the stratum number (currently 'temp'). All processing/output choices not done till later.
  outlist <- vector("list",nstrata)
  for(i in 1:nstrata){
    rcurve[[i]] <- setPatients(rcurve[[i]],Nactive=getNactive(rcurve[[i]])*stratum_probs[[i]],Ncontrol=getNcontrol(rcurve[[i]])*stratum_probs[[i]])
    outlist[[i]] <- simulate_trials(active_ecurve=active_ecurve[[i]],control_ecurve=control_ecurve[[i]],active_dcurve=active_dcurve[[i]],control_dcurve=control_dcurve[[i]],rcurve=rcurve[[i]],assess=assess,iterations=iterations,seed=seed+i,fix_events=NULL,detailed_output=TRUE,output_type="matrix")
    temp <- rep(i,nrow(outlist[[i]]))
    outlist[[i]] <- cbind(outlist[[i]],temp)
  }
  output <- do.call(rbind,outlist)
# Name stratum column, do event-fixing and output processing at this stage now data set is complete
  colnames(output)[which(colnames(output) == "temp")] <- stratum_name
  if(!is.null(fix_events)){
    output <- set_event_number(data=output,events=fix_events,output_type=output_type,detailed_output=detailed_output)
  } else {
    if(!detailed_output){
      output <- output[,c("Time","Censored","Trt",stratum_name,"Iter")]
    }
    if(output_type=="list"){
      output <- split.data.frame(as.matrix(subset(output,select=-Iter)),output[,"Iter"])
    }
  }
  return(output)
}

#Internal function to apply a new assessment time to a simulation
#Sets the assessment time, then updates all the other columns to reflect this
#Function here just to avoid repetition in simulate_trials, set_event_number and set_assess_time
update_times_assess <- function(data,assess){
    data[,"Assess"] <- assess
    data[,"RCTime"] <- data[,"Assess"] - data[,"Rec_Time"]
    data[,"Time"] <- pmin.int(data[,"ETime"],data[,"CTime"],data[,"RCTime"])
    data[,"Censored"] <- 1
    data[data[,"Time"] == data[,"ETime"],"Censored"] <- 0
    return(data)
}

#######################################################################################################
#'Adjusts simulations so that administrative censoring occurs at a fixed event number, rather than a fixed time
#'
#' Function for converting trials simulated under simulate_trials from a fixed censoring time to a fixed number of total events.\cr
#' Set up to automatically read in either matrix or list formats from simulate_trials, and only these inputs are supported.\cr
#' Note that if recruitment had not finished in the input then any increases in assessment time cannot account for the missing patients.\cr
#' It is strongly recommended to initially simulate for at least the duration of the recruitment before fixing the event number.\cr
#' This function can also be used to change format and/or slim down data for event-driven simulations.\cr
#' @param data Output file from simulate_trials. Only simulate_trials output is supported.
#' @param events Positive integer specifying the required number of events.
#' @param output_type Choice of "input" (output in same format as input),"matrix" (matrix format output) or "list" (list format output). Default="input".
#' @param detailed_output Boolean to require full details of timings of competing processes. If FALSE, the simplified data only includes the *'ed output columns - this approximately halves RAM requirements. Default=TRUE (detailed).
#' @return Returns the input simulated trial, in either matrix or list format, with modified assessment times. All columns dependent on this are also updated.
#' @author James Bell
#' @examples example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=20,iterations=5,seed=12345,detailed_output=TRUE)
#'
#' adjusted_examples <- set_event_number(data=example_sim,events=50)
#' @export
set_event_number <- function(data,events,output_type=c("input","matrix","list"),detailed_output=TRUE){

  #Required Input Check
  if(missing(data))stop("Please specify the name of the simulated data set ('data').")
  if(!exists("data"))stop("Error: specified data does not exist.")
  if(missing(events))stop("Please specify the fixed number of events ('events').")

  #Multiple choice
  output_type <- match.arg(output_type)
  #Boolean
  if(!is.logical(detailed_output))stop("Error: detailed_output argument must be boolean: default=FALSE (simplified output).")
  #Positive integer
  if(!is.numeric(events) | events%%1!=0 | events < 1 ) stop("Required event number ('events') is not a positive integer.")

  additer <- FALSE
  if(is.list(data) && !is.data.frame(data)){
    maxiter <- length(data)
    data1 <- data
    input <- "list"
    #if(output_type=="matrix"){additer <- TRUE}
    additer <- TRUE
  } else if(is.matrix(data) |(!is.list(data) & is.data.frame(data))){
    if("Iter" %in% colnames(data)){
      maxiter <- max(data[,"Iter"])
      data1 <- split.data.frame(as.matrix(data),data[,"Iter"])
      input <- "matrix"
    } else {
      maxiter <- 1
      data1 <- list(data)
      input <- "single"
    }
  } else stop("Unrecognised data format for set_event_number")
  if(output_type == "input"){
    output_type <- "matrix"
    if(input=="list"){output_type <- "list"}
  }
  if(!("ETime" %in% colnames(data1[[1]])) | !("RCTime" %in% colnames(data1[[1]]))|!("Rec_Time" %in% colnames(data1[[1]])) |!("Assess" %in% colnames(data1[[1]])) |!("CTime" %in% colnames(data1[[1]])))stop("set_event_number requires full, detailed simulation output. Please rerun simulate_trials with detailed_output=TRUE")

  results <- vector("list", maxiter)
  newassess <- rep(NA,maxiter)
  for(i in 1:maxiter){
    iterdata <- data1[[i]]
    possibles <- iterdata[iterdata[,"ETime"]<iterdata[,"CTime"],]
    sorted <- sort(possibles[,"ETime"]-possibles[,"RCTime"])
    sortE <- length(sorted)
    if(sortE <= events){
      cutoff <- sorted[sortE]+0.0001
      warning("Only",sortE,"events are possible in simulation",i,"but",events,"have been requested. Setting assessment time to just after last event.\n")
    }else{cutoff <- sorted[events]+(sorted[events+1]-sorted[events])/1000}
    newassess[i] <- mean(iterdata[,"Assess"]) + cutoff
    iterdata <- update_times_assess(data=iterdata,assess=(newassess[i]))
    if(additer){
      Iter <- rep(i,nrow(iterdata))
      iterdata <- cbind(iterdata,Iter)
    }
    results[[i]] <- iterdata
  }

  out <- rbind(do.call("rbind",results))
  out <- out[out[,"RCTime"]>=0,]
  if(!detailed_output){
    out <- subset(out,select=-c(ETime,Rec_Time,Assess,RCTime,CTime))
  }
  if(output_type=="list"){
    if(input=="single"){
      return(list(out))
    } else{
      return(split.data.frame(subset(out,select=-Iter),out[,"Iter"]))
    }
  } else {return(out)}
}

#######################################################################################################
#'Adjusts assessment time for simulations
#'
#' Function for modifying the assessment time of a simulate_trials simulation.\cr
#' Set up to automatically read in either matrix or list formats from simulate_trials, and only these inputs are supported.\cr
#' Note that if recruitment had not finished in the input then any increases in assessment time cannot account for the missing patients.\cr
#' It is strongly recommended to initially simulate for at least the duration of the recruitment before reducing the number to missing patients.\cr
#' This function can also be used to change format and/or slim down data for time-driven simulations.\cr
#' @param data Output file from simulate_trials. Only simulate_trials output is supported.
#' @param time Positive number specifying the required assessment time.
#' @param output_type Choice of "input" (output in same format as input),"matrix" (matrix format output) or "list" (list format output). Default="input".
#' @param detailed_output Boolean to require full details of timings of competing processes. If FALSE, the simplified data only includes the *'ed output columns - this approximately halves RAM requirements. Default=TRUE (detailed).
#' @return Returns the input simulated trial, in either matrix or list format, with modified assessment times. All columns dependent on this are also updated.
#' @author James Bell
#' @examples example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=20, iterations=5,seed=12345,detailed_output=TRUE)
#'
#' adjusted_example <- set_assess_time(data=example_sim,time=10)
#' @export
set_assess_time <- function(data,time,output_type=c("input","matrix","list"),detailed_output=TRUE){

  #Required Input Check
  if(missing(data))stop("Please specify the name of the simulated data set ('data').")
  if(missing(time))stop("Please specify the fixed assessment time ('time').")

  #Multiple choice
  output_type <- match.arg(output_type)
  #Boolean
  if(!is.logical(detailed_output))stop("Error: detailed_output argument must be boolean: default=FALSE (simplified output).")
  #Positive number
  if(!is.numeric(time)| time < 1 ) stop("Required event number ('time') is not a positive number.")

  if(time < 0){stop("time must be a real, positive number.")}
  output_type <- match.arg(output_type)
  # Checks for list format
  if(is.list(data) && !is.data.frame(data)){
    if(!("ETime" %in% colnames(data[[1]])) | !("RCTime" %in% colnames(data[[1]]))|!("Rec_Time" %in% colnames(data[[1]])) |!("Assess" %in% colnames(data[[1]])) |!("CTime" %in% colnames(data[[1]])))stop("set_assess_time requires full, detailed simulation output. Please rerun simulate_trials with detailed_output=TRUE")
    maxiter <- length(data)
    input <- "list"
    # If output is matrix, convert to matrix then change times in one go. May change to dplyr function later to remove loop
    if(output_type=="matrix"){
      for(i in 1:maxiter){
        Iter <- rep(i,nrow(data[[i]]))
        data[[i]] <- cbind(data[[i]],Iter)
      }
      data <- rbind(do.call("rbind",data))
      data <- update_times_assess(data=data,assess=time)
      data <- data[data[,"RCTime"]>=0,]
      if(!detailed_output){
        data <- subset(data,select=-c(ETime,Rec_Time,Assess,RCTime,CTime))
      }
      return(data)
    } else {
      # If output is list, change times in situ without converting to matrix and back (remove unwanted columns)
      if(!detailed_output){
        for(i in 1:maxiter){
          data[[i]] <- update_times_assess(data=data[[i]],assess=time)
          data[[i]] <- data[[i]][data[[i]][,"RCTime"]>=0,]
          data[[i]] <- subset(data[[i]],select=-c(ETime,Rec_Time,Assess,RCTime,CTime))
        }
      return(data)

      } else {
      # If output is list, change times in situ without converting to matrix and back
        for(i in 1:maxiter){
          data[[i]] <- update_times_assess(data=data[[i]],assess=time)
          data[[i]] <- data[[i]][data[[i]][,"RCTime"]>=0,]
        }
        return(data)
      }
    }
  } else if(is.matrix(data) |(!is.list(data) & is.data.frame(data))){
    # If input is in matrix format do everything necessary and then convert to required format at end
    if(!("ETime" %in% colnames(data)) | !("RCTime" %in% colnames(data))|!("Rec_Time" %in% colnames(data)) |!("Assess" %in% colnames(data)) |!("CTime" %in% colnames(data)))stop("set_assess_time requires full, detailed simulation output. Please rerun simulate_trials with detailed_output=TRUE")
    if("Iter" %in% colnames(data)){
      input <- "matrix"
    } else {
      input <- "single"
    }
      data <- update_times_assess(data=data,assess=time)
      data <- data[data[,"RCTime"]>=0,]
      if(!detailed_output){
        data <- subset(data,select=-c(ETime,Rec_Time,Assess,RCTime,CTime))
      }
      if(output_type=="list"){
        if(input=="single"){return(list(data))}
        return(split.data.frame(subset(data,select=-Iter),data[,"Iter"]))
      } else{
        return(data)
      }
  } else stop("Unrecognised data format for set_assess_time")
}

#######################################################################################################
#'Analyse simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for analysing simulated time-to-event trial data produced by simulate_trials.\cr
#' Automatically reads in format from simulate_trials.\cr
#' Performs log rank test and Cox regression analysis by default, but can also/instead choose RMST and/or landmark analyses.\cr
#' Option is available to perform stratified analysis using the "stratum" argument.\cr
#' If a stratum is specified, it will be included as a covariate in Cox and RMST analysis, and as a stratum in a stratified log-rank test and an inverse-precision-weighted landmark test.
#' Analysis is typically the slowest part of simulation, so parallel processing using the doParallel package is built in, enabled using the "parallel" argument.\cr
#' Use of parallel processing recommended for >10000 simulations. Ensure that the number of cores specified does not exceed number of threads provided by hardware.
#' @param data Output file from simulate_trials. Only simulate_trials output is supported, in either "list" or "matrix" format.
#' @param LR Requests log-rank test and Cox regression. Default=TRUE
#' @param RMST Requests Restricted Mean Survival Time analysis with specified (positive integer) restriction time, leave NULL for no analysis. This uses the survRM2 package. Default=NULL (no RMST analysis).
#' @param landmark Requests Landmark analysis at specified (positive integer) time, leave NULL for no analysis. Default=NULL (no landmark analysis).
#' @param stratum Specify name of column of a stratification factor and turn on stratified (LR/LM) and covariate-adjusted (Cox/RMST) analysis. By default, "", and no stratification.
#' @param parallel_cores Positive integer specifying number of cores to use. If 1 specified then no parallel processing. Default=1 (no parallel processing).
#' @return Returns a table with one row per simulation. Table contains the following columns:
#' \itemize{
#'  \item{"HR"}{ Cox Hazard Ratio (LR/Cox analysis only)}
#'  \item{"LogHR"}{ Cox Log Hazard Ratio (LR/Cox analysis only)}
#'  \item{"LogHR_SE"}{ Cox Standard Error of log Hazard Ratio (LR/Cox analysis only)}
#'  \item{"HR_Z"}{ Cox Z-Score (LR/Cox analysis only)}
#'  \item{"HR_P"}{ 1-sided Cox p-value (LR/Cox analysis only)}
#'  \item{"LR_Z"}{ Log-Rank Test Z-Score (LR/Cox analysis only)}
#'  \item{"LR_P"}{ 1-sided Log-Rank Test p-value (LR/Cox analysis only)}
#'  \item{"Events_Active"}{ Events in Active arm (LR/Cox analysis only)}
#'  \item{"Events_Control"}{ Events in Control arm (LR/Cox analysis only)}
#'  \item{"RMST_Time"}{ RMST restriction time (RMST analysis only)}
#'  \item{"RMST_Active"}{ RMST for Active arm (RMST analysis only)}
#'  \item{"RMST_Active_SE"}{ RMST Standard Error for Active arm (RMST analysis only)}
#'  \item{"RMST_Control"}{ RMST for Control arm (RMST analysis only)}
#'  \item{"RMST_Control_SE"}{ RMST Standard Error for Control arm (RMST analysis only)}
#'  \item{"RMST_Delta"}{ RMST difference between arms active-control (RMST analysis only)}
#'  \item{"RMST_Delta_SE"}{ RMST difference between arms Standard Error (RMST analysis only)}
#'  \item{"RMST_Z"}{ Z-score for RMST (RMST analysis only)}
#'  \item{"RMST_P"}{ 1-sided RMST p-value (RMST analysis only)}
#'  \item{"LM_Time"}{ Landmark time, i.e. time of survival function comparison  (Landmark analysis only)}
#'  \item{"LM_Active"}{ Survival function for active arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Active_SE"}{ Greenwood standard error for active arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Control"}{ Survival function for control arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Control_SE"}{ Greenwood standard error for control arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Delta"}{ Survival function difference between arms active-control at landmark time (Landmark analysis only)}
#'  \item{"LM_Delta_SE"}{ Greenwood standard error for difference between arms at landmark time (Landmark analysis only)}
#'  \item{"LM_Z"}{ Z-score for landmark analysis (Landmark analysis only)}
#'  \item{"LM_P"}{ 1-sided landmark analysis p-value (Landmark analysis only)}
#' }
#' @author James Bell
#' @examples example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=20,iterations=2,seed=12345,detailed_output=TRUE)
#'
#' example_analysis1 <- analyse_sim(example_sim)
#' example_analysis2 <- analyse_sim(data=example_sim,RMST=15,landmark=15)
#'
#' example_strat_sim <- simulate_trials_strata(stratum_probs=c(0.5,0.5),
#' active_ecurve=c(Weibull(250,0.8),Weibull(100,1)), control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100),assess=20,iterations=2,seed=12345)
#'
#' example_strat_analysis <- analyse_sim(data=example_strat_sim,RMST=15,landmark=15,stratum="Stratum")
#' @export
analyse_sim <- function(data, LR = TRUE, RMST = NULL, landmark=NULL, stratum="", parallel_cores=1){

# Check that the correct packages are installed for specific options.
# If not, disable options
  if(parallel_cores!= 1 && !requireNamespace("doParallel", quietly = TRUE)) {
    warning("Parallel processing requested but package \"doParallel\" needed for this. Please install it to enable parallel processing. Reverting to single-thread processing.\n")
    parallel_cores <- 1
  }
  if(!is.null(RMST) && !requireNamespace("doParallel", quietly = TRUE)) {
    warning("RMST analysis requested but package \"survRM2\" needed for this. Please install it to enable RMST analysis.\n")
    RMST <- NULL
  }

  if(!LR & is.null(RMST) & is.null(landmark)){"At least one of 'LR', 'RMST' or 'landmark' must be TRUE to specify analysis to be performed"}

 #Required Input Check
  if(missing(data))stop("Error: Please specify the name of the simulated data set ('data').")

  #Boolean
  if(!is.logical(LR))stop("Error: 'LR' argument must be boolean: default=FALSE (perform log-rank and cox analysis).")
  #Character
  if(!is.character(stratum))stop("Error: 'stratum' argument must be a string: default=FALSE (perform log-rank and cox analysis).")

  #Positive integer
  if(!is.numeric(parallel_cores) || parallel_cores%%1!=0 || parallel_cores < 1 ) stop("Error: 'parallel_cores' argument is not a positive integer.")
  if(!is.null(landmark) && (landmark%%1!=0 || landmark < 1 )) stop("Landmark analysis specified, but landmark time is not a positive integer.")
  if(!is.null(RMST) && (RMST%%1!=0 || RMST < 1 )) stop("RMST analysis specified, but restriction time is not a positive integer.")

  #Format checks and...
  #If in matrix format, split single huge data frame into one per iteration
  #Hugely improves speed of function compared to using parts of parent
  if(is.matrix(data) |(!is.list(data) & is.data.frame(data))){
    if(!("Time" %in% colnames(data)))stop("Error: 'Time' column does not exist: data not in simulate_trials format.")
    if(!("Censored" %in% colnames(data)))stop("Error: 'Censored' column does not exist: data not in simulate_trials format.")
    if(!("Trt" %in% colnames(data)))stop("Error: 'Trt' column does not exist: data not in simulate_trials format.")
    if(!("Iter" %in% colnames(data)))stop("Error: 'Iter' column does not exist: data not in simulate_trials format.")
    data <- split.data.frame(as.matrix(data),data[,"Iter"])
  } else{
    if(!is.list(data))stop("Error: Unrecognised data format.")
    if(!("Time" %in% colnames(data[[1]])))stop("Error: 'Time' column does not exist: data not in simulate_trials format.")
    if(!("Censored" %in% colnames(data[[1]])))stop("Error: 'Censored' column does not exist: data not in simulate_trials format.")
    if(!("Trt" %in% colnames(data[[1]])))stop("Error: 'Trt' column does not exist: data not in simulate_trials format.")
  }

  #Set maximum iterations
  maxiter <- length(data)
  #If parallel processing requested, load parallel library and register cores
  if(parallel_cores>1){
    registerDoParallel(parallel_cores)
  }
  # If log-rank / Cox analysis required, do it. Separate calls depending on whether parallel processing enabled or stratified analysis requested
  if(LR){
    if(parallel_cores>1){
      if(stratum==""){
        outLR <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages="survival")%dopar% LR_analysis(i,data)
      } else {
        outLR <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages="survival")%dopar% LR_cov_analysis(i,data,stratum)
      }
    } else{
      if(stratum==""){
        outLR <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages="survival")%do% LR_analysis(i,data)
      } else {
        outLR <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages="survival")%do% LR_cov_analysis(i,data,stratum)
      }
    }
    if(!is.matrix(outLR)){outLR <- t(outLR)}
    colnames(outLR) <- c("LogHR","HR","LogHR_SE","HR_Z","HR_P","LR_Z","Events_Active","Events_Control")
    outLR[,"HR_P"] <- pnorm(outLR[,"HR_Z"])
  #Note that LR_Z column at this stage is in fact the pchisq_LR, so we convert.
  #Assume that directionality of Cox is same as LR, based on LR being Score test of Cox
    outLR[,"LR_Z"] <- sign(outLR[,"LogHR"])*sqrt(outLR[,"LR_Z"])
    LR_P <- pnorm(outLR[,"LR_Z"])
    outLR <- cbind(outLR,LR_P)
    out <- outLR[,c(2,1,3,4,5,6,9,7,8)]
  }
  # If RMST analysis required, do it. Separate calls depending on whether parallel processing enabled or stratified analysis requested
  if(!is.null(RMST)){
    if(parallel_cores>1){
      if(stratum==""){
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_analysis", .packages="survRM2")%dopar% RMST_analysis(i,data,RMST)
      } else {
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_analysis", .packages="survRM2")%dopar% RMST_cov_analysis(i,data,RMST,stratum)
      }
    }else{
      if(stratum==""){
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_analysis", .packages="survRM2")%do% RMST_analysis(i,data,RMST)
      } else {
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_analysis", .packages="survRM2")%do% RMST_cov_analysis(i,data,RMST,stratum)
      }
    }

    if(!is.matrix(outRMST)){outRMST <- t(outRMST)}
    outRMST <- cbind(RMST,outRMST)
    colnames(outRMST) <- c("RMST_Restrict","RMST_Active","RMST_A_SE","RMST_Control","RMST_C_SE","RMST_Delta","RMST_D_SE","RMST_Z","RMST_P")
    out <- if(LR){cbind(out,outRMST)} else{outRMST}
  }
  # If landmark analysis required, do it. Separate calls depending on whether parallel processing enabled or stratified analysis requested
  if(!is.null(landmark)){
    if(parallel_cores>1){
      if(stratum==""){
        outLM <- foreach(i=1:maxiter,.combine="rbind",.export="LM_analysis", .packages="survival")%dopar% LM_analysis(i,data,landmark)
      } else {
        outLM <- foreach(i=1:maxiter,.combine="rbind",.export="LM_analysis", .packages="survival")%dopar% LM_cov_analysis(i,data,landmark,stratum)
      }
    }else{
      if(stratum==""){
        outLM <- foreach(i=1:maxiter,.combine="rbind",.export="LM_analysis", .packages="survival")%do% LM_analysis(i,data,landmark)
      } else {
        outLM <- foreach(i=1:maxiter,.combine="rbind",.export="LM_analysis", .packages="survival")%do% LM_cov_analysis(i,data,landmark,stratum)
      }
    }
    if(!is.matrix(outLM)){outLM <- t(outLM)}
    delta <- outLM[,1]-outLM[,3]
    delta_se <- sqrt(outLM[,2]^2+outLM[,4]^2)
    Z <- delta/delta_se
    p <- pnorm(-Z)
    outLM <- cbind(landmark,outLM,delta,delta_se,Z,p)
    colnames(outLM) <- c("LM_Time","LM_Active","LM_A_SE","LM_Control","LM_C_SE","LM_Delta","LM_D_SE","LM_Z","LM_P")
    out <- if(LR|!is.null(RMST)){cbind(out,outLM)} else{outLM}
  }
  return(out)
}

################################################################################
# LR_analysis ; Helper function for log-rank/Cox analysis in analyse_sim
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
################################################################################
LR_analysis <- function(i,data){
  iterdata <- data[[i]]
  analysis <- coxph(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"],ties="breslow")
  results <- c(coef(summary(analysis)),analysis$score,
               length(iterdata[(iterdata[,"Trt"]==2 & iterdata[,"Censored"]==0),1]),
               length(iterdata[(iterdata[,"Trt"]==1 & iterdata[,"Censored"]==0),1])
  )
}

################################################################################
# LR_cov_analysis ; Helper function for stratified log-rank/covariate-adjusted Cox analysis in analyse_sim
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  cov:     covariate
################################################################################
LR_cov_analysis <- function(i,data,cov){
  iterdata <- data[[i]]
  analysis <- coxph(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"] + iterdata[,cov],ties="breslow")
  analysisLR <- survdiff(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"] + strata(iterdata[,cov]))
  results <- c(coef(summary(analysis))[1,],analysisLR$chisq,
               length(iterdata[(iterdata[,"Trt"]==2 & iterdata[,"Censored"]==0),1]),
               length(iterdata[(iterdata[,"Trt"]==1 & iterdata[,"Censored"]==0),1])
  )
}

###############################################################################
# RMST_analysis ; Helper function for RMST analysis in analyse_sim
# Requires the 'survRM2' package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  restriction: positive value for restriction time for RMST analysis
#
# Note: This returns the SE, rather than the CI
###############################################################################
RMST_analysis <- function(i,data,restriction){
  iterdata <- data[[i]]
  results <- rep(NA,8)
  tryCatch(
    {
      analysis <- rmst2(time=iterdata[,"Time"],status=1-iterdata[,"Censored"],arm=iterdata[,"Trt"]-1,tau=restriction,alpha=0.05)
      SE <- abs(analysis$unadjusted.result[1,1]-analysis$unadjusted.result[1,2])/qnorm(0.975)
      Z <- analysis$unadjusted.result[1,1]/SE
      results <- c(
        analysis$RMST.arm1$result[1,1],analysis$RMST.arm1$result[1,2],
        analysis$RMST.arm0$result[1,1],analysis$RMST.arm0$result[1,2],
        analysis$unadjusted.result[1,1], SE, Z, pnorm(-Z)
      )
    },error=function(e){results <- rep(NA,8)}
  )
  return(results)
}

###############################################################################
# RMST_cov_analysis ; Helper function for RMST analysis in analyse_sim
# Date: 15/11/18
# Author: James Bell, james.bell.ext@boehringer-ingelheim.com
#
# Requires the 'survRM2' package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  restriction: positive value for restriction time for RMST analysis
#  stratum: name of stratum column
#
# Note: This returns the SE, rather than the CI
###############################################################################
RMST_cov_analysis <- function(i,data,restriction,stratum){
  iterdata <- data[[i]]
  results <- rep(NA,8)
  tryCatch(
    {
      analysis <- rmst2(time=iterdata[,"Time"],status=1-iterdata[,"Censored"],arm=iterdata[,"Trt"]-1,tau=restriction,alpha=0.05,covariates=iterdata[,stratum])
      SE <- abs(analysis$adjusted.result[1,1]-analysis$adjusted.result[1,2])/qnorm(0.975)
      Z <- analysis$adjusted.result[1,1]/SE
      arm1 <- analysis$adjusted.result[1,1]*analysis$adjusted.result[2,1]/(analysis$adjusted.result[2,1]-1)
      arm0 <- arm1 - analysis$adjusted.result[1,1]
      results <- c(
        arm1,NA,
        arm0,NA,
        analysis$adjusted.result[1,1], SE, Z, pnorm(-Z)
      )
    },error=function(e){results <- rep(NA,8)}
  )
  return(results)
}


###############################################################################
# LM_analysis ; Helper function for landmark analysis in analyse_sim
#
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  landmarktime: positive value for landmark analysis time
###############################################################################
LM_analysis <- function(i,data,landmarktime){
  iterdata <- data[[i]]
  max1 <- max(iterdata[iterdata[,"Trt"]== 1,"Time"])
  max2 <- max(iterdata[iterdata[,"Trt"]== 2,"Time"])
  nogo1 <- 0
  nogo2 <- 0
  if(max1<landmarktime){
    S1 <- NA
    SE1 <- NA
    nogo1 <- 1
  }
  if(max2<landmarktime){
    S2 <- NA
    SE2 <- NA
    if(nogo1 == 1){return(c(S1,SE1,S2,SE2))}
    nogo2 <- 1
  }

  survival <- summary(survfit(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"],error="greenwood"))
  survivala <- survival[c("time","surv","std.err","strata")]
  if(nogo1==0){
    survival1a <- sapply(survivala, "[",survival$strata=='iterdata[, "Trt"]=1')
    if(is.vector(survival1a)){survival1a <- t(as.matrix(survival1a))}
    if(is.na(survival1a[1,"time"]>0) | survival1a[1,"time"] > landmarktime){
      S1 <- 1
      SE1 <- 0
    }else{
        survival1a <- survival1a[survival1a[,"time"] <= landmarktime,,drop=FALSE]
        cutoff1 <- which.max(survival1a[,"time"])
        S1 <- survival1a[cutoff1,"surv"]
        SE1 <- survival1a[cutoff1,"std.err"]
    }
  }
  if(nogo2==0){
    survival2a <- sapply(survivala, "[",survival$strata=='iterdata[, "Trt"]=2')
    if(is.vector(survival2a)){survival2a <- t(as.matrix(survival2a))}
    if(is.na(survival2a[1,"time"]>0) | survival2a[1,"time"] > landmarktime){
      S2 <- 1
      SE2 <- 0
    }else{
      survival2a <- survival2a[survival2a[,"time"] <= landmarktime,,drop=FALSE]
      cutoff2 <- which.max(survival2a[,"time"])
      S2 <- survival2a[cutoff2,"surv"]
      SE2 <- survival2a[cutoff2,"std.err"]
    }
  }
  results <- c(S2,SE2,S1,SE1)
  return(results)
}

###############################################################################
# LM_cov_analysis ; Helper function for landmark analysis in analyse_sim
#
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  landmarktime: positive value for landmark analysis time
#  stratum: name of the stratum column
###############################################################################
LM_cov_analysis <- function(i,data,landmarktime,stratum){
  iterdata <- data[[i]]
  trtsplit <- split.data.frame(iterdata,iterdata[,"Trt"])
  trtsplit1 <- split.data.frame(trtsplit[[1]],trtsplit[[1]][,stratum])
  trtsplit2 <- split.data.frame(trtsplit[[2]],trtsplit[[2]][,stratum])
  length1 <- length(trtsplit1)
  length2 <- length(trtsplit2)
  max1 <- S1 <- SE1 <- rep(NA,length1)
  max2 <- S2 <- SE2 <- rep(NA,length2)
  survival1 <- vector("list",length1)
  survival2 <- vector("list",length2)

  for(j in 1:length1){
    max1[j] <- max(trtsplit1[[j]][,"Time"])
    if(max1[j] >= landmarktime){
      survival1[[j]] <- sapply(summary(survfit(Surv(Time,1-Censored)~ 1,error="greenwood",data=data.frame(trtsplit1[[j]])))[c("time","surv","std.err")],"[")
      if(is.vector(survival1[[j]])){survival1[[j]] <- t(as.matrix(survival1[[j]]))}
      if(is.na(survival1[[j]][1,"time"]>0) | survival1[[j]][1,"time"] > landmarktime){
        S1[j]  <- 1
        SE1[j] <- 0
      }else{
        temp   <- survival1[[j]][survival1[[j]][,"time"] <= landmarktime,,drop=FALSE]
        cutoff <- which.max(temp[,"time"])
        S1[j]  <- temp[cutoff,"surv"]
        SE1[j] <- temp[cutoff,"std.err"]
      }
    }
  }

  for(j in 1:length2){
    max2[j] <- max(trtsplit2[[j]][,"Time"])
    if(max2[j] >= landmarktime){
      survival2[[j]] <- sapply(summary(survfit(Surv(Time,1-Censored)~ 1,error="greenwood",data=data.frame(trtsplit2[[j]])))[c("time","surv","std.err")],"[")
      if(is.vector(survival2[[j]])){survival2[[j]] <- t(as.matrix(survival2[[j]]))}
      if(is.na(survival2[[j]][1,"time"]>0) | survival2[[j]][1,"time"] > landmarktime){
        S2[j] <- 1
        SE2[j] <- 0
      }else{
        temp <- survival2[[j]][survival2[[j]][,"time"] <= landmarktime,,drop=FALSE]
        cutoff <- which.max(temp[,"time"])
        S2[j]  <- temp[cutoff,"surv"]
        SE2[j] <- temp[cutoff,"std.err"]
      }
    }
  }

  V1   <- SE1^2
  V2   <- SE2^2
  S1T  <- sum(S1/V1)/sum(1/V1)
  S2T  <- sum(S2/V2)/sum(1/V2)
  SE1T <- sqrt(1/sum(1/V1))
  SE2T <- sqrt(1/sum(1/V2))
  results <- c(S2T,SE2T,S1T,SE1T)
  return(results)
}

#######################################################################################################
#'Summarise analyses of simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for summarising the analyses of simulated time-to-event trial data produced by analyse_sim.\cr
#' Automatically reads in format from analyse_sim; no other input format is supported.\cr
#' Automatically detects types of analysis performed and provides relevant summaries (log-rank, Cox, RMST, landmark).\cr
#' @param analysed_results Output file from analyse_sim. Only analyse_sim output is supported.
#' @param alpha1 1-sided alpha to be used for analysis. Default=0.025.
#' @return Returns a table with one row. Table contains the following columns:
#' \itemize{
#'  \item{"Simulations"}{ Number of simulations conducted}
#'  \item{"HR"}{ Exponent of Mean Cox Log Hazard Ratio (LR/Cox analysis only)}
#'  \item{"LogHR"}{ Mean Cox Log Hazard Ratio (LR/Cox analysis only)}
#'  \item{"LogHR_SE"}{ Root mean square of the Cox Standard Errors for Log Hazard Ratio (LR/Cox analysis only)}
#'  \item{"HR_Z"}{ Mean Cox Z-Score (LR/Cox analysis only)}
#'  \item{"HR_P"}{ p-value of Mean Cox Z-Score (LR/Cox analysis only)}
#'  \item{"HR_Power"}{ Simulated power of Cox-regression (LR/Cox analysis only)}
#'  \item{"HR_Failed"}{ Proportion of simulations failing to calculate a Cox HR (LR/Cox analysis only)}
#'  \item{"LR_Z"}{ Mean Log-Rank Test Z-Score (LR/Cox analysis only)}
#'  \item{"LR_P"}{ p-value of Mean Log-Rank Test Z-Score (LR/Cox analysis only)}
#'  \item{"LR_Power"}{ Simulated power of the log-rank test (LR/Cox analysis only)}
#'  \item{"LR_Failed"}{ Proportion of simulations failing to calculate a log-rank test statistic (LR/Cox analysis only)}
#'  \item{"Events_Active"}{ Mean events in active arm (LR/Cox analysis only)}
#'  \item{"Events_Control"}{ Mean events in control arm (LR/Cox analysis only)}
#'  \item{"Events_Total"}{ Mean total events(LR/Cox analysis only)}
#'  \item{"RMST_Time"}{ Restriction time for RMST analysis (RMST analysis only)}
#'  \item{"RMST_Control"}{ Mean RMST for active arm (RMST analysis only)}
#'  \item{"RMST_C_SE"}{ Root mean square of RMST Standard Errors for active arm (RMST analysis only)}
#'  \item{"RMST_Active"}{ Mean RMST for control arm (RMST analysis only)}
#'  \item{"RMST_A_SE"}{ Root mean square of RMST Standard Errors for control arm (RMST analysis only)}
#'  \item{"RMST_Delta"}{ Mean RMST difference between arms active-control (RMST analysis only)}
#'  \item{"RMST_D_SE"}{ Root mean square of RMST difference Standard Errors (RMST analysis only)}
#'  \item{"RMST_Power"}{ Simulated power of RMST (RMST analysis only)}
#'  \item{"RMST_Failed"}{ Proportion of simulations failing to calculate the RMST (RMST analysis only)}
#'  \item{"LM_Time"}{ Landmark analysis time, i.e. assessment time of Survival function (Landmark analysis only)}
#'  \item{"LM_Control"}{ Mean survival function for active arm at landmark time (Landmark analysis only)}
#'  \item{"LM_C_SE"}{ Root mean square of Greenwood standard errors for active arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Active"}{ Mean survival function for control arm at landmark time (Landmark analysis only)}
#'  \item{"LM_A_SE"}{ Root mean square of Greenwood standard errors for control arm at landmark time (Landmark analysis only)}
#'  \item{"LM_Delta"}{ Mean survival function difference active-control at landmark time (Landmark analysis only)}
#'  \item{"LM_D_SE"}{ Root mean square of Greenwood standard errors for survival differences at landmark time (Landmark analysis only)}
#'  \item{"LM_Power"}{ Power of landmark analysis (Landmark analysis only)}
#'  \item{"LM_Failed"}{ Proportion of simulations failing to calculate the survival difference at landmark time (Landmark analysis only)}
#' }
#' @author James Bell
#' @examples example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=40, iterations=5,seed=12345,detailed_output=TRUE)
#'
#' example_analysis1 <- analyse_sim(example_sim)
#' example_analysis2 <- analyse_sim(data=example_sim,RMST=30,landmark=30)
#'
#' example_summary1 <- summarise_analysis(example_analysis1)
#' example_summary2 <- summarise_analysis(example_analysis2)
#' @export
summarise_analysis <- function(analysed_results,alpha1=0.025){
  #Check argument existence
  if(missing(analysed_results))stop("Error: Please specify the name of the analysis results dataset ('analysed_results').")
  if(!is.numeric(alpha1) || alpha1 <= 0 || alpha1 > 0.5 ) stop("Please specify a positive one-sided alpha less than 0.5 using the 'alpha1' argument.")

  #Set default to not analyse anything
  LR <- FALSE
  RMST <- FALSE
  LM <- FALSE
  #Detect key columns to identify what summarising needs doing
  #Other columns required for each analysis will then be checked to ensure all is present
  if("LogHR" %in% colnames(analysed_results)){LR <- TRUE}
  if("RMST_Delta" %in% colnames(analysed_results)){RMST <- TRUE}
  if("LM_Z" %in% colnames(analysed_results)){LM <- TRUE}

  # Exit nicely if it can't detect suitable input data
  if(!(LR|RMST|LM)){stop("No suitable input analysed data identified from column headings. Please ensure input data is from analyse_sim.")}

  Simulations <- nrow(analysed_results)
  out <- cbind(Simulations)

  if(LR){
    if(!("LogHR_SE" %in% colnames(analysed_results))){stop("Analysis aborted during LR/HR stage: could not find LogHR_SE column")}
    if(!("HR_Z" %in% colnames(analysed_results))){stop("Analysis aborted during LR/HR stage: could not find HR_Z column")}
    if(!("LR_Z" %in% colnames(analysed_results))){stop("Analysis aborted during LR/HR stage: could not find LR_Z column")}
    if(!("Events_Active" %in% colnames(analysed_results))){stop("Analysis aborted during LR/HR stage: could not find Events_Active column")}
    if(!("Events_Control" %in% colnames(analysed_results))){stop("Analysis aborted during LR/HR stage: could not find Events_Control column")}

    LogHR    <- mean(analysed_results[!is.na(analysed_results[,"LogHR"]),"LogHR"])
    HR       <- round(exp(LogHR),4)
    LogHR    <- round(LogHR,4)
    LogHR_SE <- round(sqrt(sum(analysed_results[!is.na(analysed_results[,"LogHR_SE"]),"LogHR_SE"]^2)/Simulations),5)

    HR_Z <- mean(analysed_results[!is.na(analysed_results[,"HR_Z"]),"HR_Z"])
    HR_P <- round(pnorm(HR_Z,0,1),4)
    HR_Z <- round(HR_Z,4)
    LR_Z <- mean(analysed_results[,"LR_Z"])
    LR_P <- round(pnorm(LR_Z,0,1),4)
    LR_Z <- round(LR_Z,4)

    Events_Active  <- mean(analysed_results[,"Events_Active"])
    Events_Control <- mean(analysed_results[,"Events_Control"])
    Events_Total   <- round(Events_Active+Events_Control,3)
    Events_Active  <- round(Events_Active,4)
    Events_Control <- round(Events_Control,4)

    ptempLR <- analysed_results[,"LR_Z"] < qnorm(alpha1,0,1)
    ptempHR <- analysed_results[,"HR_Z"] < qnorm(alpha1,0,1)
    failuresLR <- is.na(ptempLR)
    failuresHR <- is.na(ptempHR)
    ptempLR[failuresLR] <- FALSE
    ptempHR[failuresHR] <- FALSE
    LR_Failed <- round(sum(failuresLR)/Simulations,5)
    HR_Failed <- round(sum(failuresHR)/Simulations,5)
    LR_Power <- round(sum(ptempLR)/Simulations,4)
    HR_Power <- round(sum(ptempHR)/Simulations,4)
    #Mean_Assess_Time <- round(mean(analysed_results[!is.na(analysed_results[,"Assessment_Time"]),"Assessment_Time"]),2)

    #Return the analysis summary
    out <- cbind(out,HR,LogHR,LogHR_SE,HR_Z,HR_P,HR_Power,HR_Failed,LR_Z,LR_P,LR_Power,LR_Failed,Events_Active,Events_Control,Events_Total)
  }

  if(RMST){
    if(!("RMST_Restrict" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_Restrict column")}
    if(!("RMST_Control" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_Control column")}
    if(!("RMST_C_SE" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_C_SE column")}
    if(!("RMST_Active" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_Active column")}
    if(!("RMST_A_SE" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_A_SE column")}
    if(!("RMST_Delta" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_Delta column")}
    if(!("RMST_D_SE" %in% colnames(analysed_results))){stop("Analysis aborted during RMST stage: could not find RMST_D_SE column")}

    RMST_Restrict <- analysed_results[1,"RMST_Restrict"]


    RMST_Active   <- round(mean(analysed_results[!is.na(analysed_results[,"RMST_Active"]), "RMST_Active"]),4)
    RMST_Control  <- round(mean(analysed_results[!is.na(analysed_results[,"RMST_Control"]),"RMST_Control"]),4)
    RMST_Delta    <- round(mean(analysed_results[!is.na(analysed_results[,"RMST_Delta"]),  "RMST_Delta"]),4)

    RMST_A_SE     <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"RMST_A_SE"]),"RMST_A_SE"]^2)),4)
    RMST_C_SE     <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"RMST_C_SE"]),"RMST_C_SE"]^2)),4)
    RMST_D_SE     <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"RMST_D_SE"]),"RMST_D_SE"]^2)),4)

    ptempR <- analysed_results[,"RMST_Z"] > qnorm(1-alpha1,0,1)
    failuresR <- is.na(ptempR)
    ptempR[failuresR] <- FALSE
    RMST_Failed <- round(sum(failuresR)/Simulations,4)
    RMST_Power  <- round(sum(ptempR)   /Simulations,4)

    out <- as.data.frame(cbind(out,RMST_Restrict,RMST_Active,RMST_A_SE,RMST_Control,RMST_C_SE,RMST_Delta,RMST_D_SE,RMST_Power,RMST_Failed))
  }

  if(LM){
    if(!("LM_Time" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_Time column")}
    if(!("LM_Active" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_Active column")}
    if(!("LM_A_SE" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_A_SE column")}
    if(!("LM_Control" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_Control column")}
    if(!("LM_C_SE" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_C_SE column")}
    if(!("LM_Delta" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_Delta column")}
    if(!("LM_D_SE" %in% colnames(analysed_results))){stop("Analysis aborted during Landmark stage: could not find LM_D_SE column")}

    LM_Time <- analysed_results[1,"LM_Time"]

    LM_Active  <- round(mean(analysed_results[!is.na(analysed_results[,"LM_Active"]), "LM_Active"]),4)
    LM_Control <- round(mean(analysed_results[!is.na(analysed_results[,"LM_Control"]),"LM_Control"]),4)
    LM_Delta   <- round(mean(analysed_results[!is.na(analysed_results[,"LM_Delta"]),  "LM_Delta"]),4)

    LM_A_SE <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"LM_A_SE"]),"LM_A_SE"]^2)),4)
    LM_C_SE <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"LM_C_SE"]),"LM_C_SE"]^2)),4)
    LM_D_SE <- round(sqrt(mean(analysed_results[!is.na(analysed_results[,"LM_D_SE"]),"LM_D_SE"]^2)),4)

    ptempLM <- analysed_results[,"LM_Z"] > qnorm(1-alpha1,0,1)
    failuresLM <- is.na(ptempLM)
    ptempLM[failuresLM] <- FALSE
    LM_Failed <- round(sum(failuresLM)/Simulations,4)
    LM_Power <- round(sum(ptempLM)/Simulations,4)

    #Return the analysis summary
    out <- cbind(out,LM_Time,LM_Active,LM_A_SE,LM_Control,LM_C_SE,LM_Delta,LM_D_SE,LM_Power,LM_Failed)
  }
  return(out)
}
