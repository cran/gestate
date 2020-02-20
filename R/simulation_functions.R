if(getRversion() >= "2.15.1")  utils::globalVariables(c("i","ETime","Rec_Time","Assess","RCTime","CTime","Iter"))

#######################################################################################################
#'Perform simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for simulating generalised two-arm time-to-event trial data for NPH trials with arbitrary event, censoring and recruitment distributions.\cr
#' Event and censoring distributions are specified via Curve objects, with recruitment specified through an RCurve object.
#' As it uses same architecture and similar syntax to nph_traj(), analysis results ought to be directly comparable.
#' It is designed to complement nph_traj(), either as a stochastic alternative, or as a means to validate its outputs. 
#' It can also be used to build more complex simulations by combining the outputs of multiple runs; e.g. multi-arm trials.\cr
#' Data sets created by this function are formatted so they may be automatically recognised and analysed by analyse_sim().\cr
#' @param active_ecurve Event distribution for the active arm, specified as a Curve object
#' @param control_ecurve Event distribution for the control arm, specified as a Curve object
#' @param active_dcurve Dropout/censoring distribution for the active arm, specified as a Curve object. By default, Blank(), i.e. no dropout.
#' @param control_dcurve Dropout/censoring distribution for the control arm, specified as a Curve object. By default, Blank(), i.e. no dropout.
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
  columns <- 9
  nactive <- ceiling(getNactive(rcurve))
  ncontrol <- floor(getNcontrol(rcurve))
  n <- nactive+ncontrol

  if(nactive > 0){
    asize <- nactive*iterations
    activematrix <- matrix(rep(0,asize*columns),nrow=asize)
    # Create events columns
    activematrix[,5] <- random_draw(active_ecurve,asize)
    # Create censoring columns
    activematrix[,6] <- random_draw(active_dcurve,asize)
    # Create treatment columns
    activematrix[,3] <- 2
    # Creates iterations columns
    activematrix[,4] <- rep(1:iterations, each=nactive)
  }

  if(ncontrol > 0){
    csize <- ncontrol*iterations
    controlmatrix <- matrix(rep(0,csize*columns),nrow=csize)
    # Create events columns
    controlmatrix[,5] <- random_draw(control_ecurve,csize)
    # Create censoring columns
    controlmatrix[,6] <- random_draw(control_dcurve,csize)
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
  output[,"Rec_Time"] <- random_draw(rcurve,n*iterations)
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
#' As it uses same architecture and similar syntax to nph_traj(), results in simple cases may be directly comparable to e.g. use of MixExp() or MixWei() Curves.\cr
#' Can be used to validate outputs from nph_traj().\cr
#' Data sets from this are set up to be automatically analysed with the analyse_sim function (including stratified analysis if you provide it the name of stratum column).\cr
#' @param stratum_probs Vector of probabilities that patients belong to each stratum. Must sum to 1. Its length determines the number of strata.
#' @param active_ecurve List of event distributions for the active arm, specified as a list of Curve objects. If single Curve is specified, will be used for all strata.
#' @param control_ecurve List of event distributions for the control arm, specified as a list of Curve objects. If single Curve is specified, will be used for all strata.
#' @param active_dcurve List of dropout/censoring distribution for the active arm, specified as a Curve object. If single Curve is specified, will be used for all strata. By default, Blank(), i.e. no dropout in any stratum.
#' @param control_dcurve List of dropout/censoring distribution for the control arm, specified as a Curve object. If single Curve is specified, will be used for all strata. By default, Blank(), i.e. no dropout in any stratum.
#' @param rcurve Recruitment distribution, specified as a single RCurve object.
#' @param assess Positive number for the assessment time at which administrative censoring will be performed.
#' @param fix_events Positive integer for the number of events to fix (if required), letting the assessment time vary. Alternatively, NULL for fixed time assessment with variable event numbers. Notes: Fixing event numbers overrides any specified assessment time and slows simulation considerably. Default = NULL (fixed analysis time)
#' @param stratum_name Name of the column defining the stratum. Default="Stratum".
#' @param iterations Number of simulations to perform. Depending on trial size, 10,000-20,000 is typically OK to analyse on 8GB RAM.
#' @param seed Seed number to use. Numerical, although if "Rand" is specified, a system-time-derived number will be used.
#' @param detailed_output Boolean to require full details of timings of competing processes. If FALSE, the simplified data only includes the *'ed output columns - this approximately halves RAM requirements. Default=FALSE (simplified).
#' @param output_type "matrix" or "list" specifying the type of output required. "matrix" requests a single matrix with a column "iter" to denote the simulation, while "list" creates a list with one entry per simulation. Default="matrix".
#' @import survival
#' @import foreach
#' @importFrom doParallel registerDoParallel
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
# Fix a set number of events if specified
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
#' Function for converting trials simulated by simulate_trials() or simulate_trials_strata() from a fixed censoring time to a fixed number of total events.\cr
#' It is set up to automatically read in either matrix or list formats from simulate_trials() or simulate_trials_strata(), and only these inputs are supported.\cr
#' Note that if recruitment had not finished in the input then any increases in assessment time cannot account for the missing patients.
#' It is therefore strongly recommended to initially simulate for at least the duration of the recruitment before fixing the event number.\cr
#' This function can also be used to change format and/or slim down data for event-driven simulations.\cr
#' @param data Output file from simulate_trials() or simulate_trials_strata() in either "list" or "matrix" format. Only these formats are supported.
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
      warning("Only",sortE,"events are possible in simulation",i,"but ",events," have been requested. Setting assessment time to just after last event.\n")
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
#' Function for modifying the assessment time of simulate_trials() or simulate_trials_strata() simulations.\cr
#' It is set up to automatically read in either matrix or list formats from simulate_trials() or simulate_trials_strata(), and only these inputs are supported.\cr
#' Note that if recruitment had not finished in the input then any increases in assessment time cannot account for the missing patients.
#' It is therefore strongly recommended to initially simulate for at least the duration of the recruitment before reducing the number to missing patients.\cr
#' This function can also be used to change format and/or slim down data for time-driven simulations.\cr
#' @param data Output file from simulate_trials() or simulate_trials_strata() in either "list" or "matrix" format. Only these formats are supported.
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
#' Function for analysing simulated time-to-event trial data produced by simulate_trials() or simulate_trials_strata().\cr
#' This function automatically reads in either list or matrix simulate_trials() data formats.
#' It performs log rank test and Cox regression analysis by default, but can also/instead perform RMST and/or landmark analyses.
#' Covariate adjusted (/ stratified) analysis may be selected by using the "stratum" argument.
#' If a stratum is specified, it will be included as a covariate in Cox and RMST analysis, and as a stratum in a stratified log-rank test and an inverse-precision-weighted landmark test.
#' Strata values are handled as factors, so continuous covariates are not supported.\cr
#' Analysis is typically the slowest part of simulation studies, so parallel processing using the doParallel package is built in.
#' Parallel processing is enabled by setting the number of cores in the "parallel_cores" argument.
#' Use of parallel processing is recommended for largescale (e.g. 100,000 iteration) simulations. 
#' To avoid unnecessary issues, ensure that the number of cores specified does not exceed number of threads provided by hardware.
#' @param data Output file from simulate_trials(). Only simulate_trials() or simulate_trials_strata() output is supported, in either "list" or "matrix" format.
#' @param LR Requests log-rank test and Cox regression. Default=TRUE
#' @param RMST Requests Restricted Mean Survival Time analysis with specified (positive integer) restriction time, leave NULL for no analysis. Default=NULL (no RMST analysis).
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
#' @references Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D, Takeuchi M, Uyama Y, Zhao L,
#'   Skali H, Solomon S, Jacobus S, Hughes M, Packer M, Wei LJ. Moving beyond the hazard ratio in
#'   quantifying the between-group difference in survival analysis. Journal of clinical Oncology 2014,32, 2380-2385.
#' Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subjects baseline covariates in survival analysis. 
#'   Biostatistics 2014, 15, 222-233.  
#' @examples example_sim <- simulate_trials(active_ecurve=Weibull(250,0.8),control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100), assess=20,iterations=100,seed=12345,detailed_output=TRUE)
#'
#' example_analysis1 <- analyse_sim(example_sim)
#' example_analysis2 <- analyse_sim(data=example_sim,RMST=15,landmark=15)
#'
#' example_strat_sim <- simulate_trials_strata(stratum_probs=c(0.5,0.5),
#' active_ecurve=c(Weibull(250,0.8),Weibull(100,1)), control_ecurve=Weibull(100,1),
#' rcurve=LinearR(12,100,100),assess=20,iterations=100,seed=12345)
#'
#' example_strat_analysis <- analyse_sim(data=example_strat_sim,RMST=15,landmark=15,stratum="Stratum")
#' @export
analyse_sim <- function(data, LR = TRUE, RMST = NA, landmark = NA, stratum="",parallel_cores=1){
# Check that the correct packages are installed for specific options.
# If not, disable options
  if(parallel_cores!= 1 && !requireNamespace("doParallel", quietly = TRUE)) {
    warning("Parallel processing requested but package \"doParallel\" needed for this. Please install it to enable parallel processing. Reverting to single-thread processing.\n")
    parallel_cores <- 1
  }

  if(!LR && is.na(RMST) && is.na(landmark)){"At least one of 'LR', 'RMST' or 'landmark' must be TRUE to specify analysis to be performed"}

 #Required Input Check
  if(missing(data))stop("Error: Please specify the name of the simulated data set ('data').")

  #Boolean
  if(!is.logical(LR))stop("Error: 'LR' argument must be boolean: default=FALSE (perform log-rank and cox analysis).")
  #Character
  if(!is.character(stratum))stop("Error: 'stratum' argument must be a string: default=FALSE (perform log-rank and cox analysis).")

  #Positive integer
  if(!is.numeric(parallel_cores) || parallel_cores%%1!=0 || parallel_cores < 1 ) stop("Error: 'parallel_cores' argument is not a positive integer.")
  if(is.null(landmark)){landmark <- NA}
  if(is.null(RMST)){RMST <- NA}
  if(!is.na(landmark) && (landmark%%1!=0 || landmark < 1 )) stop("Landmark analysis specified, but landmark time is not a positive integer.")
  if(!is.na(RMST) && (RMST%%1!=0 || RMST < 1 )) stop("RMST analysis specified, but restriction time is not a positive integer.")

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
        out <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages=c("survival"))%dopar% LR_analysis(i,data)
      } else {
        out <- foreach(i=1:maxiter,.combine="rbind",.export="LR_cov_analysis", .packages=c("survival"))%dopar% LR_cov_analysis(i,data,stratum)
      }
    } else{
      if(stratum==""){
        out <- foreach(i=1:maxiter,.combine="rbind",.export="LR_analysis", .packages="survival")%do% LR_analysis(i,data)
      } else {
        out <- foreach(i=1:maxiter,.combine="rbind",.export="LR_cov_analysis", .packages="survival")%do% LR_cov_analysis(i,data,stratum)
      }
    }
    if(!is.matrix(out)){out <- t(out)}
    colnames(out) <- c("HR","LogHR","LogHR_SE","HR_Z","HR_P","LR_Z","LR_P","Events_Active","Events_Control")

  #Note that at this stage: LogHR_SE column is the variance
    out[,"HR"] <- exp(out[,"LogHR"])
    out[,"LogHR_SE"] <- sqrt(out[,"LogHR_SE"])
    out[,"HR_Z"] <- out[,"LogHR"]/out[,"LogHR_SE"]
    out[,"HR_P"] <- pnorm(out[,"HR_Z"])
    out[,"LR_P"] <- pnorm(out[,"LR_Z"])
  }

  # If RMST analysis required, do it. Separate calls depending on whether parallel processing enabled or stratified analysis requested
  if(!is.na(RMST) || !is.na(landmark)){
    if(parallel_cores>1){
      if(stratum==""){
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_LM_analysis")%dopar% RMST_LM_analysis(i,data,RMST,landmark)
      } else {
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_LM_cov_analysis")%dopar% RMST_LM_cov_analysis(i,data,RMST,landmark,stratum)
      }
    }else{
      if(stratum==""){
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_LM_analysis")%do% RMST_LM_analysis(i,data,RMST,landmark)
      } else {
        outRMST <- foreach(i=1:maxiter,.combine="rbind",.export="RMST_LM_cov_analysis")%do% RMST_LM_cov_analysis(i,data,RMST,landmark,stratum)
      }
    }
    if(!is.matrix(outRMST)){outRMST <- t(outRMST)}
    colnames(outRMST) <- c("RMST_Restrict","RMST_Active","RMST_A_SE","RMST_Control","RMST_C_SE","RMST_Delta","RMST_D_SE","RMST_Z","RMST_P","LM_Time","LM_Active","LM_A_SE","LM_Control","LM_C_SE","LM_Delta","LM_D_SE","LM_Z","LM_P")
    outRMST[,"RMST_Z"] <- outRMST[,"RMST_Delta"]/outRMST[,"RMST_D_SE"]
    outRMST[,"RMST_P"] <- pnorm(-outRMST[,"RMST_Z"])
    outRMST[,"LM_Z"] <- outRMST[,"LM_Delta"]/outRMST[,"LM_D_SE"]
    outRMST[,"LM_P"] <- pnorm(-outRMST[,"LM_Z"])
    if(is.na(RMST)){outRMST <- outRMST[,10:18]}
    if(is.na(landmark)){outRMST <- outRMST[,1:9]}
    out <- if(LR){cbind(out,outRMST)} else{outRMST}
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
  fit <- coxph.fit(x=matrix(iterdata[,"Trt"],ncol=1), y=Surv(iterdata[,"Time"],1-iterdata[,"Censored"]), strata=NULL, offset=NULL, init=NULL, control=coxph.control(), method = "breslow",rownames=1:nrow(iterdata))
  output <- c(NA,fit$coefficients,fit$var,c(NA,NA),sqrt(fit$score)*sign(fit$coefficients[1]),NA,length(iterdata[(iterdata[,"Trt"]==2 & iterdata[,"Censored"]==0),1]),length(iterdata[(iterdata[,"Trt"]==1 & iterdata[,"Censored"]==0),1]))
}

################################################################################
# LR_cov_analysis ; Helper function for stratified log-rank/covariate-adjusted Cox analysis in analyse_sim
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  cov:     covariate column name
################################################################################
LR_cov_analysis <- function(i,data,cov){
  iterdata <- data[[i]]
  covs <- length(cov)
  analysis <- coxph(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"] + as.factor(iterdata[,cov]),ties="breslow")
  analysisLR <- survdiff(Surv(iterdata[,"Time"],1-iterdata[,"Censored"])~ iterdata[,"Trt"] + strata(iterdata[,cov]))
  output <- c(NA,coef(summary(analysis))[1,1],coef(summary(analysis))[1,3]^2,c(NA,NA),sqrt(analysisLR$chisq)*sign(sum(analysisLR$exp[1,])-sum(analysisLR$obs[1,])),NA,length(iterdata[(iterdata[,"Trt"]==2 & iterdata[,"Censored"]==0),1]),length(iterdata[(iterdata[,"Trt"]==1 & iterdata[,"Censored"]==0),1]))
}

###############################################################################
# RMST_LM_analysis ; Helper function for RMST & LM analysis in analyse_sim
# Requires the 'survival' default R package to be pre-loaded
# Combining LM and RMST into single function allows for considerable computational savings when both are required.
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  restriction: positive value for restriction time for RMST analysis (NA for do not calculate)
#  landmark: positive value for landmark analysis time (NA for do not calculate)
###############################################################################
RMST_LM_analysis <- function(i,data,restriction=NA,landmark=NA){


  #Inner function for calculating rmst and lm with no covariates.
  #By combining rmst and lm calculations, avoid duplicating rate-limiting tabulation step
  rmst_lm <- function (time, status, arm, restriction = NA,landmark=NA){

    #Fast inner function for RMST calculation
    rmst_fast <- function (tab, restriction){
      index     <- tab$time <= restriction
      nrisk     <- tab$n.risk[index]
      nevent    <- tab$n.event[index]
      areas     <- diff(c(0, sort(c(tab$time[index], restriction)))) * c(1, tab$surv[index])
      rmst      <- sum(areas)
      var1      <- rep(0,length(nrisk))
      var1[(nrisk - nevent) != 0] <- nevent/(nrisk * (nrisk - nevent))
      rmst_var  <- sum(cumsum(rev(areas[-1]))^2 * rev(c(var1, 0))[-1])
      rmst_se   <- sqrt(rmst_var)
      return(c(rmst, rmst_se, rmst_var))
    }

    # Fast function for KM estimate production. Rate-limiting-step is the table call.
    survfit_fast <- function (time,events){
      n      <- length(time)
      temp   <- table(factor(time),factor(events,levels=0:1))
      nevent <- as.vector(temp[, 2])
      ncens  <- as.vector(temp[, 1])
      times  <- as.numeric(row.names(temp))
      nrisk  <- n - c(0,cumsum(nevent+ncens)[-n])
      trisk  <- ifelse(nrisk == 0, 1, nrisk)
      surv   <- cumprod((trisk - nevent)/trisk)
      se     <- surv*sqrt(cumsum(nevent/(trisk * (trisk - nevent))))
      return(list(n = n, time = times, n.risk = nrisk, n.event = nevent, n.censor = ncens, surv = surv,se = se))
    }

    lm_analysis <- function(tab,landmark){
      if(is.na(tab$time[1]>0) | tab$time[1] > landmark){
        S <- 1
        SE <- 0
      }else{
        cutoff <- which.max(tab$time[tab$time <= landmark])
        S <- tab$surv[cutoff]
        SE <- tab$se[cutoff]
      }
      return(c(S,SE))
    }

    maxT0 <- max(time[arm == 0])
    maxT1 <- max(time[arm == 1])

    tab1 <- survfit_fast(time[arm == 1],status[arm == 1])
    tab0 <- survfit_fast(time[arm == 0],status[arm == 0])
 
    if(!is.na(restriction) && maxT0 >= restriction){
      wk0 <- rmst_fast(tab0, restriction)
    } else {
      wk0 <- c(NA,NA)
    }
    if(!is.na(restriction) && maxT1 >= restriction){
      wk1 <- rmst_fast(tab1, restriction)
    } else {
      wk1 <- c(NA,NA)
    }
    if(!is.na(landmark) && maxT0 >= landmark){
      lm0 <- lm_analysis(tab0,landmark)
    } else {
      lm0 <- c(NA,NA)
    }
    if(!is.na(landmark) && maxT1 >= landmark){
      lm1 <- lm_analysis(tab1,landmark)
    } else {
      lm1 <- c(NA,NA)
    }

    Z <- list()
    Z$RMST.arm1 <- wk1[1:2]
    Z$RMST.arm0 <- wk0[1:2]
    Z$RMST.diff <- c(wk1[1] - wk0[1], sqrt(wk1[3] + wk0[3]))
    Z$LM.arm1 <- lm1[1:2]
    Z$LM.arm0 <- lm0[1:2]
    Z$LM.diff <- c(lm1[1] - lm0[1], sqrt(lm1[2]^2 + lm0[2]^2))
    return(Z)
  }

  iterdata <- data[[i]]
  analysis <- rmst_lm(time=iterdata[,"Time"], status=1-iterdata[,"Censored"], arm=iterdata[,"Trt"]-1, restriction = restriction,landmark=landmark)
  return(c(restriction,analysis$RMST.arm1,analysis$RMST.arm0,analysis$RMST.diff,c(NA,NA),landmark,analysis$LM.arm1,analysis$LM.arm0,analysis$LM.diff,c(NA,NA)))
}

###############################################################################
# RMST_cov_analysis ; Helper function for RMST analysis in analyse_sim
# Requires the 'survival' default R package to be pre-loaded
# Takes as input:
#  i:   	integer of simulation to be analysed
#  data:	data set produced by the simulate_trials function
#  restriction: positive value for restriction time for RMST analysis
#  landmark: positive value for landmark analysis time (NA for do not calculate)
#  stratum: name of stratum column
###############################################################################
RMST_LM_cov_analysis <- function(i,data,restriction=NA,landmark=NA,stratum){

  #Core function for calculating RMST with covariates. Supports class covariates. 
  rmstcov <- function (time, status, arm, restriction, covariates){
    #rmstcovreg is the inner function used for performing the regression for calculating RMST with covariates
    rmstcovreg <- function (y, status, x, arm, restriction){

        xmat <- cbind(1, x)
        ncovariates <- ncol(xmat)
        yv <- pmin(y, restriction)
        ev <- status
        ev[yv == restriction] <- 1
        
        isarm1 <- arm == 1
        yv1 <- yv[isarm1]
        yord1 <- order(yv1)
        yv1 <- yv1[yord1]
        ev1 <- ev[isarm1][yord1]
        xmat1  <- xmat[isarm1,][yord1,]
        n1  <- length(ev1)
        yv1table <- table(yv1)

        isarm0 <- arm == 0
        yv0 <- yv[isarm0]
        yord0 <- order(yv0)
        yv0 <- yv0[yord0]
        ev0 <- ev[isarm0][yord0]
        xmat0  <- xmat[isarm0,][yord0,]
        n0  <- length(ev0)
        yv0table <- table(yv0)

        w1 <- ev1/rep(survfit(Surv(yv1, 1 - ev1) ~ 1)$surv, yv1table)
        w0 <- ev0/rep(survfit(Surv(yv0, 1 - ev0) ~ 1)$surv, yv0table)

        Beta <- lm(c(yv1, yv0) ~ rbind(xmat1, xmat0) - 1, weights = c(w1, w0))$coef
        score1 <- xmat1*w1*(yv1-as.vector(xmat1 %*% Beta))
        score0 <- xmat0*w0*(yv0-as.vector(xmat0 %*% Beta))

        yv1revcumsum <- (n1+1)-as.vector(rev(cumsum(rev(yv1table))))
        yv1pos   <- rep(yv1revcumsum,yv1table)
        yv1loc <- n1+1-yv1pos
        yv1forward <- rep(as.vector(cumsum(yv1table)),yv1table)

        tab1 <-  p3temp21 <- matrix(rep(NA,n1*ncovariates),ncol=ncovariates)
        for(i in 1:ncovariates){
           tab1[,i]<- rev(cumsum(rev(score1[,i])))[yv1pos]
        }

        part2_1 <- tab1*(1-ev1)/yv1loc
        p3temp11 <- part2_1/yv1loc
        for(i in 1:ncovariates){
          p3temp21[,i] <- cumsum(p3temp11[,i]) 
        }
        part3_1 <- p3temp21[yv1forward,]
        arm1 <- score1+part2_1-part3_1

        yv0revcumsum <- (n0+1)-as.vector(rev(cumsum(rev(yv0table))))
        yv0pos   <- rep(yv0revcumsum,yv0table)
        yv0loc <- n0+1-yv0pos
        yv0forward <- rep(as.vector(cumsum(yv0table)),yv0table)

        tab0 <-  p3temp20 <- matrix(rep(NA,n0*ncovariates),ncol=ncovariates)
        for(i in 1:ncovariates){
           tab0[,i]<- rev(cumsum(rev(score0[,i])))[yv0pos]
        }

        part2_0 <- tab0*(1-ev0)/yv0loc
        p3temp10 <- part2_0/yv0loc

        for(i in 1:ncovariates){
          p3temp20[,i] <- cumsum(p3temp10[,i]) 
        }
        part3_0 <- p3temp20[yv0forward,]
        arm0 <- score0+part2_0-part3_0

        gamma <- t(arm1) %*% arm1 + t(arm0) %*% arm0
        X2 <- t(xmat) %*% xmat
        varBeta <- solve(X2) %*% gamma %*% solve(X2)

        return(list(Beta,varBeta))
    }
    #function to handle class covariates, unstacking them into binary variables
    unstacker <- function(covariate){
      levs <- levels(covariate)
      width <- length(levs)-1
      if(width<1){return(NULL)}
      output <- matrix(rep(0,length(covariate)*width),ncol=width)
      for(i in 1:width){
        output[,i] <- 1*(covariate==levs[i+1])
      }
      colnames(output) <- levs[-1]
      return(output)
    }

    if(is.null(ncol(covariates))|| ncol(covariates) == 1){
      if(is.factor(covariates)){
        covariates <- unstacker(covariates)
      }
    } else{
      size <- ncol(covariates)
      tempcovs <- vector("list", length = size)
      for(i in 1:size){
        if(is.factor(covariates[,i])){
          tempcovs[[i]] <- unstacker(covariates[,i])
        } else {tempcovs[[i]] <- covariates[,i]}
      }
      covariates <- do.call("cbind",tempcovs)
    }
    aa <- rmstcovreg(time, status, as.matrix(cbind(arm, covariates)), arm, restriction)
    Z <- list()
    covariates <- as.matrix(covariates)
    props <- colSums(covariates)/nrow(covariates)
    covrange <- if(is.null(covariates)){NULL} else{3:(2+ncol(covariates))}
    Z$RMST.arm1 <- c(aa[[1]][1] + aa[[1]][2] + sum(props*aa[[1]][covrange]), sqrt(sum(outer(c(1,1,props),c(1,1,props))*aa[[2]])))
    Z$RMST.arm0 <- c(aa[[1]][1]+ sum(props*aa[[1]][covrange]), sqrt(sum(outer(c(1,props),c(1,props))*aa[[2]][-2,-2])))
    Z$RMST.diff <- c(aa[[1]][2],sqrt(aa[[2]][2,2]))
    return(Z)
  }

  #Core landmark with covariates calculation function
  lm_cov <- function(time,status,arm,stratum,landmark){
    # Fast function for KM estimate production. Rate-limiting-step is the table call.
    survfit_fast <- function (time,events){
      n      <- length(time)
      temp   <- table(factor(time),factor(events,levels=0:1))
      nevent <- as.vector(temp[, 2])
      ncens  <- as.vector(temp[, 1])
      times  <- as.numeric(row.names(temp))
      nrisk  <- n - c(0,cumsum(nevent+ncens)[-n])
      trisk  <- ifelse(nrisk == 0, 1, nrisk)
      surv   <- cumprod((trisk - nevent)/trisk)
      se     <- surv*sqrt(cumsum(nevent/(trisk * (trisk - nevent))))
      temp   <- list(n = n, time = times, n.risk = nrisk, n.event = nevent, n.censor = ncens, surv = surv,se = se)
      return(temp)
    }

    lm_analysis <- function(ft,landmark){
      if(is.na(ft$time[1]>0) | ft$time[1] > landmark){
        S  <- 1
        SE <- 0
      }else{
        cutoff <- which.max(ft$time[ft$time <= landmark])
        S  <- ft$surv[cutoff]
        SE <- ft$se[cutoff]
      }
      return(c(S,SE))
    }

    if (!is.vector(stratum) && ncol(stratum)>1){
      store  <- colSums(t(stratum)*pi^((1:ncol(stratum))-1))
      stratum <- match(store,unique(store))
    }
 
    ustrat <- unique(stratum)
    nstrat <- length(ustrat)
    if(nstrat<1)stop("Error: No strata have been specified")
    lm0 <- lm1 <- matrix(rep(NA,2*nstrat),ncol=2)
    for(i in 1: nstrat){
      relevant0 <- arm == 0 & stratum == ustrat[i]
      relevant1 <- arm == 1 & stratum == ustrat[i]
      if(max(time[relevant0]) >= landmark){
        lm0[i,] <- lm_analysis(survfit_fast(time[relevant0],status[relevant0]),landmark)
      } else {lm0[i,] <- c(NA,NA)}
      if(max(time[relevant1]) >= landmark){
        lm1[i,] <- lm_analysis(survfit_fast(time[relevant1],status[relevant1]),landmark)
      } else {lm1[i,] <- c(NA,NA)}
    }

    V0   <- lm0[,2]^2
    V1   <- lm1[,2]^2

    S0   <- sum(lm0[,1]/V0)/sum(1/V0)
    S1   <- sum(lm1[,1]/V1)/sum(1/V1)
    SE0  <- sqrt(1/sum(1/V0))
    SE1  <- sqrt(1/sum(1/V1))

    Z <- list()
    Z$LM.arm1 <- c(S1,SE1)
    Z$LM.arm0 <- c(S0,SE0)
    Z$LM.diff <- c(S1-S0, sqrt(SE1^2 + SE0^2))
    return(Z)
  }

  iterdata <- data[[i]]
  if(!is.na(restriction)){
    results <- rmstcov(time=iterdata[,"Time"],status=1-iterdata[,"Censored"],arm=iterdata[,"Trt"]-1,restriction=restriction,covariates=as.factor(iterdata[,stratum]))
    rmst <- c(restriction,results$RMST.arm1,results$RMST.arm0,results$RMST.diff,c(NA,NA))
  } else {rmst <- rep(NA,9)}
  if(!is.na(landmark)){
    analysis <- lm_cov(time=iterdata[,"Time"], status=1-iterdata[,"Censored"], arm=iterdata[,"Trt"]-1, stratum=iterdata[,stratum],landmark=landmark)
    lm <- c(landmark,analysis$LM.arm1,analysis$LM.arm0,analysis$LM.diff,c(NA,NA))
  } else {lm <- rep(NA,9)}
  return(c(rmst,lm))
}








#######################################################################################################
#'Summarise analyses of simulations of time-to-event data using arbitrary event, censoring and recruitment distributions.
#'
#' Function for summarising the analyses of simulated time-to-event trial data produced by analyse_sim().\cr
#' Automatically reads in format from analyse_sim(); no other input format is supported.\cr
#' It automatically detects types of analysis performed and provides relevant summaries (log-rank, Cox, RMST, landmark).\cr
#' @param analysed_results Output file from analyse_sim(). Only analyse_sim() output is supported.
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
