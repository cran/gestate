
##################################################################################################################################
# S4 architecture
#
##################################################################################################################################

######################################################################################################
#'Curve Class for defining distributions
#'
#' This class allows distributions to be defined. It contains all information needed to reproduce a distribution.\cr
#' References to functions that store the PDF, CDF and random number generator.\cr
#' Parameters are also stored.
#' @slot type Type of Curve (character). Typically the distribution name.
#' @slot PDF Name of the PDF function describing the Curve.
#' @slot CDF Name of the CDF function describing the Curve.
#' @slot RF Name of the random generator function describing the Curve.
#' @slot inverse Name of the inverse CDF function describing the Curve. Optional; set to NULL if unavailable.
#' @slot paramno Number of parameters required to define the distribution.
#' @slot pnames Names of parameters defining the distribution. Should be a vector of length paramno.
#' @slot pvalue Values of parameters defining the distribution. Should be a list of length paramno.
#' @author James Bell
#' @examples new("Curve", type="ExampleCurve",PDF="pdf_fn_name",CDF="CDF_fn_name",
#'   RF="random_draw_fn_name", inverse="inv_fn_name",paramno=2,pnames=c('param1','param2'),
#'   pvalue=list(1,2))
#' @export
setClass("Curve",
         slots = list(type="character", PDF = "character", CDF = "character", RF = "character", inverse="character",paramno = "numeric",
                      pnames="character",pvalue="list"))

#'RCurve Class for defining recruitment distributions
#'
#' This class extends the Curve class, adding recruitment-related quantities such as patient numbers.
#' @slot N Total number of patients recruited.
#' @slot Nactive Number of patients recruited in active arm. Nactive+Ncontrol=N.
#' @slot Ncontrol Number of patients recruited in control arm. Nactive+Ncontrol=N.
#' @slot Ratio Randomisation ratio. Nactive divided by Ncontrol = Ratio.
#' @slot Length Total length of the recruitment period.
#' @slot RF Name of the random generator function describing the Curve.
#' @slot inverse Name of the inverse CDF function describing the Curve. Optional; set to NULL if unavailable.
#' @slot paramno Number of parameters required to define the distribution.
#' @slot pnames Names of parameters defining the distribution. Should be a vector of length paramno.
#' @slot pnames Values of parameters defining the distribution. Should be a list of length paramno.
#' @slot maxF Maximum length of patient follow-up. Typically should be Inf.
#' @author James Bell
#' @examples new("RCurve", type="ExampleCurve",PDF="pdf_fn_name", CDF="CDF_fn_name",
#'   RF="random_draw_fn_name", inverse="inv_fn_name", paramno=2, pnames=c('param1','param2'),
#'   pvalue=list(1,2), 
#' N=100,Nactive=50,Ncontrol=40, Ratio=50/40, Length = 5, maxF = Inf)
#' @export
setClass("RCurve",
         slots = list(N = "numeric", Nactive = "numeric", Ncontrol = "numeric", Ratio = "numeric",Length = "numeric",maxF="numeric"),contains="Curve")

######################################################################################################
# S4 methods for the 'Curve' and 'RCurve' objects. These are mostly simple get/set methods for properties
# Individual comments are included for more interesting ones
#
# Author: James Bell, james.bell.ext@boehringer-ingelheim.com
#
######################################################################################################

######################################################################################################
#'Method for returning the PDF function for a Curve object
#'
#' This retrieves the PDF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples getPDFfunction(Weibull(100,1))
#' @export
setGeneric(name="getPDFfunction", def=function(theObject,...) {
  standardGeneric("getPDFfunction")
})

#'Method for returning the PDF function for a Curve object
#'
#' This retrieves the PDF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param x The parameter name to use in the PDF function. Default=x
#' @examples getPDFfunction(Weibull(100,1))
#' @export
setMethod(f="getPDFfunction", signature="Curve",definition=function(theObject,x="x"){
  PDF <- paste(theObject@PDF,"(x=",x)
  for(i in 1:theObject@paramno){
    PDF <- paste(PDF,",",theObject@pnames[i],"=",theObject@pvalue[i])
  }
  PDF <- paste(PDF,")")
  return(PDF)
})

######################################################################################################
#'Method for returning the CDF function for a Curve object
#'
#' This retrieves the full CDF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples getCDFfunction(Weibull(100,1))
#' @export
setGeneric(name="getCDFfunction", def=function(theObject,...) {
  standardGeneric("getCDFfunction")
})
#'Method for returning the CDF function for a Curve object
#'
#' This retrieves the CDF function of the specified Curve object as a string
#' @param theObject The name of the Curve Object
#' @param q The parameter name to use in the CDF function. Default=q
#' @examples getCDFfunction(Weibull(100,1))
#' @export
setMethod(f="getCDFfunction", signature="Curve",definition=function(theObject,q="q"){
  CDF <- paste(theObject@CDF,"(q=",q)
  for(i in 1:theObject@paramno){
    CDF <- paste(CDF,",",theObject@pnames[i],"=",theObject@pvalue[i])
  }
  CDF <- paste(CDF,")")
  return(CDF)
})

######################################################################################################
#'Method for evaluating the CDF function for a Curve object at q
#'
#' This numerically evaluates the CDF function of the Curve object at the specified q
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples evaluateCDFfunction(Weibull(100,1),10)
#' @export
setGeneric(name="evaluateCDFfunction", def=function(theObject,...) {
  standardGeneric("evaluateCDFfunction")
})
#'Method for evaluating the CDF function for a Curve object at q
#'
#' This numerically evaluates the CDF function of the Curve object at the specified q
#' @param theObject The name of the Curve Object
#' @param q The time to evaluate at
#' @examples evaluateCDFfunction(Weibull(100,1),10)
#' @export
setMethod(f="evaluateCDFfunction", signature="Curve",definition=function(theObject,q){
  return( eval(parse(text=getCDFfunction(theObject)))  )
})

######################################################################################################

#'Method for returning the RF function for a Curve object
#'
#' This retrieves the full RF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples getRFfunction(Weibull(100,1))
#' @export
setGeneric(name="getRFfunction", def=function(theObject,...) {
  standardGeneric("getRFfunction")
})

#'Method for returning the RF function for a Curve object
#'
#' This retrieves the full RF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param n The parameter name for the number of random draws to make. Default=n
#' @examples getRFfunction(Weibull(100,1))
#' @export
setMethod(f="getRFfunction", signature="Curve",definition=function(theObject,n="n"){
  RF <- paste(theObject@RF,"(n=",n)
  for(i in 1:theObject@paramno){
    RF <- paste(RF,",",theObject@pnames[i],"=",theObject@pvalue[i])
  }
  RF <- paste(RF,")")
  return(RF)
})

######################################################################################################
#'Method for creating a random draw function from a Curve object
#'
#' This creates a random draw function from the Curve object
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples createRFfunction(Weibull(100,1))
#' @export
setGeneric(name="createRFfunction", def=function(theObject,...) {
  standardGeneric("createRFfunction")
})

#'Method for creating a random draw function from a Curve object
#'
#' This creates a random draw function from the Curve object
#' @param theObject The name of the Curve Object
#' @param n The parameter name for the number of random draws to make. Default=n
#' @examples createRFfunction(Weibull(100,1))
#' @export
setMethod(f="createRFfunction", signature="Curve",definition=function(theObject,n="n"){
  return(  eval(parse(text=paste("function(",n,"){",getRFfunction(theObject,n="n"),"}")))  )
})

######################################################################################################
#'Method for taking random draws from a Curve object distribution
#'
#' This takes random draws from the Curve object distribution
#' Note that the seed should always be appropriately set before invoking this
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples random_draw(Weibull(100,1),1000)
#' @export
setGeneric(name="random_draw", def=function(theObject,...) {
  standardGeneric("random_draw")
})

#'Method for taking random draws from a Curve object distribution
#'
#' This takes random draws from the Curve object distribution
#' Note that the seed should always be appropriately set before invoking this
#' @param theObject The name of the Curve Object
#' @param n Number of random draws required
#' @examples random_draw(Weibull(100,1),1000)
#' @export
setMethod(f="random_draw", signature="Curve",definition=function(theObject,n){
  x <- createRFfunction(theObject)
  return(x(n))
})

######################################################################################################
#'Method for returning the inverse-CDF function for a Curve object
#'
#' This retrieves the full inverse-CDF function of the Curve object as a string
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples getInvfunction(Weibull(100,1))
#' @export
setGeneric(name="getInvfunction", def=function(theObject,...) {
  standardGeneric("getInvfunction")
})
#'Method for returning the inverse-CDF function for a Curve object
#'
#' This retrieves the inverse CDF function of the specified Curve object as a string
#' @param theObject The name of the Curve Object
#' @param p The probability parameter name to use in the CDF function. Default=p
#' @examples getCDFfunction(Weibull(100,1))
#' @export
setMethod(f="getInvfunction", signature="Curve",definition=function(theObject,p="p"){
  INV <- paste(theObject@inverse,"(p=",p)
  for(i in 1:theObject@paramno){
    INV <- paste(INV,",",theObject@pnames[i],"=",theObject@pvalue[i])
  }
  INV <- paste(INV,")")
  return(INV)
})

######################################################################################################
#'Method for evaluating the inverse-CDF function for a Curve object at p
#'
#' This numerically evaluates the inverse-CDF function of the Curve object at the specified p
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples evaluateInvfunction(Weibull(100,1), 0.5)
#' @export
setGeneric(name="evaluateInvfunction", def=function(theObject,...) {
  standardGeneric("evaluateInvfunction")
})
#'Method for evaluating the inverse-CDF function for a Curve object at p
#'
#' This numerically evaluates the inverse-CDF function of the Curve object at the specified p
#' @param theObject The name of the Curve Object
#' @param p The probability to evaluate the time at
#' @examples evaluateInvfunction(Weibull(100,1), 0.5)
#' @export
setMethod(f="evaluateInvfunction", signature="Curve",definition=function(theObject,p){
  return( eval(parse(text=getInvfunction(theObject)))  )
})



#######################################################################################################
#'Method for returning the CDF function for a RCurve object
#'
#' This retrieves the CDF function of the specified RCurve object as a string, writing in the assessment time parameter
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples getAssessCDFfunction(LinearR(12,100,100))
#' @export
setGeneric(name="getAssessCDFfunction", def=function(theObject,...) {
  standardGeneric("getAssessCDFfunction")
})

#'Method for returning the CDF function for a RCurve object
#'
#' This retrieves the CDF function of the specified RCurve object as a string, writing in the assessment time parameter
#' @param theObject The name of the RCurve Object
#' @param q The parameter name to use in the CDF function. Default=q
#' @examples getAssessCDFfunction(LinearR(12,100,100))
#' @export
setMethod(f="getAssessCDFfunction", signature="RCurve",definition=function(theObject,q="q"){
  CDF <- paste(theObject@CDF,"(q=",q,",assess=assess")
  for(i in 1:theObject@paramno){
    CDF <- paste(CDF,",",theObject@pnames[i],"=",theObject@pvalue[i])
  }
  CDF <- paste(CDF,")")
  return(CDF)
})

#'Method for returning a single parameter from a Curve object
#'
#' This retrieves a single parameter from a Curve object
#' @param theObject The name of the Curve Object
#' @param ... Pass-through arguments
#' @examples getParam(Weibull(100,1),1)
#' @export
setGeneric(name="getParam", def=function(theObject,...) {
  standardGeneric("getParam")
})
#'Method for returning a single parameter from a Curve object
#'
#' This retrieves a single parameter from a Curve object
#' @param theObject The name of the Curve Object
#' @param param The number of the parameter that is required
#' @examples getParam(Weibull(100,1),1)
#' @export
setMethod(f="getParam", signature="Curve",definition=function(theObject,param){
  value <- theObject@pvalue[[param]]
  names(value) <- theObject@pnames[[param]]
  return(value)
})

#'Method for returning all parameters from a Curve object as a list
#'
#' This retrieves all parameters from a Curve object as a list
#' @param theObject The name of the Curve Object
#' @examples getParams(Weibull(100,1))
#' @export
setGeneric(name="getParams", def=function(theObject) {
  standardGeneric("getParams")
})

#'Method for returning all parameters from a Curve object as a list
#'
#' This retrieves all parameters from a Curve object as a list
#' @param theObject The name of the Curve Object
#' @examples getParams(Weibull(100,1))
#' @export
setMethod(f="getParams", signature="Curve",definition=function(theObject){
  values <- theObject@pvalue
  names(values) <- theObject@pnames
  return(values)
})

# Create a method to get all parameter values as a vector
# Useful for printing them
#'Method for returning all parameters from a Curve object as a vector
#'
#' This retrieves all parameters from a Curve object as a vector
#' @param theObject The name of the Curve Object
#' @examples getParamsV(Weibull(100,1))
#' @export
setGeneric(name="getParamsV", def=function(theObject) {
  standardGeneric("getParamsV")
})

#'Method for returning all parameters from a Curve object as a vector
#'
#' This retrieves all parameters from a Curve object as a vector
#' @param theObject The name of the Curve Object
#' @examples getParamsV(Weibull(100,1))
#' @export
setMethod(f="getParamsV", signature="Curve",definition=function(theObject){
  values <- unlist(theObject@pvalue)
  size <- length(values)/length(theObject@pvalue)
  names(values) <- rep(theObject@pnames,each=size)
  return(values)
})

#'Method for returning all parameter names from a Curve object
#'
#' This retrieves all parameter names from a Curve object
#' @param theObject The name of the Curve Object
#' @examples getNames(Weibull(100,1))
#' @export
# Create a method to get all parameter names
setGeneric(name="getNames", def=function(theObject) {
  standardGeneric("getNames")
})

#'Method for returning all parameter names from a Curve object
#'
#' This retrieves all parameter names from a Curve object
#' @param theObject The name of the Curve Object
#' @examples getNames(Weibull(100,1))
#' @export
setMethod(f="getNames", signature="Curve",definition=function(theObject){
  return(theObject@pnames)
})

#'Method for returning the Curve type
#'
#' This returns the Curve object type
#' @param theObject The name of the Curve Object
#' @examples getType(Weibull(100,1))
#' @export
# Create a method to get curve type
setGeneric(name="getType", def=function(theObject) {
  standardGeneric("getType")
})

#'Method for returning the Curve type
#'
#' This returns the Curve object type
#' @param theObject The name of the Curve Object
#' @examples getType(Weibull(100,1))
#' @export
setMethod(f="getType", signature="Curve",definition=function(theObject){
  return(theObject@type)
})

#'Method for returning the total patient number from a RCurve
#'
#' This returns the RCurve total patient number
#' @param theObject The name of the RCurve Object
#' @examples getN(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get total patient number
setGeneric(name="getN", def=function(theObject) {
  standardGeneric("getN")
})

#'Method for returning the total patient number from a RCurve
#'
#' This returns the RCurve total patient number
#' @param theObject The name of the RCurve Object
#' @examples getN(LinearR(12,100,100))
#' @export
setMethod(f="getN", signature="RCurve",definition=function(theObject){
  return(theObject@N)
})

#'Method for returning the active arm patient number from a RCurve
#'
#' This returns the RCurve active arm patient number
#' @param theObject The name of the RCurve Object
#' @examples getNactive(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get active arm patient number
setGeneric(name="getNactive", def=function(theObject) {
  standardGeneric("getNactive")
})

#'Method for returning the active arm patient number from a RCurve
#'
#' This returns the RCurve active arm patient number
#' @param theObject The name of the RCurve Object
#' @examples getNactive(LinearR(12,100,100))
#' @export
setMethod(f="getNactive", signature="RCurve",definition=function(theObject){
  return(theObject@Nactive)
})

#'Method for returning the control arm patient number from a RCurve
#'
#' This returns the RCurve control arm patient number
#' @param theObject The name of the RCurve Object
#' @examples getNcontrol(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get control arm patient number
setGeneric(name="getNcontrol", def=function(theObject) {
  standardGeneric("getNcontrol")
})

#'Method for returning the control arm patient number from a RCurve
#'
#' This returns the RCurve control arm patient number
#' @param theObject The name of the RCurve Object
#' @examples getNcontrol(LinearR(12,100,100))
#' @export
setMethod(f="getNcontrol", signature="RCurve",definition=function(theObject){
  return(theObject@Ncontrol)
})

#'Method for setting N's in an RCurve
#'
#' This sets the RCurve patient numbers per arm directly and updates N and Ratio accordingly
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples setPatients(LinearR(12,100,100),200,150)
#' @export
# Create a method for RCurve to get control arm patient number
setGeneric(name="setPatients", def=function(theObject,...) {
  standardGeneric("setPatients")
})

#'Method for setting N's in an RCurve
#'
#' This sets the RCurve patient numbers per arm directly and updates N and Ratio accordingly
#' @param theObject The name of the RCurve Object
#' @param Nactive Number of patients to set in active arm
#' @param Ncontrol Number of patients to set in control arm
#' @examples setPatients(LinearR(12,100,100),200,150)
#' @export
setMethod(f="setPatients", signature="RCurve",definition=function(theObject,Nactive,Ncontrol){
  theObject@Nactive <- Nactive
  theObject@Ncontrol <- Ncontrol
  theObject@N <- Nactive+Ncontrol
  theObject@Ratio <- Nactive/Ncontrol
  return(theObject)
})

#'Method for returning the recruitment ratio from a RCurve
#'
#' This returns the RCurve recruitment ratio
#' @param theObject The name of the RCurve Object
#' @examples getRatio(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get the recruitment ratio
setGeneric(name="getRatio", def=function(theObject) {
  standardGeneric("getRatio")
})

#'Method for returning the recruitment ratio from a RCurve
#'
#' This returns the RCurve recruitment ratio
#' @param theObject The name of the RCurve Object
#' @examples getRatio(LinearR(12,100,100))
#' @export
setMethod(f="getRatio", signature="RCurve",definition=function(theObject){
  return(theObject@Ratio)
})

#'Method for returning the recruitment length from a RCurve
#'
#' This returns the RCurve recruitment length
#' @param theObject The name of the RCurve Object
#' @examples getLength(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get the recruitment ratio
setGeneric(name="getLength", def=function(theObject) {
  standardGeneric("getLength")
})

#'Method for returning the recruitment length from a RCurve
#'
#' This returns the RCurve recruitment length
#' @param theObject The name of the RCurve Object
#' @examples getLength(LinearR(12,100,100))
#' @export
setMethod(f="getLength", signature="RCurve",definition=function(theObject){
  return(theObject@Length)
})

#'Method for returning maximum duration of patient follow-up from a RCurve
#'
#' This returns the RCurve maximum patient follow-up
#' @param theObject The name of the RCurve Object
#' @examples getMaxF(LinearR(12,100,100))
#' @export
# Create a method for RCurve to get the recruitment ratio
setGeneric(name="getMaxF", def=function(theObject) {
  standardGeneric("getMaxF")
})

#'Method for returning maximum duration of patient follow-up from a RCurve
#'
#' This returns the RCurve maximum patient follow-up
#' @param theObject The name of the RCurve Object
#' @examples getMaxF(LinearR(12,100,100))
#' @export
setMethod(f="getMaxF", signature="RCurve",definition=function(theObject){
  return(theObject@maxF)
})

#'Method for calculating expected number of recruited patients at a given time from an RCurve
#'
#' This calculates the expected number of recruited patients at a given time based upon an RCurve
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples getPatients(LinearR(12,100,100),7)
#' @export
setGeneric(name="getPatients", def=function(theObject,...) {
  standardGeneric("getPatients")
})

#'Method for calculating expected number of recruited patients at a given time from an RCurve
#'
#' This calculates the expected number of recruited patients at a given time based upon an RCurve
#' @param theObject The name of the RCurve Object
#' @param x time at which to calculate expected patients recruited"
#' @examples getPatients(LinearR(12,100,100),7)
#' @export
setMethod(f="getPatients", signature="RCurve",definition=function(theObject,x){
  #Note, assess is a dummy variable used for compatibility with the function, which needs q and assess
  assess <- x+0.00001 
  q <- assess-x 
  output <- theObject@N*(1-eval(parse(text=getAssessCDFfunction(theObject))))
  return(output)
})

##################################################################################################
#'Method for plotting the CDF of a Curve object
#'
#' This plots a Curve CDF
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples plotCDF(Weibull(100,1))
#' plotCDF(Weibull(100,1),xlab="Test x label",maxT=60)
#' plotCDF(Weibull(80,0.8),overlay=TRUE,col=2,lty=2)
#' @export
# Create a method to plot CDF
setGeneric(name="plotCDF", def=function(theObject,...) {
  standardGeneric("plotCDF")
})

#'Method for plotting the CDF of a Curve object
#'
#' This plots a Curve CDF
#' @param theObject The name of the RCurve Object
#' @param overlay Boolean whether to overlay on existing one (vs start a new one). Default=FALSE
#' @param maxT Maximum time to plot up to. Default=100
#' @param increment Plotting time increment. Default=0.1
#' @param xlab X-axis label. Default="Time"
#' @param ylab Y-axis label. Default="CDF"
#' @param main title of plot. Default="CDF plot"
#' @param type type of plot (see standard graphical parameters). Default="l" (lines).
#' @param ... Standard graphical parameter arguments to be passed on to 'plot'/'lines', e.g. to change appearance of plot.
#' @examples plotCDF(Weibull(100,1))
#' plotCDF(Weibull(100,1),xlab="Test x label",maxT=60)
#' plotCDF(Weibull(80,0.8),overlay=TRUE,col=2,lty=2)
#' @export
setMethod(f="plotCDF", signature="Curve",definition=function(theObject,overlay=FALSE,maxT=100,increment=0.1,xlab="Time",ylab="CDF",main="CDF plot",type="l",...){
  q <- seq(from=0,to=maxT,by=increment)
  y <- evaluateCDFfunction(theObject,q)
  if(overlay){
    lines(q,y,type=type,...)
  } else{
    plot(q,y,xlim=c(0,maxT),ylim=c(0,1),xlab=xlab,ylab=ylab,main=main,type=type,...)
  }
})

#'Method for plotting the Survival Function of a Curve object
#'
#' This plots a Curve Survival Function
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples plotSF(Weibull(100,1))
#' plotSF(Weibull(100,1),xlab="Test x label",maxT=60)
#' plotSF(Weibull(80,0.8),overlay=TRUE,col=2,lty=2)
#' @export
# Create a method to plot SF
setGeneric(name="plotSF", def=function(theObject,...) {
  standardGeneric("plotSF")
})

#'Method for plotting the Survival Function of a Curve object
#'
#' This plots a Curve Survival Function
#' @param theObject The name of the RCurve Object
#' @param overlay Boolean whether to overlay on existing one (vs start a new one). Default=FALSE
#' @param maxT Maximum time to plot up to. Default=100
#' @param increment Plotting time increment. Default=0.1
#' @param xlab X-axis label. Default="Time"
#' @param ylab Y-axis label. Default="S(t)"
#' @param main title of plot. Default="Kaplan Meier plot"
#' @param type type of plot (see standard graphical parameters). Default="l" (lines).
#' @param ... Standard graphical parameter arguments to be passed on to 'plot'/'lines', e.g. to change appearance of plot.
#' @examples plotSF(Weibull(100,1))
#' plotSF(Weibull(100,1),xlab="Test x label",maxT=60)
#' plotSF(Weibull(80,0.8),overlay=TRUE,col=2,lty=2)
#' @export
setMethod(f="plotSF", signature="Curve",definition=function(theObject,overlay=FALSE,maxT=100,increment=0.1,xlab="Time",ylab="S(t)",main="Kaplan Meier plot",type="l",...){
  q <- seq(from=0,to=maxT,by=increment)
  y <- evaluateCDFfunction(theObject,q)
  if(overlay){
    lines(q,1-y,type=type,...)
  } else{
    plot(q,1-y,xlim=c(0,maxT),ylim=c(0,1),xlab=xlab,ylab=ylab,main=main,type=type,...)
  }
})

#'Method for plotting the Recruitment Function of a RCurve object
#'
#' This plots an RCurve Recruitment Function
#' @param theObject The name of the RCurve Object
#' @param ... Pass-through arguments
#' @examples plotRecruitment(LinearR(12,100,100))
#' plotRecruitment(LinearR(12,100,100),xlab="Test x label",maxT=60)
#' plotRecruitment(LinearR(20,90,90),overlay=TRUE,col=2,lty=2)
#' @export
# Create a method to plot Recruitment Function
setGeneric(name="plotRecruitment", def=function(theObject,...) {
  standardGeneric("plotRecruitment")
})

#'Method for plotting the Recruitment Function of a RCurve object
#'
#' This plots an RCurve Recruitment Function
#' @param theObject The name of the RCurve Object
#' @param overlay Boolean whether to overlay on existing one (vs start a new one). Default=FALSE
#' @param maxT Maximum time to plot up to. Default=100
#' @param increment Plotting time increment. Default=0.1
#' @param xlab X-axis label. Default="Time"
#' @param ylab Y-axis label. Default="Patients"
#' @param main title of plot. Default="Recruitment plot"
#' @param type type of plot (see standard graphical parameters). Default="l" (lines).
#' @param ... Standard graphical parameter arguments to be passed on to 'plot'/'lines', e.g. to change appearance of plot.
#' @examples plotRecruitment(LinearR(12,100,100))
#' plotRecruitment(LinearR(12,100,100),xlab="Test x label",maxT=60)
#' plotRecruitment(LinearR(20,90,90),overlay=TRUE,col=2,lty=2)
#' @export
setMethod(f="plotRecruitment", signature="RCurve",definition=function(theObject,overlay=FALSE,maxT=100,increment=0.1,xlab="Time",ylab="Patients",main="Recruitment plot",type="l",...){
  x <- seq(from=0,to=maxT,by=increment)
  assess <- maxT 
  q <- assess-x 
  output <- theObject@N*(1-eval(parse(text=getAssessCDFfunction(theObject))))
  
  if(overlay){
    lines(x,output,type=type,...)
  } else{
    plot(x,output,xlim=c(0,maxT),ylim=c(0,theObject@N),xlab=xlab,ylab=ylab,main=main,type=type,...)
  }
})


#################################################################################################
# Show methods to hide unnecessary content in Curve and RCurve
##################################################################################################
#'Method for displaying Curve objects neatly - replaces standard show method
#'
#' This is the standard print method for a Curve object
#' @param object The name of the Curve Object
#' @import methods
#' @examples Weibull(100,1)
#' @export
setMethod(f="show", signature="Curve",definition=function(object){
  cat("Curve Object\n")
  cat("Type:",object@type,"\n")
  cat("Distribution Parameters:\n")
  for(i in 1: object@paramno){
    cat("  ",object@pnames[i],": ",object@pvalue[[i]],"\n",sep=" ")
  }
})

#'Method for displaying RCurve objects neatly - replaces standard show method
#'
#' This is the standard print method for an RCurve Object
#' @param object The name of the RCurve Object
#' @import methods
#' @examples LinearR(12,100,100)
#' @export
setMethod(f="show", signature="RCurve",definition=function(object){
  cat("RCurve Object\n")
  cat("Type:",object@type,"\n")
  cat("N Total:",object@N,"\n")
  cat("N Active:",object@Nactive,"\n")
  cat("N Control:",object@Ncontrol,"\n")
  cat("Recruitment Ratio:",object@Ratio,"\n")
  cat("Distribution Parameters:\n")
  for(i in 1: object@paramno){
    cat("  ",object@pnames[i],": ",object@pvalue[[i]],"\n",sep=" ")
  }
})

