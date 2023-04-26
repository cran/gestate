# Testing function to check for integers
# Can do any, strictly positive, strictly negative, or positive/negative including zero
testNumber <- function(x,type=c("any","positive","positive0","negative","negative0"),integer=FALSE,is.single=FALSE){
  if(missing(x)){return(FALSE)}
  if(is.single){
    if(length(x)!= 1){return(FALSE)}
  }
  type <- match.arg(type)
  if(!all(is.numeric(x))){return(FALSE)}
  if(type=="positive"){
    if(!all(x > 0)){return(FALSE)}
  } else if(type=="positive0"){
    if(!all(x >= 0)){return(FALSE)}
  } else if(type=="negative"){
    if(!all(x < 0)){return(FALSE)}
  } else if(type=="negative0"){
    if(!all(x <= 0)){return(FALSE)}
  }
  if(integer){
    return(all(x == as.integer(x)))
  } else {return(TRUE)}
}

# Testing function to check for booleans
# Can also check for single boolean flag (is.flag=TRUE)
testBoolean <- function(x,is.flag=FALSE){
  if(missing(x)){return(FALSE)}
  if(is.flag){
    if(length(x)!= 1){return(FALSE)}
  }
  return(all(is.logical(x)))
}

# Testing function to check for probabilities
testProbability <- function(x,is.single=FALSE){
  if(missing(x)){return(FALSE)}
  if(is.single){
    if(length(x)!= 1){return(FALSE)}
  }
  if(!all(is.numeric(x))){return(FALSE)}
  return(all(x >= 0 & x <= 1))
}

# Testing function to check for strings
testString <- function(x,is.single=FALSE){
  if(missing(x)){return(FALSE)}
  if(is.single){
    if(length(x)!= 1){return(FALSE)}
  }
  if(!all(is.character(x))){return(FALSE)}
  return(TRUE)
}

# Testing function to check for Curves
# Can check for either a Curve (R=FALSE), or RCurve (R=TRUE)
testCurve <- function(x,R=FALSE,is.single=FALSE){
  if(missing(x)){return(FALSE)}
  if(is.single){
    if(length(x)!= 1){return(FALSE)}
  }
  if(R){
    return(all(class(x)[1]== "RCurve"))
  } else {
    return(all(class(x)[1]== "Curve"))    
  }
}

# Testing function to check for Covariance Matrices, making sure they are structurally valid
testCovarianceMatrix <- function(x){
   if(missing(x)){return(FALSE)}
   if(!is.matrix(x)){return(FALSE)}
   if(!all(is.numeric(x))){return(FALSE)}
   if(!all(isSymmetric(x))){return(FALSE)}
   if(!all(diag(x) >= 0)){return(FALSE)}
   #Check implied correlations are valid
   return(all(abs(t(x/diag(x))/diag(x)) <= 1))
}
