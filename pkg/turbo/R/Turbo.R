
### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### other functions
### ============================================================================


### ============================================================================
# construct initial profile from cakethickness, concentration 1
### ============================================================================

constructinitialprofile2 <- function (cakethickness, slicenumber, dx)  {
  initialprofile <- rep(0,slicenumber)
  cakeslices <- round(cakethickness / dx)
  initialprofile[1:cakeslices] <- 1
  return(initialprofile)
}

### ============================================================================
# construct initialprofile from data profile
### ============================================================================

constructinitialprofile <- function(concprofile, limits,
                                    cakethickness, slicenumber, dx) {
  initialprofile <- rep(0,slicenumber)
  cakeslices <- round(cakethickness / dx)
  thickness <- diff(limits)
  totalnumber <- sum(concprofile * thickness)
  initialprofile[1:cakeslices] <- totalnumber / cakeslices / dx
  return(initialprofile)
}

### ============================================================================
# calculates a concentration profile based on luminophore numbers
### ============================================================================

numtoconc <- function(values, limits)  {

  numtoconc <- NULL
  for(slice in 1:length(values)) {
    numtoconc <- c(numtoconc, values[slice] / (limits[slice+1] -
                                               limits [slice]))
  }
  return(numtoconc)
}


### ============================================================================
# calculates slice midpoints from slice limits
### ============================================================================

midpoints <- function(limits)  {

  midpoints <- NULL
  for(point in 1:(length(limits)-1)){
    midpoints <- c(midpoints,(limits[point]+limits[point+1])/2)
  }
  return(midpoints)
}

### ============================================================================
# calculates a (rough) profile from a model profile with regular intervals
### ============================================================================

roughprofile <- function(profile, limits, dx)  {
  roughprofile <- NULL
  mpoints <- seq(dx/2,length(profile)*dx,by=dx)
  for(slice in 1:(length(limits)-1)) {
    slicedata <- subset(profile,mpoints > limits[slice] &
                                mpoints < limits[slice+1])
                                
    roughprofile <- c(roughprofile, mean(slicedata))
  }
  return(roughprofile)
}

### ============================================================================
# estimate modelcost
### ============================================================================

modelcost <- function(profile, dataprofile, method="wssr")  {
  if (method=="ssr")
     cost <- sum((dataprofile-profile)^2)

  else if (method=="wssr")
    cost <- sum(((dataprofile-profile)^2)/(dataprofile))

  else if (method=="m1")
    cost <- sum((dataprofile-profile)^2)

  else if(method=="m2")
    cost <- sum((log(dataprofile)-log(profile))^2)

  return(cost)
}

### ============================================================================
# function modelcostprofile
### ============================================================================

modelcostprofile <- function(profile, dataprofile, method="wssr") {
  if (method=="ssr")
    cost <- (dataprofile-profile)^2

  else if (method=="wssr")
    cost <- ((dataprofile-profile)^2)/(dataprofile)

  else if (method=="m1")
    cost <- (dataprofile-profile)^2

  else if (method=="m2")
    cost <- (log(dataprofile)-log(profile))^2

  return(cost)
}






