### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### conversion of model profiles for comparison with data  
### ============================================================================


roughprofile <- function(profile, limits, dx){
  roughprofile <- NULL
  mpoints <- seq(dx/2,length(profile)*dx,by=dx)
  for(slice in 1:(length(limits)-1)){
    slicedata <- subset(profile,mpoints > limits[slice] & mpoints < limits[slice+1])
    roughprofile <- c(roughprofile, mean(slicedata))
  }
  return(roughprofile)
}
