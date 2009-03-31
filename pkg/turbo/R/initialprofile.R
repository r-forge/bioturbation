
### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### construct an initial profile based on data
### ============================================================================



initialprofile <- function(concprofile,limits,slicenumber,dx,cakethickness){
  initialprofile <- rep(0,slicenumber)
  cakeslices <- round(cakethickness / dx)
  thickness <- diff(limits)
  number <- sum(concprofile * thickness)
  initialprofile[1:cakeslices] <- number / cakeslices / dx
  return(initialprofile)
}