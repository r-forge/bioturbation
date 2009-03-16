#################################################
#
# diffusive and non local bioturbation models
# Pieter Provoost
#
# PART II: TRANSIENT TRACERS, FLUX, STEADY STATE
#
#################################################


graphics.off()
source("luminophorefunctions.R")
library(deSolve)
library(Matrix)
library(rootSolve)



#################################################
# constants
#################################################

slicenumber <- 400
dx <- 0.05
cakethickness <- 0.5

initialprofile <- rep(0,slicenumber)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

limits <- seq(0,slicenumber*dx,by=0.5)

#################################################
# parameters
#################################################

steplength <- 0.5
waitingtime <- 10
nldb <- (steplength^2)/(2*waitingtime)
db <- nldb
k <- 0.02
flux <- 0.01
fluxintroduction <- 0.1




#################################################
#################################################
#################################################



parameters <- c(steplength=steplength,waitingtime=waitingtime,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)


nonlocalJac <- function(t,profile,parameters){
with(as.list(parameters),{
slicenumber <- length(profile)
Mat <- constructtransitionmatrix(steplength,waitingtime,slicenumber,dx)
diag(Mat) <- diag(Mat) - k
return(Mat)
})
}


print(system.time(ST <- steady(runif(200),0,nonlocalmodel,parameters,
      jacfunc=nonlocalJac, jactype="fullusr")))



