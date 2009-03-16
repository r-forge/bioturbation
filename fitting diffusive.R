

#################################################
#
# luminophore profile analysis
# Pieter Provoost
#
# data considered to be luminophore numbers per
# slice, not concentrations
#
# units: days, cm
#
#################################################


graphics.off()
source("luminophorefunctions.R")
library(deSolve)
library(Matrix)

#################################################
# data and parameters
#################################################

data <- read.table("data.txt",header=TRUE)
slicenumber <- 200
db <- 0.07
dx <- 0.1
days <- 14
cakethickness <- 0.5


#################################################
# some calculations
#################################################

profiles <- unique(data$profile)
thisprofile <- subset(data,profile=="a")
limits <- unique(c(thisprofile$start,thisprofile$end))
concprofile <- numtoconc(thisprofile$lum, limits)
concmidpoints <- midpoints(limits)
initialprofile <- constructinitialprofile(concprofile,limits,cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

#################################################
# fitting diffusive model
#################################################

costfunction <- function(db,dataprofile){
parameters <- c(db=db,dx=dx,k=0)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
finalprofile <- out[dim(out)[1]-1,-1]
modelrough <- roughprofile(finalprofile,limits,dx)
cost <- sum((modelrough-dataprofile)^2)
return(cost)
}

best <- optimize(f=costfunction,interval=c(0,1),concprofile)$minimum
best


#################################################
# run model with fitted db
#################################################


parameters <- c(db=best,dx=dx,k=0)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
finalprofile <- out[dim(out)[1]-1,-1]
modelrough <- roughprofile(finalprofile,limits,dx)

windows()
maxval <- max(c(concprofile,initialprofile))
plot(initialmidpoints,initialprofile,ylim=c(0,maxval),type="l",xlab="depth",main="fitted diffusive model",ylab="luminophore concentration")
lines(concmidpoints,concprofile,col="red")
points(concmidpoints,concprofile,col="red")
lines(initialmidpoints,finalprofile,col="green")
points(concmidpoints,modelrough,col="green")

