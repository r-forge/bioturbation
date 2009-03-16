


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


#################################################
# data and parameters
#################################################



steplength <- 0.5
waitingtime <- 0.625
days <- 30

cakethickness <- 0.5
slicenumber <- 200
dx <- 0.05




#################################################
# comparison local and nonlocal models
#################################################

initialprofile <- constructinitialprofile2(cakethickness,slicenumber,dx)


db <- steplength^2/(2*waitingtime)


times <- seq(0,days,0.1)

model <- function(t,profile,parameters){
with(as.list(parameters),{
flux <- db * c(0,diff(c(profile,0))) / dx
dprofile <- diff(flux) / dx
list(c(dprofile))
})
}

parameters <- c(db=db,dx=dx)
out <- ode.band(times=times,y=initialprofile,func=model,parms=parameters,nspec=1)
finalprofile <- out[dim(out)[1]-1,-1]
modelrough <- roughprofile(finalprofile,limits,dx)

depth <- slicenumber * dx
psil <- function(l){
return((1/(steplength*sqrt(2*pi)))*exp(-(l^2)/(2*steplength^2)))
}

startint <- -depth+dx/2
endint <- depth-dx/2
intlimits <- seq(startint,endint,by=dx)
intvals <- NULL
for(i in 1:(length(intlimits)-1)){
intval <- integrate(psil,intlimits[i],intlimits[i+1])$value
intvals <- c(intvals,intval)
}
intvals <- intvals / waitingtime
intvals[slicenumber] <- intvals[slicenumber] + (1-sum(intvals))
transitionmatrix <- matrix(data=NA,nrow=slicenumber,ncol=slicenumber)
for(column in 1:slicenumber){
transitionmatrix[,column] <- intvals[(slicenumber-column+1):(2*slicenumber-column)]
rest <- rev(intvals[1:(slicenumber-column)])
transitionmatrix[(1:length(rest)),column] <- transitionmatrix[(1:length(rest)),column] + rest
}

parameters <- NULL

nlmodel <- function(t,profile,parameters){
dprofile <- transitionmatrix %*% profile - profile 
list(c(dprofile))
}

times <- seq(0,days,length=10)
library("deSolve")
nlout <- ode.band(times=times,y=initialprofile,func=nlmodel,parms=parameters,nspec=1)
nlfinalprofile <- nlout[dim(nlout)[1]-1,-1]
nlmodelrough <- roughprofile(nlfinalprofile,limits,dx)



windows()
plot(initialprofile,type="l",xlab="depth",ylab="luminophore concentration",col="lightgray")
lines(finalprofile,col="red",lwd=2)
lines(nlfinalprofile,col="blue",lwd=2)

