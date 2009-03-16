


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
binning <- 0.5

initialprofile <- rep(0,slicenumber)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

limits <- seq(0,slicenumber*dx,by=binning)
limitsmidpoints <- midpoints(limits)

#################################################
# parameters
#################################################

steplength <- 0.1
waitingtime <- 5
nldb <- (steplength^2)/(2*waitingtime)
db <- nldb
k <- 0.01
flux <- 1
fluxintroduction <- 0.5

mth <- "wssr"

#################################################
# generating profile with nonlocal model
#################################################

parameters <- c(steplength=steplength,waitingtime=waitingtime,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
nlprofile <- steady(runif(slicenumber),0,nonlocalmodel,parameters,jacfunc=nonlocaljac,jactype="fullusr")$y


#################################################
# generating profile with diffusive model
# db calculated from waiting time and step length
#################################################

parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- steady(y=initialprofile,time=0,func=diffusivemodel,parms=parameters)
diffprofile <- out$y

roughdiffprofile <- roughprofile(diffprofile,limits,dx)
roughnlprofile <- roughprofile(nlprofile,limits,dx)
nldbcost <- modelcost(roughdiffprofile,roughnlprofile,method=mth) 


#################################################
# fitting diffusive model
#################################################

costfunction <- function(db,dataprofile){
parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- steady(y=initialprofile,time=0,func=diffusivemodel,parms=parameters)
finalprofile <- out$y
roughdataprofile <- roughprofile(dataprofile,limits,dx)
roughfinalprofile <- roughprofile(finalprofile,limits,dx)
cost <- modelcost(roughfinalprofile,roughdataprofile,method=mth)
return(cost)
}

fit <- nlm(f=costfunction,1,nlprofile)
best <- fit$estimate
fitcost <- fit$minimum

fit2 <- nlm(f=costfunction,0,nlprofile)
best2 <- fit2$estimate
fitcost2 <- fit2$minimum


parameters <- c(db=best,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- steady(y=initialprofile,time=0,func=diffusivemodel,parms=parameters)
fittedprofile <- out$y

parameters <- c(db=best2,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out2 <- steady(y=initialprofile,time=0,func=diffusivemodel,parms=parameters)
fittedprofile2 <- out2$y


#################################################
# plot all profiles
#################################################

ymax <- 1000

roughnlprofile <- roughprofile(nlprofile,limits,dx)
roughfittedprofile <- roughprofile(fittedprofile,limits,dx)
roughfittedprofile2 <- roughprofile(fittedprofile2,limits,dx)


windows()
plot(initialmidpoints,initialprofile,type="l",col="gray",ylim=c(0,min(ymax,max(c(initialprofile,nlprofile)))))
lines(initialmidpoints,diffprofile,col="red")
lines(initialmidpoints,nlprofile,col="green")
lines(initialmidpoints,fittedprofile,col="blue")
points(limitsmidpoints,roughnlprofile,col="green")
points(limitsmidpoints,roughfittedprofile,col="blue")
lines(initialmidpoints,fittedprofile2,col="blue")
points(limitsmidpoints,roughfittedprofile2,col="blue")


#################################################
# numbers
#################################################

cat("\nnon local db: ",nldb,"\n",
"non local db cost: ", nldbcost,"\n",
"fitted db: ", best,"\n",
"fitted db cost: ", fitcost,"\n",
"fitted db 2: ", best2,"\n",
"fitted db 2 cost: ", fitcost2,"\n\n",sep="")


#################################################
# explore cost function
#################################################

dbs <- 10^seq(-5,3,length=100)
costs <- NULL
for(thisdb in dbs){
thiscost <- costfunction(thisdb,nlprofile)
costs <- c(costs,thiscost)
}
windows()
plot(dbs,costs,log="xy",main="cost function")


windows()
plot(limitsmidpoints,(modelcostprofile(roughfittedprofile,roughnlprofile,method=mth)),main="cost per depth")


