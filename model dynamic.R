

#################################################
#
# diffusive and non local bioturbation models
# Pieter Provoost
#
# PART I: LUMINOPHORES, NO FLUX, DYNAMIC
#
#################################################


graphics.off()
source("luminophorefunctions.R")
library(deSolve)
library(Matrix)
library(rootSolve)



################################################
# constants
#################################################

slicenumber <- 400
dx <- 0.05
cakethickness <- 0.5
binning <- 0.5

initialprofile <- constructinitialprofile2(cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

limits <- seq(0,slicenumber*dx,by=binning)
limitsmidpoints <- midpoints(limits)

#################################################
# parameters
#################################################

days <- 10
steplength <- 2
waitingtime <- 20
nldb <- (steplength^2)/(2*waitingtime)
db <- nldb
k <- 0.02
flux <- 0.1
fluxintroduction <- 0.5

mth <- "wssr"

#################################################
# generating profile with nonlocal model
#################################################

nlout <- lsoda(initialprofile,times,nonlocalmodel,parameters,jacfunc=nonlocaljac,jactype="fullusr")
nlprofile <- nlout[dim(nlout)[1],(2:dim(nlout)[2])]

#################################################
# generating profile with diffusive model
#################################################

parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=3)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
diffprofile <- out[dim(out)[1],-1]

roughdiffprofile <- roughprofile(diffprofile,limits,dx)
roughnlprofile <- roughprofile(nlprofile,limits,dx)
nldbcost <- modelcost(roughdiffprofile,roughnlprofile,method=mth)


#################################################
# fitting diffusive model
#################################################

costfunction <- function(db,dataprofile){
times <- seq(0,days,length=3)
parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
finalprofile <- out[dim(out)[1],-1]
dataprofile <- roughprofile(dataprofile,limits,dx)
finalprofile <- roughprofile(finalprofile,limits,dx)
cost <- modelcost(finalprofile,dataprofile,method=mth)
return(cost)
}

fit <- nlm(f=costfunction,1,nlprofile)
best <- fit$estimate
fitcost <- fit$minimum

fit2 <- nlm(f=costfunction,0,nlprofile)
best2 <- fit2$estimate
fitcost2 <- fit2$minimum

parameters <- c(db=best,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=3)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
fittedprofile <- out[dim(out)[1],-1]

parameters <- c(db=best2,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=3)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
fittedprofile2 <- out[dim(out)[1],-1]


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

dbs <- 10^seq(-5,0.5,length=100)
costs <- NULL
for(thisdb in dbs){
thiscost <- costfunction(thisdb,nlprofile)
costs <- c(costs,thiscost)
}
windows()
plot(dbs,costs,log="xy",main="cost function")


windows()
plot(limitsmidpoints,(modelcostprofile(roughfittedprofile,roughnlprofile,method=mth )),main="cost per depth")





