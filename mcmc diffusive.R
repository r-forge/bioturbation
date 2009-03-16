
#################################################
#
# diffusive and non local bioturbation models
# Pieter Provoost
#
# MCMC
#
#################################################


graphics.off()
source("luminophorefunctions.R")
library(deSolve)
library(Matrix)
library(rootSolve)
library(FME)







source("modcost2.R")









################################################
# settings
#################################################

slicenumber <- 400
dx <- 0.05
cakethickness <- 0.5
days <- 14
k <- 0
flux <- 0
fluxintroduction <- 0.5

datafile <- "data2.txt"
profilename <- "b"

#################################################
# data
#################################################

data <- read.table(datafile,header=TRUE)
rawprofile <- subset(data,profile==profilename)
datalimits <- unique(c(rawprofile$start,rawprofile$end))
dataprofile <- numtoconc(rawprofile$lum,datalimits)
datamidpoints <- midpoints(datalimits)
varnames <- rep("concentration",length(dataprofile))
dataprofileframe <- as.data.frame(cbind(name=varnames,depth=datamidpoints,concentration=dataprofile,err=dataprofile+0.1),stringsAsFactors=FALSE)
dataprofileframe[2:4] <- sapply(dataprofileframe[,2:4], as.numeric)
initialprofile <- constructinitialprofile(dataprofile,datalimits,cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

windows()
plot(datamidpoints,dataprofile,xlab="depth (cm)",ylab="luminophore concentration")

#################################################
# cost function
#################################################

objective <- function(db){
parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=2)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
modelprofile <- out[dim(out)[1],-1]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
modelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=roughmodelprofile))
cost <- modCost(mod=modelprofileframe,obs=dataprofileframe,x="depth",y="concentration")
return(cost)
}

wobjective <- function(db){
parameters <- c(db=db,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=2)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
modelprofile <- out[dim(out)[1],-1]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
modelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=roughmodelprofile))
cost <- modCost(mod=modelprofileframe,obs=dataprofileframe,x="depth",y="concentration",err="err")
return(cost)
}

#################################################
# fit
#################################################

fit <- modFit(p=0.01,f=objective,lower=0)
wfit <- modFit(p=0.01,f=wobjective,lower=0)

summ <- summary(fit)

#################################################
# run with fitted parameters
#################################################

times <- seq(0,days,length=2)
parameters <- c(db=fit$par[1],dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
modelprofile <- out[dim(out)[1],(2:dim(out)[2])]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)

plot(datamidpoints,dataprofile,xlab="depth (cm)",ylab="luminophore concentration",ylim=c(0,max(modelprofile)))
points(datamidpoints,roughmodelprofile,col="red")
lines(initialmidpoints,modelprofile,col="red")

times <- seq(0,days,length=2)
parameters <- c(db=wfit$par[1],dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
modelprofile <- out[dim(out)[1],(2:dim(out)[2])]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)

points(datamidpoints,roughmodelprofile,col="blue")
lines(initialmidpoints,modelprofile,col="blue")







extra <- 2.477076e-02
parameters <- c(db=extra,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
out2 <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
modelprofile2 <- out2[dim(out2)[1],(2:dim(out2)[2])]
roughmodelprofile2 <- roughprofile(modelprofile2,datalimits,dx)

#points(datamidpoints,roughmodelprofile2,col="blue")
#lines(initialmidpoints,modelprofile2,col="blue")











#################################################
# cost function
#################################################



dbs <- 10^seq(-5,1,length=100)
costs <- NULL
wcosts <- NULL
for(db in dbs){
#weights <- objective(db)$residuals$obs+0.1 
#wcosts <- c(wcosts,sum(objective(db)$residuals$res^2/weights))
wcosts <- c(wcosts,sum(wobjective(db)$residuals$res^2))
costs <- c(costs,sum(objective(db)$residuals$res^2))

}

windows()
plot(dbs,costs,log="xy",ylim=c(min(c(costs,wcosts)),max(c(costs,wcosts))),xlab="db",ylab="model cost")
points(dbs,wcosts,col="red")



#windows()
points(mcmc$pars[,1]^2/(2*mcmc$pars[,2]),mcmc$SS,col="green")


#################################################
# MCMC
#################################################

var0 <- summ$modVariance
covini <- summ$cov.scaled*2.4^2/2

mcmc <- modMCMC(p=coef(fit),f=objective,jump=covini,var0=var0,wvar0=1)
windows()
plot(mcmc)







