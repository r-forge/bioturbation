
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

datafile <- "data.txt"
profilename <- "a"

#################################################
# data
#################################################

data <- read.table(datafile,header=TRUE)
rawprofile <- subset(data,profile==profilename)
datalimits <- unique(c(rawprofile$start,rawprofile$end))
dataprofile <- numtoconc(rawprofile$lum,datalimits)
datamidpoints <- midpoints(datalimits)
dataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=dataprofile))
initialprofile <- constructinitialprofile(dataprofile,datalimits,cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))

#################################################
# cost function
#################################################

objective <- function(x){
steplength <- x[1]
waitingtime <- x[2]
parameters <- c(steplength=steplength,waitingtime=waitingtime,dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
times <- seq(0,days,length=2)
out <- lsode(initialprofile,times,nonlocalmodel,parameters,jacfunc=nonlocaljac,jactype="fullusr")
modelprofile <- out[dim(out)[1],(2:dim(out)[2])]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
modelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=roughmodelprofile))
cost <- modCost(obs=dataprofileframe,model=modelprofileframe,x="depth")
return(cost)
}

#################################################
# fit
#################################################

fit <- modFit(p=c(2,20),f=objective,lower=c(0,0))

#################################################
# run with fitted parameters
#################################################

times <- seq(0,days,length=2)
parameters <- c(steplength=fit$par[1],waitingtime=fit$par[2],dx=dx,k=k,flux=flux,fluxintroduction=fluxintroduction)
print(system.time(out <- lsode(initialprofile,times,nonlocalmodel,parameters,jacfunc=nonlocaljac,jactype="fullusr")))
modelprofile <- out[dim(out)[1],(2:dim(out)[2])]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)

summ <- summary(fit)




#################################################
# costfunction exploration
#################################################



parnumber <- 20
slist <- seq(0.5,2,length=parnumber)
wlist <- seq(5,20,length=parnumber)
costs <- matrix(data=NA,nrow=parnumber,ncol=parnumber)
dbs <- matrix(data=NA,nrow=parnumber,ncol=parnumber)


for(st in 1:length(slist)){
for(wa in 1:length(wlist)){
nldb <- slist[st]^2/(wlist[wa]*2)
x <- c(slist[st],wlist[wa])
obj <- objective(x)
cost <- sum(obj$residuals$res^2)
costs[st,wa] <- cost
dbs[st,wa] <- nldb
cat(st,wa,cost,"\n")
}
}



library(rgl)
library(RColorBrewer)



persp3d(slist,wlist,as.matrix(costs),xlab="steplength",ylab="waiting time",zlab="model cost",col="gray")
persp3d(slist,wlist,as.matrix(log(dbs)),xlab="steplength",ylab="waiting time",zlab="db",col="gray",zlim=c(0,1))

windows()
image(slist,wlist,as.matrix(costs),col=heat.colors(10),xlab="steplength",ylab="waiting time")
contour(slist,wlist,as.matrix(costs),add=TRUE,nlevels=20)

windows()
image(slist,wlist,as.matrix(dbs),col=heat.colors(10),xlab="steplength",ylab="waiting time")
contour(slist,wlist,as.matrix(dbs),add=TRUE,nlevels=20)















#################################################
# MCMC
#################################################

var0 <- summ$modVariance
covini <- summ$cov.scaled*2.4^2/2

mcmc <- modMCMC(p=coef(fit),f=objective,jump=covini,var0=var0,wvar0=1)
windows()
plot(mcmc)
windows()
pairs(mcmc)






