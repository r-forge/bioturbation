
# install package turbo from R-Forge:
# install.packages("turbo",repos="http://R-Forge.R-project.org")

library(turbo)
library(FME)
graphics.off()

# settings

slicenumber <- 200
dx <- 0.1
cakethickness <- 0.5
days <- 14

# transformation

transf <- function(prof){
return(log10(prof+1))
}

# read data

datafile <- "data.txt"
profilename <- "a"
data <- read.table(datafile,header=TRUE)
dataframe <- subset(data,profile==profilename)
datalimits <- unique(c(dataframe$start,dataframe$end))
datamidpoints <- midpoints(datalimits)
dataprofile <- numtoconc(dataframe$lum,datalimits)
dataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=dataprofile))
transfdataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=transf(dataprofile)))

# create initial profile

initprofile <- initialprofile(dataprofile,datalimits,slicenumber,dx,cakethickness)
initlimits <- seq(0,slicenumber*dx,by=dx)
initmidpoints <- midpoints(initlimits)

# cost function nonlocal model

nlobjective <- function(x){
steplength <- x[1]
waitingtime <- x[2]
times <- c(0,days)
final <- nonlocal(times,initprofile,list(waitingtime=waitingtime,steplength=steplength,dx=dx))
modelprofile <- final[2,2:(slicenumber+1)]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
modelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=roughmodelprofile))
cost <- modCost(obs=dataprofileframe,model=modelprofileframe,x="depth")
return(cost)
}

# cost function diffusive model

diffobjective <- function(x){
db <- x
times <- c(0,days)
final <- diffusive(times,initprofile,list(db=db,dx=dx))
modelprofile <- final[2,2:(slicenumber+1)]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
modelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=roughmodelprofile))
cost <- modCost(obs=dataprofileframe,model=modelprofileframe,x="depth")
return(cost)
}

# cost function diffusive model transformed

transfdiffobjective <- function(x){
db <- x
times <- c(0,days)
final <- diffusive(times,initprofile,list(db=db,dx=dx))
modelprofile <- final[2,2:(slicenumber+1)]
roughmodelprofile <- roughprofile(modelprofile,datalimits,dx)
transfmodelprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=transf(roughmodelprofile)))
cost <- modCost(obs=transfdataprofileframe,model=transfmodelprofileframe,x="depth")
return(cost)
}

# fit nonlocal model

fit <- modFit(p=c(2,20),f=nlobjective,lower=c(0,0))
fitsl <- fit$par[1]
fitwt <- fit$par[2]
fitnldb <- fitsl^2/(2*fitwt)

times <- c(0,days)
final <- nonlocal(times,initprofile,list(waitingtime=fitwt,steplength=fitsl,dx=dx))
plot(datamidpoints,dataprofile,type="n")
lines(initmidpoints,final[2,2:(slicenumber+1)],col="red",lwd=2)

# fit diffusive model

fit <- modFit(p=c(0.01),f=diffobjective,lower=c(0))
fitdb <- fit$par[1]

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitdb,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="green",lwd=2)

# plot diffusive with nonlocal db

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitnldb,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="red",lwd=1)

# fit diffusive model transformed

fit <- modFit(p=c(0.01),f=transfdiffobjective,lower=c(0))
fitdbtransf <- fit$par[1]

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitdbtransf,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="blue",lwd=2)

points(datamidpoints,dataprofile,pch=21,bg="white")

# results

cat("fitted step length: ",round(fitsl,digits=3),"\nfitted waiting time: ",round(fitwt,digits=3),"\nfitted db: ",round(fitdb,digits=3),"\nfitted db (transformed): ",round(fitdbtransf,digits=3),"\nfitted nldb: ",round(fitnldb,digits=3),"\n",sep="")







