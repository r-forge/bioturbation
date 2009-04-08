

library(turbo)
library(FME)
graphics.off()



datafile <- "z21.txt"





# settings

days <- 14
slicenumber <- 200
dx <- 0.1
cakethickness <- 0.5

# read data

data <- read.table(datafile,header=TRUE)
datalimits <- c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,7,8,9,10)
datamidpoints <- midpoints(datalimits)
dataprofile <- numtoconc(data$blue,datalimits)
dataprofileframe <- as.data.frame(cbind(depth=datamidpoints,concentration=dataprofile))

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

# plot

plot(datamidpoints,dataprofile,type="n",ylim=c(0,max(dataprofile)*1.2))

# fit nonlocal model

fit <- modFit(p=c(2,20),f=nlobjective,lower=c(0,0))
fitsl <- fit$par[1]
fitwt <- fit$par[2]
fitnldb <- fitsl^2/(2*fitwt)

times <- c(0,days)
final <- nonlocal(times,initprofile,list(waitingtime=fitwt,steplength=fitsl,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="red",lwd=1)

# fit diffusive model

fit <- modFit(p=c(0.01),f=diffobjective,lower=c(0))
fitdb <- fit$par[1]

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitdb,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="green",lwd=1)

# plot diffusive with nonlocal db

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitnldb,dx=dx))
#lines(initmidpoints,final[2,2:(slicenumber+1)],col="red",lwd=1)

# plot data

points(datamidpoints,dataprofile,pch=21,bg="white")

# results

cat("fitted step length: ",round(fitsl,digits=3),"\nfitted waiting time: ",round(fitwt,digits=3),"\nfitted db: ",round(fitdb,digits=5),"\nfitted nldb: ",round(fitnldb,digits=5),"\n",sep="")







