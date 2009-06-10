
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

# read data

datafile <- "data.txt"
profilename <- "2"
data <- read.table(datafile,header=TRUE)
dataframe <- subset(data,profile==profilename)
datalimits <- unique(c(dataframe$start,dataframe$end))
datamidpoints <- midpoints(datalimits)
dataprofile <- numtoconc(dataframe$lum,datalimits)
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

# fit nonlocal model

fit <- modFit(p=c(2,20),f=nlobjective,lower=c(0,0))
fitsl <- fit$par[1]
fitwt <- fit$par[2]
fitnldb <- fitsl^2/(2*fitwt)

times <- c(0,days)
final <- nonlocal(times,initprofile,list(waitingtime=fitwt,steplength=fitsl,dx=dx))

# plot data and nonlocal model

plot(datamidpoints,dataprofile,type="n",xlab="depth",ylab="luminophore concentration")
lines(initmidpoints,final[2,2:(slicenumber+1)],col="red",lwd=2)

# fit diffusive model

fit <- modFit(p=c(0.01),f=diffobjective,lower=c(0))
fitdb <- fit$par[1]

times <- c(0,days)
final <- diffusive(times,initprofile,list(db=fitdb,dx=dx))
lines(initmidpoints,final[2,2:(slicenumber+1)],col="green",lwd=2)

# plot data

points(datamidpoints,dataprofile,pch=21,bg="white")

# legend

legend(10,max(dataprofile)/2,pch=c(NA,NA,21),lty=c("solid","solid","blank"),lwd=c(2,2,1),col=c("red","green","black"),legend=c("nonlocal model","diffusive model","data"),bty="n")

# results

cat("\nfitted step length: ",round(fitsl,digits=3),"\nfitted waiting time: ",round(fitwt,digits=3),"\nfitted db: ",round(fitdb,digits=3),"\nfitted nldb: ",round(fitnldb,digits=3),"\n\n",sep="")






