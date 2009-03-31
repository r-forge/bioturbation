
# install package turbo from R-Forge:
# install.packages("turbo",repos="http://R-Forge.R-project.org")

library(turbo)
library(FME)

# settings

slicenumber <- 200
dx <- 0.1
cakethickness <- 0.5
days <- 14

# read data

datafile <- "data.txt"
profilename <- "a"
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

# run model once

times <- c(0,10)
final <- nonlocal(times,initprofile,list(waitingtime=10,steplength=1,dx=dx))
plot(initmidpoints,final[2,2:(slicenumber+1)])

# cost function

objective <- function(x){
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

# fit model




