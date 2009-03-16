




#################################################
# costfunction exploration
#################################################



library(rgl)
library(RColorBrewer)



steplength <- 1.4210526
waitingtime <- 12.668421

meth <- "rk4"
initialtimestep <- 3


parnumber <- 20
slist <- seq(0.5,3,length=parnumber)
wlist <- seq(0.1,20,length=parnumber)
costs <- matrix(data=NA,nrow=parnumber,ncol=parnumber)
dbs <- matrix(data=NA,nrow=parnumber,ncol=parnumber)


for(st in 1:length(slist)){
for(wa in 1:length(wlist)){

steplength <- slist[st]
waitingtime <- wlist[wa]

nldb <- steplength^2/(2*waitingtime)
slicenumber <- 200
dx <- 0.1
cakethickness <- 0.5
days <- 14
initialprofile <- constructinitialprofile(concprofile,limits,cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber * dx,by=dx))
parameters <- c(steplength=steplength,waitingtime=waitingtime,dx=dx)
times <- seq(0,days,length=2)
nlout <- rk(times=times,y=initialprofile,func=nonlocalmodel,parms=parameters,method=meth,hini=waitingtime)
nlfinalprofile <- nlout[dim(nlout)[1],-1]
nlmodelrough <- roughprofile(nlfinalprofile,limits,dx)
#windows()
maxval <- max(c(concprofile,initialprofile))
plot(initialmidpoints,initialprofile,ylim=c(0,maxval),type="l",xlab="depth",ylab="luminophore concentration",main="demo nonlocal model")
lines(concmidpoints,concprofile,col="red")
points(concmidpoints,concprofile,col="red")
lines(initialmidpoints,nlfinalprofile,col="green")
points(concmidpoints,nlmodelrough,col="green")
cost <- sum((nlmodelrough-concprofile)^2)
cat("steplength: ",steplength,"\nwaitingtime: ",waitingtime,"\ndb: ",nldb,"\ncost: ",cost,"\n\n",sep="")

costs[st,wa] <- cost
dbs[st,wa] <- nldb

}
}




#################################################
# plotting
#################################################


#costs <- read.table("costs.txt")

persp3d(slist,wlist,as.matrix(costs),xlab="steplength",ylab="waiting time",zlab="model cost",col="gray")

persp3d(slist,wlist,as.matrix(log(dbs)),xlab="steplength",ylab="waiting time",zlab="db",col="gray",zlim=c(0,1))


windows()
image(slist,wlist,as.matrix(costs),col=heat.colors(10),xlab="steplength",ylab="waiting time")
contour(slist,wlist,as.matrix(costs),add=TRUE,nlevels=20)

windows()
image(slist,wlist,as.matrix(dbs),col=heat.colors(10),xlab="steplength",ylab="waiting time")
contour(slist,wlist,as.matrix(dbs),add=TRUE,nlevels=20)







