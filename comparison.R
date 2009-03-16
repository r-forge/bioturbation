
#################################################
# comparison local and nonlocal models
#################################################



steplength <- 2
waitingtime <- 1
days <- 0.5

db <- steplength^2/(2*waitingtime)
slicenumber <- 100
dx <- 0.1
cakethickness <- 0.5

times <- seq(0,days,length=100)

initialprofile <- constructinitialprofile2(cakethickness,slicenumber,dx)
initialmidpoints <- midpoints(seq(0,slicenumber*dx,by=dx))

parameters <- c(db=db,dx=dx)
library("deSolve")
out <- ode.band(times=times,y=initialprofile,func=diffusivemodel,parms=parameters,nspec=1)
finalprofile <- out[dim(out)[1],-1]

parameters <- c(steplength=steplength,waitingtime=waitingtime,dx=dx)
nlout <- ode.band(times=times,y=initialprofile,func=nonlocalmodel,parms=parameters,nspec=1)
nlfinalprofile <- nlout[dim(out)[1],-1]

windows()
maxval <- max(initialprofile)
plot(initialmidpoints,initialprofile,ylim=c(0,maxval),main="comparison local and nonlocal model",col="gray",lty="dotted",type="l",xlab="depth",ylab="luminophore concentration")
lines(initialmidpoints,finalprofile,col="green")
lines(initialmidpoints,nlfinalprofile,col="red")
