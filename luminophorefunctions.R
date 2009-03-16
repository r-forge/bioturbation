
#################################################
# 
# luminophore profile analysis functions
# Pieter Provoost
#
#################################################


#################################################
# function nonlocaljac
#################################################


nonlocaljac <- function(t,profile,parameters){
with(as.list(parameters),{
slicenumber <- length(profile)
mat <- constructtransitionmatrix(steplength,waitingtime,slicenumber,dx)
diag(mat) <- diag(mat) - k
return(mat)
})
}



#################################################
# function diffusivemodel
#################################################


diffusivemodel <- function(t,profile,parameters){
with(as.list(parameters),{
diffusiveflux <- db * c(0,diff(c(profile,0))) / dx
dprofile <- diff(diffusiveflux) / dx - k * profile
fluxslices <- ceiling(fluxintroduction / dx)
dprofile[1:fluxslices] <- dprofile[1:fluxslices] + flux / (dx*fluxslices)
list(c(dprofile))
})
}

#################################################
# function nonlocalmodel
#################################################

nonlocalmodel <- function(t,profile,parameters){
with(as.list(parameters),{
slicenumber <- length(profile)
transitionmatrix <- constructtransitionmatrix(steplength,waitingtime,slicenumber,dx)
dprofile <- transitionmatrix %*% profile - k * profile
fluxslices <- ceiling(fluxintroduction / dx)
dprofile[1:fluxslices] <- dprofile[1:fluxslices] + flux / (dx*fluxslices)
list(c(dprofile))
})
}



#################################################
# function constructtransitionmatrix
#################################################

constructtransitionmatrix <- function(steplength,waitingtime,slicenumber,dx){
depth <- slicenumber * dx
psil <- function(l){
return((1/(steplength*sqrt(2*pi)))*exp(-(l^2)/(2*steplength^2)))
}
startint <- -depth+dx/2
endint <- depth-dx/2
intlimits <- seq(startint,endint,by=dx)
intvals <- NULL
for(i in 1:(length(intlimits)-1)){
intval <- integrate(psil,intlimits[i],intlimits[i+1])$value
intvals <- c(intvals,intval)
}
intvals <- intvals / waitingtime
intvals[slicenumber] <- -sum(intvals[-slicenumber])
transitionmatrix <- matrix(data=NA,nrow=slicenumber,ncol=slicenumber)
for(column in 1:slicenumber){
transitionmatrix[,column] <- intvals[(slicenumber-column+1):(2*slicenumber-column)]
rest <- rev(intvals[1:(slicenumber-column)])
transitionmatrix[(1:length(rest)),column] <- transitionmatrix[(1:length(rest)),column] + rest
}
return(transitionmatrix)
}

#################################################
# function constructinitialprofile2
# from cakethickness, concentration 1
#################################################

constructinitialprofile2 <- function(cakethickness,slicenumber,dx){
initialprofile <- rep(0,slicenumber)
cakeslices <- round(cakethickness / dx)
initialprofile[1:cakeslices] <- 1
return(initialprofile)
}

#################################################
# function constructinitialprofile
# from data profile
#################################################

constructinitialprofile <- function(concprofile,limits,cakethickness,slicenumber,dx){
initialprofile <- rep(0,slicenumber)
cakeslices <- round(cakethickness / dx)
thickness <- diff(limits)
totalnumber <- sum(concprofile * thickness)
initialprofile[1:cakeslices] <- totalnumber / cakeslices / dx
return(initialprofile)
}

#################################################
# function numtoconc
# calculates a concentration profile based on luminophore numbers
#################################################

numtoconc <- function(values, limits){
numtoconc <- NULL
for(slice in 1:length(values)){
numtoconc <- c(numtoconc, values[slice] / (limits[slice+1] - limits
[slice]))
}
return(numtoconc)
}


#################################################
# function midpoints
# calculates slice midpoints from slice limits
#################################################

midpoints <- function(limits){
midpoints <- NULL
for(point in 1:(length(limits)-1)){
midpoints <- c(midpoints,(limits[point]+limits[point+1])/2)
}
return(midpoints)
}




#################################################
# function roughprofile
# calculates a (rough) profile from a model profile with regular intervals
#################################################

roughprofile <- function(profile, limits, dx){
roughprofile <- NULL
mpoints <- seq(dx/2,length(profile)*dx,by=dx)
for(slice in 1:(length(limits)-1)){
slicedata <- subset(profile,mpoints > limits[slice] & mpoints < limits
[slice+1])
roughprofile <- c(roughprofile, mean(slicedata))
}
return(roughprofile)
}



#################################################
# function modelcost
#################################################

modelcost <- function(profile,dataprofile,method="wssr"){
if(method=="ssr"){
cost <- sum((dataprofile-profile)^2)
}
if(method=="wssr"){
cost <- sum(((dataprofile-profile)^2)/(dataprofile))
}
if(method=="m1"){
cost <- sum((dataprofile-profile)^2)
}
if(method=="m2"){
cost <- sum((log(dataprofile)-log(profile))^2)
}
return(cost)
}



#################################################
# function modelcostprofile
#################################################

modelcostprofile <- function(profile,dataprofile,method="wssr"){
if(method=="ssr"){
cost <- (dataprofile-profile)^2
}
if(method=="wssr"){
cost <- ((dataprofile-profile)^2)/(dataprofile)
}
if(method=="m1"){
cost <- (dataprofile-profile)^2
}
if(method=="m2"){
cost <- (log(dataprofile)-log(profile))^2
}
return(cost)
}






