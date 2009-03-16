

psil <- function(l){
return((1/(steplength*sqrt(2*pi)))*exp(-(l^2)/(2*steplength^2)))
}
steplength <- 1
windows()
curve(psil,from=-5,to=5,main="step length distribution",xlab="cm",ylab="",lwd=2,ylim=c(0,0.9))
steplength <- 0.5
curve(psil,from=-5,to=5,add=TRUE)
steplength <- 1.5
curve(psil,from=-5,to=5,add=TRUE,lty="dashed")
legend(2,0.6,legend=c(expression(paste(sigma," = 0.5 cm")),expression(paste(sigma," = 1 cm")),expression(paste(sigma," = 1.5 cm"))),lwd=c(1,2,1),bty="n",lty=c("solid","solid","dashed"))





windows()
steplength <- 1
curve(psil,from=-3,to=3,main="step length distribution",xlab="cm",ylab="")
start <- 0.3
end <- 0.5
x <- seq(start,end,length=100)
y <- psil(x)
polygon(c(start,x,end),c(0,y,0),col="#cccccc")

