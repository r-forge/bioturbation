
### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### nonlocal bioturbation
### ============================================================================

nonlocal <- function (times, initprof, parameters = list(), ...)
{

## local functions
  constructtransitionmatrix <- function(steplength, waitingtime,
                                        slicenumber, dx)  {

    depth <- slicenumber * dx

    psil <- function(l) {
      return((1 / (steplength*sqrt(2*pi))) * exp(-(l^2)/(2*steplength^2)) )
    }
    
    startint <- -depth+dx/2
    endint <- depth-dx/2
    intlimits <- seq(startint,endint,by=dx)
    intvals <- NULL
    for (i in 1:(length(intlimits)-1) ) {
      intval <- integrate(psil,intlimits[i],intlimits[i+1])$value
      intvals <- c(intvals,intval)
    }
    intvals <- intvals / waitingtime
    intvals[slicenumber] <- -sum(intvals[-slicenumber])

    transitionmatrix <- matrix(data=NA,nrow=slicenumber,ncol=slicenumber)

    for ( Col in 1:slicenumber) {
      transitionmatrix[,Col] <- intvals[(slicenumber-Col+1):(2*slicenumber-Col)]
      rest <- rev(intvals[1:(slicenumber-Col)])
      transitionmatrix[(1:length(rest)),Col] <-
            transitionmatrix[(1:length(rest)),Col] + rest
    }
    return(transitionmatrix)
  }

  nonlocaljac <- function(t, profile, parms) {
    with(as.list(parms),{
      slicenumber <- length(profile)
      mat <- constructtransitionmatrix(steplength,waitingtime,slicenumber,dx)
      diag(mat) <- diag(mat) - k
      return(mat)
    })
  }


  nonlocalmodel <- function(t, profile, parms, transitionmatrix) {
    with(as.list(parms),{
      dprofile <- transitionmatrix %*% profile - k * profile
      fluxslices <- ceiling(fluxintroduction / dx)
      dprofile[1:fluxslices] <- dprofile[1:fluxslices] + flux / (dx*fluxslices)
      list(dprofile=c(dprofile))
    })
  }

## the default parameter values

  Parms <- c(steplength=1, waitingtime=20, dx=0.05, k=0,
             flux=0, fluxintroduction=0.5)

## check parameter inputs

  if (length(parameters)) {
    nms <- names(Parms)
    Parms[(namc <- names(parameters))]<-parameters
    if (length(noNms <- namc[!namc %in% nms]) > 0)
      warning("unknown names in parameters: ", paste(noNms, collapse = ", "))
  }

## running the model
  # 1. Create transition matrix and the jacobian function
  JAC <- nonlocaljac(0,initprof,Parms)

  JacFun <- function(t,pr,par,transitionmatrix) transitionmatrix

  # 2. run dynamically
  if (length(times) > 1)
    out <- lsode (y=initprof, times=times, func=nonlocalmodel,
           parms=Parms, jacfunc=JacFun, jactype="fullusr",
           transitionmatrix=JAC, ...)
  else if ( is.null(times)) {
    out <-  nonlocalmodel(0,initprof,parms=Parms,transitionmatrix=JAC)
    warning("times is NULL; returning RATE OF CHANGE")
  }
  # or estimate steady-state
  else
    out <- steady(y=initprof, time=times, func=nonlocalmodel,
           parms=Parms,jacfunc=JacFun, jactype="fullusr",
           transitionmatrix=JAC, ...)$y


  return(out)
  
}
