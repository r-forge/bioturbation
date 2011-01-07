### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### nonlocal bioturbation
### ============================================================================

nonlocal <- function (times, initprof, parameters = list(), ...) {
	
	## transition matrix (slow)
	constructtransitionmatrix_slow <- function(steplength, waitingtime, slicenumber, dx)  {
		depth <- slicenumber * dx
		psil <- function(l) {
			return((1 / (steplength*sqrt(2*pi))) * exp(-(l^2)/(2*steplength^2)) )
		}
		transitionmatrix <- matrix(data=NA,nrow=slicenumber,ncol=slicenumber)
		for (col in 1:slicenumber) {
			for (row in 1:slicenumber) {
				transitionmatrix[row,col] <- integrate(psil,(row-col-0.5)*dx,(row-col+0.5)*dx)$value / waitingtime
			}
			# folding top
			for (row in 1:slicenumber) {
				transitionmatrix[row,col] <- transitionmatrix[row,col] + integrate(psil,((-row+1)-col-0.5)*dx,((-row+1)-col+0.5)*dx)$value / waitingtime
			}
			# folding bottom
			for (row in 1:slicenumber) {
				transitionmatrix[row,col] <- transitionmatrix[row,col] + integrate(psil,((2*slicenumber+1-row)-col-0.5)*dx,((2*slicenumber+1-row)-col+0.5)*dx)$value / waitingtime
			}
			transitionmatrix[col,col] <- - sum(transitionmatrix[,col][-col])
		}
   		return(transitionmatrix)
	}

	## transition matrix (faster)
	constructtransitionmatrix <- function(steplength, waitingtime, slicenumber, dx)  {
		depth <- slicenumber * dx
		psil <- function(l) {
			return((1 / (steplength*sqrt(2*pi))) * exp(-(l^2)/(2*steplength^2)) )
		}
		values <- NULL
		for (i in seq(-(slicenumber+0.5)*dx,(slicenumber-0.5)*dx,by=dx)) {
			values <- c(values,integrate(psil,i,i+dx)$value / waitingtime)
		}
		transitionmatrix <- matrix(data=NA,nrow=slicenumber,ncol=slicenumber)
		for (col in 1:slicenumber) {
			transitionmatrix[,col] <- values[(slicenumber+1-col+1):(2*slicenumber-1+1-col+1)]	

			# folding top
			rest <- values[1:(slicenumber+1-col+1-1)]
			transitionmatrix[1:length(rest),col] <- transitionmatrix[1:length(rest),col] + rev(rest)

			# folding bottom
			rest <- values[(2*slicenumber-1+1-col+1+1):(2*slicenumber+1)]	
			transitionmatrix[(slicenumber-length(rest)+1):slicenumber,col] <- transitionmatrix[(slicenumber-length(rest)+1):slicenumber,col] + rev(rest)

			transitionmatrix[col,col] <- - sum(transitionmatrix[,col][-col])
		}
   		return(transitionmatrix)
	}







	## jacobian
	nonlocaljac <- function(t, profile, parms) {
		with(as.list(parms),{
			slicenumber <- length(profile)
			mat <- constructtransitionmatrix(steplength,waitingtime,slicenumber,dx)
			diag(mat) <- diag(mat) - k
			return(mat)
		})
	}

	## the default parameter values
	nonlocalmodel <- function(t, profile, parms, transitionmatrix) {
		with(as.list(parms),{
			dprofile <- transitionmatrix %*% profile
			fluxslices <- ceiling(fluxintroduction / dx)
			dprofile[1:fluxslices] <- dprofile[1:fluxslices] + flux / (dx*fluxslices)
			list(dprofile=c(dprofile))
		})
	}

	## the default parameter values
	Parms <- c(steplength=1, waitingtime=20, dx=0.05, k=0, flux=0, fluxintroduction=0.5)

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
	if (length(times) > 1) {
		out <- lsode (y=initprof, times=times, func=nonlocalmodel, parms=Parms, jacfunc=JacFun, jactype="fullusr", transitionmatrix=JAC, ...)
	} else if (is.null(times)) {
		out <- nonlocalmodel(0,initprof,parms=Parms,transitionmatrix=JAC)
		warning("times is NULL; returning RATE OF CHANGE")
	}
	
	# or estimate steady-state
	else {
		out <- steady(y=initprof, time=times, func=nonlocalmodel, parms=Parms,jacfunc=JacFun, jactype="fullusr", transitionmatrix=JAC, ...)$y
	}
	
	return(out)
}
