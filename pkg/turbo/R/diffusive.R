
### ============================================================================
###
### bioturbation profile analysis functions
### Pieter Provoost - Karline Soetaert
###
### ============================================================================

### ============================================================================
### diffusive bioturbation
### ============================================================================

diffusive <- function (times, initprof, parameters = list(), ...) {

  ## local functions

  diffusivemodel <- function(t, profile, parms) {
    with (as.list(parms), {
      diffusiveflux <- db * c(0,diff(c(profile,0))) / dx
      dprofile <- diff(diffusiveflux) / dx - k * profile
      fluxslices <- ceiling(fluxintroduction / dx)
      dprofile[1:fluxslices] <- dprofile[1:fluxslices] + flux / (dx*fluxslices)
      list(c(dprofile))
    })
  }

## the default parameter values

  Parms <- c(db=0.01, dx=0.05, k=0,
             flux=0, fluxintroduction=0.5)

## check parameter inputs

  if (length(parameters)) {
    nms <- names(Parms)
    Parms[(namc <- names(parameters))]<-parameters
    if (length(noNms <- namc[!namc %in% nms]) > 0)
      warning("unknown names in parameters: ", paste(noNms, collapse = ", "))
  }

## running the model

  # 2. run dynamically
  if (length(times) > 1)
    ode.band(y=initprof, times=times, func=diffusivemodel,
        parms=Parms, nspec =1, ...)
  # or estimate steady-state
  else
    steady.band(y=initprof, times=times, func=diffusivemodel,
             parms=Parms, nspec =1, ...)$y

}

