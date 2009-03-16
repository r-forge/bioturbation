\name{diffusive}
\alias{diffusive}
\title{diffusive model}
\description{
  Runs a diffusive model dynamically
}

\usage{
diffusive(times, initprof, parameters = list(), ...)}
\arguments{
  \item{times }{times at which explicit estimates for the profile are
    desired.  The first value in \code{times} must be the initial time.
  }
  \item{initprof }{the initial (state) values for the ODE system.
  }
  \item{parameters }{vector or list of parameters used in the nonlocal function
  }
  \item{... }{additional arguments passed to the integrator.
  }
}
\value{
  A matrix with up to as many rows as elements in times and as many
  columns as elements in \code{initprof}.

  The output will have the attributes \code{istate}, and \code{rstate},
  two vectors with several useful elements. The first element of istate
  returns the conditions under which the last call to the integrator
  returned. Normal is \code{istate = 2}.  If \code{verbose = TRUE}, the
  settings of istate and rstate will be written to the screen. See the
  help for the selected integrator for details.
}
\author{
  Pieter Provoost <p.provoost@nioo.knaw.nl>
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\references{
Meysman etc..
}
\details{
This implements the diffusive type of bioturbation as in ....
}
\keyword{manip}
\examples{
times <- 0:200
yini <- rep(1,200)
out <- diffusive(times,yini)
plot(out[nrow(out),-1],ylim=c(0,1))

out2<- diffusive(times,yini, parameters = list(k=0.005))
lines(out2[nrow(out2),-1])

}

