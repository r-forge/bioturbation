\name{turbo-package}
\alias{turbo-package}
\alias{turbo}
\docType{package}
\title{
  Functions for fitting bioturbation models to tracer data
}
\description{
  Contains several bioturbation models
  ....
}

\details{
  \tabular{ll}{
    Package: \tab turbo\cr
    Type: \tab Package\cr
    Version: \tab 1.3.1\cr
    Date: \tab 2009-04-01\cr
    License: \tab GNU Public License 2 or above\cr
  }

}

\author{
  Pieter Provoost (Maintainer)

  Karline Soetaert
  
  Filip Meysman
}
\seealso{

  \code{\link{diffusive}}, diffusive model

  \code{\link{nonlocal}},  nonlocal model

  \code{\link{Luminophore}}, a test dataset.
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(diffusive)
example(nonlocal)

## run demos
demo("nonlocal") # nog te maken

## open the directory with source code of demos
browseURL(paste(system.file(package="turbo"), "/demo", sep=""))

## show package vignette - nog te maken
vignette("turbo")
edit(vignette("turbo"))
browseURL(paste(system.file(package="turbo"), "/doc", sep=""))
}
}

\keyword{ package }
