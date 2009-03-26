\name{midpoints}
\alias{midpoints}
\title{midpoints}
\description{
  Slice midpoints from slice limits
}

\usage{
midpoints(limits)}
\arguments{
  \item{limits }{A vector containing slice limits.
  }
}
\value{
  A vector containing slice midpoints.
}
\author{
  Pieter Provoost <p.provoost@nioo.knaw.nl>
}
\examples{
limits <- c(seq(0,10,by=1),seq(12,20,by=2))
mpoints <- midpoints(limits)
plot(rep(0,length(limits)),limits)
points(rep(0,length(mpoints)),mpoints,col="red")

}

