\name{numtoconc}
\alias{numtoconc}
\title{numtoconc}
\description{
  Luminophore concentrations from luminophore numbers and slice limits
}

\usage{
numtoconc(values, limits)}
\arguments{
  \item{values }{A vector containing luminophore numbers.
  }
  \item{limits }{A vector containing slice limits.
  }
}
\value{
  A vector containing luminophore concentrations.
}
\author{
  Pieter Provoost <p.provoost@nioo.knaw.nl>
}
\examples{
limits <- c(seq(0,10,by=1),seq(12,20,by=2))
numbers <- exp(seq(3,1,length.out=length(limits)))
concentrations <- numtoconc(numbers, limits)                          
plot(numbers,concentrations)

}

