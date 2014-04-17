\name{PBD-package}
\alias{PBD-package}
\alias{PBD}
\docType{package}
\title{ Protracted birth-death model of diversification }
\description{ This package computes the (maximum) likelihood of the protracted speciation model for a given set of branching times
This package is a likelihood-based statistical package to estimate parameters
under diversity-dependent diversification models.\cr \cr
   New in version 0.9\cr
   - Bug fix for phylogeny with stem age\cr
}
\details{
\tabular{ll}{
Package: \tab PBD\cr
Type: \tab Package\cr
Version: \tab 0.9\cr
Date: \tab 2014-4-4\cr
License: \tab GPL-2\cr
}
pbd_loglik computes the likelihood of the protracted birth-death model of
diversification, given a set of parameters and a data set of phylogenetic branching times.

pbd_ML finds the parameters that maximizes the likelihood computed by pbd_loglik.
}
\author{
Rampal S. Etienne
Maintainer: Rampal S. Etienne <r.s.etienne@rug.nl>
}
\references{
- Etienne, R.S. & J. Rosindell 2012. Systematic Biology 61: 204-213.\cr
- Lambert, A., H. Morlon & R.S. Etienne 2014. Journal of Mathmematical Biology. In press\cr
- Etienne, R.S., H. Morlon & A. Lambert 2014. Evolution. In press.
}
\keyword{ package }
\seealso{ \code{DDD} }
\examples{
pbd_ML(1:10)
}
