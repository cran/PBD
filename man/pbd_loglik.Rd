\name{pbd_loglik}
\alias{pbd_loglik}
\title{Loglikelihood for protracted birth-death model of diversification}
\description{
pbd_loglik computes the loglikelihood of the parameters of the protracted speciation model given a set of branching times and number of missing species
}
\usage{
pbd_loglik(
    pars1,
    pars1f = c(function(t,pars) {pars[1]},function(t,pars) {pars[2]},
function(t,pars) {pars[3]}, function(t,pars) {pars[4]}),
    pars2 = c(1,1,2,1,"lsoda",0,0),
    brts,
    missnumspec = 0
)}
\arguments{
  \item{pars1}{Vector of parameters: \cr \cr
    \code{pars1[1]} corresponds to b (= la_1 in Etienne & Rosindell R2012) = speciation initiation rate \cr
    \code{pars1[2]} corresponds to mu_1 (= mu_g in Etienne & Rosindell 2012) = extinction rate of good species \cr
    \code{pars1[3]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) = speciation completion rate \cr
    \code{pars1[4]} corresponds to mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr
When rates depend on time this time dependence should be specified in pars1f and pars1 then becomes the parameters used in pars1f \cr \cr
}
 \item{pars1f}{Vector of functions how the rates depend on time, default functions are constant functions of the parameters in pars1: \cr \cr
    \code{pars1f[1]} corresponds to time-dependence of b (= la_1 in Etienne & Rosindell R2012) = speciation initiation rate \cr
    \code{pars1f[2]} corresponds to time-dependence of mu_1 (= mu_g in Etienne & Rosindell 2012) = extinction rate of good species \cr
    \code{pars1f[3]} corresponds to tiem-dependence of la_1 (= la_2 in Etienne & Rosindell 2012) = speciation completion rate \cr
    \code{pars1f[4]} corresponds to time-dependence of mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr \cr
}
  \item{pars2}{Vector of model settings: \cr \cr
    \code{pars2[1]} {set the conditioning on non-extinction of the clade (1) or not (0) \cr \cr}
    \code{pars2[2]} {sets whether the likelihood is for the branching times (0) or the phylogeny (1) \cr \cr}
    \code{pars2[3]} {sets whether the first element of the branching times is the stem (1) or the crown (2) age \cr \cr}
    \code{pars2[4]} {sets whether the parameters and likelihood should be shown on screen (1) or not (0) \cr \cr}
    \code{pars2[5]} {sets which method should be used in the ode-solver. Default is 'lsoda'. See package deSolve for details. \cr \cr}
    \code{pars2[6]}{Sets the lower bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)\cr \cr}
    \code{pars2[7]}{Sets the upper bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)\cr \cr}  
}
  \item{brts}{A set of branching times of a phylogeny, all positive }
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny }
}
\value{The loglikelihood}
\author{ Rampal S. Etienne }
\seealso{
 \code{\link{pbd_ML}}
}
\examples{ pbd_loglik(pars1 = c(0.2,0.1,1,0.1), pars2 = c(1,1,2,0,"lsoda"),brts = 1:10) }
\keyword{ models }