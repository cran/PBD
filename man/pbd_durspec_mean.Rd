\name{pbd_durspec_mean}
\alias{pbd_durspec_mean}
\title{Mean duration of speciation under protracted birth-death model of diversification}
\description{
pbd_durspec_mean computes the mean duration of speciation under the protracted speciation model for a given set of parameters
}
\usage{
pbd_durspec_mean(
    pars
)}
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to b (= la_3 in Etienne & Rosindell R2012) = speciation initiation rate \cr
    \code{pars[2]} corresponds to la_1 (= la_2 in Etienne & Rosindell 2012) = speciation completion rate \cr
    \code{pars[3]} corresponds to mu_2 (= mu_i in ER2012) = extinction rate of incipient species \cr
}
}
\value{ The expected duration of speciation }
\author{ Rampal S. Etienne }
\seealso{
 \code{\link{pbd_durspec_density}}\cr
 \code{\link{pbd_durspec_cumdensity}}\cr
 \code{\link{pbd_durspec_mode}}\cr
 \code{\link{pbd_durspec_quantile}}\cr
 \code{\link{pbd_durspec_moment}}\cr
 \code{\link{pbd_durspec_var}}
}
\examples{ pbd_durspec_mean(pars = c(0.5,0.3,0.1)) }
\keyword{ models }