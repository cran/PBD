\name{pbd_sim}
\alias{pbd_sim}
\title{ Function to simulate the protracted speciation process }
\description{ Simulating the protracted speciation process using the Doob-Gillespie algorithm.
This function differs from pbd_sim_cpp that 1) it does not require that the speciation-initiation
rate is the same for good and incipient species, and 2) that it simulates the exact protracted
speciation process, and not the approximation made by the coalescent point process. }
\usage{ pbd_sim(
   pars,
   age
) }
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to b_1, the speciation-initiation rate of good species \cr
    \code{pars[2]} corresponds to la_1, the speciation-completion rate \cr
    \code{pars[5]} corresponds to b_2, the speciation-initiation rate of incipient species \cr 
    \code{pars[3]} corresponds to mu_1, the extinction rate of good species \cr
    \code{pars[4]} corresponds to mu_2, the extinction rate of incipient species \cr
  }
  \item{age}{ Sets the crown age for the simulation }
}
\value{
  \item{ tree }{ The tree of extant species in phylo format  }
}
\author{ Rampal S. Etienne }
\seealso{
 \code{\link{pbd_sim_cpp}}
}
\examples{ pbd_sim(c(0.2,1,0.2,0.1,0.1),15) }
\keyword{models}
