\name{pbd_sim}
\alias{pbd_sim}
\title{ Function to simulate the protracted speciation process }
\description{ Simulating the protracted speciation process using the Doob-Gillespie algorithm.
This function differs from pbd_sim_cpp that 1) it does not require that the speciation-initiation
rate is the same for good and incipient species, and 2) that it simulates the exact protracted
speciation process, and not the approximation made by the coalescent point process. }
\usage{ pbd_sim(
   pars,
   age,
   soc = 2
) }
\arguments{
  \item{pars}{Vector of parameters: \cr \cr
    \code{pars[1]} corresponds to b_1, the speciation-initiation rate of good species \cr
    \code{pars[2]} corresponds to la_1, the speciation-completion rate \cr
    \code{pars[3]} corresponds to b_2, the speciation-initiation rate of incipient species \cr 
    \code{pars[4]} corresponds to mu_1, the extinction rate of good species \cr
    \code{pars[5]} corresponds to mu_2, the extinction rate of incipient species \cr
  }
  \item{age}{ Sets the age for the simulation }
  \item{soc}{ Sets whether this age is the stem (1) or crown (2) age }
}
\value{
  \item{ out }{ A list with the following three elements:
   The first element is the tree of extant species in phylo format \cr
   The second element is a tree with one sample per species in phylo format \cr
   The third element is a matrix of all events in the simulation where \cr
   - the first column is the incipient-level label of a species \cr
   - the second column is the incipient-level label of the parent of the species \cr
   - the third column is the time at which a species is born as incipient species\cr
   - the fourth column is the time of speciation-completion of the species \cr
   If the fourth element equals -1, then the species is still incipient.   
   - the fifth column is the time of extinction of the species \cr
   If the fifth element equals -1, then the species is still extant.
   - The sixth column is the species-level label of the species \cr
  }
}
\author{ Rampal S. Etienne }
\seealso{
 \code{\link{pbd_sim_cpp}}
}
\examples{ pbd_sim(c(0.2,1,0.2,0.1,0.1),15) }
\keyword{models}
