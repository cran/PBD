\name{pbd_LR}
\alias{pbd_LR}
\title{Bootstrap likelihood ratio test of protracted birth-death model of diversification
}
\description{
   This function computes the maximum likelihood and the associated estimates of the parameters of a protracted birth-death model of diversification for a given set of phylogenetic branching times. It then performs a bootstrap likelihood ratio test of the protracted birth-death (PBD) model against the constant-rates (CR) birth-death model. Finally, it computes the power of this test.
}
\usage{
pbd_LR(
  brts,
  initparsoptPBD,
  initparsoptCR,
  missnumspec,
  outputfilename = NULL,
  seed = 42,
  endmc = 1000,
  alpha = 0.05,
  plotit = TRUE,
  parsfunc = c(function(t,pars) {pars[1]},
  function(t,pars) {pars[2]},
  function(t,pars) {pars[3]},
  function(t,pars) {pars[4]}),
  cond = 1,
  btorph = 1,
  soc = 2,
  methode = 'lsoda',
  n_low = 0,
  n_up = 0,
  tol = c(1E-6,1E-6,1E-6),
  maxiter = 2000,
  optimmethod = 'subplex',
  verbose = FALSE
)
}
\arguments{
  \item{brts}{A set of branching times of a phylogeny, all positive }
  \item{initparsoptPBD}{The initial values of the parameters that must be optimized for the protracted birth-death (PBD) model: b, mu and lambda}
  \item{initparsoptCR}{The initial values of the parameters that must be optimized for the constant-rates (CR) model: b and mu}
  \item{missnumspec}{The number of species that are in the clade but missing in the phylogeny }
  \item{outputfilename}{The name (and location) of the file where the output will be saved. Default is no save.}
  \item{seed}{The seed for the pseudo random number generator for simulating the bootstrap data}
  \item{endmc}{The number of bootstraps}
  \item{alpha}{The significance level of the test}
  \item{plotit}{Boolean to plot results or not}
  \item{parsfunc}{Specifies functions how the rates depend on time, default functions are constant functions}
  \item{cond}{ Conditioning: \cr
    cond == 0 : conditioning on stem or crown age \cr
    cond == 1 : conditioning on stem or crown age and non-extinction of the phylogeny \cr
    cond == 2 : conditioning on stem or crown age and on the total number of extant taxa (including missing species) \cr
    cond == 3 : conditioning on the total number of extant taxa (including missing species) \cr
   Note: cond == 3 assumes a uniform prior on stem age, as is the standard in constant-rate birth-death models, see
e.g. D. Aldous & L. Popovic 2004. Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J. Theor. Biol. 261: 58-66.

  }
  \item{btorph}{Sets whether the likelihood is for the branching times (0) or the phylogeny (1)}
  \item{soc}{Sets whether stem or crown age should be used (1 or 2) }
  \item{methode}{The numerical method used to solve the master equation, such as 'lsoda' or 'ode45'.}
  \item{n_low}{Sets the lower bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)}
  \item{n_up}{Sets the upper bound of the number of species on which conditioning should be done when cond = 2. Set this to 0 when conditioning should be done on precisely the number of species (default)}
  \item{tol}{Sets the tolerances in the optimization. Consists of: \cr
   reltolx = relative tolerance of parameter values in optimization \cr
   reltolf = relative tolerance of function value in optimization \cr
   abstolx = absolute tolerance of parameter values in optimization
  }
  \item{maxiter}{Sets the maximum number of iterations in the optimization}
  \item{optimmethod}{Method used in optimization of the likelihood. Current default is 'subplex'. Alternative is 'simplex'.}
  \item{verbose}{if TRUE, explanatory text will be shown}
}
\details{
The output is a list with 3 elements:
}
\value{
 \item{brtsCR}{a list of sets of branching times generated under the constant-rates model using the ML parameters under the CR model}
 \item{brtsDD}{a list of sets of branching times generated under the protracted birth-death model using the ML parameters under the PBD model}
 \item{out}{a dataframe with the parameter estimates and maximum likelihoods for protracted birth-death and constant-rates models
    \code{$model} - the model used to generate the data. 0 = unknown (for real data), 1 = CR, 2 = PBD \cr
    \code{$mc} - the simulation number for each model \cr
    \code{$b_CR} - speciation rate estimated under CR \cr
    \code{$mu_CR} - extinction rate estimated under CR \cr
    \code{$LL_CR} - maximum likelihood estimated under CR \cr
    \code{$conv_CR} - convergence code for likelihood optimization; conv = 0 means convergence \cr
    \code{$b_PBD1} - speciation-initation rate estimated under PBD for first set of initial values\cr
    \code{$mu_PB1} - extinction rate estimated under DD for first set of initial values \cr
    \code{$lambda_PB1} - speciation-completion rate estimated under PBD for first set of initial values \cr
    \code{$LL_PBD1} - maximum likelihood estimated under DD for first set of initial values \cr
    \code{$conv_PBD1} - convergence code for likelihood optimization for first set of initial values; conv = 0 means convergence \cr
    \code{$b_PBD2} - speciation-initation rate estimated under PBD for second set of initial values\cr
    \code{$mu_PB2} - extinction rate estimated under DD for second set of initial values \cr
    \code{$lambda_PB2} - speciation-completion rate estimated under PBD for second set of initial values \cr
    \code{$LL_PBD2} - maximum likelihood estimated under DD for second set of initial values \cr
    \code{$conv_PBD2} - convergence code for likelihood optimization for second set of initial values; conv = 0 means convergence \cr
    \code{$LR} - likelihood ratio between DD and CR
 }
 \item{pvalue}{p-value of the test}
 \item{LRalpha}{Likelihood ratio at the signifiance level alpha}
 \item{poweroftest}{power of the test for significance level alpha}
}
\references{
- Etienne, R.S. et al. 2016. Meth. Ecol. Evol. 7: 1092-1099, doi: 10.1111/2041-210X.12565 \cr
}
\author{ Rampal S. Etienne }

\seealso{
 \code{\link{pbd_loglik}},
 \code{\link{pbd_ML}}
}
\keyword{ models }
