\name{ordinal.mams}
\alias{ordinal.mams}

\title{
Function to design multi-arm multi-stage studies with ordinal or binary endpoints
}

\description{
The function determines (approximately) the boundaries of a multi-arm multi-stage study with ordinal or binary endpoints for a given boundary shape and finds the required number of subjects.
}

\usage{
ordinal.mams(prob=c(0.35, 0.4, 0.25), or=2, or0=1.2, K=4, J=2, alpha=0.05,
     power=0.9, r=1:2, r0=1:2, ushape="obf", lshape="fixed", ufix=NULL,
     lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, N=20, 
     parallel=TRUE, print=TRUE)
}

\arguments{
  \item{prob}{Vector of expected probabilities of falling into each category under control conditions. The elements must sum up to one (default=\code{c(0.35, 0.4, 0.25)}).}
  \item{or}{Interesting treatment effect on the scale of odds ratios (default=\code{2}).}
  \item{or0}{Uninteresting treatment effect on the scale of odds ratios (default=\code{1.2}).}
  \item{K}{Number of experimental treatments (default=\code{4}).}
  \item{J}{Number of stages (default=\code{2}).}
  \item{alpha}{One-sided familywise error rate (default=\code{0.05}).}
  \item{power}{Desired power (default=\code{0.9}).}
  \item{r}{Vector of allocation ratios (default=\code{1:2}).}
  \item{r0}{Vector ratio on control (default=\code{1:2}).}
  \item{ushape}{Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}.}
  \item{lshape}{Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"}, \code{"triangular"} and \code{"fixed"} (the default).}
  \item{ufix}{Fixed upper boundary (default=\code{NULL}). Only used if \code{shape="fixed"}.}
  \item{lfix}{Fixed lower boundary (default=\code{0}). Only used if \code{shape="fixed"}.}
  \item{nstart}{Starting point for finding the sample size (default=\code{1}).}
  \item{nstop}{Stopping point for finding the sample size (default=\code{NULL}).}
  \item{sample.size}{Logical if sample size should be found as well (default=\code{TRUE}).}
  \item{N}{Number of quadrature points per dimension in the outer integral (default=\code{20}).}
  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the computation via a user-defined strategy specified by means of the function \code{\link[future:plan]{future::plan()}}. If not set differently, the default strategy is \code{sequential}, which corresponds to a computation without parallelisation.}
  \item{print}{if \code{TRUE} (default), indicate at which stage the computation is.}
}

\details{
This function finds the (approximate) boundaries and sample size of a multi-arm multi-stage study with ordinal or binary endpoints with K active treatments plus control in which all promising treatments are continued at interim analyses as described in Magirr et al (2012). It is a wrapper around the basic \code{\link[MAMS:mams]{mams}} function to facilitate its use with ordinal and binary endpoints, following ideas of Whitehead & Jaki (2009) and Jaki & Magirr (2013). For a binary endpoint the vector \code{prob} has only two elements (success/failure, yes/no, etc.). See \code{\link{mams}} for further details on the basic methodology.\cr
}

\value{
An object of the class MAMS containing the following components:\cr 
  \item{l}{Lower boundary.}
  \item{u}{Upper boundary.}
  \item{n}{Sample size on control in stage 1.}
  \item{N}{Maximum total sample size.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{alpha}{Familywise error rate.}
  \item{alpha.star}{Cumulative familywise error rate spent by each analysis.}
  \item{power}{Power under least favorable configuration.}
  \item{rMat}{Matrix of allocation ratios. First row corresponds to control while subsequent rows are for the experimental treatments.}

}

\seealso{
  \code{\link{print.MAMS}}, \code{\link{summary.MAMS}}, \code{\link{plot.MAMS}}, \code{\link{mams}}, \code{\link{MAMS}}.
}

\references{
Jaki T., Pallmann P. and Magirr D. (2019), \emph{The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of Statistical Software}, 88(4), 1-25. Link: \doi{10.18637/jss.v088.i04}

Magirr D., Jaki T. and Whitehead J. (2012), \emph{A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection}, \bold{Biometrika}, 99(2), 494-501. Link: \doi{10.1093/biomet/ass002}

Magirr D., Stallard N. and Jaki T. (2014), \emph{Flexible sequential designs for multi-arm clinical trials}, \bold{Statistics in Medicine}, 33(19), 3269-3279. Link: \doi{10.1002/sim.6183}

Pocock S.J. (1977), \emph{Group sequential methods in the design and analysis of clinical trials}, \bold{Biometrika}, 64(2), 191-199.\cr

O'Brien P.C., Fleming T.R. (1979), \emph{A multiple testing procedure for clinical trials}, \bold{Biometrics}, 35(3), 549-556.\cr 

Whitehead J. (1997), \emph{The Design and Analysis of Sequential Clinical Trials}, \bold{Wiley}: Chichester, UK.\cr 
}

\author{Philip Pallmann}

\examples{\donttest{
## An example based on the example in Whitehead & Jaki (2009)
# 2-stage design with triangular efficacy and futility boundaries
prob <- c(0.075, 0.182, 0.319, 0.243, 0.015, 0.166)
ordinal.mams(prob=prob, or=3.06, or0=1.32, K=3, J=2, alpha=0.05,
                 power=0.9, r=1:2, r0=1:2, ushape="triangular",
                 lshape="triangular")
# same example with parallelisation via separate R sessions running in the background
future::plan(multisession)
ordinal.mams(prob=prob, or=3.06, or0=1.32, K=3, J=2, alpha=0.05,
                 power=0.9, r=1:2, r0=1:2, ushape="triangular",
                 lshape="triangular", parallel=TRUE)
future::plan("default")
}}

\keyword{design}
