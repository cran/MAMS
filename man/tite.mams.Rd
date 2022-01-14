\name{tite.mams}
\alias{tite.mams}

\title{
Function to design multi-arm multi-stage studies with time-to-event endpoints
}

\description{
The function determines (approximately) the boundaries of a multi-arm multi-stage study with time-to-event endpoints for a given boundary shape and finds the required number of events.
}

\usage{
tite.mams(hr=1.5, hr0=1.1, K=4, J=2, alpha=0.05, power=0.9,
           r=1:2, r0=1:2, ushape="obf", lshape="fixed", ufix=NULL,
           lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, N=20,
           parallel=TRUE, print=TRUE)
}

\arguments{
  \item{hr}{Interesting treatment effect on the scale of hazard ratios (default=\code{2}).}
  \item{hr0}{Uninteresting treatment effect on the scale of hazard ratios (default=\code{1.2}).}
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
  \item{N}{Number of quadrature points per dimension in the outer integral (default=20).}
  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the computation via a user-defined strategy specified by means of the function \code{\link[future:plan]{future::plan()}}. If not set differently, the default strategy is \code{sequential}, which corresponds to a computation without parallelisation.}
  \item{print}{if \code{TRUE} (default), indicate at which stage the computation is.}
  
}

\details{
This function finds the (approximate) boundaries and sample size of a multi-arm multi-stage study with time-to-event endpoints with K active treatments plus control in which all promising treatments are continued at interim analyses as described in Magirr et al (2012). It is a wrapper around the basic \code{\link[MAMS:mams]{mams}} function to facilitate its use with time-to-event endpoints, following ideas of Jaki & Magirr (2013). Note that the sample size is calculated as the required number of events, from which the total sample size can be estimated (e.g., Whitehead 2001). See \code{?mams} for further details on the basic methodology.\cr
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
Jaki T. and Magirr D. (2013), \emph{Considerations on covariates and endpoints in multi-arm multi-stage clinical trials selecting all promising treatments}, 
\bold{Statistics in Medicine}, 32(7), 1150-1163. Link: \doi{10.1002/sim.5669}

Jaki T., Pallmann P. and Magirr D. (2019), \emph{The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of Statistical Software}, 88(4), 1-25. Link: \doi{10.18637/jss.v088.i04}

Magirr D., Jaki T. and Whitehead J. (2012), \emph{A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection}, \bold{Biometrika}, 99(2), 494-501. Link: \doi{10.1093/biomet/ass002}

Whitehead J. (2001), \emph{Predicting the duration of sequential survival studies}, \bold{Drug Information Journal}, 35(4), 1387-1400.
}

\author{Philip Pallmann, Dominic Magirr}

\examples{\donttest{
## An example 2-stage design with triangular efficacy and futility boundaries
tite.mams(hr=2, hr0=1.5, K=3, J=2, alpha=0.05, power=0.9,
          r=1:2, r0=1:2, ushape="triangular", lshape="triangular")
}}

\keyword{design}
