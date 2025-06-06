% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simultaneous.R
\name{tite.mams}
\alias{tite.mams}
\title{Function to design multi-arm multi-stage studies with time-to-event
endpoints}
\usage{
tite.mams(
  hr = 1.5,
  hr0 = 1.1,
  K = 4,
  J = 2,
  alpha = 0.05,
  power = 0.9,
  r = 1:2,
  r0 = 1:2,
  ushape = "obf",
  lshape = "fixed",
  ufix = NULL,
  lfix = 0,
  nstart = 1,
  nstop = NULL,
  sample.size = TRUE,
  Q = 20,
  parallel = TRUE,
  print = TRUE
)
}
\arguments{
\item{hr}{Interesting treatment effect on the scale of hazard ratios
(default=2).}

\item{hr0}{Uninteresting treatment effect on the scale of hazard ratios
(default=1.2).}

\item{K}{Number of experimental treatments (default=4).}

\item{J}{Number of stages (default=2).}

\item{alpha}{One-sided familywise error rate (default=0.05).}

\item{power}{Desired power (default=0.9).}

\item{r}{Vector of allocation ratios (default=1:2).}

\item{r0}{Vector ratio on control (default=1:2).}

\item{ushape}{Shape of upper boundary. Either a function specifying the shape
or one of "pocock", "obf" (the default), "triangular" and "fixed".}

\item{lshape}{Shape of lower boundary. Either a function specifying the shape
or one of "pocock", "obf", "triangular" and "fixed" (the default).}

\item{ufix}{Fixed upper boundary (default=NULL). Only used if shape="fixed".}

\item{lfix}{Fixed lower boundary (default=0). Only used if shape="fixed".}

\item{nstart}{Starting point for finding the sample size (default=1).}

\item{nstop}{Stopping point for finding the sample size (default=NULL).}

\item{sample.size}{Logical if sample size should be found as well
(default=TRUE).}

\item{Q}{Number of quadrature points per dimension in the outer integral
(default=20).}

\item{parallel}{if TRUE (default), allows parallelisation of the computation
via a user-defined strategy specified by means of the function
future::plan(). If not set differently, the default strategy is sequential,
which corresponds to a computation without parallelisation.}

\item{print}{if TRUE (default), indicate at which stage the computation is.}
}
\value{
An object of the class MAMS containing the following components:
\item{l}{Lower boundary.}
\item{u}{Upper boundary.}
\item{n}{Sample size on control in stage 1.}
\item{N}{Maximum total sample size.}
\item{K}{Number of experimental treatments.}
\item{J}{Number of stages in the trial.}
\item{alpha}{Familywise error rate.}
\item{alpha.star}{Cumulative familywise error rate spent by each analysis.}
\item{power}{Power under least favorable configuration.}
\item{rMat}{Matrix of allocation ratios. First row corresponds to control
while subsequent rows are for the experimental treatments.}
}
\description{
The function determines (approximately) the boundaries of a
multi-arm multi-stage study with time-to-event endpoints for a given boundary
shape and finds the required number of events.
}
\details{
This function finds the (approximate) boundaries and sample size of
a multi-arm multi-stage study with time-to-event endpoints with K active
treatments plus control in which all promising treatments are continued at
interim analyses as described in Magirr et al (2012). It is a wrapper around
the basic mams function to facilitate its use with time-to-event endpoints,
following ideas of Jaki & Magirr (2013). Note that the sample size is
calculated as the required number of events, from which the total sample
size can be estimated (e.g., Whitehead 2001). See ?mams for further details
on the basic methodology.
}
\examples{
\donttest{
## An example 2-stage design with triangular efficacy and futility boundaries
tite.mams(hr=2, hr0=1.5, K=3, J=2, alpha=0.05, power=0.9,
          r=1:2, r0=1:2, ushape="triangular", lshape="triangular")
}
}
\author{
Philip Pallmann, Dominic Magirr
}
