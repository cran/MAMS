\name{new.bounds}
\alias{new.bounds}

\title{
Function to update boundaries based on observed sample sizes
}

\description{
The function determines updated boundaries of a multi-arm multi-stage study based on observed number of observations per arm.
}

\usage{
new.bounds(K=3, J=2, alpha=0.05, nMat=matrix(c(10, 20), nrow=2, ncol=4), 
           u=NULL, l=NULL, ushape="obf", lshape="fixed", ufix=NULL,
           lfix=0, N=20, parallel=TRUE, print=TRUE)
}

\arguments{
  \item{K}{Number of experimental treatments (default=\code{3}).}
  \item{J}{Number of stages (default=\code{2}).}
  \item{alpha}{One-sided familywise error rate (default=\code{0.05}).}
  \item{nMat}{Jx(K+1) dimensional matrix of observed/expected sample sizes. Rows correspond to stages and columns to arms. First column is control (default: 2x4 matrix with 10 subjects per stage and arm).}
  \item{u}{Vector of previously used upper boundaries (default=\code{NULL}).}
  \item{l}{Vector of previously used upper boundaries (default=\code{NULL}).}
  \item{ushape}{Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{lshape}{Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"}, \code{"triangular"} and \code{"fixed"} (the default). See details.}
  \item{ufix}{Fixed upper boundary (default=\code{NULL}). Only used if \code{shape="fixed"}.}
  \item{lfix}{Fixed lower boundary (default=\code{0}). Only used if \code{shape="fixed"}.}
  \item{N}{Number of quadrature points per dimension in the outer integral (default=20).}
  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the computation via a user-defined strategy specified by means of the function \code{\link[future:plan]{future::plan()}}. If not set differently, the default strategy is \code{sequential}, which corresponds to a computation without parallelisation.}
  \item{print}{if \code{TRUE} (default), indicate at which stage the computation is.}
}

\details{
This function finds the boundaries for a given matrix of sample sizes in multi-arm multi-stage study with K active treatments plus control. The vectors \code{u} and \code{l} are the boundaries used so far while \code{u.shape} and \code{l.shape} specify the shape to the boundaries for the remaining analysis. By specifying \code{u} and \code{l} as \code{NULL}, a design using only the shapes given by \code{ushape} and \code{lshape} can be found for any sample sizes per stage and arm.\cr

The shape of the boundaries (\code{ushape}, \code{lshape}) are either using the predefined shapes following Pocock (1977), O'Brien & Fleming (1979) or the triangular Test (Whitehead, 1997) using options \code{"pocock"}, \code{"obf"} or \code{"triangular"} respectively, are constant (option \code{"fixed"}) or supplied in as a function. If a function is passed it should require exactly one argument specifying the number of stages and return a vector of the same length. The lower boundary shape is required to be non-decreasing while the upper boundary shape needs to be non-increasing. If a fixed lower boundary is used, \code{lfix} must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2} to ensure that it is smaller than the upper boundary.\cr
}

\value{
An object of the class MAMS containing the following components: \cr 
  \item{l}{Lower boundary.} 
  \item{u}{Upper boundary.}
  \item{n}{Sample size on control in stage 1.}
  \item{N}{Maximum total sample size.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{alpha}{Familywise error rate.}
  \item{power}{Power under least favorable configuration.}
  \item{rMat}{Matrix of allocation ratios. First row corresponds to control and second row to experimental treatments.}

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

\author{Thomas Jaki, Dominic Magirr and Dominique-Laurent Couturier}

\examples{\donttest{
# Note that some of these examples may take a few minutes to run

# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary with 
# equal sample size per arm and stage. Results are equivalent to using 
# mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, ushape="obf", lshape="fixed", 
#      lfix=0, sample.size=FALSE)

new.bounds(K=4, J=2, alpha=0.05, nMat=matrix(c(10, 20), nrow=2, ncol=5), u=NULL, l=NULL, 
           ushape="obf", lshape="fixed", lfix=0)

# A 2-stage design that was designed to use an O'Brien & Fleming efficacy and 
# zero futility boundary with equal sample size per arm and stage (n=14).
# The observed sample size after stage one are 10, 10, 18, 10, 13 for each arm 
# while the original upper bounds used are (3.068, 2.169) for stage 1. 
# The updated bounds are (3.068, 2.167).

new.bounds(K=4, J=2, alpha=0.05, 
           nMat=matrix(c(10, 28, 10, 28, 18, 28, 10, 28, 13, 28), nrow=2, ncol=5), 
           u=3.068, l=0, ushape="obf", lshape="fixed", lfix=0)
           
# same using parallelisation via separate R sessions running in the background
future::plan(multisession)
new.bounds(K=4, J=2, alpha=0.05, 
           nMat=matrix(c(10, 28, 10, 28, 18, 28, 10, 28, 13, 28), nrow=2, ncol=5), 
           u=3.068, l=0, ushape="obf", lshape="fixed", lfix=0)
future::plan("default")
           
}}

\keyword{design}
