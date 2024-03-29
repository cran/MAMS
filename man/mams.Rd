\name{mams}
\alias{mams}

\title{
Function to design multi-arm multi-stage studies with normal endpoints
}

\description{
The function determines the boundaries of a multi-arm multi-stage study for a given boundary shape and finds the required number of subjects.
}

\usage{
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.75, p0=0.5, 
     delta=NULL, delta0=NULL, sd=NULL, ushape="obf", lshape="fixed",
     ufix=NULL, lfix=0, nstart=1, nstop=NULL, sample.size=TRUE, N=20,
     type="normal", parallel=TRUE, print=TRUE)
}

\arguments{
  \item{K}{Number of experimental treatments (default=\code{4}).}
  \item{J}{Number of stages (default=\code{2}).}
  \item{alpha}{One-sided familywise error rate (default=\code{0.05}).}
  \item{power}{Desired power (default=\code{0.9}).}
  \item{r}{Vector of allocation ratios (default=\code{1:2}).}
  \item{r0}{Vector ratio on control (default=\code{1:2}).}
  \item{p}{Interesting treatment effect on the probability scale. See Details (default=\code{0.75}).}
  \item{p0}{Uninteresting treatment effect on the probability scale. See Details (default=\code{0.5}).}
  \item{delta}{Interesting treatment effect on the traditional scale. See Details (default=\code{NULL}).}
  \item{delta0}{Uninteresting treatment effect on the traditional scale. See Details (default=\code{NULL}).}
  \item{sd}{Standard deviation, assumed to be known. See Details (default=\code{NULL}).}
  \item{ushape}{Shape of upper boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"} (the default), \code{"triangular"} and \code{"fixed"}. See details.}
  \item{lshape}{Shape of lower boundary. Either a function specifying the shape or one of \code{"pocock"}, \code{"obf"}, \code{"triangular"} and \code{"fixed"} (the default). See details.}
  \item{ufix}{Fixed upper boundary (default=\code{NULL}). Only used if \code{shape="fixed"}.}
  \item{lfix}{Fixed lower boundary (default=\code{0}). Only used if \code{shape="fixed"}.}
  \item{nstart}{Starting point for finding the sample size (default=\code{1}).}
  \item{nstop}{Stopping point for finding the sample size (default=\code{NULL}).}
  \item{sample.size}{Logical if sample size should be found as well (default=\code{TRUE}).}
  \item{N}{Number of quadrature points per dimension in the outer integral (default=\code{20}).}
  \item{type}{Will be changed automatically by the wrappers \code{\link[MAMS:tite.mams]{tite.mams()}} (to \code{"tite"}) and \code{\link[MAMS:ordinal.mams]{ordinal.mams()}} (to \code{"ordinal"}) to customise the output.}
  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the computation via a user-defined strategy specified by means of the function \code{\link[future:plan]{future::plan()}}. If not set differently, the default strategy is \code{sequential}, which corresponds to a computation without parallelisation.}
  \item{print}{if \code{TRUE} (default), indicate at which stage the computation is.}
  
}

\details{
This function finds the boundaries and sample size of a multi-arm multi-stage study with K active treatments plus control in which all promising treatments are continued at interim analyses as described in Magirr et al (2012). At each interim analysis the test statistics are compared to the lower (futility) bound and any treatment whose corresponding test statistic falls below that bound is discontinued. Similarly if any test statistic exceeds the upper (efficacy) bound the null hypothesis corresponding to that treatment can be rejected and superiority of that treatment over control claimed. At the same time the study is stopped. If at least one test statistic exceeds the lower bound and none exceeds the upper bound the study is continued and further patients are recruited to all remaining experimental treatments plus control.\cr

The design is found under the least favorable configuration, which requires an interesting treatment effect \code{p} that if present we would like to find with high probability and an uninteresting effect \code{p0}. Both \code{p} and \code{p0} are parameterized as \eqn{P(X_k > X_0 ) = p}{P(X_k > X_0 ) = p}, that is the probability of a randomly selected person on treatment k observing a better outcome than a random person on control. For \code{p=0.5} the experimental treatment and control perform equally well. The advantage of this paramterization is that no knowledge about the variance is required. To convert traditional effect sizes, \eqn{\delta}{delta} to this format use \eqn{p=\Phi(\frac{\delta}{\sqrt{2}\sigma})}{Phi(delta/(2^0.5*sigma))}. Alternatively, the interesting and uninteresting effect size can also be specified directly on the traditional scale of \code{delta} and \code{delta} with an additional specification of the standard deviation \code{sd} assumed to be known.\cr

The shape of the boundaries (\code{ushape}, \code{lshape}) are either using the predefined shapes following Pocock (1977), O'Brien & Fleming (1979) or the triangular Test (Whitehead, 1997) using options \code{"pocock"}, \code{"obf"} or \code{"triangular"} respectively, are constant (option \code{"fixed"}) or supplied in as a function. If a function is passed it should require exactly one argument specifying the number of stages and return a vector of the same length. The lower boundary shape is required to be non-decreasing while the upper boundary shape needs to be non-increasing. If a fixed lower boundary is used, \code{lfix} must be smaller than \eqn{\Phi^{-1}(1-\alpha)/2}{Phi(1-alpha)/2} to ensure that it is smaller than the upper boundary.\cr

The default starting point for finding the sample size is \code{nstart=1}, and the default point where the search is stopped (when \code{nstop=NULL}) is 3 times the sample size of the corresponding fixed single-stage design.\cr

Computation of designs with more than four stages are very time consuming and not advised. The parameter \code{sample.size} controls whether the required sample size is computed as well. Setting to \code{FALSE} approximately halves the compuation time. 

For designs with more than 2 stages, parallelisation of the computation by means of the packages \code{future} and \code{future.apply} lead to decreased computation times when choosing a parallelisation strategy like, for example, \code{multicore} (using separate forked \pkg{R} processes, available to unix/osx users) or \code{multisession} (using separate \pkg{R} sessions, available to all users) (refer to \code{\link[future:plan]{future::plan()}} for detail). \cr 
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
  \code{\link{print.MAMS}}, \code{\link{summary.MAMS}}, \code{\link{plot.MAMS}}, \code{\link{new.bounds}}, \code{\link{ordinal.mams}}, \code{\link{tite.mams}}, \code{\link{MAMS}}.
}

\references{
Jaki T., Pallmann P. and Magirr D. (2019), \emph{The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of Statistical Software}, 88(4), 1-25. Link: \doi{10.18637/jss.v088.i04}

Magirr D., Jaki T. and Whitehead J. (2012), \emph{A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection}, \bold{Biometrika}, 99(2), 494-501. Link: \doi{10.1093/biomet/ass002}

Pocock S.J. (1977), \emph{Group sequential methods in the design and analysis of clinical trials}, \bold{Biometrika}, 64(2), 191-199.\cr

O'Brien P.C., Fleming T.R. (1979), \emph{A multiple testing procedure for clinical trials}, \bold{Biometrics}, 35(3), 549-556.\cr 

Whitehead J. (1997), \emph{The Design and Analysis of Sequential Clinical Trials}, \bold{Wiley}: Chichester, UK.\cr 
}

\author{Thomas Jaki, Dominic Magirr and Dominique-Laurent Couturier}

\examples{\donttest{
## A fixed sample (single stage) design specified on the p scale
mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=0.65, p0=0.55)

## The same design specified on the delta scale
mams(K=4, J=1, alpha=0.05, power=0.9, r=1, r0=1, p=NULL, p0=NULL,
     delta=0.545, delta0=0.178, sd=1)

## An example in Table 1 of Magirr et al (2012)
# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.65, p0=0.55, 
     ushape="obf", lshape="fixed", lfix=0, nstart=40)

# Note that these examples may take a few minutes to run

## 3-stage design with Triangular efficacy and futility boundary
mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65, p0=0.55, 
     ushape="triangular", lshape="triangular", nstart=30)

## Example of a custom boundary function without sample size evaluation
mams(K=6, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.7, p0=0.5, 
     ushape=function(x)return(x:1), lshape="fixed", lfix=0, 
     sample.size=FALSE)

## Different allocation ratios between control and experimental treatments. 
## Twice as many patients are randomized to control at each stage.
mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=c(2, 4), p=0.65, p0=0.55, 
     ushape="obf", lshape="fixed", lfix=0, nstart=30)

## 
## example considering different parallelisation strategies 
## 


# parallel = FALSE (future framework not used)
set.seed(1)
system.time(
print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65, p0=0.55, 
     ushape="triangular", lshape="triangular", nstart=30, parallel = FALSE))
)
# parallel = TRUE (default) with default strategy (sequential computation)
plan(sequential)
set.seed(1)
system.time(
print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65, p0=0.55, 
     ushape="triangular", lshape="triangular", nstart=30))
)
# parallel = TRUE(default) with multisession strategy (parallel computation)
plan(multisession)
set.seed(1)
system.time(
print(mams(K=4, J=3, alpha=0.05, power=0.9, r=1:3, r0=1:3, p=0.65, p0=0.55, 
     ushape="triangular", lshape="triangular", nstart=30))
)
plan("default")       
     
}}

\keyword{ design }
