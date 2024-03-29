\name{mams.sim}
\alias{mams.sim}

\title{
Simulating multi-arm multi-stage designs
}

\description{
The function simulates multi-arm multi-stage designs and estimates power and expected sample size.
}

\usage{
mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5),
         u=c(3.068, 2.169), l=c(0.000, 2.169),
         pv=rep(0.5, 4), deltav=NULL, sd=NULL, ptest=1, parallel=TRUE)
}

\arguments{
  \item{nsim}{Number of simulations (default=\code{10000}).}
  \item{nMat}{Jx(K+1) dimensional matrix of observed/expected sample sizes. Rows correspond to stages and columns to arms. First column is control (default: 2x5 matrix with 44 subjects per stage and arm).}
  \item{u}{Vector of previously used upper boundaries (default=\code{NULL}).}
  \item{l}{Vector of previously used upper boundaries (default=\code{NULL}).}
  \item{pv}{Vector of size K of true treatment effects on the probability scale. See Details (default=\code{rep(0.5, 4)}).}
  \item{deltav}{Vector of size K of true treatment effects on the traditional scale. See Details (default=\code{NULL}).}
  \item{sd}{Standard deviation. See Details (default=\code{NULL}).}
  \item{ptest}{Vector of treatment numbers for determining power. For example, c(1, 2) will count rejections of one or both hypotheses for testing treatments 1 and 2 against control.}
  \item{parallel}{if \code{TRUE} (default), allows parallelisation of the computation via a user-defined strategy specified by means of the function \code{\link[future:plan]{future::plan()}}. If not set differently, the default strategy is \code{sequential}, which corresponds to a computation without parallelisation.}
}

\details{
This function simulates multi-arm multi-stage studies for a given matrix of sample sizes and boundaries given by the vectors \code{u} and \code{l}. The effect difference between each experimental treatment and control is given by \code{pv} and is parameterized as \eqn{P(X_k > X_0 ) = p}{P(X_k > X_0 ) = p}. That is the probability of a randomly selected person on treatment k observing a better outcome than a random person on control. For \code{pv=rep(0.5,4} the experimental treatments and control perform equally well (i.e. the global null hypothesis is true). The advantage of this paramterization is that no knowledge about the variance is required. To convert traditional effect sizes, \eqn{\delta}{delta} to this format use \eqn{p=\Phi(\frac{\delta}{\sqrt{2}\sigma})}{Phi(delta/(2^0.5*sigma))}. Alternatively, the effect size can also be specified directly on the traditional scale of \code{deltav} with an additional specification of the standard deviation \code{sd}.\cr

The function returns the probability of rejecting any hypothesis (\code{typeI}), the power to reject the first hypothesis when the first treatment has the largest estimated effect, the proportion of rejections of the hypothesis specified by \code{ptest} (\code{prop.rej}) as well as the expected sample size. \cr
}

\value{
An object of the class MAMS.sim containing the following components: \cr 

  \preformatted{res$typeI <- mean(unlist(reps["rej",]))}
  \preformatted{res$power <- mean(unlist(reps["pow",]))}
  \preformatted{res$prop.rej <- rej/nsim}
  \preformatted{res$exss <- mean(unlist(reps["ess",]))}

  \item{l}{Lower boundary.} 
  \item{u}{Upper boundary.}
  \item{n}{Sample size on control in stage 1.}
  \item{N}{Maximum total sample size.}
  \item{K}{Number of experimental treatments.}
  \item{J}{Number of stages in the trial.}
  \item{rMat}{Matrix of allocation ratios. First row corresponds to control and second row to experimental treatments.}
  \item{nsim}{Number of simulation runs.}
  \item{typeI}{The proportion any hypothesis is rejected.}
  \item{power}{The proportion the first hypothesis is rejected and the corresponding test statistic is largest.}
  \item{ptest}{The vector \code{ptest}.}
  \item{prop.rej}{The proportion of times at least one of the hypothesis specified by \code{ptest} is rejected.}
  \item{exss}{The expected sample size.}

}

\seealso{
  \code{\link{print.MAMS.sim}}, \code{\link{summary.MAMS.sim}}, \code{\link{mams}}, \code{\link{MAMS}}.
}

\references{
Jaki T., Pallmann P. and Magirr D. (2019), \emph{The R Package MAMS for Designing Multi-Arm Multi-Stage Clinical Trials}, \bold{Journal of Statistical Software}, 88(4), 1-25. Link: \doi{10.18637/jss.v088.i04}

Magirr D., Jaki T. and Whitehead J. (2012), \emph{A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection}, \bold{Biometrika}, 99(2), 494-501. Link: \doi{10.1093/biomet/ass002}
}

\author{Thomas Jaki, Dominic Magirr and Dominique-Laurent Couturier}

\examples{\donttest{
# Note that some of these examples may take a few minutes to run

# 2-stage design with O'Brien & Fleming efficacy and zero futility boundary with 
# equal sample size per arm and stage. Design can be found using
# mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, ushape="obf", lshape="fixed", 
#      lfix=0, p=0.65, p0=0.55)

# under global null hypothesis (using the pv scale)
mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=rep(0.5, 4), ptest=1)

# under global null hypothesis (using the deltav scale)
mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=NULL, deltav=rep(0, 4), sd=1, ptest=1)

# under LFC
mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2)

# when all treatments doing similarly well
mams.sim(nsim=10000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=c(0.63, 0.62, 0.60, 0.61), ptest=4)

## 
## example considering different parallelisation strategies 
## 

# parallel = FALSE (future framework not used)
set.seed(1)
system.time(
print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2, parallel=FALSE))
)
# parallel = TRUE (default) with default strategy (sequential computation)
plan(sequential)
set.seed(1)
system.time(
print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2))
)
# parallel = TRUE (default) with multisession strategy (parallel computation)
plan(multisession)
set.seed(1)
system.time(
print(mams.sim(nsim=25000, nMat=matrix(c(44, 88), nrow=2, ncol=5), u=c(3.068, 2.169),
         l=c(0.000, 2.169), pv=c(0.65, 0.55, 0.55, 0.55), ptest=1:2))
)
plan("default")       

}}


\keyword{design}
