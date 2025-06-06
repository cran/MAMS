% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{summary.MAMS}
\alias{summary.MAMS}
\title{Generic summary function for class MAMS.}
\usage{
\method{summary}{MAMS}(
  object,
  digits = max(3, getOption("digits") - 4),
  extended = FALSE,
  ...
)
}
\arguments{
\item{object}{An output object of class MAMS}

\item{digits}{Number of significant digits to be printed.}

\item{extended}{TRUE or FALSE}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Text output.
}
\description{
Produces a detailed summary of an object from class MAMS
}
\examples{
\donttest{
# 2-stage design with triangular boundaries
res <- mams(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2,
             p=0.65, p0=0.55,
             ushape="triangular", lshape="triangular", nstart=30)

summary(res)
}
}
\author{
Dominique-Laurent Couturier
}
