\name{actuar-deprecated}
\alias{actuar-deprecated}
\alias{minvGauss-deprecated}
\alias{levinvGauss-deprecated}
\alias{mgfinvGauss-deprecated}
\title{Deprecated Functions in Package \pkg{actuar}}
\description{
  These functions are provided for compatibility with older versions
  of \pkg{actuar} only, and may be defunct as soon as the next release.
}
\usage{
# Deprecated in \pkg{actuar} 2.0-0
minvGauss(order, nu, lambda)
levinvGauss(limit, nu, lambda, order = 1)
mgfinvGauss(x, nu, lambda, log= FALSE)
}
\details{
  Functions \code{[m,lev,mgf]invGauss} complemented functions
  \code{[d,p,q,r]invGauss} of package \pkg{SuppDists}. From version
  2.0-0, \pkg{actuar} has its own, complete, set of support functions
  for the inverse gaussian distribution; see \code{\link{dinvgauss}}.
}
\keyword{misc}
