\name{PhaseType}
\alias{dphasetype}
\alias{pphasetype}
\alias{qphasetype}
\alias{rphasetype}
\alias{mphasetype}
\alias{mgfphasetype}
\title{The Phase-type Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and moment generating function for the Phase-type 
  distribution with parameters \code{pi}, \code{T} and
  \code{m}.
}
\usage{
  dphasetype(x, pi, T, m, log = FALSE)
  pphasetype(q, pi, T, m, lower.tail = TRUE, log.p = FALSE)
  qphasetype(p, pi, T, m, lower.tail = TRUE, log.p = FALSE)
  rphasetype(n, pi, T, m)
  mphasetype(order, pi, T, m)
  mgfphasetype(x, pi, T, m, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{pi, T, m}{parameters of the Phase-type distribution. Must be a vector, a square matrix and an integer respectively}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
}
\details{
  The Phase-type distribution with parameters \code{pi}, \code{T} and \code{m}
  has density:
\deqn{f(x) = pi * exp(Tx) * t_0}
  and a distribution function 
\deqn{F(x) = 1-pi * exp(Tx) * 1_m}
where \code{pi} is the initial probability vector (of length \code{m}),
\code{T} the \eqn{m\times m} subintensity matrix and \eqn{1_m} 
the one vector of \eqn{R^m}, \eqn{t_0} the exit rates vector
defined as \eqn{-T * 1_m}. Warning in the previous expression, 
\eqn{exp(Tx)} denotes the matrix exponential defined as the Taylor series
\eqn{\sum_{i=0}^{+\infty} (Tx)^i/i!}.
More details can be found in Bladt(2005).
}
\value{
  \code{dphasetype} gives the density,
  \code{pphasetype} gives the distribution function,
  \code{qphasetype} gives the quantile function,
  \code{rphasetype} generates random deviates,
  \code{mphasetype} gives the \eqn{k}th raw moment, and
  \code{mgfphasetype} gives the moment generating function.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}

\references{
  Bladt, M. (2005),
  \emph{A review on phase-type distribution and their use
in risk theory}, ASTIN Bul.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Christophe Dutang
}
\examples{
#erlang distribution Er(3,2)
 T <- cbind(c(-2,0,0),c(2,-2,0),c(0,2,-2))
 m <- 3
 pi <- c(1,0,0)
 x<- 0:10
#density
 dphasetype(x,pi,T,m)
#comparison with dgamma
 dgamma(x,3,2)
#cumulative distribution
 pphasetype(x,pi,T,m)
#comparison with pgamma
 pgamma(x,3,2)
#random generator
 rphasetype(10, pi, T, m) 
#comparison with rgamma
 rgamma(10,3,2) 
#moment function
 mphasetype(0:10, pi, T, m)
#comparison with mgamma
 mgamma(0:10, 3, 2)
#moment generating function
 x <- -10:2
 mgfphasetype(x, pi, T, m)
#comparison with mgamma
 mgfgamma(x, 3, 2)
}
\keyword{distribution}
