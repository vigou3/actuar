\name{ZeroModifiedPoisson}
\alias{ZeroModifiedPoisson}
\alias{ZMpoisson}
\alias{dzmpois}
\alias{pzmpois}
\alias{qzmpois}
\alias{rzmpois}
\title{The Zero-Modified Poisson Distribution}
\description{
  Density function, distribution function, quantile function, random
  generation for the Zero-Modified Poisson distribution with parameter
  \code{lambda} and arbitrary probability at zero \code{p0}.
}
\usage{
dzmpois(x, lambda, p0, log = FALSE)
pzmpois(q, lambda, p0, lower.tail = TRUE, log.p = FALSE)
qzmpois(p, lambda, p0, lower.tail = TRUE, log.p = FALSE)
rzmpois(n, lambda, p0)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of values to return.}
  \item{lambda}{vector of (non negative) means.}
  \item{p0}{probability mass at zero. \code{0 <= p0 <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-modified Poisson distribution is a discrete mixture between a
  degenerate distribution at zero and a (standard) Poisson. The
  probability mass function is \eqn{p(0) = p_0}{p(0) = p0} and
  \deqn{%
    p(x) = \frac{(1-p_0)}{(1-e^{-\lambda})} f(x)}{%
    p(x) = (1-p0)/(1-exp(-lambda)) f(x)}
  for \eqn{x = 1, 2, ...}, \eqn{\lambda > 0} and \eqn{0 \le
  p_0 \le 1}{0 \le p0 \le 1}, where \eqn{f(x)} is the probability mass
  function of the Poisson.
  The cumulative distribution function is
  \deqn{P(x) = p_0 + (1 - p_0) \left(\frac{F(x) - F(0)}{1 - F(0)}\right).}{%
    P(x) = p0 + (1 - p0) [F(x) - F(0)]/[1 - F(0)].}

  The mean is \eqn{(1-p_0) \mu}{(1-p0)m} and the variance is
  \eqn{(1-p_0) \sigma^2 + p_0(1-p_0) \mu^2}{(1-p0)v + p0(1-p0)m^2},
  where \eqn{\mu}{m} and \eqn{\sigma^2}{v} are the mean and variance of
  the zero-truncated Poisson.

  In the terminology of Klugman et al. (2012), the zero-modified
  Poisson is a member of the \eqn{(a, b, 1)} class of distributions
  with \eqn{a = 0} and \eqn{b = \lambda}.

  The special case \code{p0 == 0} is the zero-truncated Poisson.

  If an element of \code{x} is not integer, the result of
  \code{dzmpois} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dzmpois} gives the (log) probability mass function,
  \code{pzmpois} gives the (log) distribution function,
  \code{qzmpois} gives the quantile function, and
  \code{rzmpois} generates random deviates.

  Invalid \code{lambda} or \code{p0} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rzmpois}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}zmpois} use \code{\{d,p,q\}pois} for all
  but the trivial input values and \eqn{p(0)}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dpois}} for the standard Poisson distribution.

  \code{\link{dztpois}} for the zero-truncated Poisson distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dzmpois(0:5, lambda = 1, p0 = 0.2)
(1-0.2) * dpois(0:5, lambda = 1)/ppois(0, 1, lower = FALSE) # same

## simple relation between survival functions
pzmpois(0:5, 1, p0 = 0.2, lower = FALSE)
(1-0.2) * ppois(0:5, 1, lower = FALSE) /
    ppois(0, 1, lower = FALSE) # same

qzmpois(pzmpois(0:10, 1, p0 = 0.7), 1, p0 = 0.7)
}
\keyword{distribution}
