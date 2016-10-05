\name{ZeroModifiedLogarithmic}
\alias{ZeroModifiedLogarithmic}
\alias{ZMLogarithmic}
\alias{dzmlogarithmic}
\alias{pzmlogarithmic}
\alias{qzmlogarithmic}
\alias{rzmlogarithmic}
\title{The Zero-Modified Logarithmic Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Modified Logarithmic (or log-series)
  distribution with parameter \code{prob} and arbitrary probability at
  zero \code{p0}.
}
\usage{
dzmlogarithmic(x, prob, p0, log = FALSE)
pzmlogarithmic(q, prob, p0, lower.tail = TRUE, log.p = FALSE)
qzmlogarithmic(p, prob, p0, lower.tail = TRUE, log.p = FALSE)
rzmlogarithmic(n, prob, p0)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{prob}{parameter. \code{0 <= prob < 1}.}
  \item{p0}{probability mass at zero. \code{0 <= p0 <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-modified logarithmic distribution with \code{prob} \eqn{= p}
  and \code{p0} \eqn{= p_0}{= p0} is a discrete mixture between a
  degenerate distribution at zero and a (standard) logarithmic. The
  probability mass function is \eqn{p(0) = p_0}{p(0) = p0} and
  \deqn{%
    p(x) = (1-p_0) f(x)}{p(x) = (1-p0) f(x)}
  for \eqn{x = 1, 2, \ldots}, \eqn{0 < p < 1} and \eqn{0 \le
  p_0 \le 1}{0 \le p0 \le 1}, where \eqn{f(x)} is the probability mass
  function of the logarithmic.
  The cumulative distribution function is
  \deqn{P(x) = p_0 + (1 - p_0) F(x)}{P(x) = p0 + (1 - p0) F(x).}

  The special case \code{p0 == 0} is the standard logarithmic.

  The zero-modified logarithmic distribution is the limiting case of the
  zero-modified negative binomial distribution with \code{size}
  parameter equal to \eqn{0}. Note that in this context, parameter
  \code{prob} generally corresponds to the probability of \emph{failure}
  of the zero-truncated negative binomial.

  If an element of \code{x} is not integer, the result of
  \code{dzmlogarithmic} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
  \code{dzmlogarithmic} gives the probability mass function,
  \code{pzmlogarithmic} gives the distribution function,
  \code{qzmlogarithmic} gives the quantile function, and
  \code{rzmlogarithmic} generates random deviates.

  Invalid \code{prob} or \code{p0} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rzmlogarithmic}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}zmlogarithmic} use
  \code{\{d,p,q\}logarithmic} for all but the trivial input values and
  \eqn{p(0)}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dlogarithmic}} for the logarithmic distribution.

  \code{\link{dztnbinom}} for the zero modified negative binomial
  distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
p <- 1/(1 + 0.5)
dzmlogarithmic(1:5, prob = p, p0 = 0.6)
(1-0.6) * dlogarithmic(1:5, p)/plogarithmic(0, p, lower = FALSE) # same

## simple relation between survival functions
pzmlogarithmic(0:5, p, p0 = 0.2, lower = FALSE)
(1-0.2) * plogarithmic(0:5, p, lower = FALSE)/plogarithmic(0, p, lower = FALSE) # same

qzmlogarithmic(pzmlogarithmic(0:10, 0.3, p0 = 0.6), 0.3, p0 = 0.6)
}
\keyword{distribution}
