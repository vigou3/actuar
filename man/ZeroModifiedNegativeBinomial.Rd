\name{ZeroModifiedNegativeBinomial}
\alias{ZeroModifiedNegativeBinomial}
\alias{ZMNegativeBinomial}
\alias{ZMNegBinomial}
\alias{dzmnbinom}
\alias{pzmnbinom}
\alias{qzmnbinom}
\alias{rzmnbinom}
\title{The Zero-Modified Negative Binomial Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Modified Negative Binomial distribution with
  parameters \code{size} and \code{prob}, and arbitrary probability at
  zero \code{p0}.
}
\usage{
dzmnbinom(x, size, prob, p0, log = FALSE)
pzmnbinom(q, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
qzmnbinom(p, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
rzmnbinom(n, size, prob, p0)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{size}{target for number of successful trials, or dispersion
    parameter. Must be positive, need not be integer.}
  \item{prob}{parameter. \code{0 < prob <= 1}.}
  \item{p0}{probability mass at zero. \code{0 <= p0 <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-modified negative binomial distribution with \code{size}
  \eqn{= r}, \code{prob} \eqn{= p} and \code{p0} \eqn{= p_0}{= p0} is a
  discrete mixture between a degenerate distribution at zero and a
  (standard) negative binomial. The probability mass function is
  \eqn{p(0) = p_0}{p(0) = p0} and
  \deqn{%
    p(x) = \frac{(1-p_0)}{(1-p^r)} f(x)}{%
    p(x) = (1-p0)/(1-p^r) f(x)}
  for \eqn{x = 1, 2, \ldots}, \eqn{r \ge 0}, \eqn{0 < p < 1} and \eqn{0 \le
  p_0 \le 1}{0 \le p0 \le 1}, where \eqn{f(x)} is the probability mass
  function of the negative binomial.
  The cumulative distribution function is
  \deqn{P(x) = p_0 + (1 - p_0) \left(\frac{F(x) - F(0)}{1 - F(0)}\right)}{%
    P(x) = p0 + (1 - p0) [F(x) - F(0)]/[1 - F(0)].}

  The mean is \eqn{(1-p_0) \mu}{(1-p0)m} and the variance is
  \eqn{(1-p_0) \sigma^2 + p_0(1-p_0) \mu^2}{(1-p0)v + p0(1-p0)m^2},
  where \eqn{\mu}{m} and \eqn{\sigma^2}{v} are the mean and variance of
  the zero-truncated negative binomial.

  In the terminology of Klugman et al. (2012), the zero-modified
  negative binomial is a member of the \eqn{(a, b, 1)} class of
  distributions with \eqn{a = 1-p} and \eqn{b = (r-1)(1-p)}.

  The special case \code{p0 == 0} is the zero-truncated negative
  binomial.

  The limiting case \code{size == 0} is the zero-modified logarithmic
  distribution with parameters \code{1 - prob} and \code{p0}.

  Unlike the standard negative binomial functions, parametrization
  through the mean \code{mu} is not supported to avoid ambiguity as
  to whether \code{mu} is the mean of the underlying negative binomial
  or the mean of the zero-modified distribution.

  If an element of \code{x} is not integer, the result of
  \code{dzmnbinom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dzmnbinom} gives the (log) probability mass function,
  \code{pzmnbinom} gives the (log) distribution function,
  \code{qzmnbinom} gives the quantile function, and
  \code{rzmnbinom} generates random deviates.

  Invalid \code{size}, \code{prob} or \code{p0} will result in return
  value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rzmnbinom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}zmnbinom} use \code{\{d,p,q\}nbinom} for all
  but the trivial input values and \eqn{p(0)}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dnbinom}} for the negative binomial distribution.

  \code{\link{dztnbinom}} for the zero-truncated negative binomial
  distribution.

  \code{\link{dzmgeom}} for the zero-modified geometric and
  \code{\link{dzmlogarithmic}} for the zero-modified logarithmic, which
  are special cases of the zero-modified negative binomial.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Example 6.3 of Klugman et al. (2012)
p <- 1/(1 + 0.5)
dzmnbinom(1:5, size = 2.5, prob = p, p0 = 0.6)
(1-0.6) * dnbinom(1:5, 2.5, p)/pnbinom(0, 2.5, p, lower = FALSE) # same

## simple relation between survival functions
pzmnbinom(0:5, 2.5, p, p0 = 0.2, lower = FALSE)
(1-0.2) * pnbinom(0:5, 2.5, p, lower = FALSE) /
    pnbinom(0, 2.5, p, lower = FALSE) # same

qzmnbinom(pzmnbinom(0:10, 2.5, 0.3, p0 = 0.1), 2.5, 0.3, p0 = 0.1)
}
\keyword{distribution}
