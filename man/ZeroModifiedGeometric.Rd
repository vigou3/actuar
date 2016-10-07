\name{ZeroModifiedGeometric}
\alias{ZeroModifiedGeometric}
\alias{Zmgeometric}
\alias{dzmgeom}
\alias{pzmgeom}
\alias{qzmgeom}
\alias{rzmgeom}
\title{The Zero-Modified Geometric Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Modified Geometric distribution with
  parameter \code{prob} and arbitrary probability at zero \code{p0}.
}
\usage{
dzmgeom(x, prob, p0, log = FALSE)
pzmgeom(q, prob, p0, lower.tail = TRUE, log.p = FALSE)
qzmgeom(p, prob, p0, lower.tail = TRUE, log.p = FALSE)
rzmgeom(n, prob, p0)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{prob}{parameter. \code{0 < prob <= 1}.}
  \item{p0}{probability mass at zero. \code{0 <= p0 <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-modified geometric distribution with \code{prob} \eqn{= p}
  and \code{p0} \eqn{= p_0}{= p0} is a discrete mixture between a
  degenerate distribution at zero and a (standard) geometric. The
  probability mass function is \eqn{p(0) = p_0}{p(0) = p0} and
  \deqn{%
    p(x) = \frac{(1-p_0)}{(1-p)} f(x)}{%
    p(x) = (1-p0)/(1-p) f(x)}
  for \eqn{x = 1, 2, \ldots}, \eqn{0 < p < 1} and \eqn{0 \le
  p_0 \le 1}{0 \le p0 \le 1}, where \eqn{f(x)} is the probability mass
  function of the geometric.
  The cumulative distribution function is
  \deqn{P(x) = p_0 + (1 - p_0) \left(\frac{F(x) - F(0)}{1 - F(0)}\right)}{%
    P(x) = p0 + (1 - p0) [F(x) - F(0)]/[1 - F(0)].}

  The mean is \eqn{(1-p_0) \mu}{(1-p0)m} and the variance is
  \eqn{(1-p_0) \sigma^2 + p_0(1-p_0) \mu^2}{(1-p0)v + p0(1-p0)m^2},
  where \eqn{\mu}{m} and \eqn{\sigma^2}{v} are the mean and variance of
  the zero-truncated geometric.

  In the terminology of Klugman et al. (2012), the zero-modified
  geometric is a member of the \eqn{(a, b, 1)} class of
  distributions with \eqn{a = 1-p} and \eqn{b = 0}.

  The special case \code{p0 == 0} is the zero-truncated geometric.

  If an element of \code{x} is not integer, the result of
  \code{dzmgeom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dzmgeom} gives the (log) probability mass function,
  \code{pzmgeom} gives the (log) distribution function,
  \code{qzmgeom} gives the quantile function, and
  \code{rzmgeom} generates random deviates.

  Invalid \code{prob} or \code{p0} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rzmgeom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}zmgeom} use \code{\{d,p,q\}geom} for all but
  the trivial input values and \eqn{p(0)}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dgeom}} for the geometric distribution.

  \code{\link{dztgeom}} for the zero-truncated geometric distribution.

  \code{\link{dzmnbinom}} for the zero-modified negative binomial, of
  which the zero-modified geometric is a special case.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
p <- 1/(1 + 0.5)
dzmgeom(1:5, prob = p, p0 = 0.6)
(1-0.6) * dgeom(1:5, p)/pgeom(0, p, lower = FALSE) # same

## simple relation between survival functions
pzmgeom(0:5, p, p0 = 0.2, lower = FALSE)
(1-0.2) * pgeom(0:5, p, lower = FALSE)/pgeom(0, p, lower = FALSE) # same

qzmgeom(pzmgeom(0:10, 0.3, p0 = 0.6), 0.3, p0 = 0.6)
}
\keyword{distribution}
