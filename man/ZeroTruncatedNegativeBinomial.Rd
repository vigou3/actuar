\name{ZeroTruncatedNegativeBinomial}
\alias{ZeroTruncatedNegativeBinomial}
\alias{ZTNegativeBinomial}
\alias{ZTNegBinomial}
\alias{dztnbinom}
\alias{pztnbinom}
\alias{qztnbinom}
\alias{rztnbinom}
\title{The Zero Truncated Negative Binomial Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the zero truncated negative binomial distribution with
  parameters \code{size} and \code{prob}.
}
\usage{
dztnbinom(x, size, prob, log = FALSE)
pztnbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
qztnbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
rztnbinom(n, size, prob)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{size}{target for number of successful trials, or dispersion
    parameter. Must be positive, need not be integer.}
  \item{prob}{parameter. \code{0 < prob <= 1}}.
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The random variable \eqn{Z} with probability mass function
  \deqn{P[Z = x] = \frac{P[X = x]}{1 - P[X = 0]}}{%
    P[Z = x] = P[X = x]/(1 - P[X = 0])}
  and distribution function
  \deqn{P[Z \le x] = \frac{P[X \le x] - P[X = 0]}{1 - P[X = 0]}}{%
    P[Z <= x] = (P[X <= x] - P[X = 0])/(1 - P[X = 0])},
  for \eqn{x = 1, 2, \ldots}, has a zero truncated negative binomial
  distribution when the random variable \eqn{X} has a
  \link[=NegBinomial]{negative binomial} distribution.

  The zero truncated negative binomial distribution with \code{size}
  \eqn{= r} and \code{prob} \eqn{= p} has probability mass function
  \deqn{%
    P[X = x] = \frac{\Gamma(x + r) p^r (1 - p)^x}{\Gamma(r) x! (1 - p^r)}}{%
    P[X = x] = [\Gamma(x + r) p^r (1 - p)^x]/[\Gamma(n) x! (1 - p^r)]}
  for \eqn{x = 1, 2, \ldots}, \eqn{r \ge 0} and \eqn{0 < p < 1}, and
  \eqn{P[X = 1] = 1} when \eqn{p = 1}.

  The mean is \eqn{r(1 - p)/(p(1 - p^r))} and the variance is
  \eqn{[r(1 - p)(1 - p^r(1 + r(1 - p))]/[p(1 - p^r)^2]}.

  In the terminology of Klugman et al. (2012), the zero truncated
  negative binomial is a member of the \eqn{(a, b, 1)} family of
  distributions with \eqn{a = 1 - p} and \eqn{b = (r - 1)(1 - p)}.

  The limiting case \code{size == 0} is the
  \link[=Logarithmic]{logarithmic} distribution with parameter \code{1 -
  prob}.

  Unlike the standard negative binomial functions, parametrization
  through the mean \code{mu} is not supported to avoid ambiguity as
  to whether \code{mu} is the mean of the underlying negative binomial
  or the mean of the zero truncated distribution.

  If an element of \code{x} is not integer, the result of
  \code{dztnbinom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
  \code{dztnbinom} gives the probability mass function,
  \code{pztnbinom} gives the distribution function,
  \code{qztnbinom} gives the quantile function, and
  \code{rztnbinom} generates random deviates.

  Invalid \code{size} or \code{prob} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rztnbinom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}ztnbinom} use their equivalent
  \code{\{d,p,q\}nbinom} function for all but the trivial input values.

  \code{rztnbinom} uses the simple inversion algorithm suggested by
  Peter Dalgaard on the r-help mailing list on 1 May 2005 %
  (\url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}).
}
\references{
  Johnson, N. L., Kemp, A. W. and Kotz, S. (2005), \emph{Univariate
  Discrete Distributions, Third Edition}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dnbinom}} for the negative binomial distribution.

  \code{\link{dztgeom}} for the zero truncated geometric and
  \code{\link{dlogarithmic}} for the logarithmic, which are a special
  case of the zero truncated negative binomial.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Example 6.3 of Klugman et al. (2012)
p <- 1/(1 + 0.5)
dztnbinom(c(1, 2, 3), size = 2.5, prob = p)
dnbinom(c(1, 2, 3), 2.5, p)/pnbinom(0, 2.5, p, lower = FALSE) # same

## Limiting cases
pztnbinom(1, 2, prob = 1)        # point mass at 1
dztnbinom(2, size = 1, 0.25)     # == dztgeom(2, 0.25)
dztnbinom(2, size = 0, 0.25)     # == dlogarithmic(2, 0.75)

qztnbinom(pztnbinom(1:10, 2.5, 0.3), 2.5, 0.3)

x <- rztnbinom(1000, size = 2.5, prob = 0.4)
y <- sort(unique(x))
plot(y, table(x)/length(x), type = "h", lwd = 2,
     pch = 19, col = "black", xlab = "x", ylab = "p(x)",
     main = "Empirical vs theoretical probabilities")
points(y, dztnbinom(y, size = 2.5, prob = 0.4),
       pch = 19, col = "red")
legend("topright", c("empirical", "theoretical"),
       lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
}
\keyword{distribution}
