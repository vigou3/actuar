\name{ZeroTruncatedNegativeBinomial}
\alias{ZeroTruncatedNegativeBinomial}
\alias{ZTNegativeBinomial}
\alias{ZTNegBinomial}
\alias{dztnbinom}
\alias{pztnbinom}
\alias{qztnbinom}
\alias{rztnbinom}
\title{The Zero-Truncated Negative Binomial Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Truncated Negative Binomial distribution with
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
  \item{prob}{parameter. \code{0 < prob <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-truncated negative binomial distribution with \code{size}
  \eqn{= r} and \code{prob} \eqn{= p} has probability mass function
  \deqn{%
    p(x) = \frac{\Gamma(x + r) p^r (1 - p)^x}{\Gamma(r) x! (1 - p^r)}}{%
    p(x) = [\Gamma(x+r) p^r (1-p)^x]/[\Gamma(n) x! (1-p^r)]}
  for \eqn{x = 1, 2, \ldots}, \eqn{r \ge 0} and \eqn{0 < p < 1}, and
  \eqn{p(1) = 1} when \eqn{p = 1}.
  The cumulative distribution function is
  \deqn{P(x) = \frac{F(x) - F(0)}{1 - F(0)},}{%
    P(x) = [F(x) - F(0)]/[1 - F(0)],}
  where \eqn{F(x)} is the distribution function of the standard negative
  binomial.

  The mean is \eqn{r(1-p)/(p(1-p^r))} and the variance is
  \eqn{[r(1-p)(1 - (1 + r(1-p))p^r)]/[p(1-p^r)]^2}.

  In the terminology of Klugman et al. (2012), the zero-truncated
  negative binomial is a member of the \eqn{(a, b, 1)} class of
  distributions with \eqn{a = 1-p} and \eqn{b = (r-1)(1-p)}.

  The limiting case \code{size == 0} is the
  \link[=Logarithmic]{logarithmic} distribution with parameter \code{1 -
  prob}.

  Unlike the standard negative binomial functions, parametrization
  through the mean \code{mu} is not supported to avoid ambiguity as
  to whether \code{mu} is the mean of the underlying negative binomial
  or the mean of the zero-truncated distribution.

  If an element of \code{x} is not integer, the result of
  \code{dztnbinom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dztnbinom} gives the (log) probability mass function,
  \code{pztnbinom} gives the (log) distribution function,
  \code{qztnbinom} gives the quantile function, and
  \code{rztnbinom} generates random deviates.

  Invalid \code{size} or \code{prob} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rztnbinom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}ztnbinom} use \code{\{d,p,q\}nbinom} for all
  but the trivial input values and \eqn{p(0)}.

  \code{rztnbinom} uses the simple inversion algorithm suggested by
  Peter Dalgaard on the r-help mailing list on 1 May 2005 %
  (\url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}).
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dnbinom}} for the negative binomial distribution.

  \code{\link{dztgeom}} for the zero-truncated geometric and
  \code{\link{dlogarithmic}} for the logarithmic, which are special
  cases of the zero-truncated negative binomial.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Example 6.3 of Klugman et al. (2012)
p <- 1/(1 + 0.5)
dztnbinom(c(1, 2, 3), size = 2.5, prob = p)
dnbinom(c(1, 2, 3), 2.5, p)/pnbinom(0, 2.5, p, lower = FALSE) # same

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
       lty = c(1, NA), lwd = 2, pch = c(NA, 19), col = c("black", "red"))
}
\keyword{distribution}
