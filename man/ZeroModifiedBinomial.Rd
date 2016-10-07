\name{ZeroModifiedBinomial}
\alias{ZeroModifiedBinomial}
\alias{ZMBinomial}
\alias{dzmbinom}
\alias{pzmbinom}
\alias{qzmbinom}
\alias{rzmbinom}
\title{The Zero-Modified Binomial Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Modified Binomial distribution with
  parameters \code{size} and \code{prob}, and probability at zero
  \code{p0}.
}
\usage{
dzmbinom(x, size, prob, p0, log = FALSE)
pzmbinom(q, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
qzmbinom(p, size, prob, p0, lower.tail = TRUE, log.p = FALSE)
rzmbinom(n, size, prob, p0)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{size}{number of trials (strictly positive integer).}
  \item{prob}{probability of success on each trial. \code{0 <= prob <= 1}.}
  \item{p0}{probability mass at zero. \code{0 <= p0 <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-modified binomial distribution with \code{size} \eqn{= n},
  \code{prob} \eqn{= p} and \code{p0} \eqn{= p_0}{= p0} is a discrete
  mixture between a degenerate distribution at zero and a (standard)
  binomial. The probability mass function is \eqn{p(0) = p_0}{p(0) = p0}
  and
  \deqn{%
    p(x) = \frac{(1-p_0)}{(1 - (1-p)^n)} f(x)}{%
    p(x) = (1-p0)/[1 - (1-p)^n] f(x)}
  for \eqn{x = 1, \ldots, n}, \eqn{0 < p \le 1} and \eqn{0 \le
  p_0 \le 1}{0 \le p0 \le 1}, where \eqn{f(x)} is the probability mass
  function of the binomial.
  The cumulative distribution function is
  \deqn{P(x) = p_0 + (1 - p_0) \left(\frac{F(x) - F(0)}{1 - F(0)}\right)}{%
    P(x) = p0 + (1 - p0) [F(x) - F(0)]/[1 - F(0)].}

  The mean is \eqn{(1-p_0) \mu}{(1-p0)m} and the variance is
  \eqn{(1-p_0) \sigma^2 + p_0(1-p_0) \mu^2}{(1-p0)v + p0(1-p0)m^2},
  where \eqn{\mu}{m} and \eqn{\sigma^2}{v} are the mean and variance of
  the zero-truncated binomial.

  In the terminology of Klugman et al. (2012), the zero-modified
  binomial is a member of the \eqn{(a, b, 1)} class of
  distributions with \eqn{a = -p/(1-p)} and \eqn{b = (n+1)p/(1-p)}.

  The special case \code{p0 == 0} is the zero-truncated binomial.

  If an element of \code{x} is not integer, the result of
  \code{dzmbinom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dzmbinom} gives the probability mass function,
  \code{pzmbinom} gives the distribution function,
  \code{qzmbinom} gives the quantile function, and
  \code{rzmbinom} generates random deviates.

  Invalid \code{size}, \code{prob} or \code{p0} will result in return
  value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rzmbinom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}zmbinom} use \code{\{d,p,q\}binom} for all
  but the trivial input values and \eqn{p(0)}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dbinom}} for the binomial distribution.

  \code{\link{dztbinom}} for the zero-truncated binomial distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dzmbinom(1:5, size = 5, prob = 0.4, p0 = 0.2)
(1-0.2) * dbinom(1:5, 5, 0.4)/pbinom(0, 5, 0.4, lower = FALSE) # same

## simple relation between survival functions
pzmbinom(0:5, 5, 0.4, p0 = 0.2, lower = FALSE)
(1-0.2) * pbinom(0:5, 5, 0.4, lower = FALSE) /
    pbinom(0, 5, 0.4, lower = FALSE) # same

qzmbinom(pzmbinom(1:10, 10, 0.6, p0 = 0.1), 10, 0.6, p0 = 0.1)

n <- 8; p <- 0.3; p0 <- 0.025
x <- 0:n
title <- paste("ZM Binomial(", n, ", ", p, ", p0 = ", p0,
               ") and Binomial(", n, ", ", p,") PDF",
               sep = "")
plot(x, dzmbinom(x, n, p, p0), type = "h", lwd = 2, ylab = "p(x)",
     main = title)
points(x, dbinom(x, n, p), pch = 19, col = "red")
legend("topright", c("ZT binomial probabilities", "Binomial probabilities"),
       col = c("black", "red"), lty = c(1, 0), lwd = 2, pch = c(NA, 19))
}
\keyword{distribution}
