\name{ZeroTruncatedPoisson}
\alias{ZeroTruncatedPoisson}
\alias{dztpois}
\alias{pztpois}
\alias{qztpois}
\alias{rztpois}
\title{The Zero Truncated Poisson Distribution}
\description{
  Density function, distribution function, quantile function, random generation
  for the zero truncated Poisson distribution with parameter \code{lambda}.
}
\usage{
dztpois(x, lambda, log = FALSE)
pztpois(q, lambda, lower.tail = TRUE, log.p = FALSE)
qztpois(p, lambda, lower.tail = TRUE, log.p = FALSE)
rztpois(n, lambda)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of values to return.}
  \item{lambda}{vector of (non negative) means.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The Zero Truncated Poisson distribution has probability mass function
  \deqn{p(x) = \frac{e^{-/lambda} \lambda^x}{x! (1 - e^{-\lambda})}
      = \frac{\lambda^x}{x! (e^{\lambda} - 1)}}{%
      p(x) = lambda^x exp(-lambda)/[x! (1 - exp(-lambda))]
      = lambda^x/[x! (e^lambda - 1)]}
  for \eqn{x = 1, 2, ...}. If \eqn{P(x)} is the cumulative distribution
  function of the Zero Truncated Poisson and \eqn{F(x)} is that of the
  regular Poisson distribution, then we have:
  \deqn{P(x) = \frac{F(x) - F(0)}{1 - F(0)}}{%
    P(x) = [F(x) - F(0)][1 - F(0)]}
  for \eqn{x = 1, 2, ...} and \eqn{F(0) = e^{-\lambda}}{exp(-lambda)}.

  The Zero Truncated Poisson is a member of the \eqn{(a, b, 1)} family
  of distributions with \eqn{a = 0}, \eqn{b = \lambda}{b = lambda} and
  \deqn{p(1) = \frac{\lambda}{e^\lambda - 1}.}{%
    p(1) = lambda/(exp(lambda) - 1).}

  Most computations are done using the corresponding function for the
  Poisson distribution, so one can expect the same behaviour.
}
\value{
  \code{dztpois} gives the (log) density,
  \code{pztpois} gives the (log) distribution function,
  \code{qztpois} gives the quantile function, and
  \code{rztpois} generates random deviates.

  Invalid \code{lambda} will result in return value \code{NaN}, with a warning.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2008),
  \emph{Loss Models, From Data to Decisions, Third Edition}, Wiley.
}
\seealso{
  \code{\link{dpois}} for the standard Poisson distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dztpois(0:10, lambda = 1)

qztpois(pztpois(1:8, 1), 1)

x <- seq(0, 8)
plot(x, dztpois(x, 2), type = "h", ylab = "p(x)", main = "Zero Truncated Poisson(2) and Poisson(2) PDF")
points(x, dpois(x, 2), col = "red")
legend("topright", c("ZT Poisson probabilities", "Poisson probabilities"),
       col = c("black", "red"), lty = c(1, 0), pch = c(NA, 1))
}
\keyword{distribution}
