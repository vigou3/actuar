\name{ZeroTruncatedPoisson}
\alias{ZeroTruncatedPoisson}
\alias{ZTPoisson}
\alias{dztpois}
\alias{pztpois}
\alias{qztpois}
\alias{rztpois}
\title{The Zero-Truncated Poisson Distribution}
\description{
  Density function, distribution function, quantile function, random
  generation for the Zero-Truncated Poisson distribution with parameter
  \code{lambda}.
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
  The zero-truncated Poisson distribution has probability mass function
  \deqn{%
    p(x) = \frac{e^{-/lambda} \lambda^x}{x! (1 - e^{-\lambda})}
    = \frac{\lambda^x}{x! (e^{\lambda} - 1)}}{%
    p(x) = lambda^x exp(-lambda)/[x! (1 - exp(-lambda))]
         = lambda^x/[x! (e^lambda - 1)]}
  for \eqn{x = 1, 2, ...}, and \eqn{p(1) = 1} when \eqn{\lambda = 0}.
  The cumulative distribution function is
  \deqn{P(x) = \frac{F(x) - F(0)}{1 - F(0)},}{%
    P(x) = [F(x) - F(0)]/[1 - F(0)],}
  where \eqn{F(x)} is the distribution function of the standard Poisson.

  The mean is \eqn{\lambda/(1 - e^{-\lambda})^2}{\lambda/(1 -
    exp(-\lambda))} and the variance is
  \eqn{\lambda[1 - (\lambda+1)e^{-\lambda}]/(1 - e^{-\lambda})^2}{%
    \lambda[1 - (\lambda+1)exp(-\lambda)]/(1 - exp(-\lambda))^2}.

  In the terminology of Klugman et al. (2012), the zero-truncated
  Poisson is a member of the \eqn{(a, b, 1)} class of distributions
  with \eqn{a = 0} and \eqn{b = \lambda}.

  If an element of \code{x} is not integer, the result of
  \code{dztpois} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dztpois} gives the (log) probability mass function,
  \code{pztpois} gives the (log) distribution function,
  \code{qztpois} gives the quantile function, and
  \code{rztpois} generates random deviates.

  Invalid \code{lambda} will result in return value \code{NaN}, with a
  warning.

  The length of the result is determined by \code{n} for
  \code{rztpois}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}ztpois} use \code{\{d,p,q\}pois} for all
  but the trivial input values and \eqn{p(0)}.

  \code{rztpois} uses the simple inversion algorithm suggested by
  Peter Dalgaard on the r-help mailing list on 1 May 2005 %
  (\url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}).
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dpois}} for the standard Poisson distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dztpois(1:5, lambda = 1)
dpois(1:5, lambda = 1)/ppois(0, 1, lower = FALSE) # same

pztpois(1, lambda = 0)          # point mass at 1

qztpois(pztpois(1:10, 1), 1)

x <- seq(0, 8)
plot(x, dztpois(x, 2), type = "h", lwd = 2, ylab = "p(x)",
     main = "Zero-Truncated Poisson(2) and Poisson(2) PDF")
points(x, dpois(x, 2), pch = 19, col = "red")
legend("topright", c("ZT Poisson probabilities", "Poisson probabilities"),
       col = c("black", "red"), lty = c(1, 0), lwd = 2, pch = c(NA, 19))
}
\keyword{distribution}
