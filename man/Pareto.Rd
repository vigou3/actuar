\name{Pareto}
\alias{Pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\alias{mpareto}
\alias{levpareto}
\alias{pareto2}
\alias{dpareto2}
\alias{ppareto2}
\alias{qpareto2}
\alias{rpareto2}
\alias{mpareto2}
\alias{levpareto2}
\title{The Pareto Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Pareto distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
dpareto(x, shape, scale, log = FALSE)
ppareto(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
qpareto(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
rpareto(n, shape, scale)
mpareto(order, shape, scale)
levpareto(limit, shape, scale, order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Pareto distribution with parameters \code{shape} \eqn{=
    \alpha}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\alpha \theta^\alpha}{(x + \theta)^{\alpha + 1}}}{%
    f(x) = a s^a / (x + s)^(a + 1)}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\theta}{s > 0}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}, \eqn{-1 < k < \alpha}{-1 < k < shape}.

  The \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]},
  \eqn{k > -1} and \eqn{\alpha - k}{shape - k} not a
  negative integer.
}
\value{
  \code{dpareto} gives the density,
  \code{ppareto} gives the distribution function,
  \code{qpareto} gives the quantile function,
  \code{rpareto} generates random deviates,
  \code{mpareto} gives the \eqn{k}th raw moment, and
  \code{levpareto} gives the \eqn{k}th moment of the limited loss variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levpareto} computes the limited expected value using
  \code{\link{betaint}}.

  Distribution also known as the Pareto Type II or Lomax distribution.
  See also Kleiber and Kotz (2003) for alternative names and
  parametrizations.
}
\references{
  Kleiber, C. and Kotz, S. (2003), \emph{Statistical Size Distributions
  in Economics and Actuarial Sciences}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto(qpareto(p, 2, 3), 2, 3)

## variance
mpareto(2, 4, 1) - mpareto(1, 4, 1)^2

## case with shape - order > 0
levpareto(10, 3, scale = 1, order = 2)

## case with shape - order < 0
levpareto(10, 1.5, scale = 1, order = 2)
}
\keyword{distribution}
