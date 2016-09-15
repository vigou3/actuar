\name{Paralogistic}
\alias{Paralogistic}
\alias{dparalogis}
\alias{pparalogis}
\alias{qparalogis}
\alias{rparalogis}
\alias{mparalogis}
\alias{levparalogis}
\title{The Paralogistic Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Paralogistic distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
dparalogis(x, shape, rate = 1, scale = 1/rate, log = FALSE)
pparalogis(q, shape, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
qparalogis(p, shape, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
rparalogis(n, shape, rate = 1, scale = 1/rate)
mparalogis(order, shape, rate = 1, scale = 1/rate)
levparalogis(limit, shape, rate = 1, scale = 1/rate,
             order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The paralogistic distribution with parameters \code{shape} \eqn{=
    \alpha}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\alpha^2 (x/\theta)^\alpha}{%
      x [1 + (x/\theta)^\alpha)^{\alpha + 1}}}{%
    f(x) = a^2 (x/s)^a / (x [1 + (x/s)^a]^(a + 1))}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\theta > 0}{b > 0}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, \eqn{-\alpha < k < \alpha^2}{-shape < k <
  shape^2}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]}, \eqn{k > -\alpha}{k > -shape}
  and \eqn{\alpha - k/\alpha}{shape - k/shape} not a negative integer.
}
\value{
  \code{dparalogis} gives the density,
  \code{pparalogis} gives the distribution function,
  \code{qparalogis} gives the quantile function,
  \code{rparalogis} generates random deviates,
  \code{mparalogis} gives the \eqn{k}th raw moment, and
  \code{levparalogis} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levparalogis} computes the limited expected value using
  \code{\link{betaint}}.

  See Kleiber and Kotz (2003) for alternative names and
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
exp(dparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pparalogis(qparalogis(p, 2, 3), 2, 3)

## variance
mparalogis(2, 2, 3) - mparalogis(1, 2, 3)^2

## case with shape - order/shape > 0
levparalogis(10, 2, 3, order = 2)

## case with shape - order/shape < 0
levparalogis(10, 1.25, 3, order = 2)
}
\keyword{distribution}
