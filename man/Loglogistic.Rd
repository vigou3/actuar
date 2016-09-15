\name{Loglogistic}
\alias{Loglogistic}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{rllogis}
\alias{mllogis}
\alias{levllogis}
\title{The Loglogistic Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Loglogistic distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
dllogis(x, shape, rate = 1, scale = 1/rate, log = FALSE)
pllogis(q, shape, rate = 1, scale = 1/rate,
        lower.tail = TRUE, log.p = FALSE)
qllogis(p, shape, rate = 1, scale = 1/rate,
        lower.tail = TRUE, log.p = FALSE)
rllogis(n, shape, rate = 1, scale = 1/rate)
mllogis(order, shape, rate = 1, scale = 1/rate)
levllogis(limit, shape, rate = 1, scale = 1/rate,
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
  The loglogistic distribution with parameters \code{shape} \eqn{=
    \gamma}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\gamma (x/\theta)^\gamma}{%
      x [1 + (x/\theta)^\gamma]^2}}{%
    f(x) = a (x/s)^a / (x [1 + (x/s)^a]^2)}
  for \eqn{x > 0}, \eqn{\gamma > 0}{a > 0} and \eqn{\theta > 0}{b > 0}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}, \eqn{-\gamma < k < \gamma}{-shape < k < shape}.

  The \eqn{k}th limited moment at some limit \eqn{d} is \eqn{E[\min(X,
  d)^k]}{E[min(X, d)^k]}, \eqn{k > -\gamma}{k > -shape}
  and \eqn{1 - k/\gamma}{1 - k/shape} not a negative integer.
}
\value{
  \code{dllogis} gives the density,
  \code{pllogis} gives the distribution function,
  \code{qllogis} gives the quantile function,
  \code{rllogis} generates random deviates,
  \code{mllogis} gives the \eqn{k}th raw moment, and
  \code{levllogis} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levllogis} computes the limited expected value using
  \code{\link{betaint}}.

  Also known as the Fisk distribution. See also Kleiber and Kotz (2003)
  for alternative names and parametrizations.
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
exp(dllogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pllogis(qllogis(p, 2, 3), 2, 3)

## mean
mllogis(1, 2, 3)

## case with 1 - order/shape > 0
levllogis(10, 2, 3, order = 1)

## case with 1 - order/shape < 0
levllogis(10, 2/3, 3, order = 1)
}
\keyword{distribution}
