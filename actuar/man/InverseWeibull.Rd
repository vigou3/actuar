\name{InverseWeibull}
\alias{InverseWeibull}
\alias{dinvweibull}
\alias{pinvweibull}
\alias{qinvweibull}
\alias{rinvweibull}
\alias{minvweibull}
\alias{levinvweibull}
\alias{dlgompertz}
\alias{plgompertz}
\alias{qlgompertz}
\alias{rlgompertz}
\alias{mlgompertz}
\alias{levlgompertz}
\title{The Inverse Weibull Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Inverse Weibull distribution
  with parameters \code{shape} and \code{scale}.
}
\usage{
  dinvweibull(x, shape, rate = 1, scale = 1/rate, log = FALSE)
  pinvweibull(q, shape, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
  qinvweibull(p, shape, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
  rinvweibull(n, shape, rate = 1, scale = 1/rate)
  minvweibull(order, shape, rate = 1, scale = 1/rate)
  levinvweibull(limit, shape, rate = 1, scale = 1/rate,
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
  The Inverse Weibull distribution with parameters \code{shape} \eqn{=
    \tau}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\tau (\theta/x)^\tau e^{-(\theta/x)^\tau}}{x}}{%
    f(x) = a (s/x)^a exp(-(s/x)^a)/x}
  for \eqn{x > 0}, \eqn{\tau > 0}{a > 0} and \eqn{\theta > 0}{s > 0}.

  The special case \code{shape == 1} is an
  \link[=dinvexp]{Inverse Exponential} distribution.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]}.
}
\value{
  \code{dinvweibull} gives the density,
  \code{pinvweibull} gives the distribution function,
  \code{qinvweibull} gives the quantile function,
  \code{rinvweibull} generates random deviates,
  \code{minvweibull} gives the \eqn{k}th raw moment, and
  \code{levinvweibull} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Distribution also knonw as the log-Gompertz.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dinvweibull(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvweibull(qinvweibull(p, 2, 3), 2, 3)
mlgompertz(-1, 3, 3)
levinvweibull(10, 2, 3, order = 2)
}
\keyword{distribution}
