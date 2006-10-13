\name{InverseExponential}
\alias{InverseExponential}
\alias{dinvexp}
\alias{pinvexp}
\alias{qinvexp}
\alias{rinvexp}
\alias{minvexp}
\alias{levinvexp}
\title{The Inverse Exponential Distribution}
\description{
  Density, distribution function, quantile function, random generation
  raw moments and limited moments for the Inverse Exponential
  distribution with parameter \code{scale}.
}
\usage{
  dinvexp(x, rate = 1, scale = 1/rate, log = FALSE)
  pinvexp(q, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
  qinvexp(p, rate = 1, scale = 1/rate, lower.tail = TRUE, log.p = FALSE)
  rinvexp(n, rate = 1, scale = 1/rate)
  minvexp(order, rate = 1, scale = 1/rate)
  levinvexp(limit, rate = 1, scale = 1/rate, order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{scale}{parameter. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Inverse Exponential distribution with parameter \code{scale}
  \eqn{\theta}{= s} has density:
  \deqn{f(x) = \frac{s e^{-s/x}}{x^2}}{f(x) = s exp(-s/x)/x^2}
  for \eqn{x > 0} and \eqn{\theta > 0}{s > 0}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
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
exp(dinvexp(2, 2, log = TRUE))
p <- (1:10)/10
pinvexp(qinvexp(p, 2), 2)
minvexp(-2, 3)
}
\keyword{distribution}
