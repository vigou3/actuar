\name{Invexp}
\alias{invexp}
\alias{print.invexp}
\alias{dinvexp}
\alias{pinvexp}
\alias{qinvexp}
\alias{rinvexp}
\alias{minvexp}
\alias{levinvexp}
\title{The Inverse Exponential Distribution}
\description{
  Density, distribution function, quantile function, random generation
  raw moments, and limited moments for the Inverse Exponential distribution with parameter \code{scale}. 
}
\usage{
  dinvexp(x, scale, log = FALSE)
  pinvexp(q, scale, lower.tail = TRUE, log.p = FALSE)
  qinvexp(p, scale, lower.tail = TRUE, log.p = FALSE)
  rinvexp(n, scale)
  minvexp(k, scale, log = FALSE)
  levinvexp(d, scale, order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
   \item{d}{limit of the loss variable.}
  \item{scale}{parameter. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k, order}{order of the moment.}
}
\details{
  The Inverse Exponential distribution with parameter \code{scale} = a has
  density: \deqn{p(x) = \frac{a e^{\frac{-a}{x}}}{x^2}}{p(x) = a exp(-a / x) /
    x^2} for \eqn{x > 0} and \eqn{a > 0}.
}
\value{
\code{dinvexp} gives the density, \code{pinvexp} gives the distribution function,
\code{qinvexp} gives the quantile function, \code{rinvexp} generates random
deviates, \code{minvexp} gives the kth raw moment, and \code{levinvexp} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\references{
  Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dinvexp(2, 2, log = TRUE))
p <- (1:10)/10
pinvexp(qinvexp(p, 2), 2)
minvexp(2, 3) - minvexp(1, 3) ^ 2
levinvexp(10, 3, 2)
}
\keyword{distribution}
