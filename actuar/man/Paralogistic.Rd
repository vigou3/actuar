
\name{Paralogis}
\alias{paralogis}
\alias{print.paralogis}
\alias{dparalogis}
\alias{pparalogis}
\alias{qparalogis}
\alias{rparalogis}
\alias{mparalogis}
\alias{levparalogis}
\title{The Paralogistic Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Paralogistic distribution with parameters \code{shape} and \code{scale}.
}
\usage{
  dparalogis(x, shape, scale, log = FALSE)
  pparalogis(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
  qparalogis(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
  rparalogis(n, shape, scale)
  mparalogis(k, shape, scale, log = FALSE)
  levparalogis(d, shape, scale, order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{d}{limit of the loss variable.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k, order}{order of the moment.}
}
\details{
  The Paralogistic distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{a^2 (\frac{x}{b})^a}{x(1 + (\frac{x}{b})^a)^{(a +
      1)}}}{p(x) = a^2 (x/b)^a / (x(1 + (x/b)^a)^(a+1))}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}.
}
\value{
\code{dparalogis} gives the density, \code{pparalogis} gives the distribution function,
\code{qparalogis} gives the quantile function, \code{rparalogis} generates random
deviates, \code{mparalogis} gives the kth raw moment, and \code{levparalogis} gives
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
exp(dparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pparalogis(qparalogis(p, 2, 3), 2, 3)
mparalogis(2, 2, 3) - mparalogis(1, 2, 3) ^ 2
levparalogis(10, 2, 3, 2)
}
\keyword{distribution}
