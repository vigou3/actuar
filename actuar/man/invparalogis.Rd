
\name{Invparalogis}
\alias{invparalogis}
\alias{print.invparalogis}
\alias{dinvparalogis}
\alias{pinvparalogis}
\alias{qinvparalogis}
\alias{rinvparalogis}
\alias{minvparalogis}
\alias{levinvparalogis}
\title{The Inverse Paralogistic Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Inverse Paralogistic distribution with parameters \code{shape} and \code{scale}.
}
\usage{
  dinvparalogis(x, shape, scale, log = FALSE)
  pinvparalogis(q, shape, scale, lower.tail = TRUE,
           log.p = FALSE)
  qinvparalogis(p, shape, scale, lower.tail = TRUE,
           log.p = FALSE)
  rinvparalogis(n, shape, scale)
  minvparalogis(k, shape, scale, log = FALSE)
  levinvparalogis(d, shape, scale, order = 1, log = FALSE)
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
  The Inverse Paralogistic distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{a^2 (\frac{x}{b})^{a^2}}{x(1 + (\frac{x}{b})^a)^{(a
	+ 1)}}}{p(x) = a^2 (x / b)^(a^2) / (x(1 + (x/b)^a)^(a + 1))}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
\code{dinvparalogis} gives the density, \code{pinvparalogis} gives the distribution function,
\code{qinvparalogis} gives the quantile function, \code{rinvparalogis} generates random
deviates, \code{minvparalogis} gives the kth raw moment, and \code{levinvparalogis} gives
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
exp(dinvparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvparalogis(qinvparalogis(p, 2, 3), 2, 3)
minvparalogis(-1, 2, 2)
levinvparalogis(10, 1, 2, 2)
}
\keyword{distribution}
