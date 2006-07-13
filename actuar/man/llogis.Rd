
\name{Llogis}
\alias{llogis}
\alias{print.llogis}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{rllogis}
\alias{mllogis}
\alias{levllogis}
\title{The Loglogistic Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Loglogistic distribution with parameters \code{shape} and \code{scale}.
}
\usage{
  dllogis(x, shape, scale, log = FALSE)
  pllogis(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
  qllogis(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
  rllogis(n, shape, scale)
  mllogis(k, shape, scale, log = FALSE)
  levllogis(d, shape, scale, order = 1, log = FALSE)
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
  The Loglogistic distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{a (\frac{x}{b})^a}{x(1 + (\frac{x}{b})^a)^2}}{p(x)
    = a (x/b)^a /(x(1 + (x/b)^a)^2)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}.
}
\value{
\code{dllogis} gives the density, \code{pllogis} gives the distribution function,
\code{qllogis} gives the quantile function, \code{rllogis} generates random
deviates, \code{mllogis} gives the kth raw moment, and \code{levllogis} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\note{
Also known as Fisk distribution.
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
exp(dllogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pllogis(qllogis(p, 2, 3), 2, 3)
mllogis(2, 2, 3) - mllogis(1, 2, 3) ^ 2
levllogis(10, 2, 3, 2)
}
\keyword{distribution}
