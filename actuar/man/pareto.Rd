
\name{Pareto}
\alias{pareto}
\alias{print.pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\alias{mpareto}
\alias{levpareto}
\alias{pareto2}
\alias{print.pareto2}
\alias{dpareto2}
\alias{ppareto2}
\alias{qpareto2}
\alias{rpareto2}
\alias{mpareto2}
\alias{levpareto2}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Pareto distribution with parameters \code{shape} and \code{scale}.
}
\usage{
  dpareto(x, shape, scale, log = FALSE)
  ppareto(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
  qpareto(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
  rpareto(n, shape, scale)
  mpareto(k, shape, scale, log = FALSE)
  levpareto(d, shape, scale, order = 1, log = FALSE)
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
  The Pareto distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{a b^a}{(x + b)^{a + 1}}}{p(x) = a b^a / (x + b)^(a
    + 1)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
\code{dpareto} gives the density, \code{ppareto} gives the distribution function,
\code{qpareto} gives the quantile function, \code{rpareto} generates random
deviates, \code{mpareto} gives the kth raw moment, and \code{levpareto} gives
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
exp(dpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto(qpareto(p, 2, 3), 2, 3)
mpareto(1, 2, 3)
levpareto(10, 2, 3, 2)
}
\keyword{distribution}
