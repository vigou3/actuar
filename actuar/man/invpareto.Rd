
\name{Invpareto}
\alias{invpareto}
\alias{print.invpareto}
\alias{dinvpareto}
\alias{pinvpareto}
\alias{qinvpareto}
\alias{rinvpareto}
\alias{minvpareto}
\title{The Inverse Pareto Distribution}
\description{
  Density, distribution function, quantile function, random generation
  raw moments and limited moments for the Inverse Pareto distribution with parameters \code{shape} and \code{scale}.
}
\usage{
  dinvpareto(x, shape, scale, log = FALSE)
  pinvpareto(q, shape, scale, lower.tail = TRUE,
           log.p = FALSE)
  qinvpareto(p, shape, scale, lower.tail = TRUE,
           log.p = FALSE)
  rinvpareto(n, shape, scale)
  minvpareto(k, shape, scale, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k}{order of the moment.}
}
\details{
  The Inverse Pareto distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{a b x^{(a - 1)}}{(x + b)^{(a + 1)}}}{p(x) = a b x^(a - 1) /
    (x + b)^(a + 1)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
\code{dinvpareto} gives the density, \code{pinvpareto} gives the distribution function,
\code{qinvpareto} gives the quantile function, \code{rinvpareto} generates random
deviates, and \code{minvpareto} gives the kth raw moment.

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
exp(dinvpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvpareto(qinvpareto(p, 2, 3), 2, 3)
minvpareto(-1, 2, 2)
}
\keyword{distribution}
