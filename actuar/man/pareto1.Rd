
\name{Pareto1}
\alias{pareto1}
\alias{print.pareto1}
\alias{dpareto1}
\alias{ppareto1}
\alias{qpareto1}
\alias{rpareto1}
\title{The Single-parameter Pareto Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Single-parameter Pareto distribution with parameters \code{shape} and \code{min}.
}
\usage{
  dpareto1(x, shape, min, log = FALSE)
  ppareto1(q, shape, min, lower.tail = TRUE, log.p = FALSE)
  qpareto1(p, shape, min, lower.tail = TRUE, log.p = FALSE)
  rpareto1(n, shape, min)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, min}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
}
\details{
  The Single-parameter Pareto distribution with parameters \code{shape} = a
  and \code{min} = b has density:
  \deqn{p(x) = \frac{a b^a}{x^{(a + 1)}}}{p(x) = a b^a / x^(a + 1)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}.
  
  Although there appears to be two parameters, only \code{shape} is a true
  parameter. The value of \code{min} must be set in advance.
}
\value{
\code{dpareto1} gives the density, \code{ppareto1} gives the distribution function,
\code{qpareto1} gives the quantile function, and \code{rpareto1} generates random
deviates.

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
exp(dpareto1(5, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto1(qpareto1(p, 2, 3), 2, 3)
}
\keyword{distribution}
