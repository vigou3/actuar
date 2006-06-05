
\name{Pareto}
\alias{pareto}
\alias{print.pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\alias{pareto2}
\alias{print.pareto2}
\alias{dpareto2}
\alias{ppareto2}
\alias{qpareto2}
\alias{rpareto2}
\title{The Pareto Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Pareto distribution with parameters 'shape' and 'scale'.
  
}
\usage{
dpareto(x, shape, scale, log = FALSE)
ppareto(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
qpareto(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
rpareto(n, shape, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  The Pareto distribution with parameters 'shape' = a
  and 'scale' = b has density:
  \deqn{p(x) = \frac{a b^a}{(x + b)^{a + 1}}}{p(x) = a b^a / (x + b)^(a
    + 1)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
'dpareto' gives the density, 'ppareto' gives the distribution function,
'qpareto' gives the quantile function, and 'rpareto' generates random
deviates.

Invalid arguments will results in return warning.
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
}
\keyword{distribution}
