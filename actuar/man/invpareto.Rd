
\name{Invpareto}
\alias{invpareto}
\alias{print.invpareto}
\alias{dinvpareto}
\alias{pinvpareto}
\alias{qinvpareto}
\alias{rinvpareto}
\title{The Inverse Pareto Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Pareto distribution with parameters 'shape' and 'scale'.
  
}
\usage{
dinvpareto(x, shape, scale, log = FALSE)
pinvpareto(q, shape, scale, lower.tail = TRUE,
         log.p = FALSE)
qinvpareto(p, shape, scale, lower.tail = TRUE,
         log.p = FALSE)
rinvpareto(n, shape, scale)
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
  The Inverse Pareto distribution with parameters 'shape' = a
  and 'scale' = b has density:
  \deqn{p(x) = \frac{a b x^{(a - 1)}}{(x + b)^(a + 1)}}{p(x) = a b x^(a - 1) /
    (x + b)^(a + 1)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
'dinvpareto' gives the density, 'pinvpareto' gives the distribution function,
'qinvpareto' gives the quantile function, and 'rinvpareto' generates random
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
exp(dinvpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvpareto(qinvpareto(p, 2, 3), 2, 3)
}
\keyword{distribution}
