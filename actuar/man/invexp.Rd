\name{Invexp}
\alias{invexp}
\alias{print.invexp}
\alias{dinvexp}
\alias{pinvexp}
\alias{qinvexp}
\alias{rinvexp}
\title{The Inverse Exponential Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Exponential distribution with parameter 'scale'.
  
}
\usage{
dinvexp(x, shape1, scale, shape2, log = FALSE)
pinvexp(q, scale, lower.tail = TRUE, log.p = FALSE)
qinvexp(p, scale, lower.tail = TRUE, log.p = FALSE)
rinvexp(n, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{scale}{parameter. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  The Inverse Exponential distribution with parameter 'scale' = a has
  density: \deqn{p(x) = \frac{a e^{\frac{-a}{x}}}{x^2}}{a exp(-a / x) /
    x^2} for \eqn{x > 0} and \eqn{a > 0}.
}
\value{
'dinvexp' gives the density, 'pinvexp' gives the distribution function,
'qinvexp' gives the quantile function, and 'rinvexp' generates random
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
exp(dinvexp(2,log = TRUE))
p <- (1:10)/10
pinvexp(qinvexp(p,2),2)
}
\keyword{distribution}
