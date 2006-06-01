
\name{Invweibull}
\alias{invweibull}
\alias{print.invweibull}
\alias{dinvweibull}
\alias{pinvweibull}
\alias{qinvweibull}
\alias{rinvweibull}
\alias{lgompertz}
\alias{print.lgompertz}
\alias{dlgompertz}
\alias{plgompertz}
\alias{qlgompertz}
\alias{rlgompertz}
\title{The Inverse Weibull Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Weibull distribution with parameters 'scale' and 'shape'.
  
}
\usage{
dinvweibull(x, scale, shape, log = FALSE)
pinvweibull(q, scale, shape, lower.tail = TRUE,
         log.p = FALSE)
qinvweibull(p, scale, shape, lower.tail = TRUE,
         log.p = FALSE)
rinvweibull(n, scale, shape)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{scale, shape}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}

  If 'shape' = 1, it is an Inverse Exponential distribution.
  
}
\details{
  The Inverse Weibull distribution with parameters 'scale' = a
  and 'shape' = b has density:
  \deqn{p(x) = \frac{b (\frac{a}{x})^b e^{(\frac{-a}{x})^b}}{x}}{p(x) = b
  (a/x)^b exp((-a/x)^b) / x}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. 
}
\value{
'dinvweibull' gives the density, 'pinvweibull' gives the distribution function,
'qinvweibull' gives the quantile function, and 'rinvweibull' generates random
deviates.

Invalid arguments will results in return warning.
}
\note{
Also knonw as log-Gompertz distribution.
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
exp(dinvweibull(2,3,4,log = TRUE))
p <- (1:10)/10
pinvweibull(qinvweibull(p,2,3),2,3)
}
\keyword{distribution}
