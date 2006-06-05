
\name{Invtrgamma}
\alias{invtrgamma}
\alias{print.invtrgamma}
\alias{dinvtrgamma}
\alias{pinvtrgamma}
\alias{qinvtrgamma}
\alias{rinvtrgamma}
\title{The Inverse Transformed Gamma Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Transformed Gamma distribution with parameters 'shape1',
  'scale' and 'shape2'.
  
}
\usage{
dinvtrgamma(x, shape1, scale, shape2, log = FALSE)
pinvtrgamma(q, shape1, scale, shape2, lower.tail = TRUE,
         log.p = FALSE)
qinvtrgamma(p, shape1, scale, shape2, lower.tail = TRUE,
         log.p = FALSE)
rinvtrgamma(n, shape1, scale, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  The Inverse Transformed Gamma distribution with parameters 'shape1' = a,
  'scale' = c and shape2 = b has density:
  \deqn{p(x) = \frac{b (\frac{c}{x})^b e^{(\frac{-c}{x})^b}}{x
      \Gamma(a)}}{p(x) = b (c/x)^b exp((-c/x)^b) / x Gamma(x)}
  for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's
  'gamma()' and defined in its help.

  If 'shape2' = 1, it is an Inverse Gamma distribution. If 'shape2' = 'shape1' =
  1. it is an Inverse Exponential distribution.

  If 'shape1' = 1, it is an Inverse Weibull distribution.
}
\value{
'dinvtrgamma' gives the density, 'pinvtrgamma' gives the distribution function,
'qinvtrgamma' gives the quantile function, and 'rinvtrgamma' generates random
deviates.

Invalid arguments will results in return warning.
}
\note{
Also known as inverse generalized gamma distribution.
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
exp(dinvtrgamma(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pinvtrgamma(qinvtrgamma(p, 2, 3, 4), 2, 3, 4)
}
\keyword{distribution}
