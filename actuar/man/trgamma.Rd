
\name{Trgamma}
\alias{trgamma}
\alias{print.trgamma}
\alias{dtrgamma}
\alias{ptrgamma}
\alias{qtrgamma}
\alias{rtrgamma}
\title{The Transformed Gamma Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Transformed Gamma distribution with parameters 'shape1',
  'scale' and 'shape2'.
  
}
\usage{
dtrgamma(x, shape1, scale, shape2, log = FALSE)
ptrgamma(q, shape1, scale, shape2, lower.tail = TRUE,
         log.p = FALSE)
qtrgamma(p, shape1, scale, shape2, lower.tail = TRUE,
         log.p = FALSE)
rtrgamma(n, shape1, scale, shape2)
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
  The Transformed Gamma distribution with parameters 'shape1' = a,
  'scale' = c and shape2 = b has density:
  \deqn{p(x) = \frac{b (\frac{x}{c})^{(a b)} e^{-(\frac{x}{c})^b}}{x
      \Gamma(a)}}{p(x) = b (x / c)^(a b) exp(-(x / c)^b) / (x gamma(a))}
  for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's
  'gamma()' and defined in its help.

  Let Y has the Gamma probability density function with parameters
  'shape1' and 'scale'. \deqn{Y^{\frac{1}{shape2}}}{Y^(1/shape2)}
  has the Transformed Gamma probability density function with parameters
  'shape1', 'scale' and 'shape2'.

  If 'shape2' = 1, it is a Gamma distribution. If 'shape2' = 'shape1' =
  1. it is an Exponential distribution.

  If 'shape1' = 1, it is a Weibull distribution.
}
\value{
'dtrgamma' gives the density, 'ptrgamma' gives the distribution function,
'qtrgamma' gives the quantile function, and 'rtrgamma' generates random
deviates.

Invalid arguments will results in return warning.
}
\note{
Also known as generalized gamma distribution.
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
exp(dtrgamma(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrgamma(qtrgamma(p, 2, 3, 4), 2, 3, 4)
}
\keyword{distribution}
