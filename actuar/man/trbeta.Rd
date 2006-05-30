
\name{Trbeta}
\alias{trbeta}
\alias{print.trbeta}
\alias{dtrbeta}
\alias{ptrbeta}
\alias{qtrbeta}
\alias{rtrbeta}
\title{The Transformed Beta Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Transformed Beta distribution with parameters 'shape1',
  'scale', 'shape2' and 'shape3'.
  
}
\usage{
dtrbeta(x, shape1, scale, shape2, shape3, log = FALSE)
ptrbeta(q, shape1, scale, shape2, shape3, lower.tail = TRUE,
        log.p = FALSE)
qtrbeta(p, shape1, scale, shape2, shape3, lower.tail = TRUE,
        log.p = FALSE)
rtrbeta(n, shape1, scale, shape2, shape3)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{shape1, scale, shape2, shape3}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  
  The Transformed Beta distribution with parameters 'shape1' = a,
  'scale' = c, 'shape2' = b and 'shape3' = d has density:
  \deqn{p(x) = \frac{\Gamma(a + d)}{\Gamma(a)\Gamma(d)}
    \frac{b(frac{x}{c})^{(b d)}}{x(1 + (\frac{x}{c})^b)^{(a + d)}}}{%
    (gamma(a + d) / gamma(a)*gamma(d))((b (x / c)^(a b)) / (x
    (1 + (x / c)^b)^(a + d)))} for \eqn{x > 0, a > 0, b > 0, c > 0} and
  \eqn{d > 0}. Here Gamma(a) is the function implemented by R's
  'gamma()' and defined in its help.

  Let Y has the beta probability density function with parameters
  'shape1' and shape3'.  \deqn{\text{scale}(\frac{Y}{1 - Y})^(\frac{1}{ text{shape2}})}{%
    scale (Y / (1 - Y))^(1 / shape2)} has the
  Transformed Beta probability density function with parameters
  'shape1', 'scale', 'shape2' and 'shape3'.

  If 'shape3' = 1, it is a Burr distribution. If 'shape3' = 'shape2' =
  1, it is a Pareto distribution. If 'shape3' = 'shape1' = 1, it is a
  Loglogistic distribution. If 'shape3' = 1 and 'shape1' = 'shape2', it
  is a Paralogistic distribution.

  If 'shape2' = 1, it is a Generalized Pareto distribution. If 'shape2'
  = 'shape1' = 1, it is an Inverse Pareto distribution.

  If 'shape1' = 1, it is an Inverse Burr distribution. If 'shape1' = 1
  and 'shape2' = 'shape3', it is an Inverse Paralogistic distribution.
}
\value{
'dtrbeta' gives the density, 'ptrbeta' gives the distribution function,
'qtrbeta' gives the quantile function, and 'rtrbeta' generates random
deviates.

Invalid arguments will results in return warning.
 
}
\note{
  Also known as generalized beta of the second kind and Pearson Type VI.
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
exp(dtrbeta(2, 2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrbeta(qtrbeta(p, 2, 3, 4, 5), 2, 3, 4, 5)
qtrbeta(0.3, 2, 3, 4, 5, lower.tail = FALSE)
}
\keyword{distribution}
