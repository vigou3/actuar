
\name{Invgamma}
\alias{invgamma}
\alias{print.invgamma}
\alias{dinvgamma}
\alias{pinvgamma}
\alias{qinvgamma}
\alias{rinvgamma}
\title{The Inverse Gamma Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Gamma distribution with parameters 'shape' and 'scale'.
  
}
\usage{
dinvgamma(x, shape, scale, log = FALSE)
pinvgamma(q, shape, scale, lower.tail = TRUE,
         log.p = FALSE)
qinvgamma(p, shape, scale, lower.tail = TRUE,
         log.p = FALSE)
rinvgamma(n, shape, scale)
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
  The Inverse Gamma distribution with parameters 'shape' = a
  and 'scale' = b has density:
  \deqn{p(x) = \frac{(\frac{b}{x})^a e^{\frac{-b}{x}}}{x
      \Gamma(a)}}{(b/x)^a exp(-b/x) / (x gamma(a))}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. Here Gamma(a) is the function implemented by R's
  'gamma()' and defined in its help.

  If 'shape' = 1, it is an Inverse Exponential distribution. 
}
\value{
'dinvgamma' gives the density, 'pinvgamma' gives the distribution function,
'qinvgamma' gives the quantile function, and 'rinvgamma' generates random
deviates.

Invalid arguments will results in return warning.

\note{
Also known as Vinci distribution.
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
exp(dinvgamma(2,3,4,log = TRUE))
p <- (1:10)/10
pinvgamma(qinvgamma(p,2,3),2,3)
}
\keyword{distribution}
