
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
  for the Inverse Transformed Gamma distribution with parameters \code{shape1},
  \code{scale} and \code{shape2}. 
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
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}  
}
\details{
  The Inverse Transformed Gamma distribution with parameters \code{shape1} = a,
  \code{scale} = b and \code{shape2} = c has density:
  \deqn{p(x) = \frac{c (\frac{b}{x})^c e^{(\frac{-b}{x})^c}}{x
      \Gamma(a)}}{p(x) = c (b/x)^c exp((-b/x)^c) / x Gamma(x)}
  for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's
  \code{gamma()} and defined in its help.

  The Inverse Transformed Gamma distribution defines a family of
  distributions with the following special cases:
  \itemize{
    \item \code{shape2} == 1 is an Inverse Gamma distribution,
    \item \code{shape2} == \code{shape1} == 1 is an Inverse Exponential
    distribution, and
    \item \code{shape1} == 1 is an Inverse Weibull distribution.
  }
}
\value{
\code{dinvtrgamma} gives the density, \code{pinvtrgamma} gives the distribution function,
\code{qinvtrgamma} gives the quantile function, and \code{rinvtrgamma} generates random
deviates.

Invalid arguments will return warning.
}
\note{
Also known as the inverse generalized gamma distribution.
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
