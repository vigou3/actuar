
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
  for the Transformed Gamma distribution with parameters \code{shape1},
  \code{scale} and \code{shape2}.
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
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
}
\details{
  The Transformed Gamma distribution with parameters \code{shape1} = a,
  \code{scale} = b and \code{shape2} = c has density:
  \deqn{p(x) = \frac{c (\frac{x}{b})^{(a c)} e^{-(\frac{x}{b})^c}}{x
      \Gamma(a)}}{p(x) = c (x / b)^(a c) exp(-(x / b)^c) / (x gamma(a))}
  for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's
  \code{gamma()} and defined in its help.

  Let Y be a Gamma random variable with parameters
  \code{shape1} and \code{scale}. \deqn{Y^{\frac{1}{shape2}}}{Y^(1/shape2)}
  has the Transformed Gamma probability density function with parameters
  \code{shape1}, \code{scale} and \code{shape2}.

  The Transformed Gamma probability distribution defines a family of
  distributions with the following special cases:
  \itemize{
    \item \code{shape2} == 1 is a Gamma distribution,
    \item \code{shape2} == \code{shape1} == 1 is an Exponential
    distribution, and
    \item \code{shape1} == 1 is a Weibull distribution.
  }
}
\value{
\code{dtrgamma} gives the density, \code{ptrgamma} gives the distribution function,
\code{qtrgamma} gives the quantile function, and \code{rtrgamma} generates random
deviates.

Invalid arguments will return warning.
}
\note{
Also known as the generalized gamma distribution.
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
