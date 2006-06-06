
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
  for the Inverse Gamma distribution with parameters \code{shape} and \code{scale}.
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
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
}
\details{
  The Inverse Gamma distribution with parameters \code{shape} = a
  and \code{scale} = b has density:
  \deqn{p(x) = \frac{(\frac{b}{x})^a e^{\frac{-b}{x}}}{x
      \Gamma(a)}}{p(x) = (b/x)^a exp(-b/x) / (x gamma(a))}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}. Here Gamma(a) is the function implemented by R's
  \code{gamma()} and defined in its help.

  The Inverse Gamma distribution has the following special case:
  \itemize{
    \item \code{shape1} == 1 is an Inverse Exponential distribution.
  }
}
\value{
\code{dinvgamma} gives the density, \code{pinvgamma} gives the distribution function,
\code{qinvgamma} gives the quantile function, and \code{rinvgamma} generates random
deviates.

Invalid arguments will return warning.
}
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
exp(dinvgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvgamma(qinvgamma(p, 2, 3), 2, 3)
}
\keyword{distribution}
