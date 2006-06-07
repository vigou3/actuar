
\name{Invgauss}
\alias{invgauss}
\alias{print.invgauss}
\alias{dinvgauss}
\alias{pinvgauss}
\alias{qinvgauss}
\alias{rinvgauss}
\title{The Inverse Gauss Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Inverse Gauss distribution with parameters \code{mean} and \code{scale}.
}
\usage{
  dinvgauss(x, mean, scale, log = FALSE)
  pinvgauss(q, mean, scale, lower.tail = TRUE,
           log.p = FALSE)
  qinvgauss(p, mean, scale, lower.tail = TRUE,
           log.p = FALSE)
  rinvgauss(n, mean, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mean, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
}
\details{
  The Inverse Gauss distribution with parameters \code{mean} = a
  and \code{scale} = b has density:
  \deqn{p(x) = (\frac{scale}{2 \pi x^3})^{0.5} e^{\frac{-scale z^2}{2
	x}}}{p(x) = (scale / (2 pi x^3))^(1/2) exp(-scale z^2 / (2 x))}
  with \eqn{z = \frac{x - mean}{mean}}{(x - mean) / mean} and for \eqn{x
    > 0, a > 0} and \eqn{b > 0}.
}
\value{
\code{dinvgauss} gives the density, \code{pinvgauss} gives the distribution function,
\code{qinvgauss} gives the quantile function, and \code{rinvgauss} generates random
deviates.

Invalid arguments will return warning.
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
exp(dinvgauss(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvgauss(qinvgauss(p, 2, 3), 2, 3)
}
\keyword{distribution}
