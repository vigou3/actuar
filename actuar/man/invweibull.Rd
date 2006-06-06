
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
  for the Inverse Weibull distribution with parameters \code{scale} and \code{shape}.  
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
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{scale, shape}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
}
\details{
  The Inverse Weibull distribution with parameters \code{scale} = a
  and \code{shape} = b has density:
  \deqn{p(x) = \frac{b (\frac{a}{x})^b e^{(\frac{-a}{x})^b}}{x}}{p(x) = b
  (a/x)^b exp((-a/x)^b) / x}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}.

  The Inverse Weibull distribution has the following special case:
  \itemize{
    \item \code{shape} == 1 is an Inverse Exponential distribution.
  }
}
\value{
\code{dinvweibull} gives the density, \code{pinvweibull} gives the distribution function,
\code{qinvweibull} gives the quantile function, and \code{rinvweibull} generates random
deviates.

Invalid arguments will return warning.
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
exp(dinvweibull(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvweibull(qinvweibull(p, 2, 3), 2, 3)
}
\keyword{distribution}
