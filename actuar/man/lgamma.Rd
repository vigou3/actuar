
\name{Lgamma}
\alias{lgamma}
\alias{print.lgamma}
\alias{dlgamma}
\alias{plgamma}
\alias{qlgamma}
\alias{rlgamma}
\title{The Loggamma Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Loggamma distribution with parameters \code{shapelog} and \code{ratelog}.
}
\usage{
  dlgamma(x, shapelog, ratelog, log = FALSE)
  plgamma(q, shapelog, ratelog, lower.tail = TRUE,
           log.p = FALSE)
  qlgamma(p, shapelog, ratelog, lower.tail = TRUE,
           log.p = FALSE)
  rlgamma(n, shapelog, ratelog)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shapelog, ratelog}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
}
\details{
  The Loggamma distribution with parameters \code{shapelog} = a
  and \code{ratelog} = b has density:
  \deqn{p(x) = dgamma(\log(x), a,\frac{1}{b})}{p(x) = dgamma(log(x),
    shapelog, 1 / ratelog)}
  for \eqn{x > 1, a > 0} and \eqn{b > 0}. 

  Let X be a Gamma random variable with parameters \code{shapelog} and
  \code{ratelog}. Then, \deqn{Y = e^X}{Y = exp(X)} has the Loggamma
  probability density function with parameters \code{shapelog} and \code{ratelog}. 
}
\value{
\code{dlgamma} gives the density, \code{plgamma} gives the distribution function,
\code{qlgamma} gives the quantile function, and \code{rlgamma} generates random
deviates.

Invalid arguments will return warning.
}
\references{
  Appendix of Hogg & Klugman, Loss Distributions, Wiley, 1984.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dlgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
plgamma(qlgamma(p, 2, 3), 2, 3)
}
\keyword{distribution}
