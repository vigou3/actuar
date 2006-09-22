\name{coverage}
\alias{coverage}
\alias{Coverage}
\title{Density and Cumulative Distribution Function for Modified Data}
\description{
  Compute probability density or cumulative distribution function of the
  payment per payment or payment per loss random variable under any
  combination of the following coverage modifications: deductible, limit,
  coinsurance, inflation.
}
\usage{
coverage(dist, deductible = 0, franchise = FALSE,
         limit = Inf, coinsurance = 1, inflation = 0,
         per.loss = FALSE, cdf = FALSE)
}
\arguments{
  \item{dist}{character string of the name of a probability law.}
  \item{deductible}{a unique positive value.}
  \item{franchise}{logical; \code{TRUE} for a franchise deductible,
    \code{FALSE} (default) for an ordinary deductible.}
  \item{limit}{a unique positive value larger than \code{deductible}.}
  \item{coinsurance}{a unique value between 0 and 1; the proportion of
    coinsurance.}
  \item{inflation}{a unique value between 0 and 1; the rate of
    inflation.}
  \item{per.loss}{logical; \code{TRUE} for the per loss distribution,
    \code{FALSE} (default) for the per payment distribution.}
  \item{cdf}{logical; \code{TRUE} for the cumulative distribution
    function, \code{FALSE} (default) for the probability distribution
    function.}
}
\details{
  \code{coverage()} returns the pdf or cdf of distribution \code{dist}
  under coverage modifications. Functions \code{pdist} and, if \code{cdf
    = FALSE}, \code{ddist} must exist.

  See \code{vignette("coverage")} for the exact definitions of the
  per payment and per loss random variables under an ordinary or
  franchise deductible.
}
\value{
  An object of mode \code{"function"} with the same arguments as
  \code{pdist} or \code{ddist}, except \code{"lower.tail"},
  \code{"log.p"} or \code{"log"}, which are not supported.
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
## Default case: pdf of the per payment random variable with
## an ordinary deductible
coverage("gamma", deductible = 1)

## Add a limit
f <- coverage("gamma", deductible = 1, limit = 7)
f(0, shape = 3, rate = 1)
f(2, shape = 3, rate = 1)
f(6, shape = 3, rate = 1)
f(8, shape = 3, rate = 1)
curve(dgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 0.3))    # original
curve(f(x, 3, 1), xlim = c(0.01, 5.99), lty = 2, add = TRUE) # modified
points(6, f(6, 3, 1), type = "h", pch = 16)
points(6, f(6, 3, 1), pch = 16)

## With no deductible, all distributions below are identical
coverage("weibull", limit = 5)
coverage("weibull", per.loss = TRUE, limit = 5)
coverage("weibull", franchise = TRUE, limit = 5)
coverage("weibull", per.loss = TRUE, franchise = TRUE, limit = 5)
}
\keyword{models}
