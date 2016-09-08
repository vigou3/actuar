\name{coverage}
\alias{coverage}
\alias{Coverage}
\title{Density and Cumulative Distribution Function for Modified Data}
\description{
  Compute probability density function or cumulative distribution
  function of the payment per payment or payment per loss random
  variable under any combination of the following coverage
  modifications: deductible, limit, coinsurance, inflation.
}
\usage{
coverage(pdf, cdf, deductible = 0, franchise = FALSE,
         limit = Inf, coinsurance = 1, inflation = 0,
         per.loss = FALSE)
}
\arguments{
  \item{pdf, cdf}{function object or character string naming a function
    to compute, respectively, the probability density function and cumulative
    distribution function of a probability law.}
  \item{deductible}{a unique positive numeric value.}
  \item{franchise}{logical; \code{TRUE} for a franchise deductible,
    \code{FALSE} (default) for an ordinary deductible.}
  \item{limit}{a unique positive numeric value larger than
    \code{deductible}.}
  \item{coinsurance}{a unique value between 0 and 1; the proportion of
    coinsurance.}
  \item{inflation}{a unique value between 0 and 1; the rate of
    inflation.}
  \item{per.loss}{logical; \code{TRUE} for the per loss distribution,
    \code{FALSE} (default) for the per payment distribution.}
}
\details{
  \code{coverage} returns a function to compute the probability
  density function (pdf) or the cumulative distribution function (cdf)
  of the distribution of losses under coverage modifications. The pdf
  and cdf of unmodified losses are \code{pdf} and \code{cdf},
  respectively.

  If \code{pdf} is specified, the pdf is returned; if \code{pdf} is
  missing or \code{NULL}, the cdf is returned. Note that \code{cdf} is
  needed if there is a deductible or a limit.
}
\value{
  An object of mode \code{"function"} with the same arguments as
  \code{pdf} or \code{cdf}, except \code{"lower.tail"},
  \code{"log.p"} and \code{"log"}, which are not supported.
}
\note{
  Setting arguments of the function returned by \code{coverage} using
  \code{\link{formals}} may very well not work as expected.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\seealso{
  \code{vignette("coverage")} for the exact definitions of the
  per payment and per loss random variables under an ordinary or
  franchise deductible.
}
\examples{
## Default case: pdf of the per payment random variable with
## an ordinary deductible
coverage(dgamma, pgamma, deductible = 1)

## Add a limit
f <- coverage(dgamma, pgamma, deductible = 1, limit = 7)
f <- coverage("dgamma", "pgamma", deductible = 1, limit = 7) # same
f(0, shape = 3, rate = 1)
f(2, shape = 3, rate = 1)
f(6, shape = 3, rate = 1)
f(8, shape = 3, rate = 1)
curve(dgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 0.3))    # original
curve(f(x, 3, 1), xlim = c(0.01, 5.99), col = 4, add = TRUE) # modified
points(6, f(6, 3, 1), pch = 21, bg = 4)

## Cumulative distribution function
F <- coverage(cdf = pgamma, deductible = 1, limit = 7)
F(0, shape = 3, rate = 1)
F(2, shape = 3, rate = 1)
F(6, shape = 3, rate = 1)
F(8, shape = 3, rate = 1)
curve(pgamma(x, 3, 1), xlim = c(0, 10), ylim = c(0, 1))    # original
curve(F(x, 3, 1), xlim = c(0, 5.99), col = 4, add = TRUE)  # modified
curve(F(x, 3, 1), xlim = c(6, 10), col = 4, add = TRUE)    # modified

## With no deductible, all distributions below are identical
coverage(dweibull, pweibull, limit = 5)
coverage(dweibull, pweibull, per.loss = TRUE, limit = 5)
coverage(dweibull, pweibull, franchise = TRUE, limit = 5)
coverage(dweibull, pweibull, per.loss = TRUE, franchise = TRUE,
         limit = 5)

## Coinsurance alone; only case that does not require the cdf
coverage(dgamma, coinsurance = 0.8)
}
\keyword{models}
