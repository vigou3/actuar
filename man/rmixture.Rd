\name{rmixture}
\alias{rmixture}
\title{Simulation from Discrete Mixtures}
\description{
  Generate random variates from a discrete mixture of distributions.
}
\usage{
rmixture(n, probs, models)
}
\arguments{
  \item{n}{number of random variates to generate. If \code{length(n) >
      1}, only the value is used.}
  \item{probs}{numeric non-negative vector specifying the probability
    for each model; is internally normalized to sum 1. Infinite
    and missing values are not allowed. Values are recycled as necessary
    to match the length of \code{models}.}
  \item{models}{vector of expressions specifying the simulation models
    with the number of variates omitted (see details). Models are
    recycled as necessary to match the length of \code{probs}.}
}
\details{
  \code{rmixture} generates variates from a discrete mixture, that is
  random variable with a probability density function of the form
  \deqn{f(x) = p_1 f_1(x) + ... + p_n f_n(x),}
  where \eqn{f_1, \dots, f_n} are densities and \eqn{\sum_{i = 1}^n p_i
    = 1}{p_1 + \dots + p_n = 1}.

  The values in \code{probs} will be internally normalized to be
  used as probabilities \eqn{p_1 + \dots + p_n}.

  The specification of simulation models uses the syntax of
  \code{\link{rcomphierarc}}. Models \eqn{f_1, \dots, f_n} are expressed in a
  semi-symbolic fashion using an object of mode
  \code{\link[base]{expression}} where each element is a complete call
  to a random number generation function, with the number of variates
  omitted.

  The argument of the random number generation functions for the number
  of variates to simulate \strong{must} be named \code{n}.
}
\value{
  A vector of random deviates from the mixture with density \eqn{f(x)}.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{rcompound}} to simulate from compound models.

  \code{\link{rcomphierarc}} to simulate from compound hierarchical models.
}
\examples{
## Mixture of two exponentials (with means 1/3 and 1/7) with equal
## probability.
rmixture(10, 0.5, expression(rexp(3), rexp(7)))
rmixture(10, 42, expression(rexp(3), rexp(7))) # same

## Mixture of two lognormals with different probabilities.
rmixture(10, probs = c(0.554, 0.446),
         models = expression(rlnorm(3.575, 0.637),
                             rlnorm(4.555, 0.265)))
}
\keyword{datagen}
