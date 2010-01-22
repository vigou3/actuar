\name{CTE}
\alias{CTE}
\alias{CTE.aggregateDist}
\title{Conditional Tail Expectation}
\description{
  Conditional Tail Expectation.
}
\usage{
CTE(x, \dots)

\method{CTE}{aggregateDist}(x, conf.level = c(0.9, 0.95, 0.99),
         names = TRUE, \dots)
}  
\arguments{
  \item{x}{an \R object.}
  \item{conf.level}{numeric vector of probabilities with values
    in \eqn{[0, 1)}}.
  \item{names}{logical; if true, the result has a \code{names}
    attribute.  Set to \code{FALSE} for speedup with many \code{probs}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The Conditional Tail Expectation measures the average of losses above
  the Value at Risk for some given confidence level, that is \eqn{E[X|X
  > VaR(X)]} where \eqn{X} is the loss random variable.

  \code{CTE} is a generic function with, currently, only a method for
  objects of class \code{"aggregateDist"}.

  For the recursive, convolution and simulation methods of
  \code{\link{aggregateDist}}, the CTE is computed from the definition
  using the empirical cdf.

  For the normal approximation method, an explicit formula exists:
  \deqn{\mu + \frac{\sigma}{(1 - \alpha)} \sqrt{2 \pi}
    e^{-\mathrm{VaR(X)}^2/2},}{%
    m + s exp(-VaR(X)^2/2)/((1 - a) * sqrt(2 pi)),}
  where \eqn{\mu}{m} is the mean, \eqn{\sigma}{s} the standard
  deviation and \eqn{\alpha}{a} the confidence level.

  For the Normal Power approximation, the CTE is computed from the
  definition using \code{\link[stats]{integrate}}.
}
\value{
  A numeric vector, named if \code{names} is \code{TRUE}.
}
\seealso{
  \code{\link{aggregateDist}}; \code{\link{VaR}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
  Tommy Ouellet
}
\examples{
model.freq <- expression(data = rpois(7))
model.sev <- expression(data = rnorm(9, 2))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
CTE(Fs)
}
\keyword{univar}
