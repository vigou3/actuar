\name{CTE}
\alias{CTE}
\alias{CTE.aggregateDist}
\title{Conditional Tail Expectation}
\description{
  Conditional Tail Expectation (CTE) for objects of class \code{"aggregateDist"}.
}
\usage{
\method{CTE}{aggregateDist}(x, conf.level = 0.99, \dots)
}  
\arguments{
  \item{x}{an object of class \code{"aggregateDist"}.}
  \item{conf.level}{level of confidence, with values from 0 to 1.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The Conditional Tail Expectation measures the mean loss given the loss
  is greater than the Value at Risk. Mathematically, it can be described
  as \eqn{E[x|x>VaR(x)]} for a given level of confidence.
}
\value{
  An object of class \code{numeric}.
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
CTE(Fs, conf.level = 0.95)
}
\keyword{univar}
