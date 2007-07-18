\name{VaR}
\alias{VaR}
\alias{VaR.aggregateDist}
\title{Value at Risk}
\description{
  Value at Risk (VaR) for objects of class \code{"aggregateDist"}.
}
\usage{
\method{VaR}{aggregateDist}(x, conf.level = 0.99, \dots)
}  
\arguments{
  \item{x}{an object of class \code{"aggregateDist"}.}
  \item{conf.level}{level of confidence, with values from 0 to 1.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The Value at Risk is taken directly from the quantile metehod for
  \code{"aggregateDist"} objects, which uses the cumulative distribution
  function defined in \code{x}. Linear interpolation is available for
  step functions.
}
\value{
  An object of class \code{numeric}.
}
\seealso{
  \code{\link{aggregateDist}}; \code{\link{quantile.aggregateDist}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
  Tommy Ouellet
}
\examples{
model.freq <- expression(data = rpois(7))
model.sev <- expression(data = rnorm(9, 2))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
VaR(Fs, conf.level = 0.95)
}
\keyword{univar}
