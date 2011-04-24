\name{quantile.aggregateDist}
\alias{quantile.aggregateDist}
\alias{VaR.aggregateDist}
\title{Quantiles of Aggregate Claim Amount Distribution}
\description{
  Quantile and Value-at-Risk methods for objects of class
  \code{"aggregateDist"}.
}
\usage{
\method{quantile}{aggregateDist}(x, 
         probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
         smooth = FALSE, names = TRUE, \dots)

\method{VaR}{aggregateDist}(x, conf.level = c(0.9, 0.95, 0.99),
         smooth = FALSE, names = TRUE, \dots)
}  
\arguments{
  \item{x}{an object of class \code{"aggregateDist"}.}
  \item{probs, conf.level}{numeric vector of probabilities with values
    in \eqn{[0, 1)}.}
  \item{smooth}{logical; when \code{TRUE} and \code{x} is a step
    function, quantiles are linearly interpolated between knots.}
  \item{names}{logical; if true, the result has a \code{names}
    attribute.  Set to \code{FALSE} for speedup with many \code{probs}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The quantiles are taken directly from the cumulative distribution
  function defined in \code{x}. Linear interpolation is available for
  step functions.
}
\value{
  A numeric vector, named if \code{names} is \code{TRUE}.
}
\seealso{
  \code{\link{aggregateDist}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
  Louis-Philippe Pouliot
}
\examples{
model.freq <- expression(data = rpois(3))
model.sev <- expression(data = rlnorm(10, 1.5))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
quantile(Fs, probs = c(0.25, 0.5, 0.75))
VaR(Fs)
}
\keyword{univar}
