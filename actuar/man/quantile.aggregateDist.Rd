\name{quantile.aggregateDist}
\alias{quantile.aggregateDist}
\title{Quantiles of Aggregate Claim Amount Distribution}
\description{
  Quantile method for objects of class \code{"aggregateDist"}.
}
\usage{
\method{quantile}{aggregateDist}(x, 
         probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
         approx.lin = FALSE, names = TRUE, \dots)
}  
\arguments{
  \item{x}{an \code{\link{aggregateDist}} object.}
  \item{approx.lin}{logical; when \code{TRUE} and \code{x} is a step
    function, quantiles are linearly interpolated between knots.}
  \item{probs}{numeric vector of probabilities with values in [0, 1)}
  \item{names}{logical; if true, the result has a 'names' attribute.  Set to
    'FALSE' for speedup with many 'probs'.}
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
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
  Louis-Philippe Pouliot
}
\examples{
model.freq <- list(dist1 = "pois", par1 = c(lambda = 4))
model.sev <- list(dist1 = "lnorm", par1 = c(meanlog = 10, sdlog = 1.5))
Fs <- aggregateDist("simulation", model.freq = model.freq,
                    model.sev = model.sev, nb.simul = 1000)
quantile(Fs, probs = c(0.25, 0.5, 0.75))
}
\keyword{internal}
