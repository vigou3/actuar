\name{quantile.aggregateDist}
\alias{quantile.aggregateDist}
\title{Quantiles for Objects of Class aggregateDist }
\description{
  The \code{quantile} method for class \code{aggregateDist} produces
  quantiles for the given distribution.  
}
\usage{
\method{quantile}{aggregateDist}(x, 
         probs = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995),
         approx.lin = FALSE, names = TRUE, \dots)
}  
\arguments{
  \item{x}{an \code{\link{aggregateDist}} object.}
  \item{approx.lin}{ logical; whether to approximate linearly or not
    between points. Constant approximation by default. Ignored if object
  is issued of methods other than the Normal or Normal Power approximations.}
  \item{probs}{numeric vector of probabilities with values in [0,1)}
  \item{names}{logical; if true, the result has a 'names' attribute.  Set to
    'FALSE' for speedup with many 'probs'.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  This function is a method for the generic function quantile() for class
  \code{aggregateDist}.
  The quantiles are not calculated from a sample distribution as they
  are in the default method. The task is easier since
  \code{aggregateDist} objects contain a discrete empirical
  CDF. Quantiles are consequently taken directly using the CDF, and are
  subject to a constant approximation by default. The linear
  approximation is less appropriate in that context since the CDF is
  a step function.

  If the \code{aggregateDist} object was created with the \code{normal}
  or the \code{np2} (Normal Power) method, the object is consenquently a continuous
  distribution, the result of a normal standardization. Therefore a call to
  qnorm() will be made and no approximation is required.

  \deqn{v(X) = \frac{(X - \mu_x)}{\sigma_X}}{
    v(X) = (X - mu_x)/sigma_x}
  
}
\value{
  A numeric vector, named if \code{names} is \code{TRUE}.}

\author{
	Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
	Louis-Philippe Pouliot
}
\examples{
model.freq <- list(dist = "pois", par = c(lambda = 4))
model.sev <- list(dist = "lnorm", par = c(meanlog = 10, sdlog = 1.5))
res <- aggregateDist("simulation", model.freq = model.freq, model.sev = model.sev, n = 1000)
quantile(res, probs = c(0.25, 0.5, 0.75))
}
\keyword{internal}
