\name{CTE}
\alias{CTE}
\alias{TVaR}
\alias{CTE.aggregateDist}
\title{Conditional Tail Expectation}
\description{
  Conditional Tail Expectation, also called Tail Value-at-Risk.

  \code{TVaR} is an alias for \code{CTE}.
}
\usage{
CTE(x, \dots)

\method{CTE}{aggregateDist}(x, conf.level = c(0.9, 0.95, 0.99),
         names = TRUE, \dots)

TVaR(x, \dots)
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
  The Conditional Tail Expectation (or Tail Value-at-Risk) measures the
  average of losses above the Value at Risk for some given confidence
  level, that is \eqn{E[X|X > \mathrm{VaR}(X)]} where \eqn{X} is the loss random
  variable.

  \code{CTE} is a generic function with, currently, only a method for
  objects of class \code{"aggregateDist"}.

  For the recursive, convolution and simulation methods of
  \code{\link{aggregateDist}}, the CTE is computed from the definition
  using the empirical cdf.

  For the normal approximation method, an explicit formula exists:
  \deqn{\mu + \frac{\sigma}{(1 - \alpha) \sqrt{2 \pi}}
    e^{-\mathrm{VaR}(X)^2/2},}{%
    m + s exp(-VaR(X)^2/2)/((1 - a) * sqrt(2 pi)),}
  where \eqn{\mu}{m} is the mean, \eqn{\sigma}{s} the standard
  deviation and \eqn{\alpha}{a} the confidence level.

  For the Normal Power approximation, the explicit formula given in
  Castañer et al. (2013) is
  \deqn{\mu + \frac{\sigma}{(1 - \alpha) \sqrt{2 \pi}}
    e^{-\mathrm{VaR}(X)^2/2}
    \left( 1 + \frac{\gamma}{6} \mathrm{VaR}(X) \right),}{%
    m + s exp(-VaR(X)^2/2)/((1 - a) * sqrt(2 pi)) (1 + g * VaR(X)/6),}
  where, as above, \eqn{\mu}{m} is the mean, \eqn{\sigma}{s} the standard
  deviation, \eqn{\alpha}{a} the confidence level and \eqn{\gamma}{g} is
  the skewness.

  the CTE is computed from the
  definition using \code{\link[stats]{integrate}}.
}
\value{
  A numeric vector, named if \code{names} is \code{TRUE}.
}
\seealso{
  \code{\link{aggregateDist}}; \code{\link{VaR}}
}
\references{
  Castañer, A. and Claramunt, M.M. and Mármol, M. (2013), Tail value at
  risk. An analysis with the Normal-Power approximation. In \emph{Statistical
  and Soft Computing Approaches in Insurance Problems}, pp. 87-112. Nova
  Science Publishers, 2013. ISBN 978-1-62618-506-7.
  \url{https://www.novapublishers.com/catalog/product_info.php?products_id=43105}
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
