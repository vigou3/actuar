\name{mde}
\alias{Mde}
\alias{mde}
\title{Minimum Distance Estimation}
\description{
  Minimum distance estimators for individual and grouped data.
}
\usage{
mde(x, fun, start, measure = c("CvM", "chi-square", "LAS"),
    weights = NULL, ...)
}
\arguments{
  \item{x}{a vector or an object of class \code{"grouped data"}.}
  \item{fun}{either a character string or a function returning a
    cumulative distribution evaluated at its first argument.}
  \item{start}{a named list giving the parameters to be optimized with
    initial values}
  \item{measure}{either \code{"CvM"} for Cramer-von Mises' method,
    \code{"chi-square"} for modified chi-square's method, or \code{"LAS"}
    for layer average severity's method.}
  \item{weights}{weights. If arguments \code{weights} is missing,
    all weights are equal to 1 for \code{measure = "CvM"} and
    \code{measure = "LAS"}. For \code{measure = "chi-square"}, weights
    are equal to \code{1/nj}.}
  \item{\dots}{Additional parameters, either for \code{fun} or for
    \code{optim}. In particular, it can be used to specify bounds via
    \code{lower} or \code{upper} or both.  If arguments of \code{fun}
    are included they will be held fixed.}
}
\details{
  The Cramer-von Mises method (\code{"CvM"}) minimizes the distance
  between the theoretical cdf and the empirical cdf (for individual
  data) or the ogive (for grouped data).

  The modified chi-square method (\code{"chi-square"}) minimizes the
  chi-square statistic for grouped data, that is the difference between
  the expected and observed frequency within each group.

  The layer average severity method (\code{"LAS"}) minimizes the
  difference between the theoretical and empirical limited expected
  value within each group for grouped data.
  
  Optimization is performed using \code{\link{optim}}. For
  one-dimensional problems the Nelder-Mead method is used and for
  multi-dimensional problems the BFGS method, unless arguments named
  \code{lower} or \code{upper} are supplied when \code{L-BFGS-B} is used
  or \code{method} is supplied explicitly.
}
\value{
  An object of class \code{"mde"}, a list with two components:
  
  \item{estimate}{the parameter estimates, and}  
  \item{distance}{the distance.}
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
\dontrun{mde(t, pgamma, d, "chi-square")}
}
\keyword{distribution}
\keyword{htest}
