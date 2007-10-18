\name{mde}
\alias{Mde}
\alias{mde}
\title{Minimum Distance Estimation}
\description{
  Minimum distance fitting of univariate distributions, allowing
  parameters to be held fixed if desired.
}
\usage{
mde(x, fun, start, measure = c("CvM", "chi-square", "LAS"),
    weights = NULL, ...)
}
\arguments{
  \item{x}{a vector or an object of class \code{"grouped data"} (in
    which case only the first column of frequencies is used).}
  \item{fun}{function returning a cumulative distribution (for
    \code{measure = "CvM"} and \code{measure = "chi-square"}) or a
    limited expected value (for \code{measure = "LAS"}) evaluated at
    its first argument.}
  \item{start}{a named list giving the parameters to be optimized with
    initial values}
  \item{measure}{either \code{"CvM"} for the Cramer-von Mises method,
    \code{"chi-square"} for the modified chi-square method, or \code{"LAS"}
    for the layer average severity method.}
  \item{weights}{weights; see details.}
  \item{\dots}{Additional parameters, either for \code{fun} or for
    \code{optim}. In particular, it can be used to specify bounds via
    \code{lower} or \code{upper} or both.  If arguments of \code{fun}
    are included they will be held fixed.}
}
\details{
  The Cramer-von Mises method (\code{"CvM"}) minimizes the squared
  difference between the theoretical cdf and the empirical cdf at the
  data points (for individual data) or the ogive at the knots (for
  grouped data).

  The modified chi-square method (\code{"chi-square"}) minimizes the
  modified chi-square statistic for grouped data, that is the squared
  difference between the expected and observed frequency within each
  group.

  The layer average severity method (\code{"LAS"}) minimizes the
  squared difference between the theoretical and empirical limited
  expected value within each group for grouped data.

  All sum of squares can be weighted. If arguments \code{weights} is
  missing, weights default to 1 for \code{measure = "CvM"} and
  \code{measure = "LAS"}; for \code{measure = "chi-square"}, weights
  default to \eqn{1/n_j}{1/n[j]}, where \eqn{n_j}{n[j]} is the frequency
  in group \eqn{j = 1, \dots, r}.
  
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
## Individual data example
data(dental)
mde(dental, pexp, start = list(rate = 1/200), measure = "CvM")

## Example 2.21 of Klugman et al. (1998)
data(gdental)
mde(gdental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200), measure = "LAS")

## Two-parameter distribution example
try(mde(gdental, ppareto, start = list(shape = 3, scale = 600),
        measure = "CvM")) # no convergence

## Working in log scale often solves the problem
pparetolog <- function(x, shape, scale)
    ppareto(x, exp(shape), exp(scale))

( p <- mde(gdental, pparetolog, start = list(shape = log(3),
           scale = log(600)), measure = "CvM") )
exp(p$estimate)
}
\keyword{distribution}
\keyword{htest}
