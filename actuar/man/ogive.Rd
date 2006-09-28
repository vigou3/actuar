\name{ogive}
\alias{ogive}
\alias{Ogive}
\alias{print.ogive}
\alias{summary.ogive}
\alias{knots.ogive}
\alias{plot.ogive}
\title{Ogive for grouped data}
\description{
  Compute a smoothed empirical distribution function for grouped data.
}
\usage{
ogive(x, y = NULL, \dots)
\method{print}{ogive}(x, digits = getOption("digits") - 2, \dots)
\method{summary}{ogive}(object, \dots)
\method{knots}{ogive}(Fn, \dots)
\method{plot}{ogive}(x, \dots, main = NULL, xlab = "x", ylab = "F(x)")
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"} or a vector of group
    boundaries in \code{ogive}; for the methods, an object of class
    \code{"ogive"}, typically.}
  \item{y}{a vector of group frequencies; used only if \code{x} does not
    inherit from class \code{"grouped.data"}.}
  \item{digits}{number of significant digits to use, see
    \code{\link{print}}.}
  \item{Fn, object}{an \R object inheriting from \code{"ogive"}.}
  \item{main}{main title.}
  \item{xlab, ylab}{labels of x and y axis.}
  \item{\dots}{arguments to be passed to subsequent methods.}
}
\details{
  The ogive of a group data set links the values of the empirical
  cumulative distribution known at group boundaries by straight line
  segments, resulting in an approximation of the empirical cdf.
  
  The equation of the ogive is
  \deqn{F_n(x) = \frac{(c_j - x) F_n(c_{j-1}) + (x - c_{j-1})
      F_n(c_j)}{c_j - c_{j - 1}}
    {Fn(x) = ((c_j - x) Fn(c_{j-1})) + (x - c_{j-1})
      Fn(c_j))/(c_j - c_{j-1})}
  for \eqn{c_{j-1} < x \leq c_j}.
}
\value{
  For \code{ogive}, a function of class \code{"ogive"}, inheriting from the
  \code{"\link{function}"} class.
}
\references{
  Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{approxfun}}, which is used to compute the ogive;
  \code{\link{stepfun}} for useful documentation (although the ogive is
  not a step function).
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
data(gdental)
Fn <- ogive(gdental)
Fn
summary(Fn)
knots(Fn)            # the group boundaries

Fn(knots(Fn))        # true values of the empirical cdf
Fn(c(80, 200, 2000)) # linear interpolations

plot(Fn)
}
\keyword{dplot}
\keyword{hplot}
