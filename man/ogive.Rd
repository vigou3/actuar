\name{ogive}
\alias{ogive}
\alias{print.ogive}
\alias{summary.ogive}
\alias{knots.ogive}
\alias{plot.ogive}
\title{Ogive for Grouped Data}
\description{
  Compute a smoothed empirical distribution function for grouped data.
}
\usage{
ogive(x, y = NULL, \dots)

\method{print}{ogive}(x, digits = getOption("digits") - 2, \dots)

\method{summary}{ogive}(object, \dots)

\method{knots}{ogive}(Fn, \dots)

\method{plot}{ogive}(x, main = NULL, xlab = "x", ylab = "F(x)", \dots)
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
  The ogive is a linear interpolation of the empirical cumulative
  distribution function.

  The equation of the ogive is
  \deqn{G_n(x) = \frac{(c_j - x) F_n(c_{j - 1}) +
      (x - c_{j - 1}) F_n(c_j)}{c_j - c_{j - 1}}}{%
    Gn(x) = ((c[j] - x) Fn(c[j-1]) +
    (x - c[j-1]) Fn(c[j]))/(c[j] - c[j-1])}
  for \eqn{c_{j-1} < x \leq c_j}{c[j-1] < x <= c[j]} and where
  \eqn{c_0, \dots, c_r}{c[0], \dots, c[r]} are the \eqn{r + 1} group
  boundaries and \eqn{F_n}{Fn} is the empirical distribution function of
  the sample.
}
\value{
  For \code{ogive}, a function of class \code{"ogive"}, inheriting from the
  \code{"\link{function}"} class.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{quantile.grouped.data}} for the inverse function;
  \code{\link{approxfun}}, which is used to compute the ogive;
  \code{\link{stepfun}} for related documentation (even though the ogive
  is not a step function).
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
