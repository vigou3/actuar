\name{ogive}
\alias{ogive}
\alias{print.ogive}
\alias{summary.ogive}
\alias{knots.ogive}
\alias{plot.ogive}
\title{Ogive for Grouped Data}
\description{
  Compute a smoothed empirical distribution function from either a
  grouped.data object, a provided set of group boundaries and group
  frequencies, or from individual data using automatic or specified
  breakpoints.
}
\usage{
ogive(\dots, breaks = "Sturges", nclass = NULL, group = FALSE)

\method{print}{ogive}(x, digits = getOption("digits") - 2, \dots)

\method{summary}{ogive}(object, \dots)

\method{knots}{ogive}(Fn, \dots)

\method{plot}{ogive}(x, main = NULL, xlab = "x", ylab = "F(x)", \dots)
}
\arguments{
  \item{\dots}{in \code{ogive}: a single object inheriting from class
    \code{"grouped.data"}, or one vector of group boundaries and one
    vector of group frequencies, or a single vector of individual data.
    In the methods: arguments to be passed to subsequent methods.}
  \item{breaks, nclass, group}{arguments passed to \code{\link{hist}}
    for the grouping of individual data; see also
    \code{\link{grouped.data}}.}
  \item{x}{an object of class \code{"ogive"}.}
  \item{digits}{number of significant digits to use, see
    \code{\link{print}}.}
  \item{Fn, object}{an \R object inheriting from \code{"ogive"}.}
  \item{main}{main title.}
  \item{xlab, ylab}{labels of x and y axis.}
}
\details{
  The ogive is a linear interpolation of the empirical cumulative
  distribution function.

  The equation of the ogive is
  \deqn{G_n(x) = \frac{(c_j - x) F_n(c_{j - 1}) +
      (x - c_{j - 1}) F_n(c_j)}{c_j - c_{j - 1}}}{%
    Gn(x) = 1/(c[j] - c[j-1]) * [(c[j] - x) Fn(c[j-1])
            + (x - c[j-1]) Fn(c[j])]}
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
## Most common usage: create ogive from grouped data object.
Fn <- ogive(gdental)
Fn
summary(Fn)
knots(Fn)                      # the group boundaries

Fn(knots(Fn))                  # true values of the empirical cdf
Fn(c(80, 200, 2000))           # linear interpolations

plot(Fn)                       # graphical representation

## Alternative: create ogive directly from individual data
## without first creating a grouped data object.
ogive(dental)                  # automatic class boundaries
ogive(dental, breaks = c(0, 50, 200, 500, 1500, 2000))
}
\keyword{dplot}
\keyword{hplot}
