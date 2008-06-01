\name{elev}
\alias{elev}
\alias{elev.default}
\alias{elev.grouped.data}
\alias{print.elev}
\alias{summary.elev}
\alias{knots.elev}
\alias{plot.elev}
\title{Empirical Limited Expected Value}
\description{
  Compute the empirical limited expected value for individual or grouped
  data.
}
\usage{
elev(x, ...)

\method{elev}{default}(x, \dots)

\method{elev}{grouped.data}(x, \dots)

\method{print}{elev}(x, digits = getOption("digits") - 2, \dots)

\method{summary}{elev}(object, \dots)

\method{knots}{elev}(Fn, \dots)

\method{plot}{elev}(x, \dots, main = NULL, xlab = "x", ylab = "Empirical LEV")
}
\arguments{
  \item{x}{a vector or an object of class \code{"grouped.data"} (in
    which case only the first column of frequencies is used); for the
    methods, an object of class \code{"elev"}, typically.}
  \item{digits}{number of significant digits to use, see
    \code{\link{print}}.}
  \item{Fn, object}{an \R object inheriting from \code{"ogive"}.}
  \item{main}{main title.}
  \item{xlab, ylab}{labels of x and y axis.}
  \item{\dots}{arguments to be passed to subsequent methods.}
}
\details{
  The limited expected value (LEV) at \eqn{u} of a random variable
  \eqn{X} is
  \eqn{E[X \wedge u] = E[\min(X, u)]}{E[X ^ u] = E[min(X, u)]}. For
  individual data \eqn{x_1, \dots, x_n}{x[1], \dots, x[n]}, the
  empirical LEV \eqn{E_n[X \wedge u]}{En[X ^ u]} is thus
  \deqn{E_n[X \wedge u] = \frac{1}{n}
    \left( \sum_{x_j < u} x_j + \sum_{x_j \geq u} u \right).}{%
    En[X ^ u] = (sum(x[j] < u; 1) + sum(x[j] >= u; u))/n.}

  Methods of \code{elev} exist for individual data or for grouped data
  created with \code{\link{grouped.data}}. The formula in this case is
  too long to show here. See the reference for details.
}
\value{
  For \code{elev}, a function of class \code{"elev"}, inheriting from the
  \code{"\link{function}"} class.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{stepfun}} for related documentation (even though the
  empirical LEV is not a step function).
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
data(gdental)
lev <- elev(gdental)
lev
summary(lev)
knots(lev)            # the group boundaries

lev(knots(lev))       # empirical lev at boundaries
lev(c(80, 200, 2000)) # and at other limits

plot(lev, type = "o", pch = 16)
}
\keyword{dplot}
\keyword{hplot}
