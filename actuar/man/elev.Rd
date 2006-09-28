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
\method{plot}{elev}(x, \dots, main = NULL, xlab = "x",
   ylab = "Empirical LEV")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a vector or an object of class \code{"grouped.data"}.}
  \item{digits}{number of significant digits to use, see
    \code{\link{print}}.}
  \item{Fn, object}{an \R object inheriting from \code{"ogive"}.}
  \item{main}{main title.}
  \item{xlab, ylab}{labels of x and y axis.}
  \item{\dots}{arguments to be passed to subsequent methods.}
}
\details{
  The formula of the ELEV is...
}
\value{
  For \code{elev}, a function of class \code{"elev"}, inheriting from the
  \code{"\link{function}"} class.
}
\references{
    Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{stepfun}} for useful documentation (although the empirical
  LEV is not a step function).
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

lev(knots(lev))       # empirical lev at boudaries
lev(c(80, 200, 2000)) # and at other limits

plot(lev, type = "o", pch = 16)
}
\keyword{dplot}
\keyword{hplot}
