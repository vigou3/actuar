\name{grouped.data}
\alias{grouped.data}
\alias{Grouped.data}
\alias{print.grouped.data}
\title{Grouped data}
\description{
   This function creates grouped data, tightly coupled collections of
     variables used as the fundamental data structure by some of R's
     modeling software.
}
\usage{
grouped.data(x, y = NULL, ...)
\method{print}{grouped.data}(x, ...)
}
\arguments{
  \item{x, y}{the \code{x} and \code{y} arguments provide the x and y coordinates for
    the plot. \code{Length(x)} must be equal to \code{length(y) + 1}. If
    \code{length(x)} = \code{length(y)}, \code{y[1]} is ignored.
    Also, \code{x} may contains an object of class \code{data.frame}}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  A grouped data is a list of two variables (cj and nj) with
  \code{length(cj)} = \code{length(nj) - 1}. First data of nj is unused. cj contains group boundaries and nj contains the number of
  data in each groups.
}
\value{
  Returns an object of \code{class} \code{"grouped.data"}. An object of class \code{"grouped.data"} is a list with the following
  components:
    \item{cj}{boundaries.}
    \item{nj}{number of data in each groups}
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\note{
  See also \code{\link{ogive}}. 

}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
cj <- c(0, 10, 20, 30, 40)
nj <- c(4, 12, 11, 2)
data <- grouped.data(cj, nj)
data
}
\keyword{hplot}
