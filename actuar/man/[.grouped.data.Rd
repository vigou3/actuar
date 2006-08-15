\name{[.grouped.data}
\alias{[.grouped.data}
\title{Extract parts of an Object}
\description{
  Operators acting on grouped data to extract parts.
}
\usage{
\method{[}{grouped.data}(x, i, j)
}
\arguments{
  \item{x}{object from which to extract element(s).}
  \item{i, j}{name: indices specifying elements to extract}
}
\details{
     These operators are generic.  You can write methods to handle
     indexing of specific classes of objects.
}
\value{
  Returns an object of \code{class} \code{"grouped.data"} or an object
of \code{class} \code{"numeric"} (if only nj are asked). 
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\note{
  See also \code{\link{grouped.data}}. 
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
cj <- c(0, 10, 20, 30, 40)
nj <- c(4, 12, 11, 2)
data <- grouped.data(cj, nj)
data[1,]
data[1,2]
data[,2]
}
\keyword{hplot}
