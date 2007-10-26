\name{VaR}
\alias{VaR}
\title{Value at Risk}
\description{
  Value at Risk.
}
\usage{
VaR(x, \dots)
}  
\arguments{
  \item{x}{an \R object.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  This is a generic function with, currently, only a method for objects
  of class \code{"aggregateDist"}.
}
\value{
  An object of class \code{numeric}.
}
\seealso{
  \code{\link{VaR.aggregateDist}}, \code{\link{aggregateDist}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
  Tommy Ouellet
}
\keyword{univar}
