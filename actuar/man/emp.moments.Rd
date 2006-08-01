\name{emp.moments}
\alias{empmoments}
\alias{emp.moments}
\alias{emp.moments.default}
\title{Empirical Moments}
\description{
  Function for calculating raw empirical moments for individual and grouped data.
}
\usage{
## default S3 method:
emp.moments(x, order = 1, ...)
}
\arguments{
  \item{x}{individual or grouped data.}
  \item{order}{order of the moment.}
  \item{\dots}{potential further arguments required by generic.}
}
\value{
  \code{emp.moments} gives the kth raw empirical moment of the variable.
  
  Infinite value and negative \code{order} will return warning.
}
\details{
  See also \code{emp.moments.grouped.data}.
}
\references{
  Klugman, Panjer & Willmot, Loss Models, Second Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
t<-c(1,12,3,6,5,3,0)
emp.moments(t, k = 2)
}
\keyword{models}
