\name{rearrangepf}
\alias{rearrangepf}
\title{Reorganization of a portfolio}
\description{
  \code{rearrangepf} reorganizes a portfolio of claim amounts stored in a two
  dimension list.  
}
\usage{
rearrangepf(pf)
}
\arguments{
  \item{pf}{a two dimension list containing individual claim amounts}
}
\value{
  A list made of:
  \item{aggregate}{matrix of aggregate losses for each contract and
    each year}
  \item{frequencies}{matrix of the number of claims per contract and
    per year}
  \item{severities}{a list of two elements:
    \describe{
      \item{claims}{matrix of the individual claim amounts for the first
	n - 1 years of observations}
      \item{claims.last}{matrix of the individual claim amounts for the
	last year of observations}}
  }
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Sébastien Auclair
}
\seealso{\code{\link{simpf}} to create a portfolio of data}
\examples{
modelfreq <- list(dist1 = "pois",
                  par1 = list(lambda = quote(Lambda * weights)),
                  dist2 = "gamma",
                  par2 = c(shape = 2, rate = 1))
modelsev<-list(dist1 = "lnorm",
               par1 = list(meanlog = quote(Theta), sdlog = 1),
               dist2 = "norm",
               par2 = c(mean = 5, sd = 1))
( x <- simpf(25, 4, modelfreq, modelsev) )
rearrangepf(x$data)
}
\keyword{file}
