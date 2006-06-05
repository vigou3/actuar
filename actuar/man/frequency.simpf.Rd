\name{frequency.simpf}
\alias{frequency.simpf}
\title{Compute the Frequency of Claims in Portfolios }
\description{
  \code{'frequency'} method for class \code{"simpf"}. Computes the
  frequency of claims for subsets of portfolio data. 
}
\usage{
\method{frequency}{simpf}(x, by.contract = TRUE, by.year = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{simpf} object, created by \code{\link{simpf}}}
  \item{by.contract}{boolean, whether to aggregate data by contract or not}
  \item{by.year}{boolean, whether to aggregate data by year or not}
  \item{\dots}{potential further arguments required by generic.}

}
\details{
  This function is a method for the generic function frequency() of
  package \code{stats}. It computes the frequency of the aggregate
  losses in a simulated portfolio, created by function simpf(). 
}
\value{
  The dimensions of the matrix returned depend on the way data is
  aggregated with respect to the call. The same results are obtained
  using \code{length} as a \code{FUN} argument  in function
  \code{\link{aggregate}}. See help on this function for precisions
  concerning the effects of arguments \code{by.contract} and \code{by.year}.  
     
}
\seealso{
  \code{\link{simpf}} \code{\link{aggregate}}
  \code{\link[stats]{frequency}}
}

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
modelfreq <- list(dist1 = "pois",
                  par1 = list(lambda = quote(Lambda * weights)),
                  dist2 = "gamma",
                  par2 = c(shape = 2, rate = 1))
modelsev<-list(dist1 = "lnorm",
               par1 = list(meanlog = quote(Theta), sdlog = 1),
               dist2 = "norm",
               par2 = c(mean = 5, sd = 1))
data(hachemeister)
weights <- hachemeister$weights/mean(hachemeister$weights)
x <- simpf(5, 12, modelfreq, modelsev, weights)
frequency(x, by.contract = FALSE)
}

\keyword{models}
