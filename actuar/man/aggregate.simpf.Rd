\name{aggregate.simpf}
\alias{aggregate.simpf}
\alias{frequency.simpf}
\title{Compute Summary Statistics of a Portfolio}
\description{
  \code{aggregate} method for class \code{"simpf"}. Splits portfolio
  data into subsets, computes summary statistics for each.

  \code{frequency} method for class \code{"simpf"}. Computes the
  frequency of claims for subsets of portfolio data.
}
\usage{
\method{aggregate}{simpf}(x, by = c("contract", "year"),
          FUN = sum, \dots)

\method{frequency}{simpf}(x, by = c("contract", "year"), \dots)
}
\arguments{
  \item{x}{a \code{simpf} object, typically created by \code{\link{simpf}}}
  \item{by}{vector, indicates how data should be aggregated}
  \item{FUN}{the function to be applied to data subsets}
  \item{\dots}{potential further arguments required by generic.}

}
\details{
  \code{aggregate.simpf} and \code{frequency.simpf} are respectively
  methods for functions \code{aggregate} and \code{frequency}, both from
  package \code{stats}. Both methods offer the possibility to aggregate
  the claims with respect to the contracts, the years, both or neither
  through argument \code{by}.
  Using \code{frequency.simpf} is equivalent to using
  \code{aggregate.simpf} with argument \code{FUN} equal to
  \code{length}.   
}
\value{
  The object returned depends on the way data is aggregated with respect
  to the call:

  If data is aggregated by both contract and year, a matrix will
  be returned with FUN applied to all losses vectors in all years and in
  all contracts. Therefore, the dimensions of this resulting matrix will
  be the same as the dimensions of the data matrix listed in the
  \code{simpf} object, as long as order is respected in the
  \sQuote{by} vector. Otherwise, the resulting matrix will be transposed.   

  If the argument to \code{by} is of length 1, function \code{FUN} will
  be applied to either each contract for all years or to each year for
  all contracts. Therefore, the resulting object will be a vector of
  length equal to either the number of years or the number of contracts.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\seealso{
  \code{\link{aggregate}} \code{\link{simpf}} \code{\link{severity}}
  \code{\link{frequency}} \code{\link{frequency.simpf}}
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
aggregate(x, by = "contract", FUN = mean)
frequency(x)
}
\keyword{models}
