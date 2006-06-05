\name{aggregate.simpf}
\alias{aggregate.simpf}
\title{Compute Summary Statistics of a Portfolio}
\description{
  \code{'aggregate'} method for class \code{"simpf"}. Splits portfolio
  data into subsets, computes summary statistics for each. Results are
  most commonly returned in a \code{matrix}. 
}
\usage{
\method{aggregate}{simpf}(x, by.contract = TRUE, by.year = TRUE, FUN =
sum, \dots)
}
\arguments{
  \item{x}{a \code{"simpf"} object, typically created by \code{\link{simpf}}}
  \item{by.contract}{boolean, whether to aggregate data by contract or not}
  \item{by.year}{boolean, whether to aggregate data by year or not}
  \item{TOL}{maximum relative error in the iterative procedure}
  \item{FUN}{a function to apply to data subsets}
  \item{\dots}{potential further arguments required by generic.}

}
\details{
  This function is a method for the generic function aggregate() of
  package \code{stats}. It computes summary statistics for of the
  aggregate losses of a simulated portfolio, created by function
  simpf(). It is theoretically possible to use any scalar function as a
  FUN argument, but only functions for summary statistics make sense in
  the case of a \code{simpf} object. 
}
\value{
  The object returned depends on the way data is aggregated with respect
  to the call:
  If both \code{by.contract} and \code{by.year} are TRUE, a matrix will
  be returned with FUN applied to all losses vectors in all years and in
  all contracts. Therefore, the dimensions of this resulting matrix will
  be the same as the dimensions of the data matrix listed in the
  \code{simpf} object.
  If \code{by.contract} only is TRUE, the summary statistic is applied
  to each contract, for all years. If \code{by.year} only is TRUE, the
  summary statistic is applied to each year, for all
  contracts. Therefore, if either \code{by.contract} emph{only} or
  \code{by.year} emph{only} is TRUE, the object returned will be
  respectively a column vector or a row vector, both of class \code{matrix}.   
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\seealso{
  \code{\link{aggregate}} \code{\link{simpf}} \code{\link{severity}}
  \code{\link{frequency}} \code{\link{frequency.simpf}}
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
aggregate(x, by.contract = FALSE, FUN = mean)
}
\keyword{models}
