\name{aggregate.simpf}
\alias{aggregate.simpf}
\alias{frequency.simpf}
\title{Summary Statistics of a Portfolio}
\description{
  \code{aggregate} method for \link[base]{class} \code{"simpf"}. Splits
  portfolio data into subsets, computes summary statistics for each.

  \code{frequency} method for \link[base]{class}
  \code{"simpf"}. Computes the frequency of claims for subsets of
  portfolio data.
}
\usage{
\method{aggregate}{simpf}(x, by = names(x$nodes), FUN = sum, \dots)

\method{frequency}{simpf}(x, by = names(x$nodes))
}
\arguments{
  \item{x}{an object of class \code{"simpf"}, typically created
    with \code{\link{simpf}}.}
  \item{by}{character vector of grouping elements using the level names
    of the portfolio in \code{x}. The names can be abbreviated. Defaults
    to summaries by elements of the matrix in \code{x$data}.}
  \item{FUN}{the function to be applied to data subsets.}
  \item{\dots}{optional arguments to \code{FUN}.}
}
\details{
  \code{aggregate.simpf} and \code{frequency.simpf} are respectively
  methods for functions \code{\link[stats]{aggregate}} and
  \code{\link[stats]{frequency}}, both from package \code{stats}.

  By default, \code{aggregate.simpf} computes the aggregate claim amount
  for the grouping specified in \code{by}. Any other statistic based on
  the individual claim amounts can be used through argument \code{FUN}.

  \code{frequency.simpf} is equivalent to using \code{aggregate.simpf}
  with argument \code{FUN} equal to \code{if (identical(x, NA)) NA else
  length(x)}.
}
\value{

  A matrix obtained by binding the appropriate node identifiers
  extracted from \code{x$classification} and the summaries per
  grouping. If the name of the last level is in \code{by}, the column
  names of \code{x$data} are retained; otherwise, the column name is
  replaced by the deparsed name of \code{FUN}.

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
