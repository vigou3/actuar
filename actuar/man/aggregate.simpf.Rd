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

  A matrix or vector depending on the groupings specified in \code{by}.

  If at least one level other than the last one is used for grouping,
  the result is a matrix obtained by binding the appropriate node
  identifiers extracted from \code{x$classification} and the summaries
  per grouping. If the last level is used, the column names of
  \code{x$data} are retained; otherwise, the column name is replaced by
  the deparsed name of \code{FUN}.

  If \emph{only} the last level is used (column summaries), a named
  vector is returned.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\seealso{
  \code{\link{simpf}} to generate data;
  \code{\link{severity.simpf}} to obtain individual claim amounts.
}
\examples{
nodes <- list(sector = 3, unit = c(3, 4),
              employer = c(3, 4, 3, 4, 2, 3, 4), year = 5)
model.freq <- expression(sector = rexp(1),
                         unit = rexp(sector),
                         employer = rgamma(unit, 1),
                         year = rpois(contract))
model.sev <- expression(sector = rnorm(6, 0.1),
                        unit = rnorm(sector, 1),
                        employer = rnorm(unit, 1),
                        year = rlnorm(contract, 1))
pf <- simpf(nodes, model.freq, model.sev)

aggregate(pf)        # aggregate claim amount by employer and year
aggregate(pf, by = "sector")	     # by sector 
aggregate(pf, by = "y")		     # by year
aggregate(pf, by = c("s", "u"), var) # variance by sector and unit
frequency(pf)			     # number of claims
}
\keyword{models}
