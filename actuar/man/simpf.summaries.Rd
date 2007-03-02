\name{simpf.summaries}
\alias{aggregate.simpf}
\alias{frequency.simpf}
\alias{severity.simpf}
\alias{weights.simpf}
\title{Summary Statistics of a Portfolio}
\description{
  Methods for \link[base]{class} \code{"simpf"} objects.

  \code{aggregate} splits portfolio data into subsets and computes
  summary statistics for each.

  \code{frequency} computes the frequency of claims for subsets of
  portfolio data.
  
  \code{severity} extracts the individual claim amounts.

  \code{weights} extracts the matrix of weights.
}
\usage{
\method{aggregate}{simpf}(x, by = names(x$nodes), FUN = sum, \dots)

\method{frequency}{simpf}(x, by = names(x$nodes), \dots)

\method{severity}{simpf}(x, by = head(names(x$node), -1), splitcol = NULL, \dots)

\method{weights}{simpf}(object, \dots)
}
\arguments{
  \item{x, object}{an object of class \code{"simpf"}, typically
    created with \code{\link{simpf}}.}
  \item{by}{character vector of grouping elements using the level names
    of the portfolio in \code{x}. The names can be abbreviated.}
  \item{FUN}{the function to be applied to data subsets.}
  \item{splitcol}{columns of the data matrix to extract separately;
    usual matrix indexing methods are supported.}
  \item{\dots}{optional arguments to \code{FUN}.}
}
\details{
  By default, \code{aggregate.simpf} computes the aggregate claim amount
  for the grouping specified in \code{by}. Any other statistic based on
  the individual claim amounts can be used through argument \code{FUN}.

  \code{frequency.simpf} is equivalent to using \code{aggregate.simpf}
  with argument \code{FUN} equal to \code{if (identical(x, NA)) NA else
    length(x)}.

  \code{severity.simpf} extracts individual claim amounts of a portfolio
  by groupings using the default method of
  \code{\link{severity}}. Argument \code{splitcol} allows to get the
  individual claim amounts of specific columns separately.

  \code{weights.simpf} extracts the weight matrix of a portfolio.
}
\value{
  A matrix or vector depending on the groupings specified in \code{by}.

  For the \code{aggregate} and \code{frequency} methods: if at least
  one level other than the last one is used for grouping, the result
  is a matrix obtained by binding the appropriate node identifiers
  extracted from \code{x$classification} and the summaries per
  grouping. If the last level is used, the column names of
  \code{x$data} are retained; otherwise, the column name is replaced
  by the deparsed name of \code{FUN}.  If only the last level is used
  (column summaries), a named vector is returned.

  For the \code{severity} method: a list of two elements:
  \item{first}{\code{NULL} or a matrix of claim amounts for the
    columns not specified in \code{splitcol} with the appropriate node
    identifiers extracted from \code{x$classification};}
  \item{last}{same as above, but for the columns specified in
    \code{splitcol}.}

  For the \code{weights} method: the weight matrix of the portfolio with
  node identifiers.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\seealso{
  \code{\link{simpf}}
}
\examples{
nodes <- list(sector = 3, unit = c(3, 4),
              employer = c(3, 4, 3, 4, 2, 3, 4), year = 5)
model.freq <- expression(sector = rexp(1),
                         unit = rexp(sector),
                         employer = rgamma(unit, 1),
                         year = rpois(employer))
model.sev <- expression(sector = rnorm(6, 0.1),
                        unit = rnorm(sector, 1),
                        employer = rnorm(unit, 1),
                        year = rlnorm(employer, 1))
pf <- simpf(nodes, model.freq, model.sev)

aggregate(pf)        # aggregate claim amount by employer and year
aggregate(pf, by = "sector")	      # by sector 
aggregate(pf, by = "y")		      # by year
aggregate(pf, by = c("s", "u"), mean) # average claim amount

frequency(pf)			      # number of claims

severity(pf)			      # claim amounts by row
severity(pf, by = "year")	      # by column
severity(pf, by = c("s", "u"))        # by unit
severity(pf, splitcol = "year.5")     # last year separate
severity(pf, splitcol = 5)            # same
severity(pf, splitcol = c(FALSE, FALSE, FALSE, FALSE, TRUE)) # same

weights(pf)
}
\keyword{models}
\keyword{methods}
