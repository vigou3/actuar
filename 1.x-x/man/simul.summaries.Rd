\name{simul.summaries}
\alias{simul.summaries}
\alias{aggregate.portfolio}
\alias{frequency.portfolio}
\alias{severity.portfolio}
\alias{weights.portfolio}
\title{Summary Statistics of a Portfolio}
\description{
  Methods for \link[base]{class} \code{"portfolio"} objects.

  \code{aggregate} splits portfolio data into subsets and computes
  summary statistics for each.

  \code{frequency} computes the frequency of claims for subsets of
  portfolio data.

  \code{severity} extracts the individual claim amounts.

  \code{weights} extracts the matrix of weights.
}
\usage{
\method{aggregate}{portfolio}(x, by = names(x$nodes), FUN = sum,
        classification = TRUE, prefix = NULL, \dots)

\method{frequency}{portfolio}(x, by = names(x$nodes),
        classification = TRUE, prefix = NULL, \dots)

\method{severity}{portfolio}(x, by = head(names(x$node), -1), splitcol = NULL,
        classification = TRUE, prefix = NULL, \dots)

\method{weights}{portfolio}(object, classification = TRUE, prefix = NULL, \dots)
}
\arguments{
  \item{x, object}{an object of class \code{"portfolio"}, typically
    created with \code{\link{simul}}.}
  \item{by}{character vector of grouping elements using the level names
    of the portfolio in \code{x}. The names can be abbreviated.}
  \item{FUN}{the function to be applied to data subsets.}
  \item{classification}{boolean; if \code{TRUE}, the node identifier
    columns are included in the output.}
  \item{prefix}{characters to prefix column names with; if \code{NULL},
        sensible defaults are used when appropriate.}
  \item{splitcol}{columns of the data matrix to extract separately;
    usual matrix indexing methods are supported.}
  \item{\dots}{optional arguments to \code{FUN}, or passed to or from
    other methods.}
}
\details{
  By default, \code{aggregate.portfolio} computes the aggregate claim amounts
  for the grouping specified in \code{by}. Any other statistic based on
  the individual claim amounts can be used through argument \code{FUN}.

  \code{frequency.portfolio} is equivalent to using \code{aggregate.portfolio}
  with argument \code{FUN} equal to \code{if (identical(x, NA)) NA else
    length(x)}.

  \code{severity.portfolio} extracts individual claim amounts of a portfolio
  by groupings using the default method of \code{\link{severity}}.
  Argument \code{splitcol} allows to get the individual claim amounts of
  specific columns separately.

  \code{weights.portfolio} extracts the weight matrix of a portfolio.
}
\value{
  A matrix or vector depending on the groupings specified in \code{by}.

  For the \code{aggregate} and \code{frequency} methods: if at least one
  level other than the last one is used for grouping, the result is a
  matrix obtained by binding the appropriate node identifiers extracted
  from \code{x$classification} if \code{classification = TRUE}, and the
  summaries per grouping. If the last level is used for grouping, the
  column names of \code{x$data} are retained; if the last level is not
  used for grouping, the column name is replaced by the deparsed name of
  \code{FUN}. If only the last level is used (column summaries), a named
  vector is returned.

  For the \code{severity} method: a list of two elements:
  \item{main}{\code{NULL} or a matrix of claim amounts for the
    columns not specified in \code{splitcol}, with the appropriate node
    identifiers extracted from \code{x$classification} if
    \code{classification = TRUE};}
  \item{split}{same as above, but for the columns specified in
    \code{splitcol}.}

  For the \code{weights} method: the weight matrix of the portfolio with
  node identifiers if \code{classification = TRUE}.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Louis-Philippe Pouliot.
}
\seealso{
  \code{\link{simul}}
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
pf <- simul(nodes, model.freq, model.sev)

aggregate(pf)            # aggregate claim amount by employer and year
aggregate(pf, classification = FALSE) # same, without node identifiers
aggregate(pf, by = "sector")	      # by sector
aggregate(pf, by = "y")		      # by year
aggregate(pf, by = c("s", "u"), mean) # average claim amount

frequency(pf)			      # number of claims
frequency(pf, prefix = "freq.")       # more explicit column names

severity(pf)			      # claim amounts by row
severity(pf, by = "year")	      # by column
severity(pf, by = c("s", "u"))        # by unit
severity(pf, splitcol = "year.5")     # last year separate
severity(pf, splitcol = 5)            # same
severity(pf, splitcol = c(FALSE, FALSE, FALSE, FALSE, TRUE)) # same

weights(pf)

## For portfolios with weights, the following computes loss ratios.
\dontrun{aggregate(pf, classif = FALSE) / weights(pf, classif = FALSE)}
}
\keyword{models}
\keyword{methods}
