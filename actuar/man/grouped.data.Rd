\name{grouped.data}
\alias{grouped.data}
\title{Grouped data}
\description{
  Creation of grouped data objects, allowing for consistent
  representation and manipulation of data presented in a frequency per
  group form.
}
\usage{
grouped.data(..., row.names = NULL, check.rows = FALSE,
             check.names = TRUE)
}
\arguments{
  \item{\dots}{these arguments are either of the form \code{value} or
    \code{tag = value}. See Details.}
  \item{row.names, check.rows, check.names}{arguments identical to those
    of \code{data.frame}.}
}
\details{
  A grouped data object is a special form of data frame consisting of:
  \enumerate{
    \item one column of contiguous group boundaries;
    \item one or more columns of frequencies within each class.
  }

  The first argument to the function will be taken as the vector of
  group boundaries. This vector must be exactly one element longer than
  the other arguments, which will be taken as vectors of group
  frequencies. All arguments are coerced to data frames.

  Extraction and replacement methods exist for \code{grouped.data}
  objects, but working on non adjacent groups will most likely yield
  useless results.  
}
\value{
  An object of \code{class} \code{c("grouped.data", "data.frame")} with
  an environment containing the vector \code{cj} of group boundaries.
}
\references{
  Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
}
\seealso{
  \code{\link{[.grouped.data}} for subsetting methods,
  \code{\link{data.frame}} for usual data frame creation and
  manipulation.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Mathieu Pigeon and Louis-Philippe Pouliot
}
\examples{
## Most common usage
cj <- c(0, 25, 50, 100, 250, 500, 1000)
nj <- c(30, 31, 57, 42, 45, 10)
(x <- grouped.data(Group = cj, Frequency = nj))
class(x)

## Multiple frequency columns are supported
x <- sample(1:100, 9)
y <- sample(1:100, 9)
grouped.data(cj = 1:10, nj.1 = x, nj.2 = y)

x[, 1] # group boundaries
x[, 2] # group frequencies
}
\keyword{classes}
\keyword{methods}
