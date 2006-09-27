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
  \item{\dots}{these arguments are of either the form \code{value} or
    \code{tag = value}. See Details.}
  \item{row.names, check.rows, check.names}{arguments identical to those
    of \code{data.frame}.}
}
\details{
  A grouped data object is a special form of data frame consisting of
  \enumerate{
    \item one column of contiguous group boundaries;
    \item one or more columns of frequencies within each class.
  }

  The first argument to the function will be taken as the vector of
  group boundaries. This vector must be exactly one element longer than
  the other arguments, which will be taken as vectors of group
  frequencies. All arguments are coerced to data frames.

  Extraction and replacement methods exist for \code{grouped.data}
  objects, but working on non adjacent groups can yield to
  unpredictable results.  
}
\value{
  An object of \code{class} \code{c("grouped.data", "data.frame")} with
  an environment containing the vector \code{cj} of group boundaries.
}
\references{Klugman, Panjer & Willmot, Loss Models, Wiley, 1998.
}
\seealso{\code{\link{data.frame}} for usual data frame creation and
  manipulation, \code{\link{[.grouped.data}} for subsetting methods.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Mathieu Pigeon and Louis-Philippe Pouliot
}
\examples{
## Grouped dental data set of Klugman et al. (1998)
cj <- c(0, 25, 50, 100, 150, 250, 500, 1000, 1500, 2500, 4000)
nj <- c(30, 31, 57, 42, 65, 84, 45, 10, 11,  3)
(x <- grouped.data("amount" = cj, "number" = nj))

x[, 1] # group boundaries
x[, 2] # group frequencies
}
\keyword{hplot}
