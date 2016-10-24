\name{grouped.data}
\alias{grouped.data}
\title{Grouped data}
\description{
  Creation of grouped data objects, either from a provided set of group
  boundaries and group frequencies, or from individual data using
  automatic or specified breakpoints.
}
\usage{
grouped.data(\dots, breaks = "Sturges", include.lowest = TRUE,
             right = TRUE, row.names = NULL,
             nclass = NULL, group = FALSE,
             check.rows = FALSE, check.names = TRUE)
}
\arguments{
  \item{\dots}{these arguments are either of the form \code{value} or
    \code{tag = value}. See Details.}
  \item{breaks}{same as for \code{\link{hist}}, namely one of:
    \itemize{
      \item{a vector giving the breakpoints between groups;}
      \item{a function to compute the vector of breakpoints;}
      \item{a single number giving the number of groups;}
      \item{a character string naming an algorithm to compute the
      number of groups (see \code{\link{hist}});}
      \item{a function to compute the number of groups.}
    }
    In the last three cases the number is a suggestion only; the
    breakpoints will be set to \code{\link{pretty}} values.  If
    \code{breaks} is a function, the first element in \code{\dots}
    is supplied to it as the only argument.
  }
  \item{include.lowest}{logical; if \code{TRUE}, a data point equal to
    the \code{breaks} value will be included in the first (or last, for
    \code{right = FALSE}) group. Used only for individual data; see
    Details.}
  \item{right}{logical; indicating if the intervals should be closed on
    the right (and open on the left) or vice versa.}
  \item{nclass}{numeric (integer); equivalent to \code{breaks} for a
    scalar or character argument.}
  \item{group}{logical; an alternative way to force grouping of
    individual data.}
  \item{row.names, check.rows, check.names}{arguments identical to those
    of \code{\link{data.frame}}.}
}
\details{
  A grouped data object is a special form of data frame consisting of
  one column of contiguous group boundaries and one or more columns of
  frequencies within each group.

  The function can create a grouped data object from two types of
  arguments.
  \enumerate{
    \item{Group boundaries and frequencies. This is the default mode of
      operation if the call has at least two elements in \code{\dots}.

      The first argument will then be taken as the vector of group
      boundaries. This vector must be exactly one element longer than
      the other arguments, which will be taken as vectors of group
      frequencies. All arguments are coerced to data frames.}
    \item{Individual data. This mode of operation is active if there
      is a single argument in \code{\dots}, or if either \code{breaks}
      or \code{nclass} is specified or \code{group} is \code{TRUE}.

      Arguments of \code{\dots} are first grouped using
      \code{\link{hist}}. If needed, breakpoints are set using the first
      argument.}
  }

  Missing (\code{NA}) frequencies are replaced by zeros, with a
  warning.

  Extraction and replacement methods exist for \code{grouped.data}
  objects, but working on non adjacent groups will most likely yield
  useless results.
}
\value{
  An object of \code{class} \code{c("grouped.data", "data.frame")} with
  an environment containing the vector \code{cj} of group boundaries.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\seealso{
  \code{\link{[.grouped.data}} for extraction and replacement methods.

  \code{\link{data.frame}} for usual data frame creation and
  manipulation.

  \code{\link{hist}} for details on the calculation of breakpoints.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Mathieu Pigeon and Louis-Philippe Pouliot
}
\examples{
## Most common usage using a predetermined set of group
## boundaries and group frequencies.
cj <- c(0, 25, 50, 100, 250, 500, 1000)
nj <- c(30, 31, 57, 42, 45, 10)
(x <- grouped.data(Group = cj, Frequency = nj))
class(x)

x[, 1] # group boundaries
x[, 2] # group frequencies

## Multiple frequency columns are supported
x <- sample(1:100, 9)
y <- sample(1:100, 9)
grouped.data(cj = 1:10, nj.1 = x, nj.2 = y)

## Alternative usage with grouping of individual data.
grouped.data(x)                         # automatic breakpoints
grouped.data(x, breaks = 7)             # forced number of groups
grouped.data(x, breaks = c(0,25,75,100))    # specified groups
grouped.data(x, y, breaks = c(0,25,75,100)) # multiple data sets

\dontrun{## Providing two or more data sets and automatic breakpoints is
## very error-prone since the range of the first data set has to
## include the ranges of all the other data sets.
range(x)
range(y)
grouped.data(x, y, group = TRUE)}
}
\keyword{classes}
\keyword{methods}
