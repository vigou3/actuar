\name{Extract.grouped.data}
\alias{Extract.grouped.data}
\alias{[.grouped.data}
\alias{[<-.grouped.data}
\title{Extract or Replace Parts of a Grouped Data Object}
\description{
  Extract or replace subsets of grouped data objects.
}
\usage{
\method{[}{grouped.data}(x, i, j)
\method{[}{grouped.data}(x, i, j) <- value
}
\arguments{
  \item{x}{an object of class \code{grouped.data}.}

  \item{i, j}{elements to extract or replace. \code{i, j} are
    \code{numeric} or \code{character} or, for \code{[} only, empty.
    Numeric values are coerced to integer as if by
    \code{\link[base]{as.integer}}. For replacement by \code{[}, a logical
    matrix is allowed, but not replacement in the group boundaries and
    group frequencies simultaneously.}

  \item{value}{a suitable replacement value.}
}
\details{
  Objects of class \code{"grouped.data"} can mostly be indexed like data
  frames, with the following restrictions:
  \enumerate{
    \item For \code{[}, the extracted object must keep a group
      boundaries column and at least one group frequencies column to
      remain of class \code{"grouped.data"};
    \item For \code{[<-}, it is not possible to replace group boundaries
      and group frequencies simultaneously;
    \item When replacing group boundaries, \code{length(value) ==
      length(i) + 1}.
  }

  \code{x[, 1]} will return the plain vector of group boundaries.
  
  Replacement of non adjacent group boundaries is not possible for
  obvious reasons.
  
  Otherwise, extraction and replacement should work just like for data
  frames.
}
\value{
  For \code{[} an object of class \code{"grouped.data"}, a data frame or a
  vector.

  For \code{[<-} an object of class \code{"grouped.data"}.
}
\note{
  Currently \code{[[}, \code{[[<-}, \code{$} and \code{$<-} are not 
  specifically supported, but should work as usual on group frequency
  columns.
}
\seealso{
  \code{\link[base]{[.data.frame}} for extraction and replacement
  methods of data frames, \code{\link{grouped.data}} to create grouped
  data objects.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(gdental)

(x <- gdental[1])         # select column 1
class(x)                  # no longer a grouped.data object
class(gdental[2])         # same
gdental[, 1]              # group boundaries
gdental[, 2]              # group frequencies

gdental[1:4,]             # a subset
gdental[c(1, 3, 5),]      # avoid this

gdental[1:2, 1] <- c(0, 30, 60) # modified boundaries
gdental[, 2] <- 10              # modified frequencies
\dontrun{gdental[1, ] <- 2}   # not allowed
}
\keyword{manip}
\keyword{array}
