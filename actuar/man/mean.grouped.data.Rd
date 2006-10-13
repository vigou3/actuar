\name{mean.grouped.data}
\alias{mean.grouped.data}
\title{Extract or Replace Parts of a Grouped Data Object}
\description{
  Mean of grouped data objects.
}
\usage{
\method{mean}{grouped.data}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{grouped.data}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The mean of grouped data with \eqn{n_j}{nj} values in group \eqn{(c_{j
      - 1}, c_j]}{(c{j-1}, c{j}]}, \eqn{j = 1, \dots, r} is
  \deqn{\sum_{j = 1}^r \frac{c_{j - 1} + c_j}{2}\, n_j.}{%
    sum ((c{j-1} + c{j})/2 * nj).}
}
\value{
  A named vector of means.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(gdental)
mean(gdental)
}
\keyword{univar}
