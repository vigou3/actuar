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
  \item{x}{an object of class \code{"grouped.data"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The mean of grouped data with group boundaries \eqn{c_1, \dots,
    c_r}{cj[1], \dots, cj[r]} and group frequencies \eqn{n_1, \dots,
    n_r}{nj[1], \dots, nj[r]} is
  \deqn{\sum_{j = 1}^r \frac{c_{j - 1} + c_j}{2}\, n_j.}{%
    sum((cj[j-1] + cj[j])/2 * nj).}
}
\value{
  A named vector of means.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{emm}} to compute higher moments.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(gdental)
mean(gdental)
}
\keyword{univar}
