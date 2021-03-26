\name{var.grouped.data}
\alias{var.grouped.data}
\title{Sample Variance}
\description{
  Sample variance of grouped data objects.
}
\usage{
var.grouped.data(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The sample variance of grouped data with group boundaries \eqn{c_1, 
    \dots, c_r}{c[1], \dots, c[r]} and group frequencies \eqn{n_1, 
    \dots, n_r}{n[1], \dots, n[r]} is
  \deqn{\sum_{j = 1}^r \frac{n_j ( n (c_{j-1} + c_j) - \sum_{k=1}^r {(c_{k - 1} + c_k} ))}{4 n^2 (n-1)},}{%
        sum(j; (n[j] * ( n * (c[j-1] + c[j]) -  sum(k; (c[k - 1] +
	c[k]))/(4 * n^2 * (n-1)),}
  where \eqn{n = \sum_{j = 1}^r n_j}{n = sum(j; n[j])}.
}
\value{
  A named vector of sample variances.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{sd.grouped.data}} to compute the sample standard deviation;
  \code{\link{emm}} to compute higher moments.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Walter Garcia-Fontes \email{walter.garcia@upf.edu}
}
\examples{
data(gdental)
var.grouped.data(gdental)
}
\keyword{univar}
