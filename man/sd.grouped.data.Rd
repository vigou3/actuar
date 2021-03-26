\name{sd.grouped.data}
\alias{sd.grouped.data}
\title{Sample Standard Deviation}
\description{
  Sample standard deviation of grouped data objects.
}
\usage{
sd.grouped.data(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The sample sample standard deviation  of grouped data with group
  boundaries \eqn{c_1, \dots, c_r}{c[1], \dots, c[r]} and group
  frequencies \eqn{n_1, \dots, n_r}{n[1], \dots, n[r]} is the squared
  root of the sample variance.
}
\value{
  A named vector of means.
}
\seealso{
  \code{\link{grouped.data}} to create grouped data objects;
  \code{\link{var.grouped.data}} to compute the sample variance;
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
sd.grouped.data(gdental)
}
\keyword{univar}
