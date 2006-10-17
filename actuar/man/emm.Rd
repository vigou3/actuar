\name{emm}
\alias{emm}
\alias{emm.default}
\alias{emm.grouped.data}
\title{Empirical Moments}
\description{
  Raw empirical moments for individual and grouped data.
}
\usage{
emm(x, order = 1, \dots)
\method{emm}{default}(x, order = 1, \dots)
\method{emm}{grouped.data}(x, order = 1, \dots)
}
\arguments{
  \item{x}{a vector or matrix of individual data, or an object of class
    \code{"grouped data"}.}
  \item{order}{order of the moment. Must be positive.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  Arguments \code{\dots} are passed to \code{\link{colMeans}};
  \code{na.rm = TRUE} may be useful for individual data with missing
  values. 

  For individual data, the \eqn{k}th empirical moment is
  \deqn{\sum_{j = 1}^n x_j^k.}{sum(x[j]^k).}

  For grouped data with group boundaries \eqn{c_1, \dots, c_r}{cj[1],
    \dots, cj[r]} and group frequencies \eqn{n_1, \dots, n_r}{nj[1],
    \dots, nj[r]}, the \eqn{k}th empirical moment is
  \deqn{\sum_{j = 1}^r \frac{n_j (c_j^k - c_{j - 1}^k)}{%
      n (k + 1) (c_j - c_{j - 1})},}{%
    sum((nj[j] * (cj[j]^k - cj[j-1]^k)/(n * (k+1) * (cj[j] -
    cj[j-1]))),}
  where \eqn{n = \sum_{j = 1}^r n_j}{n = sum(nj)}.
}
\value{
  A vector or matrix of moments.
}
\seealso{
  \code{\link{mean}} and \code{\link{mean.grouped.data}} for simpler
  access to the first moment.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
## Individual data
data(dental)
emm(dental, order = 1:3)

## Grouped data
data(gdental)
emm(gdental)
x <- grouped.data(cj = gdental[, 1],
                  nj1 = sample(1:100, nrow(gdental)),
                  nj2 = sample(1:100, nrow(gdental)))
emm(x) # same as mean(x)
}
\keyword{univar}
