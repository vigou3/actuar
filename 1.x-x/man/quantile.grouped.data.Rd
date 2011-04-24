\name{quantile.grouped.data}
\alias{quantile.grouped.data}
\title{Quantiles of Grouped Data}
\description{
  Sample quantiles corresponding to the given probabilities for objects
  of class \code{"grouped.data"}.
}
\usage{
\method{quantile}{grouped.data}(x, probs = seq(0, 1, 0.25),
         names = TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"grouped.data"}.}
  \item{probs}{numeric vector of probabilities with values
    in \eqn{[0, 1]}.}
  \item{names}{logical; if true, the result has a \code{names}
    attribute.  Set to \code{FALSE} for speedup with many \code{probs}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The quantile function is the inverse of the ogive, that is a linear
  interpolation of the empirical quantile function.

  The equation of the quantile function is
  \deqn{x = \frac{c_j (F_n(c_{j - 1}) - q) +
      c_{j - 1} (q - F_n(c_j)}{F_n(c_j) - F_n(c_{j - 1})}}{%
    x = (c[j] (Fn(c[j-1]) - q) +
    c[j-1] (q - Fn(c[j])))/(Fn(c[j]) - Fn(c[j-1]))}
  for \eqn{0 \leq q \leq c_j}{0 <= q <= 1} and where \eqn{c_0, \dots,
  c_r}{c[0], \dots, c[r]} are the \eqn{r + 1} group
  boundaries and \eqn{F_n}{Fn} is the empirical distribution function of
  the sample.
}
\value{
  A numeric vector, named if \code{names} is \code{TRUE}.
}
\seealso{
  \code{\link{ogive}} for the smoothed empirical distribution of which
  \code{quantile.grouped.data} is an inverse;
  \code{\link{grouped.data}} to create grouped data objects.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
data(gdental)
quantile(gdental)
Fn <- ogive(gdental)
Fn(quantile(gdental))		# inverse function
}
\keyword{univar}
