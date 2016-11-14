\name{UniformSupp}
\alias{UniformSupp}
\alias{munif}
\alias{levunif}
\alias{mgfunif}
\title{Moments and Moment Generating Function of the Uniform Distribution}
\description{
  Raw moments, limited moments and moment generating function for the
  Uniform distribution from \code{min} to \code{max}.
}
\usage{
munif(order, min = 0, max = 1)
levunif(limit, min = 0, max =1, order = 1)
mgfunif(t, min = 0, max = 1, log = FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{min, max}{lower and upper limits of the distribution. Must be finite.}
  \item{limit}{limit of the random variable.}
  \item{t}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]} and the moment
  generating function is \eqn{E[e^{tX}]}.
}
\value{
  \code{munif} gives the \eqn{k}th raw moment,
  \code{levunif} gives the \eqn{k}th moment of the limited random
  variable, and
  \code{mgfunif} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Uniform}}.
}
\references{
  \url{http://en.wikipedia.org/wiki/Uniform_distribution_\%28continuous\%29}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
munif(-1)
munif(1:5)
levunif(3, order=1:5)
levunif(3, 2, 4)
mgfunif(1, 1, 2)
}
\keyword{distribution}
