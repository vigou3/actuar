\name{NormalSupp}
\alias{NormalSupp}
\alias{mnorm}
\alias{mgfnorm}
\title{Moments and Moment generating function of the Normal Distribution}
\description{
  Raw moments and moment generating function for the normal distribution with
  mean equal to \code{mean} and standard deviation equal to \code{sd}.
}
\usage{
mnorm(order, mean = 0, sd = 1)
mgfnorm(t, mean = 0, sd = 1, log = FALSE)
}
\arguments{
  \item{order}{vector of integers; order of the moment.}
  \item{mean}{vector of means.}
  \item{sd}{vector of standard deviations.}
  \item{t}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the moment generating function is
  \eqn{E[e^{tX}]}.

  Only integer moments are supported.
}
\value{
  \code{mnorm} gives the \eqn{k}th raw moment and
  \code{mgfnorm} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Normal}}
}
\references{
  Johnson, N. L. and Kotz, S. (1970), \emph{Continuous Univariate
    Distributions, Volume 1}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
mgfnorm(0:4,1,2)
mnorm(3)
}
\keyword{distribution}
