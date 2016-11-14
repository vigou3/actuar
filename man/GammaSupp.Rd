\name{GammaSupp}
\alias{GammaSupp}
\alias{mgamma}
\alias{levgamma}
\alias{mgfgamma}
\title{Moments and Moment Generating Function of the Gamma Distribution}
\description{
  Raw moments, limited moments and moment generating function for the
  Gamma distribution with parameters \code{shape} and \code{scale}.
}
\usage{
mgamma(order, shape, rate = 1, scale = 1/rate)
levgamma(limit, shape, rate = 1, scale = 1/rate, order = 1)
mgfgamma(t, shape, rate = 1, scale = 1/rate, log = FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{rate}{an alternative way to specify the scale.}
  \item{shape, scale}{shape and scale parameters. Must be strictly
    positive.}
  \item{t}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]} and the moment
  generating function is \eqn{E[e^{tX}]}, \eqn{k >
  -\alpha}{k > -shape}.
}
\value{
  \code{mgamma} gives the \eqn{k}th raw moment,
  \code{levgamma} gives the \eqn{k}th moment of the limited loss
  variable, and
  \code{mgfgamma} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link[stats]{GammaDist}}
}
\references{
  Johnson, N. L. and Kotz, S. (1970), \emph{Continuous Univariate
    Distributions, Volume 1}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang and Mathieu Pigeon
}
\examples{
mgamma(2, 3, 4) - mgamma(1, 3, 4)^2
levgamma(10, 3, 4, order = 2)
mgfgamma(1,3,2)
}
\keyword{distribution}
