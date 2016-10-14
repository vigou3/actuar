\name{InvGaussSupp}
\alias{InvGaussSupp}
\alias{minvGauss-deprecated}
\alias{levinvGauss}
\alias{mgfinvGauss}
\alias{minvgauss}
\alias{levinvgauss}
\alias{mgfinvgauss}
\title{Moments and Moment Generating Function of the Inverse Gaussian Distribution}
\description{
  Raw moments, limited moments and moment generating function for the
  Inverse Gaussian distribution with parameters \code{nu} and
  \code{lambda}.
}
\usage{
minvGauss(order, nu, lambda)
levinvGauss(limit, nu, lambda, order = 1)
mgfinvGauss(x, nu, lambda, log= FALSE)

minvgauss(order, nu, lambda)
levinvgauss(limit, nu, lambda, order = 1)
mgfinvgauss(x, nu, lambda, log= FALSE)
}
\arguments{
  \item{order}{order of the moment. Only \code{order = 1} is
    supported by \code{levinvGauss}.}
  \item{limit}{limit of the loss variable.}
  \item{nu, lambda}{parameters. Must be strictly positive.}
  \item{x}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]} and the moment generating function is
  \eqn{E[e^{xX}]}, \eqn{k = 1, 2, \dots}.

  The moment generating function is defined for \code{x <= lambda/(2 *
    nu^2)}.
}
\value{
  \code{minvGauss} gives the \eqn{k}th raw moment,
  \code{levinvGauss} gives the \eqn{k}th moment of the limited loss
  variable, and
  \code{mgfinvGauss} gives the moment generating function in \code{x}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link[SuppDists]{invGauss}} in package \pkg{SuppDists} for the
  density function, distribution function, quantile function and random
  number generator.
}
\references{
  Chhikara, R. S. and Folk, T. L. (1989), \emph{The Inverse Gaussian
    Distribution: Theory, Methodology and Applications}, Decker.

  Seshadri, D. N. (1989), \emph{The Inverse Gaussian Distribution:
  Statistical Theory and Applications}, Springer.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
minvGauss(2, 3, 4)
levinvGauss(10, 3, 4)
mgfinvGauss(0.2, 3, 4)
}
\keyword{distribution}
