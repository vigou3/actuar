\name{InvGaussSupp}
\alias{InvGaussSupp}
\alias{minvGauss}
\alias{levinvGauss}
\alias{mgfinvGauss}
\title{Raw, Limited Moments and Moment generating function of the inverse Gaussian Distribution}
\description{
  Raw moments, limited moments, and moment generating function for the inverse Gaussian distribution with parameters \code{nu} and \code{lambda}.
}
\usage{
  minvGauss(order, nu, lambda )
  levinvGauss(limit, nu, lambda, order = 1)
  mgfinvGauss(t, nu, lambda, log= FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{nu}{the mean of the inverse Gaussian distribution. Must be strictly positive}
  \item{lambda}{the scale paramter. Must be strictly positive}
  \item{t}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}. 
  The moment generating function is given by \eqn{E[e^{tX}]}.
  
}
\value{
  \code{minvGauss} gives the \eqn{k}th raw moment and
  \code{levinvGauss} gives the \eqn{k}th moment of the limited loss
  variable. Only \eqn{k=1} is supported.
  \code{mgfinvGauss} gives the moment generating function in \code{t}
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  Details on the inverse Gaussian distribution can be found in the 2 references below. For the density, the distribution fuction, the quantile and the random generator of the inverse Gaussian distribution, use the package SuppDists, and the functions \code{{d,p,q,r}invGauss}.
}
\references{
  Chhikara & Folk (1989), 
 \emph{The inverse Gaussian distribution: theory, methodology and application}, Decker

  Seshadri (1989), 
  \emph{The inverse Gaussian distribution: statistical theory and applications}, Springer

}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
minvGauss(2, 3, 4) 
levinvGauss(10, 3, 4)
mgfinvGauss(1,3,2)
}
\keyword{distribution}
