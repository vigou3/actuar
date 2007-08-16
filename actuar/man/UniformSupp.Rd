\name{UniformSupp}
\alias{UniformSupp}
\alias{munif}
\alias{levunif}
\alias{mgfunif}
\title{Moments, Limited Moments and Moment generating function of the Uniform Distribution}
\description{
  Raw moments, limited moments and moment generating function for the Uniform distribution with parameters \code{min} and \code{max}.
}
\usage{
  munif(order, min = 0, max = 1)
  levunif(limit, min = 0, max =1, order = 1)
  mgfunif(t, min = 0, max = 1, log = FALSE)
}
\arguments{
  \item{min}{first parameter of the uniform distribution.}
  \item{max}{second parameter of the uniform distribution.}
  \item{t}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the random variable.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}. 
  The moment generating function is given by \eqn{E[e^{tX}]}.
}
\value{
  \code{munif} gives the \eqn{k}th raw moment. 
  \code{levunif} gives the \eqn{k}th moment of the limited random
  variable.
  \code{mgfunif} gives the moment generating function in \code{t}.
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Uniform}} for details on the Uniform distribution and
  functions \code{{d,p,q,r}unif}.
}
\references{
  Wikipedia, The uniform distribution (continuous)
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
munif(-1)
munif(1:5)
levunif(3, order=1:5)
levunif(3, 2,4)
mgfunif(1,1,2)
}
\keyword{distribution}
