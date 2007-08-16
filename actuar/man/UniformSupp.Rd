\name{UniformSupp}
\alias{UniformSupp}
\alias{munif}
\alias{mgfunif}
\title{Moments and Moment generating function of the Uniform Distribution}
\description{
  moment generating function for the Uniform distribution with
  parameters \code{min} and \code{max}.
}
\usage{
  munif(order, min = 0, max = 1)
  mgfunif(t, min = 0, max = 1, log = FALSE)
}
\arguments{
  \item{min}{first parameter of the uniform distribution.}
  \item{max}{second parameter of the uniform distribution.}
  \item{t}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned.}
  \item{order}{order of the moment.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}.
  The moment generating function is given by \eqn{E[e^{tX}]}.
}
\value{
  \code{munif} gives the \eqn{k}th raw moment. 
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
mgfunif(1,1,2)
}
\keyword{distribution}
