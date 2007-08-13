\name{UniformSupp}
\alias{UniformSupp}
\alias{mgfunif}
\title{Moment generating function of the Uniform Distribution}
\description{
  moment generating function for the Uniform distribution with
  parameters \code{min} and \code{max}.
}
\usage{
  mgfunif(x, min = 0, max = 1, log = FALSE)
}
\arguments{
  \item{min}{first parameter of the uniform distribution}
  \item{max}{second parameter of the uniform distribution}
  \item{x}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned}
}
\details{
  The moment generating function is given by \eqn{E[e^{xX}]}.
}
\value{
  \code{mgfunif} gives the moment generating function in \code{x}
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Uniform}} for details on the Uniform distribution and
  functions \code{{d,p,q,r}unif}.
}
\references{
  Johnson & Kotz, Loss Distributions, Wiley, 1970
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
mgfunif(1,1,2)
}
\keyword{distribution}
