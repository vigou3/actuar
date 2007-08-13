\name{NormalSupp}
\alias{NormalSupp}
\alias{mgfnorm}
\title{Moment generating function of the Normal Distribution}
\description{
  moment generating function for the Uniform distribution with
  parameters \code{min} and \code{max}.
}
\usage{
  mgfnorm(x, mean = 0, sd = 1, log = FALSE)
}
\arguments{
  \item{mean}{mean of the Normal distribution}
  \item{sd}{standard deviation of the Normal distribution}
  \item{x}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned}
}
\details{
  The moment generating function is given by \eqn{E[e^{xX}]}.
}
\value{
  \code{mgfnorm} gives the moment generating function in \code{x}
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Normal}} for details on the Normal distribution and
  functions \code{{d,p,q,r}norm}.
}
\references{
  Johnson & Kotz, Loss Distributions, Wiley, 1970
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
mgfnorm(0:4,1,2)
}
\keyword{distribution}
