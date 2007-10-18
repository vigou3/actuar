\name{NormalSupp}
\alias{NormalSupp}
\alias{mnorm}
\alias{mgfnorm}
\title{Moments and Moment generating function of the Normal Distribution}
\description{
  moments and the moment generating function for the Normal distribution with
  parameters \code{mean} and \code{sd}.
}
\usage{
  mnorm(order, mean = 0, sd = 1)
  mgfnorm(x, mean = 0, sd = 1, log = FALSE)
}
\arguments{
  \item{mean}{mean of the Normal distribution}
  \item{sd}{standard deviation of the Normal distribution}
  \item{x}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned}
  \item{order}{integer, the order of the moment}
}
\details{
  The moment is given by \eqn{E[X^k]}, where \eqn{k} denotes the order.
  The moment generating function is given by \eqn{E[e^{xX}]}.
  Only integer order is supported for the \eqn{k}th moment.
}
\value{
  \code{mnorm} gives the moment.
  \code{mgfnorm} gives the moment generating function in \code{x}.
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
mnorm(3)
}
\keyword{distribution}
