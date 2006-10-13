\name{ExponentialMoments}
\alias{ExponentialMoments}
\alias{mexp}
\alias{levexp}
\title{Raw and Limited Moments of the Exponential Distribution}
\description{
  Raw moments and limited moments for the exponential distribution with
  rate \code{rate} (i.e., mean \code{1/rate}).
}
\usage{
  mexp(order, rate = 1)
  levexp(limit, rate = 1, order = 1)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{rate}{vector of rates.}
}
\value{
  \code{mexp} gives the \eqn{k}th raw moment and
  \code{levexp} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Exponential}} for details on the exponential distribution and
  functions \code{{d,p,q,r}exp}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
mexp(2, 3) - mexp(1, 3)^2
levexp(10, 3, 2)
}
\keyword{distribution}
