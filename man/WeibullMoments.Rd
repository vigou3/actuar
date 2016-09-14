\name{WeibullMoments}
\alias{WeibullMoments}
\alias{mweibull}
\alias{levweibull}
\title{Raw and Limited Moments of the Weibull Distribution}
\description{
  Raw moments and limited moments for the Weibull distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
mweibull(order, shape, scale = 1)
levweibull(limit, shape, scale = 1, order = 1)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{shape, scale}{shape and scale parameters, the latter defaulting
    to 1.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}, \eqn{k >
  -\tau}{k > -shape}.
}
\value{
  \code{mweibull} gives the \eqn{k}th raw moment and
  \code{levweibull} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Weibull}} for details on the Weibull distribution and
  functions \code{[dpqr]weibull}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
mweibull(2, 3, 4) - mweibull(1, 3, 4)^2
levweibull(10, 3, 4, order = 2)
}
\keyword{distribution}
