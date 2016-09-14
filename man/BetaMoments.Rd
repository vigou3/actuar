\name{BetaMoments}
\alias{BetaMoments}
\alias{mbeta}
\alias{levbeta}
\title{Raw and Limited Moments of the Beta Distribution}
\description{
  Raw moments and limited moments for the (central) Beta distribution
  with parameters \code{shape1} and \code{shape2}.
}
\usage{
mbeta(order, shape1, shape2)
levbeta(limit, shape1, shape2, order = 1)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{shape1, shape2}{positive parameters of the Beta distribution.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}, \eqn{k > -\alpha}{k
    > -shape1}.

  The noncentral beta distribution is not supported.
}
\value{
  \code{mbeta} gives the \eqn{k}th raw moment and
  \code{levbeta} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a
  warning.
}
\seealso{
  \code{\link[stats]{Beta}} for details on the beta distribution and
  functions \code{[dpqr]beta}.
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
mbeta(2, 3, 4) - mbeta(1, 3, 4)^2
levbeta(10, 3, 4, order = 2)
}
\keyword{distribution}
