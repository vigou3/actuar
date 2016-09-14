\name{LognormalMoments}
\alias{LognormalMoments}
\alias{mlnorm}
\alias{levlnorm}
\title{Raw and Limited Moments of the Lognormal Distribution}
\description{
  Raw moments and limited moments for the Lognormal distribution whose
  logarithm has mean equal to \code{meanlog} and standard deviation
  equal to \code{sdlog}.
}
\usage{
mlnorm(order, meanlog = 0, sdlog = 1)
levlnorm(limit, meanlog = 0, sdlog = 1, order = 1)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{meanlog, sdlog}{mean and standard deviation of the distribution
    on the log scale with default values of \code{0} and \code{1}
    respectively.}
}
\value{
  \code{mlnorm} gives the \eqn{k}th raw moment and
  \code{levlnorm} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Lognormal}} for details on the lognormal distribution and
  functions \code{[dpqr]lnorm}.
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
mlnorm(2, 3, 4) - mlnorm(1, 3, 4)^2
levlnorm(10, 3, 4, order = 2)
}
\keyword{distribution}
