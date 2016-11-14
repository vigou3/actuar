\name{ChisqSupp}
\alias{ChisqSupp}
\alias{mchisq}
\alias{levchisq}
\alias{mgfchisq}
\title{Moments and Moment Generating Function of the (non-central) Chi-Squared Distribution}
\description{
  Raw moments, limited moments and moment generating function for the
  chi-squared (\eqn{\chi^2}{chi^2}) distribution with \code{df} degrees
  of freedom and optional non-centrality parameter \code{ncp}.
}
\usage{
mchisq(order, df, ncp = 0)
levchisq(limit, df, ncp = 0, order = 1)
mgfchisq(t, df, ncp = 0, log= FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{df}{degrees of freedom (non-negative, but can be non-integer).}
  \item{ncp}{non-centrality parameter (non-negative).}
  \item{t}{numeric vector.}
  \item{log}{logical; if \code{TRUE}, the cumulant generating function
    is returned.}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]} and the moment generating
  function is \eqn{E[e^{tX}]}.

  Only integer moments are supported for the non central Chi-square
  distribution (\code{ncp > 0}).

  The limited expected value is supported for the centered Chi-square
  distribution (\code{ncp = 0}).
}
\value{
  \code{mchisq} gives the \eqn{k}th raw moment,
  \code{levchisq} gives the \eqn{k}th moment of the limited loss
  variable, and
  \code{mgfchisq} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link[stats]{Chisquare}}
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.

  Johnson, N. L. and Kotz, S. (1970), \emph{Continuous Univariate
    Distributions, Volume 1}, Wiley.
}
\author{
  Christophe Dutang, Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
mchisq(2, 3, 4)
levchisq(10, 3, order = 2)
mgfchisq(0.25, 3, 2)
}
\keyword{distribution}
