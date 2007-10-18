\name{ChisqSupp}
\alias{ChisqSupp}
\alias{mchisq}
\alias{levchisq}
\alias{mgfchisq}
\title{Raw, Limited Moments and Moment generating function of the non central Chi-square Distribution}
\description{
  Raw moments, limited moments, and moment generating function for the non central Chi-square distribution with
  parameters \code{shape} and \code{scale}.
}
\usage{
  mchisq (order, df, ncp = 0)
  levchisq(limit, df, ncp = 0, order = 1)
  mgfchisq(x, df, ncp = 0, log= FALSE)
}
\arguments{
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{df}{degree of freedom of the non central Chi-square distribution. Must be strictly positive}
  \item{ncp}{the non central paramter. Must be positive or zero}
  \item{x}{a real}
  \item{log}{logical, if \code{TRUE}, the log of the m.g.f. is returned}
}
\details{
  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]}. 
  The moment generating function is given by \eqn{E[e^{xX}]}.
  Note that the limited expected value is not supported for the non central Chi-square distribution (\code{ncp}>0), and the order of the moment must be an integer for the non central Chi-square distribution(\code{ncp}>0).
}
\value{
  \code{mchisq} gives the \eqn{k}th raw moment and
  \code{levchisq} gives the \eqn{k}th moment of the limited loss
  variable.
  \code{mgfchisq} gives the moment generating function in \code{x}
  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\seealso{
  \code{\link{Chisquare}} for details on the non central Chi-square distribution and
  functions \code{{d,p,q,r}chisq}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
 
  Johnson & Kotz, Loss Distributions, Wiley, 1970
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Christophe Dutang
}
\examples{
mchisq(2, 3, 4) 
levchisq(10, 3, order = 2)
mgfchisq(1,3,2)
}
\keyword{distribution}
