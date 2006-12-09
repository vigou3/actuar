\name{discretize}
\alias{discretize}
\alias{discretise}
\title{Discretization of a Continuous Distribution}
\description{
  Compute a discrete probability mass function from a continuous
  cumulatve distribution function with various methods.

  \code{discretise} is an alias for \code{discretize}.
}
\usage{
discretize(expr, from, to, step = 1, method = c("upper", "lower"),
           by = step, xlim = NULL)

discretise(expr, from, to, step = 1, method = c("upper", "lower"),
           by = step, xlim = NULL)
}
\arguments{
  \item{expr}{an expression written as a function of \code{x}, or
    alternatively the name of a function, giving the cdf to discretize.}
  \item{from, to}{the range over which the function will be discretized.}
  \item{step}{numeric; the discretization step (or span, or lag).}
  \item{method}{discretization method to use.}
  \item{by}{an alias for \code{step}.}
  \item{xlim}{numeric of length 2; if specified, it serves as default
    for \code{c(from, to)}.}
}
\details{
  Usage is similar to \code{\link{curve}}.

  \code{discretize} returns the probability mass function (pmf) of the
  random variable corresponding to the discretized cdf speficied in
  \code{expr}.

  Method \code{"upper"} is the forward difference of the cdf \eqn{F}:
  \deqn{f(x) = F(x + h) - F(x).}

  Method \code{"lower"} is the backward difference of the cdf \eqn{F}:
  \deqn{f(x) = F(x) - F(x - h).}
}
\value{
  A numeric vector of probabilities suitable for use in
  \code{\link{aggregateDist}}.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{aggregateDist}}
}
\examples{
x <- seq(from = 0, to = qgamma(1 - 1E-6, 2, 1), by = 0.5)
fx <- discretize(pgamma(x, 2, 1), method = "upper",
           from = 0, to = tail(x, 1), step = 0.5)
plot(stepfun(head(x, -1), diffinv(fx)))
curve(pgamma(x, 2, 1), add = TRUE)
}
\keyword{distribution}
\keyword{models}
