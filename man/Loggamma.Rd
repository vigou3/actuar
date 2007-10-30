\name{Loggamma}
\alias{Loggamma}
\alias{dlgamma}
\alias{plgamma}
\alias{qlgamma}
\alias{rlgamma}
\alias{mlgamma}
\alias{levlgamma}
\title{The Loggamma Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Loggamma distribution with
  parameters \code{shapelog} and \code{ratelog}.
}
\usage{
dlgamma(x, shapelog, ratelog, log = FALSE)
plgamma(q, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
qlgamma(p, shapelog, ratelog, lower.tail = TRUE, log.p = FALSE)
rlgamma(n, shapelog, ratelog)
mlgamma(order, shapelog, ratelog)
levlgamma(limit, shapelog, ratelog, order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shapelog, ratelog}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Loggamma distribution with parameters \code{shapelog} \eqn{=
    \alpha}{= a} and \code{ratelog} \eqn{= \lambda}{= b} has density:
  \deqn{f(x) = \frac{\lambda^\alpha}{\Gamma(\alpha)}%
    \frac{(\log x)^{\alpha - 1}}{x^{\lambda + 1}}}{%
    f(x) = (b^a (log(x))^(a - 1))/(Gamma(a) * x^(b + 1))}
  for \eqn{x > 1}, \eqn{\alpha > 0}{a > 0} and \eqn{\lambda > 0}{b >
    0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by R's \code{\link{gamma}()} and defined in its help.)

  The Loggamma is the distribution of the random variable
  \eqn{e^X}{exp(X)},  where \eqn{X} has a Gamma distribution with
  shape parameter \eqn{alpha}{a} and scale parameter
  \eqn{1/\lambda}{1/b}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}.
}
\value{
  \code{dlgamma} gives the density,
  \code{plgamma} gives the distribution function,
  \code{qlgamma} gives the quantile function,
  \code{rlgamma} generates random deviates,
  \code{mlgamma} gives the \eqn{k}th raw moment, and
  \code{levlgamma} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\references{
  Hogg, R. V. and Klugman, S. A. (1984), \emph{Loss Distributions},
  Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dlgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
plgamma(qlgamma(p, 2, 3), 2, 3)
mlgamma(2, 3, 4) - mlgamma(1, 3, 4)^2
levlgamma(10, 3, 4, order = 2)
}
\keyword{distribution}
