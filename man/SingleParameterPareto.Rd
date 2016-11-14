\name{SingleParameterPareto}
\alias{SingleParameterPareto}
\alias{dpareto1}
\alias{ppareto1}
\alias{qpareto1}
\alias{rpareto1}
\alias{mpareto1}
\alias{levpareto1}
\title{The Single-parameter Pareto Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Single-parameter Pareto
  distribution with parameter \code{shape}.
}
\usage{
dpareto1(x, shape, min, log = FALSE)
ppareto1(q, shape, min, lower.tail = TRUE, log.p = FALSE)
qpareto1(p, shape, min, lower.tail = TRUE, log.p = FALSE)
rpareto1(n, shape, min)
mpareto1(order, shape, min)
levpareto1(limit, shape, min, order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{shape}{parameter. Must be strictly positive.}
  \item{min}{lower bound of the support of the distribution.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The single-parameter Pareto distribution with parameter \code{shape}
  \eqn{= \alpha}{= a} has density:
  \deqn{f(x) = \frac{\alpha \theta^\alpha}{x^{\alpha + 1}}}{%
    f(x) = a b^a/x^(a + 1)}
  for \eqn{x > \theta}{x > b}, \eqn{\alpha > 0}{a > 0} and \eqn{\theta >
    0}{b > 0}.

  Although there appears to be two parameters, only \code{shape} is a true
  parameter. The value of \code{min} \eqn{= \theta}{= b} must be set in
  advance.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, \eqn{k < \alpha}{k < shape} and the \eqn{k}th
  limited moment at some limit \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X,
  d)^k]}, \eqn{x \ge \theta}{x \ge min}.
}
\value{
  \code{dpareto1} gives the density,
  \code{ppareto1} gives the distribution function,
  \code{qpareto1} gives the quantile function,
  \code{rpareto1} generates random deviates,
  \code{mpareto1} gives the \eqn{k}th raw moment, and
  \code{levpareto1} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
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
exp(dpareto1(5, 3, 4, log = TRUE))
p <- (1:10)/10
ppareto1(qpareto1(p, 2, 3), 2, 3)
mpareto1(2, 3, 4) - mpareto(1, 3, 4) ^ 2
levpareto(10, 3, 4, order = 2)
}
\keyword{distribution}
