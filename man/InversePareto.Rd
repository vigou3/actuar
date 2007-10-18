\name{InversePareto}
\alias{InversePareto}
\alias{dinvpareto}
\alias{pinvpareto}
\alias{qinvpareto}
\alias{rinvpareto}
\alias{minvpareto}
\alias{levinvpareto}
\title{The Inverse Pareto Distribution}
\description{
  Density, distribution function, quantile function, random generation
  raw moments and limited moments for the Inverse Pareto distribution
  with parameters \code{shape} and \code{scale}.
}
\usage{
  dinvpareto(x, shape, scale, log = FALSE)
  pinvpareto(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
  qinvpareto(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
  rinvpareto(n, shape, scale)
  minvpareto(order, shape, scale)
  levinvpareto(limit, shape, scale, order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Inverse Pareto distribution with parameters \code{shape} \eqn{=
    \tau}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\tau \theta x^{\tau - 1}}{%
      (x + \theta)^{\tau + 1}}}{%
    f(x) = a s x^(a - 1)/(x + s)^(a + 1)}
  for \eqn{x > 0}, \eqn{\tau > 0}{a > 0} and \eqn{\theta > 0}{s > 0}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}.
}
\value{
  \code{dinvpareto} gives the density,
  \code{pinvpareto} gives the distribution function,
  \code{qinvpareto} gives the quantile function,
  \code{rinvpareto} generates random deviates,
  \code{minvpareto} gives the \eqn{k}th raw moment, and
  \code{levinvpareto} calculates the \eqn{k}th limited moment.

  Invalid arguments will result in return value \code{NaN}, with a warning.
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
exp(dinvpareto(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvpareto(qinvpareto(p, 2, 3), 2, 3)
minvpareto(-1, 2, order = 2)
}
\keyword{distribution}
