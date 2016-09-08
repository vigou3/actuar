\name{InverseParalogistic}
\alias{InverseParalogistic}
\alias{dinvparalogis}
\alias{pinvparalogis}
\alias{qinvparalogis}
\alias{rinvparalogis}
\alias{minvparalogis}
\alias{levinvparalogis}
\title{The Inverse Paralogistic Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Inverse Paralogistic
  distribution with parameters \code{shape} and \code{scale}.
}
\usage{
dinvparalogis(x, shape, rate = 1, scale = 1/rate, log = FALSE)
pinvparalogis(q, shape, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
qinvparalogis(p, shape, rate = 1, scale = 1/rate,
              lower.tail = TRUE, log.p = FALSE)
rinvparalogis(n, shape, rate = 1, scale = 1/rate)
minvparalogis(order, shape, rate = 1, scale = 1/rate)
levinvparalogis(limit, shape, rate = 1, scale = 1/rate,
                order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Inverse Paralogistic distribution with parameters \code{shape}
  \eqn{= \tau}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\tau^2 (x/\theta)^{\tau^2}}{%
      x [1 + (x/\theta)^\tau]^{\tau + 1}}}{%
    f(x) = a^2 (x/s)^(a^2)/(x [1 + (x/s)^a]^(a + 1))}
  for \eqn{x > 0}, \eqn{\tau > 0}{a > 0} and \eqn{\theta > 0}{b > 0}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X, d)^k]}.
}
\value{
  \code{dinvparalogis} gives the density,
  \code{pinvparalogis} gives the distribution function,
  \code{qinvparalogis} gives the quantile function,
  \code{rinvparalogis} generates random deviates,
  \code{minvparalogis} gives the \eqn{k}th raw moment, and
  \code{levinvparalogis} gives the \eqn{k}th moment of the limited loss
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
exp(dinvparalogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvparalogis(qinvparalogis(p, 2, 3), 2, 3)
minvparalogis(-1, 2, 2)
levinvparalogis(10, 2, 2, order = 1)
}
\keyword{distribution}
