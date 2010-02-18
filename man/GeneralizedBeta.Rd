\name{GeneralizedBeta}
\alias{GeneralizedBeta}
\alias{dgenbeta}
\alias{pgenbeta}
\alias{qgenbeta}
\alias{rgenbeta}
\alias{mgenbeta}
\alias{levgenbeta}
\title{The Generalized Beta Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments and limited moments for the Generalized Beta distribution
  with parameters \code{shape1}, \code{shape2}, \code{shape3} and
  \code{scale}.
}
\usage{
dgenbeta(x, shape1, shape2, shape3, rate = 1, scale = 1/rate,
         log = FALSE)
pgenbeta(q, shape1, shape2, shape3, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
qgenbeta(p, shape1, shape2, shape3, rate = 1, scale = 1/rate,
         lower.tail = TRUE, log.p = FALSE)
rgenbeta(n, shape1, shape2, shape3, rate = 1, scale = 1/rate)
mgenbeta(order, shape1, shape2, shape3, rate = 1, scale = 1/rate)
levgenbeta(limit, shape1, shape2, shape3, rate = 1, scale = 1/rate,
           order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, shape2, shape3, scale}{parameters. Must be strictly
    positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{ 
  The Generalized Beta distribution with parameters \code{shape1} \eqn{=
    \alpha}{= a}, \code{shape2} \eqn{= \beta}{= b}, \code{shape3}
  \eqn{= \tau}{= c} and \code{scale} \eqn{= \theta}{= s}, has
  density:
  \deqn{f(x) = \frac{\Gamma(\alpha + \beta)}{\Gamma(\alpha)\Gamma(\beta)}
    (x/\theta)^{\alpha \tau} (1 - (x/\theta)^\tau)^{\beta - 1}
    \frac{\tau}{x}}{%
      f(x) = Gamma(a + b)/(Gamma(a) * Gamma(b))
      (c (x/s)^(ac) [1 - (x/s)^c]^(b - 1))/x}
  for \eqn{0 < x < \theta}{0 < x < s}, \eqn{\alpha > 0}{a > 0},
  \eqn{\beta > 0}{b > 0}, \eqn{\tau > 0}{c > 0} and \eqn{\theta > 0}{s >
    0}. (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by \R's \code{\link{gamma}()} and defined in its help.)

  The Generalized Beta is the distribution of the random variable
  \deqn{\theta X^{1/\tau},}{s X^(1/c),}
  where \eqn{X} has a Beta distribution with parameters \eqn{\alpha}{a}
  and \eqn{\beta}{b}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k}th limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]}.
}
\value{
  \code{dgenbeta} gives the density,
  \code{pgenbeta} gives the distribution function,
  \code{qgenbeta} gives the quantile function,
  \code{rgenbeta} generates random deviates,
  \code{mgenbeta} gives the \eqn{k}th raw moment, and
  \code{levgenbeta} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  This is \emph{not} the generalized three-parameter beta distribution
  defined on page 251 of Johnson et al, 1995.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2008),
  \emph{Loss Models, From Data to Decisions, Third Edition}, Wiley.

  Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995)
  \emph{Continuous Univariate Distributions}, volume 2, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
exp(dgenbeta(2, 2, 3, 4, 0.2, log = TRUE))
p <- (1:10)/10
pgenbeta(qgenbeta(p, 2, 3, 4, 0.2), 2, 3, 4, 0.2)
mgenbeta(2, 1, 2, 3, 0.25) - mgenbeta(1, 1, 2, 3, 0.25) ^ 2
levgenbeta(10, 1, 2, 3, 0.25, order = 2)
}
\keyword{distribution}
