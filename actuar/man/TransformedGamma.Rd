\name{TransformedGamma}
\alias{TransformedGamma}
\alias{dtrgamma}
\alias{ptrgamma}
\alias{qtrgamma}
\alias{rtrgamma}
\alias{mtrgamma}
\alias{levtrgamma}
\title{The Transformed Gamma Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Transformed Gamma distribution
  with parameters \code{shape1}, \code{shape2} and \code{scale}.
}
\usage{
  dtrgamma(x, shape1, shape2, rate = 1, scale = 1/rate,
           log = FALSE)
  ptrgamma(q, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
  qtrgamma(p, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
  rtrgamma(n, shape1, shape2, rate = 1, scale = 1/rate)
  mtrgamma(order, shape1, shape2, rate = 1, scale = 1/rate,
           log = FALSE)
  levtrgamma(limit, shape1, shape2, rate = 1, scale = 1/rate,
             order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, shape2, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Transformed Gamma distribution with parameters \code{shape1} \eqn{=
    \alpha}{= a}, \code{shape2} \eqn{= \tau}{= b} and \code{scale}
  \eqn{= \theta}{= s}, has density:
  \deqn{f(x) = \frac{\tau u^\alpha e^{-u}}{x \Gamma(\alpha)}, %
      \quad u = (x/\theta)^\tau}{%
    f(x) = b u^a exp(-u) / (x Gamma(a)), u = (x/s)^b}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0}, \eqn{\tau > 0}{b > 0}
  and \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by R's \code{\link{gamma}()} and defined in its help.)

  The Transformed Gamma is the distribution of the random variable
  \eqn{\theta X^{1/\tau},}{s X^(1/b),}
  where \eqn{X} has a Gamma distribution with shape parameter
  \eqn{\alpha}{a} and scale parameter \eqn{1} or, equivalently, of the
  random variable
  \eqn{Y^{1/\tau}}{Y^(1/b)}
  with \eqn{Y} a Gamma distribution with shape parameter \eqn{\alpha}{a}
  and scale parameter \eqn{\theta^\tau}{s^b}.

  The Transformed Gamma probability distribution defines a family of
  distributions with the following special cases:
  \itemize{
    \item A \link[=dgamma]{Gamma} distribution when \code{shape2 == 1};
    \item A \link[=dweibull]{Weibull} distribution when \code{shape1 ==
      1};
    \item An \link[=dexp]{Exponential} distribution when \code{shape2 ==
      shape1 == 1}.
  }
}
\value{
  \code{dtrgamma} gives the density,
  \code{ptrgamma} gives the distribution function,
  \code{qtrgamma} gives the quantile function,
  \code{rtrgamma} generates random deviates,
  \code{mtrgamma} gives the \eqn{k}th raw moment, and
  \code{levtrgamma} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Distribution also known as the Generalized Gamma.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dtrgamma(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrgamma(qtrgamma(p, 2, 3, 4), 2, 3, 4)
mtrgamma(2, 3, 4, 5) - mtrgamma(1, 3, 4, 5) ^ 2
levtrgamma(10, 3, 4, 5, 2)
}
\keyword{distribution}
