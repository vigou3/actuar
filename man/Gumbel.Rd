\name{Gumbel}
\alias{Gumbel}
\alias{dgumbel}
\alias{pgumbel}
\alias{qgumbel}
\alias{rgumbel}
\alias{mgumbel}
\alias{mgfgumbel}
\title{The Gumbel Distribution}
\description{
  Density function, distribution function, quantile function, random
  generation and raw moments for the Gumbel extreme value distribution
  with parameters \code{alpha} and \code{scale}.
}
\usage{
dgumbel(x, alpha, scale, log = FALSE)
pgumbel(q, alpha, scale, lower.tail = TRUE, log.p = FALSE)
qgumbel(p, alpha, scale, lower.tail = TRUE, log.p = FALSE)
rgumbel(n, alpha, scale)
mgumbel(order, alpha, scale)
mgfgumbel(t, alpha, scale, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{alpha}{location parameter.}
  \item{scale}{parameter. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment. Only values \eqn{1} and \eqn{2} are
    supported.}
  \item{t}{numeric vector.}
}
\details{
  The Gumbel distribution with parameters \code{alpha} \eqn{=
    \alpha}{= a} and \code{scale} \eqn{= \theta}{= s} has distribution
  function:
  \deqn{F(x) = \exp[-\exp(-(x - \alpha)/\theta)]}{%
    F(x) = exp[-exp(-(x - a)/s)],}
  for \eqn{-\infty < x < \infty}{-Inf < x < Inf}, \eqn{-\infty < a <
    \infty}{-Inf < a < Inf} and \eqn{\theta > 0}{s > 0}.

  The mode of the distribution is in \eqn{\alpha}{a}, the mean is
  \eqn{\alpha + \gamma\theta}{a + g * s}, where \eqn{\gamma}{g} \eqn{=
  0.57721566} is the Euler-Mascheroni constant, and the variance is
  \eqn{\pi^2 \theta^2/6}{(pi * s)^2/6}.
}
\value{
  \code{dgumbel} gives the density,
  \code{pgumbel} gives the distribution function,
  \code{qgumbel} gives the quantile function,
  \code{rgumbel} generates random deviates,
  \code{mgumbel} gives the \eqn{k}th raw moment, \eqn{k = 1, 2}, and
  \code{mgfgamma} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Distribution also knonw as the generalized extreme value distribution
  Type-I.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dgumbel(c(-5, 0, 10, 20), 0.5, 2)

p <- (1:10)/10
pgumbel(qgumbel(p, 2, 3), 2, 3)

curve(pgumbel(x, 0.5, 2), from = -5, to = 20, col = "red")
curve(pgumbel(x, 1.0, 2), add = TRUE, col = "green")
curve(pgumbel(x, 1.5, 3), add = TRUE, col = "blue")
curve(pgumbel(x, 3.0, 4), add = TRUE, col = "cyan")

a <- 3; s <- 4
mgumbel(1, a, s)                        # mean
a - s * digamma(1)                      # same

mgumbel(2, a, s) - mgumbel(1, a, s)^2   # variance
(pi * s)^2/6                            # same
}
\keyword{distribution}
