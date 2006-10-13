\name{InverseBurr}
\alias{InverseBurr}
\alias{dinvburr}
\alias{pinvburr}
\alias{qinvburr}
\alias{rinvburr}
\alias{minvburr}
\alias{levinvburr}
\title{The Inverse Burr Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Inverse Burr distribution with
  parameters \code{shape1}, \code{shape2} and \code{scale}.
}
\usage{
  dinvburr(x, shape1, shape2, rate = 1, scale = 1/rate,
           log = FALSE)
  pinvburr(q, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
  qinvburr(p, shape1, shape2, rate = 1, scale = 1/rate,
           lower.tail = TRUE, log.p = FALSE)
  rinvburr(n, shape1, shape2, rate = 1, scale = 1/rate)
  minvburr(order, shape1, shape2, rate = 1, scale = 1/rate)
  levinvburr(limit, shape1, shape2, rate = 1, scale = 1/rate,
             order = 1)
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
  The Inverse Burr distribution with parameters \code{shape1} \eqn{=
    \tau}{= a}, \code{shape2} \eqn{= \gamma}{= b} and \code{scale}
  \eqn{= \theta}{= s}, has density:
  \deqn{f(x) = \frac{\tau \gamma (x/\theta)^{\gamma \tau}}{%
      x [1 + (x/\theta)^\gamma]^{\tau + 1}}}{%
    f(x) = a b (x/s)^(ba)/(x [1 + (x/s)^b]^(a + 1))}
  for \eqn{x > 0}, \eqn{\tau > 0}{a > 0}, \eqn{\gamma > 0}{b > 0} and
  \eqn{\theta > 0}{s > 0}.

  The Inverse Burr is the distribution of the random variable
  \deqn{\theta \left(\frac{X}{1 - X}\right)^{1/\gamma},}{%
    s (X/(1 - X))^(1/b),}
  where \eqn{X} has a Beta distribution with parameters \eqn{\tau}{a}
  and \eqn{1}.

  The Inverse Burr distribution has the following special cases:
  \itemize{
    \item A \link[=dllogis]{Loglogistic} distribution when \code{shape1
      == 1};
    \item An \link[=dinvpareto]{Inverse Pareto} distribution when
    \code{shape2 == 1};
    \item An \link[=dinvparalogis]{Inverse Paralogistic} distribution
    when \code{shape1 == shape2}.
  }

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the \eqn{k} limited moment at some limit
  \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]}.
}
\value{
  \code{dinvburr} gives the density,
  \code{invburr} gives the distribution function,
  \code{qinvburr} gives the quantile function,
  \code{rinvburr} generates random deviates,
  \code{minvburr} gives the \eqn{k}th raw moment, and
  \code{levinvburr} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Also known as the Dagum distribution.
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
exp(dinvburr(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pinvburr(qinvburr(p, 2, 3, 1), 2, 3, 1)
minvburr(2, 1, 2, 3) - minvburr(1, 1, 2, 3) ^ 2
levinvburr(10, 1, 2, 3, 2)
}
\keyword{distribution}
