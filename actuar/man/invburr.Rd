\name{Invburr}
\alias{invburr}
\alias{print.invburr}
\alias{dinvburr}
\alias{pinvburr}
\alias{qinvburr}
\alias{rinvburr}
\alias{minvburr}
\alias{levinvburr}
\title{The Inverse Burr Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited, moments for the Inverse Burr distribution with parameters \code{shape1}, \code{scale} and \code{shape2}.
}
\usage{
  dinvburr(x, shape1, scale, shape2, log = FALSE)
  pinvburr(q, shape1, scale, shape2, lower.tail = TRUE,
        log.p = FALSE)
  qinvburr(p, shape1, scale, shape2, lower.tail = TRUE,
        log.p = FALSE)
  rinvburr(n, shape1, scale, shape2)
  minvburr(k, shape1, scale, shape2, log = FALSE)
  levinvburr(d, shape1, scale, shape2, order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{d}{limit of the loss variable.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k, order}{order of the moment.}
}
\details{
  The Inverse Burr distribution with parameters \code{shape1} = a, \code{scale} = b and
  \code{shape2} = c has density: \deqn{p(x) = \frac{a c (\frac{x}{b})^{(a
	c)}}{x (1 + (\frac{x}{b})^c)^{(a + 1)}}}{p(x) = a c (x/b)^(a c) / (x (1
    + (x/b)^c)^(a + 1))} for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}.

  The Inverse Burr distribution has the following special cases:
  \itemize{
    \item \code{shape2} == 1 is an Inverse Pareto distribution,
    \item \code{shape2} == \code{shape1} is an Inverse Paralogistic
    distribution, and
    \item \code{shape1} == 1 is a Loglogistic distribution.
  }
}
\value{
\code{dinvburr} gives the density, \code{invburr} gives the distribution function,
\code{qinvburr} gives the quantile function, \code{rinvburr} generates random
deviates, code{minvburr} gives the kth raw moment, and \code{levinvburr} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\note{
  Also known as Dagum distribution.
}
\references{
  See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
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
