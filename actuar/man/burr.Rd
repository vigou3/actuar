\name{Burr}
\alias{burr}
\alias{print.burr}
\alias{dburr}
\alias{pburr}
\alias{qburr}
\alias{rburr}
\alias{mburr}
\alias{levburr}
\title{The Burr Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Burr distribution with parameters \code{shape1}, \code{scale},
  and \code{shape2}.
}
\usage{
  dburr(x, shape1, scale, shape2, log = FALSE)
  pburr(q, shape1, scale, shape2, lower.tail = TRUE,
        log.p = FALSE)
  qburr(p, shape1, scale, shape2, lower.tail = TRUE,
        log.p = FALSE)
  rburr(n, shape1, scale, shape2)
  mburr(k, shape1, scale, shape2, log = FALSE)
  levburr(d, shape1, scale, shape2, order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{d}{limit of the loss variable.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p). Arguments unused for \code{mburr} and \code{levburr} functions.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k, order}{order of the moment.}
}
\details{
  The Burr distribution with parameters \code{shape1} = a, \code{scale} = b and
  \code{shape2} = c has density: \deqn{p(x) = \frac{a c (\frac{x}{b})^c}{x(1
      + (\frac{x}{b})^c)^{(a + 1)}}}{p(x) = (a c (x / b)^c) / (x [1 + (x /
    b)^c]^(a + 1))} for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}.

  The Burr distribution has the following special cases:
  \itemize{
    \item \code{shape2} == 1 is a Pareto distribution, 
    \item \code{shape2} == \code{shape1} is a Paralogistic distribution,
    and
    \item \code{shape1} == 1 is a Loglogistic distribution.
  }    
}
\value{
\code{dburr} gives the density, \code{pburr} gives the distribution function,
\code{qburr} gives the quantile function, \code{rburr} generates random
deviates, \code{mburr} gives the kth raw moment, and \code{levburr} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\note{
  Also known as Burr Type XII and Singh-Maddala.
}
\references{
  Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dburr(2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
pburr(qburr(p, 2, 3, 1), 2, 3, 1)
mburr(2, 1, 2, 3) - mburr(1, 1, 2, 3) ^ 2
levburr(10, 1, 2, 3, 2)
}
\keyword{distribution}
