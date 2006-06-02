\name{Burr}
\alias{burr}
\alias{print.burr}
\alias{dburr}
\alias{pburr}
\alias{qburr}
\alias{rburr}
\title{The Burr Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Burr distribution with parameters 'shape1', 'scale' and 'shape2'.
  
}
\usage{
dburr(x, shape1, scale, shape2, log = FALSE)
pburr(q, shape1, scale, shape2, lower.tail = TRUE,
      log.p = FALSE)
qburr(p, shape1, scale, shape2, lower.tail = TRUE,
      log.p = FALSE)
rburr(n, shape1, scale, shape2)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{shape1, scale, shape2}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  The Burr distribution with parameters 'shape1' = a, 'scale' = c and
  'shape2' = b has density: \deqn{p(x) = \frac{a b (\frac{x}{c})^b}{x(1
      + (\frac{x}{c})^b)^(a + 1)}}{p(x) = (a b (x / c)^b) / (x [1 + (x /
    c)^b]^(a + 1))} for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}.

  If 'shape2' = 1, it is a Pareto distribution. If 'shape2' = 'shape1',
  it is a Paralogistic distribution.

  If 'shape1' = 1, it is a Loglogistic distribution.
}
\value{
'dburr' gives the density, 'pburr' gives the distribution function,
'qgamma' gives the quantile function, and 'rburr' generates random
deviates.

Invalid arguments will results in return warning.
 
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
}
\keyword{distribution}
