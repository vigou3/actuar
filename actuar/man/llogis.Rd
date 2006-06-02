
\name{Llogis}
\alias{llogis}
\alias{print.llogis}
\alias{dllogis}
\alias{pllogis}
\alias{qllogis}
\alias{rllogis}
\title{The Loglogistic Distribution}
\description{
  Density, distribution function, quantile function and random generation
  for the Loglogistic distribution with parameters 'shape' and 'scale'.
  
}
\usage{
dllogis(x, shape, scale, log = FALSE)
pllogis(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
qllogis(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
rllogis(n, shape, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If 'length(n) > 1', the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if 'TRUE', probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if 'TRUE' (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].} 
  
}
\details{
  The Loglogistic distribution with parameters 'shape' = a
  and 'scale' = b has density:
  \deqn{p(x) = \frac{a (\frac{x}{b})^a}{x(1 + (\frac{x}{b})^a)^2}}{p(x)
    = a (x/b)^a /(x(1 + (x/b)^a)^2)}
  for \eqn{x > 0, a > 0} and \eqn{b > 0}.
  
}
\value{
'dllogis' gives the density, 'pllogis' gives the distribution function,
'qllogis' gives the quantile function, and 'rllogis' generates random
deviates.

Invalid arguments will results in return warning.
}
\note{
Also known as Fisk distribution.
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
exp(dllogis(2, 3, 4, log = TRUE))
p <- (1:10)/10
pllogis(qllogis(p, 2, 3), 2, 3)
}
\keyword{distribution}
