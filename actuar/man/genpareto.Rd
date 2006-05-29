
\name{Genpareto}
\alias{genpareto}
\alias{print.genpareto}
\alias{dgenpareto}
\alias{pgenpareto}
\alias{qgenpareto}
\alias{rgenpareto}
\title{The Generalized Pareto Distribution}
\description{
  Density, distribution function, quantile function and random function
  for the Generalized Pareto distribution with parameters 'shape1', 'scale' and 'shape2'.
  
}
\usage{
dgenpareto(x, shape1, scale, shape2, log = FALSE)
pgenpareto(q, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
qgenpareto(p, shape1, scale, shape2, lower.tail = TRUE, log.p = FALSE)
rgenpareto(n, shape1, scale, shape2)
}
\arguments{
  \item{q}{vector of quantiles.}
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
  
  The Generalized Pareto distribution with parameters 'shape1' = a, 'scale' = c and
  'shape2' = b has density: \deqn{p(x) = (gamma(a + b) /
    gamma(a)*gamma(b))((c^a x^(c-1)) / (x + c)^(a + b))} for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's 'gamma()' and defined in its help.
}
\value{
'dgenpareto' gives the density, 'pgenpareto' gives the distribution function,
'qgenpareto' gives the quantile function, and 'rgenpareto' generates random
deviates.

Invalid arguments will results in return warning.
 
}
\note{
  It' a beta of the second kind.
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
exp(dgenpareto(3,3,4,4,log = TRUE))
p <- (1:10)/10
pgenpareto(qgenpareto(p,3,3,1),3,3,1)
qgenpareto(3,3,4,4,lower.tail = FALSE)
}
\keyword{distribution}
