
\name{Genpareto}
\alias{genpareto}
\alias{print.genpareto}
\alias{dgenpareto}
\alias{pgenpareto}
\alias{qgenpareto}
\alias{rgenpareto}
\alias{mgenpareto}
\alias{levgenpareto}
\title{The Generalized Pareto Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Generalized Pareto distribution with parameters \code{shape1}, \code{scale}, and \code{shape2}.
}
\usage{
  dgenpareto(x, shape1, scale, shape2, log = FALSE)
  pgenpareto(q, shape1, scale, shape2, lower.tail = TRUE,
             log.p = FALSE)
  qgenpareto(p, shape1, scale, shape2, lower.tail = TRUE,
             log.p = FALSE)
  rgenpareto(n, shape1, scale, shape2)
  mgenpareto(k, shape1, scale, shape2, log = FALSE)
  levgenpareto(d, shape1, scale, shape2, order = 1, log = FALSE)
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
  The Generalized Pareto distribution with parameters \code{shape1} = a, \code{scale} = b and
  \code{shape2} = c has density: \deqn{p(x) = \frac{\Gamma(a +
      c)}{\Gamma(a)\Gamma(c)}\frac{b^a x^{(c - 1)}}{(x + b)^{(a +
      c)}}}{p(x) = (gamma(a + c) / gamma(a)*gamma(c))((b^a x^(b-1)) / (x + b)^(a + c))} for \eqn{x > 0, a > 0, b > 0} and
  \eqn{c > 0}. Here Gamma(a) is the function implemented by R's \code{gamma()} and defined in its help.

  Let Y be a Beta random variable with parameters
  \code{shape1} and \code{shape2}. Then, \deqn{scale(\frac{Y}{1 - Y})}{
    scale (Y / (1 - Y))} has the Generalized Pareto probability density function with parameters
  \code{shape1}, \code{scale} and \code{shape2}.

  The Generalized Pareto distribution has the following special cases:
  \itemize{
    \item \code{shape2} == 1 is a Pareto distribution, and
    \item \code{shape2} == \code{shape1} == 1 is an Inverse Pareto
    distribution.
  }
}
\value{
\code{dgenpareto} gives the density, \code{pgenpareto} gives the distribution function,
\code{qgenpareto} gives the quantile function, \code{rgenpareto} generates random
deviates, \code{mgenpareto} gives the kth raw moment, and \code{levgenpareto} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\note{
  Also known as the Beta of the second kind.
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
exp(dgenpareto(3, 3, 4, 4, log = TRUE))
p <- (1:10)/10
pgenpareto(qgenpareto(p, 3, 3, 1), 3, 3, 1)
qgenpareto(.3, 3, 4, 4, lower.tail = FALSE)
mgenpareto(1, 3, 2, 1) ^ 2
levgenpareto(10, 1, 2, 3, 2)
}
\keyword{distribution}
