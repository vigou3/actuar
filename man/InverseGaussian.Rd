\name{InverseGaussian}
\alias{InverseGaussian}
\alias{dinvgauss}
\alias{pinvgauss}
\alias{qinvgauss}
\alias{rinvgauss}
\alias{minvgauss}
\alias{levinvgauss}
\alias{mgfinvgauss}
\title{The Inverse Gaussian Distribution}
\description{
  Density function, distribution function, quantile function, random
  generation, raw moments, limited moments and moment generating
  function for the Inverse Gaussian distribution with parameters
  \code{mean} and \code{shape}.
}
\usage{
dinvgauss(x, mean, shape = 1, dispersion = 1/shape,
          log = FALSE)
pinvgauss(q, mean, shape = 1, dispersion = 1/shape,
          lower.tail = TRUE, log.p = FALSE)
qinvgauss(p, mean, shape = 1, dispersion = 1/shape,
          lower.tail = TRUE, log.p = FALSE,
          tol = 1e-14, maxit = 100, echo = FALSE, trace = echo)
rinvgauss(n, mean, shape = 1, dispersion = 1/shape)
minvgauss(order, mean, shape = 1, dispersion = 1/shape)
levinvgauss(limit, mean, shape = 1, dispersion = 1/shape, order = 1)
mgfinvgauss(t, mean, shape = 1, dispersion = 1/shape, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mean, shape}{parameters. Must be strictly positive. Infinite
    values are supported.}
  \item{dispersion}{an alternative way to specify the shape.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment. Only \code{order = 1} is
    supported by \code{levinvgauss}.}
  \item{limit}{limit of the loss variable.}
  \item{tol}{small positive value. Tolerance to assess convergence in
    the Newton computation of quantiles.}
  \item{maxit}{positive integer; maximum number of recursions in the
    Newton computation of quantiles.}
  \item{echo, trace}{logical; echo the recursions to screen in the
    Newton computation of quantiles.}
  \item{t}{numeric vector.}
}
\details{
  The inverse gaussian distribution with parameters \code{mean} \eqn{=
    \mu} and \code{dispersion} \eqn{= \phi} has density:
  \deqn{f(x) = \left( \frac{1}{2 \pi \phi x^3} \right)^{1/2}
    \exp\left( -\frac{(x - \mu)^2}{2 \mu^2 \phi x} \right),}{%
    f(x) = sqrt(/(2 \pi \phi x^3)) *
           exp(-((x - \mu)^2)/(2 \mu^2 \phi x)),}
  for \eqn{x \ge 0}, \eqn{\mu > 0} and \eqn{\phi > 0}.

  The limiting case \eqn{\mu = \infty}{\mu = Inf} has a well defined
  density function and distribution function, but has no finite moments.

  The limiting case \eqn{\phi = 0} is an infinite spike in \eqn{x = 0}.

  If the random variable \eqn{X} is IG\eqn{(\mu, \phi)}, then
  \eqn{X/\mu} is IG\eqn{(1, \phi \mu)}{(1, \phi * \mu)}.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, \eqn{k = 1, 2, \dots}, the limited expected
  value at some limit \eqn{d} is \eqn{E[\min(X, d)]}{E[min(X, d)]} and
  the moment generating function is \eqn{E[e^{tX}]}.

  The moment generating function of the inverse guassian is defined for
    \code{t <= 1/(2 * mean^2 * phi)}.
}
\value{
  \code{dinvgauss} gives the density,
  \code{pinvgauss} gives the distribution function,
  \code{qinvgauss} gives the quantile function,
  \code{rinvgauss} generates random deviates,
  \code{minvgauss} gives the \eqn{k}th raw moment,
  \code{levinvgauss} gives the limited expected value, and
  \code{mgfinvgauss} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Functions \code{dinvgauss}, \code{pinvgauss} and \code{qinvgauss} are
  C implementations of functions of the same name in package
  \pkg{statmod}; see Giner and Smyth (2016).

  Devroye (1986, chapter 4) provides a nice presentation of the
  algorithm to generate random variates from an inverse gaussian
  distribution.
}
\references{
  Giner, G. and Smyth, G. K. (2016), \dQuote{\pkg{statmod}: Probability
    Calculations for the Inverse Gaussian Distribution}, \emph{R
    Journal}, vol. 8, no 1, p. 339-351.
  \url{https://journal.r-project.org/archive/2016-1/giner-smyth.pdf}

  Chhikara, R. S. and Folk, T. L. (1989), \emph{The Inverse Gaussian
    Distribution: Theory, Methodology and Applications}, Decker.

  Devroye, L. (1986), \emph{Non-Uniform Random Variate Generation},
  Springer-Verlag. \url{http://luc.devroye.org/rnbookindex.html}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dinvgauss(c(-1, 0, 1, 2, Inf), mean = 1.5, dis = 0.7)
dinvgauss(c(-1, 0, 1, 2, Inf), mean = Inf, dis = 0.7)
dinvgauss(c(-1, 0, 1, 2, Inf), mean = 1.5, dis = Inf) # spike at zero

## Typical graphical representations of the inverse gaussian
## distribution. First fixed mean and varying shape; second
## varying mean and fixed shape.
col = c("red", "blue", "green", "cyan", "yellow", "black")
par = c(0.125, 0.5, 1, 2, 8, 32)
curve(dinvgauss(x, 1, par[1]), from = 0, to = 2, col = col[1])
for (i in 2:6)
    curve(dinvgauss(x, 1, par[i]), add = TRUE, col = col[i])

curve(dinvgauss(x, par[1], 1), from = 0, to = 2, col = col[1])
for (i in 2:6)
    curve(dinvgauss(x, par[i], 1), add = TRUE, col = col[i])

pinvgauss(qinvgauss((1:10)/10, 1.5, shape = 2), 1.5, 2)

minvgauss(1:4, 1.5, 2)

levinvgauss(c(0, 0.5, 1, 1.2, 10, Inf), 1.5, 2)
}
\keyword{distribution}
