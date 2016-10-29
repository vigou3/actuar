\name{PoissonInverseGaussian}
\alias{PoissonInverseGaussian}
\alias{PIG}
\alias{dpoisinvgauss}
\alias{ppoisinvgauss}
\alias{qpoisinvgauss}
\alias{rpoisinvgauss}
\alias{dpig}
\alias{ppig}
\alias{qpig}
\alias{rpig}
\title{The Poisson-Inverse Gaussian Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Poisson-inverse Gaussian discrete distribution with
  parameters \code{mean} and \code{shape}.
}
\usage{
dpoisinvgauss(x, mean, shape = 1, dispersion = 1/shape,
              log = FALSE)
ppoisinvgauss(q, mean, shape = 1, dispersion = 1/shape,
              lower.tail = TRUE, log.p = FALSE)
qpoisinvgauss(p, mean, shape = 1, dispersion = 1/shape,
              lower.tail = TRUE, log.p = FALSE)
rpoisinvgauss(n, mean, shape = 1, dispersion = 1/shape)
}
\arguments{
  \item{x}{vector of (positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{mean, shape}{parameters. Must be strictly positive. Infinite
    values are supported.}
  \item{dispersion}{an alternative way to specify the shape.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The Poisson-inverse Gaussian distribution is the result of the
  continuous mixture between a Poisson distribution and an inverse
  Gaussian, that is, the distribution with probability mass function
  \deqn{%
    p(x) = \int_0^\infty \frac{\lambda^x e^{-\lambda}}{x!}\,
           g(\lambda; \mu, \phi)\, d\lambda,}{%
    p(x) = int_0^Inf (y^x exp(-y))/x! g(y; \mu, \phi) dy,}
  where \eqn{g(\lambda; \mu, \phi)}{g(y; \mu, \phi)} is the density
  function of the inverse Gaussian distribution with parameters
  \code{mean} \eqn{= \mu} and \code{dispersion} \eqn{= \phi} (see
  \code{\link{dinvgauss}}).

  The resulting probability mass function is
  \deqn{%
    p(x) = \sqrt{\frac{2}{\pi \phi}}
    \frac{e^{(\phi\mu)^{-1}}}{x!}
    \left(
      \sqrt{2\phi\left(1 + \frac{1}{2\phi\mu^2}\right)}
    \right)^{-(x - \frac{1}{2})}
    K_{x - \frac{1}{2}}
    \left(
      \sqrt{\frac{2}{\phi}\left(1 + \frac{1}{2\phi\mu^2}\right)}
    \right),}{%
    p(x) = sqrt(2/(\pi \phi)) exp(1/(\phi \mu))/x!
           * [\sqrt(2 \phi (1 + 1/(2 \phi \mu^2)))]^(-(x-1/2))
           * K(\sqrt((2/\phi) (1 + 1/(2 \phi \mu^2))); x-1/2),}
  for \eqn{x = 0, 1, \dots}, \eqn{\mu > 0}, \eqn{\phi > 0} and where
  \eqn{K_\nu(x)}{K(x; \nu)} is the modified Bessel function of the third
  kind implemented by \R's \code{\link{besselK}()} and defined in its
  help.

  The limiting case \eqn{\mu = \infty}{\mu = Inf} has well defined
  probability mass and distribution functions, but has no finite
  strictly positive, integer moments. The pmf in this case reduces to
  \deqn{%
    p(x) = \sqrt{\frac{2}{\pi \phi}}
    \frac{1}{x!}
    (\sqrt{2\phi})^{-(x - \frac{1}{2})}
    K_{x - \frac{1}{2}}(\sqrt{2/\phi}).}{%
    p(x) = sqrt(2/(\pi \phi)) 1/x! [\sqrt(2 \phi)]^(-(x-1/2))
           * K(\sqrt(2/\phi); x-1/2).}

  The limiting case \eqn{\phi = 0} is a degenerate distribution in
  \eqn{x = 0}.

  If an element of \code{x} is not integer, the result of
  \code{dpoisinvgauss} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
  \code{dpoisinvgauss} gives the probability mass function,
  \code{ppoisinvgauss} gives the distribution function,
  \code{qpoisinvgauss} gives the quantile function, and
  \code{rpoisinvgauss} generates random deviates.

  Invalid arguments will result in return value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rpoisinvgauss}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  \code{[dpqr]pig} are aliases for \code{[dpqr]poisinvgauss}.

  \code{qpoisinvgauss} is based on \code{qbinom} et al.; it uses the
  Cornish--Fisher Expansion to include a skewness correction to a normal
  approximation, followed by a search.
}
\references{
  Holla, M. S. (1966), \dQuote{On a Poisson-Inverse Gaussian
  Distribution}, \emph{Metrika}, vol. 15, p. 377-384.

  Johnson, N. L., Kemp, A. W. and Kotz, S. (2005), \emph{Univariate
  Discrete Distributions, Third Edition}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.

  Shaban, S. A., (1981) \dQuote{Computation of the poisson-inverse
    gaussian distribution}, \emph{Communications in Statistics - Theory
    and Methods}, vol. 10, no. 14, p. 1389-1399.
}
\seealso{
  \code{\link{dpois}} for the Poisson distribution,
  \code{\link{dinvgauss}} for the inverse Gaussian distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Tables I and II of Shaban (1981)
x <- 0:2
sapply(c(0.4, 0.8, 1), dpoisinvgauss, x = x, mean = 0.1)
sapply(c(40, 80, 100, 130), dpoisinvgauss, x = x, mean = 1)

qpoisinvgauss(ppoisinvgauss(0:10, 1, dis = 2.5), 1, dis = 2.5)

x <- rpoisinvgauss(1000, 1, dis = 2.5)
y <- sort(unique(x))
plot(y, table(x)/length(x), type = "h", lwd = 2,
     pch = 19, col = "black", xlab = "x", ylab = "p(x)",
     main = "Empirical vs theoretical probabilities")
points(y, dpoisinvgauss(y, 1, dis = 2.5),
       pch = 19, col = "red")
legend("topright", c("empirical", "theoretical"),
       lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
}
\keyword{distribution}
