\name{Logarithmic}
\alias{Logarithmic}
\alias{dlogarithmic}
\alias{plogarithmic}
\alias{qlogarithmic}
\alias{rlogarithmic}
\alias{log-series}
\title{The Logarithmic Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Logarithmic (or log-series) distribution with parameter
  \code{prob}.
}
\usage{
dlogarithmic(x, prob, log = FALSE)
plogarithmic(q, prob, lower.tail = TRUE, log.p = FALSE)
qlogarithmic(p, prob, lower.tail = TRUE, log.p = FALSE)
rlogarithmic(n, prob)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{prob}{parameter. \code{0 <= prob < 1}}.
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The Logarithmic (or log-series) distribution with parameter
  \code{prob} \eqn{= \theta}{= p} has probability mass function
  \deqn{%
    P[X = x] = \frac{a \theta^x}{x}}{%
    P[X = x] = a p^x / x},
  with \eqn{a = -1/\log(1 - \theta)}{a = -1/log(1 - p)} and for \eqn{x =
    1, 2, \ldots}, \eqn{0 \le p < 1}.

  The logarithmic distribution is the limiting case of the
  \link[=ZTNegBinomial]{zero truncated negative binomial distribution}
  with \code{size} parameter equal to \eqn{0}. Note that in this
  context, parameter \code{prob} generally corresponds to the
  probability of \emph{failure} of the zero truncated negative binomial.

  If an element of \code{x} is not integer, the result of
  \code{dlogarithmic} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
  \code{dlogarithmic} gives the probability mass function,
  \code{plogarithmic} gives the distribution function,
  \code{qlogarithmic} gives the quantile function, and
  \code{rlogarithmic} generates random deviates.

  Invalid \code{prob} will result in return value \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rlogarithmic}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  \code{qlogarithmic} is based on \code{qbinom} et al.; it uses the
  Cornish--Fisher Expansion to include a skewness correction to a normal
  approximation, followed by a search.

  \code{rlogarithmic} is an implementation of the LS and LK algorithms
  of Kemp (1981) with automatic selection. As suggested by Devroye
  (1986), the LS algorithm is used when \code{prob < 0.95}, and the LK
  algorithm otherwise.
}
\references{
  Johnson, N. L., Kemp, A. W. and Kotz, S. (2005), \emph{Univariate
  Discrete Distributions, Third Edition}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.

  Kemp, A. W. (1981), \dQuote{Efficient Generation of Logarithmically
  Distributed Pseudo-Random Variables}, \emph{Journal of the Royal
  Statistical Society, Series C}, vol. 30, p. 249-253.
  \url{http://www.jstor.org/stable/2346348}

  Devroye, L. (1986), \emph{Non-Uniform Random Variate Generation},
  Springer-Verlag. \url{http://luc.devroye.org/rnbookindex.html}
}
\seealso{
  \code{\link{dztnbinom}} for the zero truncated negative binomial
  distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
## Table 1 of Kemp (1981) [also found in Johnson et al. (2005), chapter 7]
p <- c(0.1, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9999)
round(rbind(dlogarithmic(1, p),
            dlogarithmic(2, p),
            plogarithmic(9, p, lower.tail = FALSE),
            -p/((1 - p) * log(1 - p))), 2)

qlogarithmic(plogarithmic(1:10, 0.9), 0.9)

x <- rlogarithmic(1000, 0.8)
y <- sort(unique(x))
plot(y, table(x)/length(x), type = "h", lwd = 2,
     pch = 19, col = "black", xlab = "x", ylab = "p(x)",
     main = "Empirical vs theoretical probabilities")
points(y, dlogarithmic(y, prob = 0.8),
       pch = 19, col = "red")
legend("topright", c("empirical", "theoretical"),
       lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
}
\keyword{distribution}
