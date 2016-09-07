\name{ZeroTruncatedBinomial}
\alias{ZeroTruncatedBinomial}
\alias{ZTBinomial}
\alias{dztbinom}
\alias{pztbinom}
\alias{qztbinom}
\alias{rztbinom}
\title{The Zero Truncated Binomial Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the zero truncated binomial distribution with
  parameters \code{size} and \code{prob}.
}
\usage{
dztbinom(x, size, prob, log = FALSE)
pztbinom(q, size, prob, lower.tail = TRUE, log.p = FALSE)
qztbinom(p, size, prob, lower.tail = TRUE, log.p = FALSE)
rztbinom(n, size, prob)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{size}{number of trials (strictly positive integer).}
  \item{prob}{probability of success on each trial. \code{0 <= prob <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero truncated binomial distribution with \code{size}
  \eqn{= n} and \code{prob} \eqn{= p} has probability mass function
  \deqn{%
    p(x) = {n \choose x} \frac{p^x (1 - p)^{n-x}}{1 - (1 - p)^n}}{%
    p(x) = choose(n, x) [p^x (1-p)^(n-x)]/[1 - (1-p)^n]}
  for \eqn{x = 1, \ldots, n} and \eqn{0 < p \le 1}, and \eqn{p(1) = 1}
  when \eqn{p = 0}.

  The mean is \eqn{np/(1 - (1-p)^n)} and the variance is
  \eqn{np[(1-p) - (1-p+np)(1-p)^n]/[1 - (1-p)^n]^2}.

  In the terminology of Klugman et al. (2012), the zero truncated
  binomial is a member of the \eqn{(a, b, 1)} family of
  distributions with \eqn{a = -p/(1-p)} and \eqn{b = (n+1)q/(1-p)}.

  If an element of \code{x} is not integer, the result of
  \code{dztbinom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{F(x) \ge p}, where \eqn{F} is the distribution function.
}
\value{
  \code{dztbinom} gives the probability mass function,
  \code{pztbinom} gives the distribution function,
  \code{qztbinom} gives the quantile function, and
  \code{rztbinom} generates random deviates.

  Invalid \code{size} or \code{prob} will result in return value
  \code{NaN}, with a warning.

  The length of the result is determined by \code{n} for
  \code{rztbinom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}ztbinom} use \code{\{d,p,q\}binom} for all
  but the trivial input values and \eqn{p(0)}.

  \code{rztbinom} uses the simple inversion algorithm suggested by
  Peter Dalgaard on the r-help mailing list on 1 May 2005 %
  (\url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}).
}
\references{
  Johnson, N. L., Kemp, A. W. and Kotz, S. (2005), \emph{Univariate
  Discrete Distributions, Third Edition}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dbinom}} for the binomial distribution.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
dztbinom(1:5, size = 5, prob = 0.4)
dbinom(1:5, 5, 0.4)/pbinom(0, 5, 0.4, lower = FALSE) # same

pztbinom(1, 2, prob = 0)        # point mass at 1

qztbinom(pztbinom(1:10, 10, 0.6), 10, 0.6)

x <- rztbinom(1000, size = 10, prob = 0.4)
y <- sort(unique(x))
plot(y, table(x)/length(x), type = "h", lwd = 2,
     pch = 19, col = "black", xlab = "x", ylab = "p(x)",
     main = "Empirical vs theoretical probabilities")
points(y, dztbinom(y, size = 10, prob = 0.4),
       pch = 19, col = "red")
legend("topright", c("empirical", "theoretical"),
       lty = c(1, NA), pch = c(NA, 19), col = c("black", "red"))
}
\keyword{distribution}
