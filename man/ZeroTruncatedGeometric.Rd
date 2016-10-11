\name{ZeroTruncatedGeometric}
\alias{ZeroTruncatedGeometric}
\alias{ZTGeometric}
\alias{dztgeom}
\alias{pztgeom}
\alias{qztgeom}
\alias{rztgeom}
\title{The Zero-Truncated Geometric Distribution}
\description{
  Density function, distribution function, quantile function and random
  generation for the Zero-Truncated Geometric distribution with
  parameter \code{prob}.
}
\usage{
dztgeom(x, prob, log = FALSE)
pztgeom(q, prob, lower.tail = TRUE, log.p = FALSE)
qztgeom(p, prob, lower.tail = TRUE, log.p = FALSE)
rztgeom(n, prob)
}
\arguments{
  \item{x}{vector of (strictly positive integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{prob}{parameter. \code{0 < prob <= 1}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The zero-truncated geometric distribution with \code{prob} \eqn{= p}
  has probability mass function
  \deqn{%
    p(x) = p (1-p)^{x - 1}}{%
    p(x) = p (1-p)^(x-1)}
  for \eqn{x = 1, 2, \ldots} and \eqn{0 < p < 1}, and
  \eqn{p(1) = 1} when \eqn{p = 1}.
  The cumulative distribution function is
  \deqn{P(x) = \frac{F(x) - F(0)}{1 - F(0)},}{%
    P(x) = [F(x) - F(0)]/[1 - F(0)],}
  where \eqn{F(x)} is the distribution function of the standard geometric.

  The mean is \eqn{1/p} and the variance is \eqn{(1-p)/p^2}.

  In the terminology of Klugman et al. (2012), the zero-truncated
  geometric is a member of the \eqn{(a, b, 1)} class of
  distributions with \eqn{a = 1-p} and \eqn{b = 0}.

  If an element of \code{x} is not integer, the result of
  \code{dztgeom} is zero, with a warning.

  The quantile is defined as the smallest value \eqn{x} such that
  \eqn{P(x) \ge p}, where \eqn{P} is the distribution function.
}
\value{
  \code{dztgeom} gives the (log) probability mass function,
  \code{pztgeom} gives the (log) distribution function,
  \code{qztgeom} gives the quantile function, and
  \code{rztgeom} generates random deviates.

  Invalid \code{prob} will result in return value \code{NaN}, with a
  warning.

  The length of the result is determined by \code{n} for
  \code{rztgeom}, and is the maximum of the lengths of the
  numerical arguments for the other functions.
}
\note{
  Functions \code{\{d,p,q\}ztgeom} use \code{\{d,p,q\}geom} for all but
  the trivial input values and \eqn{p(0)}.

  \code{rztgeom} uses the simple inversion algorithm suggested by
  Peter Dalgaard on the r-help mailing list on 1 May 2005 %
  (\url{https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html}).
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\seealso{
  \code{\link{dgeom}} for the geometric distribution.

  \code{\link{dztnbinom}} for the zero-truncated negative binomial, of
  which the zero-truncated geometric is a special case.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
p <- 1/(1 + 0.5)
dztgeom(c(1, 2, 3), prob = p)
dgeom(c(1, 2, 3), p)/pgeom(0, p, lower = FALSE) # same
dgeom(c(1, 2, 3) - 1, p)                        # same

pztgeom(1, prob = 1)        # point mass at 1

qztgeom(pztgeom(1:10, 0.3), 0.3)
}
\keyword{distribution}
