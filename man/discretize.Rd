\name{discretize}
\alias{discretize}
\alias{discretise}
\title{Discretization of a Continuous Distribution}
\description{
  Compute a discrete probability mass function from a continuous
  cumulative distribution function (cdf) with various methods.

  \code{discretise} is an alias for \code{discretize}.
}
\usage{
discretize(cdf, from, to, step = 1,
           method = c("upper", "lower", "rounding", "unbiased"),
           lev, by = step, xlim = NULL)

discretise(cdf, from, to, step = 1,
           method = c("upper", "lower", "rounding", "unbiased"),
           lev, by = step, xlim = NULL)
}
\arguments{
  \item{cdf}{an expression written as a function of \code{x}, or
    alternatively the name of a function, giving the cdf to discretize.}
  \item{from, to}{the range over which the function will be discretized.}
  \item{step}{numeric; the discretization step (or span, or lag).}
  \item{method}{discretization method to use.}
  \item{lev}{an expression written as a function of \code{x}, or
    alternatively the name of a function, to compute the limited
    expected value of the distribution corresponding to
    \code{cdf}. Used only with the \code{"unbiased"} method.}
  \item{by}{an alias for \code{step}.}
  \item{xlim}{numeric of length 2; if specified, it serves as default
    for \code{c(from, to)}.}
}
\details{
  Usage is similar to \code{\link{curve}}.

  \code{discretize} returns the probability mass function (pmf) of the
  random variable obtained by discretization of the cdf specified in
  \code{cdf}.

  Let \eqn{F(x)} denote the cdf, \eqn{E[\min(X, x)]}{E[min(X, x)]]} the
  limited expected value at \eqn{x}, \eqn{h} the step, \eqn{p_x}{p[x]}
  the probability mass at \eqn{x} in the discretized distribution and
  set \eqn{a =} \code{from} and \eqn{b =} \code{to}.

  Method \code{"upper"} is the forward difference of the cdf \eqn{F}:
  \deqn{p_x = F(x + h) - F(x)}{p[x] = F(x + h) - F(x)}
  for \eqn{x = a, a + h, \dots, b - step}.

  Method \code{"lower"} is the backward difference of the cdf \eqn{F}:
  \deqn{p_x = F(x) - F(x - h)}{p[x] = F(x) - F(x - h)} for \eqn{x = a +
    h, \dots, b} and \eqn{p_a = F(a)}{p[a] = F(a)}.

  Method \code{"rounding"} has the true cdf pass through the
  midpoints of the intervals \eqn{[x - h/2, x + h/2)}:
  \deqn{p_x = F(x + h/2) - F(x - h/2)}{p[x] = F(x + h/2) - F(x - h/2)}
  for \eqn{x = a + h, \dots, b - step} and \eqn{p_a = F(a + h/2)}{p[a] =
  F(a + h/2)}. The function assumes the cdf is continuous. Any
  adjusment necessary for discrete distributions can be done via
  \code{cdf}.

  Method \code{"unbiased"} matches the first moment of the discretized
  and the true distributions. The probabilities are as follows:
  \deqn{p_a = \frac{E[\min(X, a)] - E[\min(X, a + h)]}{h} + 1 - F(a)}{%
    p[a] = (E[min(X, a)] - E[min(X, a + h)])/h + 1 - F(a)}
  \deqn{p_x = \frac{2 E[\min(X, x)] - E[\min(X, x - h)] - E[\min(X, x +
      h)]}{h}, \quad a < x < b}{%
    p[x] = (2 E[min(X, x)] - E[min(X, x - h)] - E[min(X, x + h)])/h, a < x < b}
  \deqn{p_b = \frac{E[\min(X, b)] - E[\min(X, b - h)]}{h} - 1 + F(b),}{%
    p[b] = (E[min(X, b)] - E[min(X, b - h)])/h - 1 + F(b).}
}
\value{
  A numeric vector of probabilities suitable for use in
  \code{\link{aggregateDist}}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{aggregateDist}}
}
\examples{
x <- seq(0, 5, 0.5)

op <- par(mfrow = c(1, 1), col = "black")

## Upper and lower discretization
fu <- discretize(pgamma(x, 1), method = "upper",
                 from = 0, to = 5, step = 0.5)
fl <- discretize(pgamma(x, 1), method = "lower",
                 from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fu)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(x, diffinv(fl)), pch = 19, add = TRUE)
par(col = "black")

## Rounding (or midpoint) discretization
fr <- discretize(pgamma(x, 1), method = "rounding",
                 from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fr)), pch = 19, add = TRUE)
par(col = "black")

## First moment matching
fb <- discretize(pgamma(x, 1), method = "unbiased",
                 lev = levgamma(x, 1), from = 0, to = 5, step = 0.5)
curve(pgamma(x, 1), xlim = c(0, 5))
par(col = "blue")
plot(stepfun(x, diffinv(fb)), pch = 19, add = TRUE)

par(op)
}
\keyword{distribution}
\keyword{models}
