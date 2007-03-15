\name{discretize}
\alias{discretize}
\alias{discretise}
\title{Discretization of a Continuous Distribution}
\description{
  Compute a discrete probability mass function from a continuous
  cumulatve distribution function (cdf) with various methods.

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
  the probability mass at \eqn{x = a, a + h, \dots, b} in the
  discretized distribution and set \eqn{a =} \code{from} and \eqn{b =}
  \code{to}.
  
  Method \code{"upper"} is the forward difference of the cdf \eqn{F}:
  \deqn{p_x = F(x + h) - F(x).}{p[x] = F(x + h) - F(x).}

  Method \code{"lower"} is the backward difference of the cdf \eqn{F}:
  \deqn{p_x = F(x) - F(x - h).}{p[x] = F(x) - F(x - h).}

  Method \code{"rounding"} has the discretized cdf pass through the
  midpoints of the intervals \eqn{[x - h/2, x + h/2)}:
  \deqn{p_x = F(x + h/2) - F(x - h/2),}{p[x] = F(x + h/2) - F(x - h/2),}
  with \eqn{a \equiv a - h/2}{a := a - h/2} and \eqn{b \equiv b - h/2}{b
    := b + h/2}. The function assumes the cdf is continuous. Any
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
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\seealso{
  \code{\link{aggregateDist}}
}
\examples{
x <- seq(from = 0, to = qgamma(1 - 1E-6, 2, 1), by = 0.5)

## Upper and lower discretization
fxu <- discretize(pgamma(x, 2, 1), method = "upper",
            from = 0, to = tail(x, 1), step = 0.5)
fxl <- discretize(pgamma(x, 2, 1), method = "lower",
            from = 0, to = tail(x, 1), step = 0.5)
plot(stepfun(head(x, -1), diffinv(fxu)), do.points = FALSE)
plot(stepfun(x, diffinv(fxl)), do.points = FALSE, add = TRUE)
curve(pgamma(x, 2, 1), add = TRUE, col = "blue")

## Rounding (or midpoint) discretization
fx <- discretize(pgamma(x, 2, 1), method = "rounding",
           from = 0, to = tail(x, 1), step = 0.5)
plot(stepfun(x, diffinv(fx)), do.points = FALSE)
curve(pgamma(x, 2, 1), add = TRUE, col = "blue")

## First moment matching
fx <- discretize(pgamma(x, 2, 1), method = "unbiased",
           lev = levgamma(x, 2, 1), from = 0, to = tail(x, 1),
           step = 0.5)
plot(stepfun(x, diffinv(fx)), do.points = FALSE)
curve(pgamma(x, 2, 1), add = TRUE, col = "blue")
drop(crossprod(x, fx))		# expected value preserved
}
\keyword{distribution}
\keyword{models}
